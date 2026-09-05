#include "fabm_driver.h"

!-----------------------------------------------------------------------
! Seagrass (eelgrass, Zostera marina) module for ERSEM
!
! Named `seagrass`, not `macrophyte`: macrophyte conventionally spans
! seagrasses (rooted vascular plants, porewater uptake through roots)
! and macroalgae (holdfast, no roots, no porewater pathway), which are
! always modelled separately. This module is a rooted vascular plant.
! Renamed 2026-08-19; see nippon-steel/docs/21.
!
! Design: nippon-steel/docs/10-macrophyte-module-design.md and docs/21 (2026-08-14).
! Bottom-attached plant with above-ground (AG) and below-ground (BG)
! structural biomass plus a rhizome non-structural-carbohydrate reserve
! (NSC). Replaces the suspended primary-producer proxy for seagrass.
!
! v1 process set:
!   - gross production: p_max * CTMI(T) * f_light(canopy) * f_quota * AGc
!   - canopy self-shading on near-bottom PAR (LAI from AG carbon)
!   - leaf uptake of NO3/NH4/PO4 from the water, root uptake from the
!     porewater pools (layers 1-2), NH4 preferred, Droop-type quota
!   - alkalinity bookkeeping: +1 eq per mol NO3, -1 per mol NH4, +1 per
!     mol PO4 taken up; leaf terms on pelagic TA, root terms on the
!     benthic alkalinity pool (benTA -> G5 layer 1)
!   - AG <-> NSC translocation (storage of excess production, dark/low
!     light remobilisation); BG structural growth paid from NSC
!   - respiration: AG (+NSC maintenance) against pelagic O2 with a
!     Monod guard; BG against benthic layer-1 O2 (G2) with a Monod
!     guard and extra anoxic mortality
!   - mortality/sloughing: AG split between pelagic POM (R6) and the
!     plant detritus layer (seagrass_Q6); BG to the detritus layer
!
! v1 simplifications (documented in the design doc, section 6):
!   - all respired CO2 is returned to pelagic DIC (no porewater DIC)
!   - no epiphytes, no shoot demography, no root oxygen loss
!   - PAR received is the host-provided bottom PAR; in the 0-D box this
!     is the box PAR after background extinction
!
! Deactivation contract: if no instance of this model appears in
! fabm.yaml, nothing is registered and results are identical to a build
! without this module (same convention as benthic_cao with iswCaO=0).
!-----------------------------------------------------------------------

module ersem_seagrass

   use fabm_types
   use ersem_shared

   implicit none

   private

   type,extends(type_base_model),public :: type_ersem_seagrass
      ! Own bottom state variables
      type (type_bottom_state_variable_id) :: id_AGc, id_AGn, id_AGp
      type (type_bottom_state_variable_id) :: id_BGc, id_BGn, id_BGp
      type (type_bottom_state_variable_id) :: id_NSCc

      ! Pelagic dependencies (bottom cell)
      type (type_state_variable_id) :: id_O2o, id_O3c, id_TA
      type (type_state_variable_id) :: id_N1p, id_N3n, id_N4n
      type (type_state_variable_id) :: id_R6c, id_R6n, id_R6p

      ! Benthic dependencies
      type (type_bottom_state_variable_id) :: id_K1p1, id_K1p2
      type (type_bottom_state_variable_id) :: id_K3n1, id_K3n2
      type (type_bottom_state_variable_id) :: id_K4n1, id_K4n2
      type (type_bottom_state_variable_id) :: id_G2o
      type (type_bottom_state_variable_id) :: id_benTA
      type (type_bottom_state_variable_id) :: id_Q6c, id_Q6n, id_Q6p

      ! External grazers (optional grazing closure, docs/13 2026-08-15)
      type (type_horizontal_dependency_id) :: id_gr1c, id_gr2c

      ! Environment
      type (type_dependency_id) :: id_ETW, id_par

      ! Diagnostics
      type (type_horizontal_diagnostic_variable_id) :: id_gpp, id_npp
      type (type_horizontal_diagnostic_variable_id) :: id_resp, id_fT, id_fI
      type (type_horizontal_diagnostic_variable_id) :: id_graz
      ! Nitrogen exchange terms, so a nitrogen budget can be CLOSED from
      ! output alone (nippon-steel docs/113 s17): the mat's net stock
      ! change is only a floor on its gross uptake, because it loses
      ! nitrogen to mortality and grazing at the same time, and the
      ! water-column nitrate budget could not be closed without these.
      type (type_horizontal_diagnostic_variable_id) :: id_uN3l, id_uN4l
      type (type_horizontal_diagnostic_variable_id) :: id_relN4, id_uN3r, id_uN4r

      ! Parameters
      real(rk) :: p_max, alpha, a_lai, k_can
      real(rk) :: Tmin, Topt, Tmax, ctmi_a, ctmi_b
      real(rk) :: qn_min, qn_max, qp_min, qp_max
      real(rk) :: qn_bg, qp_bg
      real(rk) :: vmax_n, vmax_p, hN3, hN4, hP, hKn, hKp, f_root
      ! docs/21 s10.2 + review 2026-09-04 s1.3. isw_nupt = 0 keeps the
      ! legacy AMMONIUM-PRIORITY rule (nitrate receives only what ammonium
      ! leaves), which forces nitrate uptake UP where ammonium is scarce.
      ! isw_nupt = 1 computes substrate-limited POTENTIALS for all four
      ! source/form combinations and caps their sum by the one quota
      ! demand, after ersem/primary_producer's own reviewed structure.
      integer  :: isw_nupt
      real(rk) :: psiN4, psiN3
      real(rk) :: srs_ag, srs_bg, r_nsc, pu_ra, hO2, hG2o
      real(rk) :: pq, rq_o2c
      real(rk) :: tau_store, tau_mob, k_bg
      real(rk) :: rs_target, k_alloc, q_nsc
      real(rk) :: sd_ag, sd_bg, sd_anx, nsc_starve, sd_starve, f_pel
      real(rk) :: g_max, h_ag, pe_gr
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self, configunit)
      class (type_ersem_seagrass), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      real(rk) :: a, b, ab2

      ! Set time unit to d-1 (ERSEM convention)
      self%dt = 86400._rk

      ! --- Parameters (defaults follow docs/10 section 5) -----------------
      call self%get_parameter(self%p_max, 'p_max', '1/d', &
         'maximum gross production at Topt', default=0.12_rk, minimum=0.0_rk)
      call self%get_parameter(self%alpha, 'alpha', '(W/m^2)^-1', &
         'initial slope of the P-I curve (tanh)', default=0.06_rk, minimum=0.0_rk)
      call self%get_parameter(self%a_lai, 'a_lai', 'm^2/mg C', &
         'leaf area per unit AG carbon', default=4.0e-5_rk, minimum=0.0_rk)
      call self%get_parameter(self%k_can, 'k_can', '-', &
         'canopy attenuation coefficient per unit LAI', default=0.7_rk, minimum=0.0_rk)

      call self%get_parameter(self%Tmin, 'Tmin', 'degrees_Celsius', &
         'CTMI minimum temperature for growth', default=2.0_rk)
      call self%get_parameter(self%Topt, 'Topt', 'degrees_Celsius', &
         'CTMI optimal temperature for growth', default=18.0_rk)
      call self%get_parameter(self%Tmax, 'Tmax', 'degrees_Celsius', &
         'CTMI maximum temperature for growth', default=30.0_rk)
      if (self%Tmin >= self%Topt) call self%fatal_error('initialize','CTMI requires Tmin < Topt')
      if (self%Topt >= self%Tmax) call self%fatal_error('initialize','CTMI requires Topt < Tmax')
      a = self%Topt - self%Tmin
      b = self%Topt - self%Tmax
      ab2 = (a * b)**2
      self%ctmi_a = -(a + b) / ab2
      self%ctmi_b = (a * b + (a + b) * self%Topt) / ab2

      call self%get_parameter(self%qn_min, 'qn_min', 'mmol N/mg C', &
         'minimum AG nitrogen quota', default=0.003_rk, minimum=0.0_rk)
      call self%get_parameter(self%qn_max, 'qn_max', 'mmol N/mg C', &
         'maximum AG nitrogen quota', default=0.008_rk, minimum=0.0_rk)
      call self%get_parameter(self%qp_min, 'qp_min', 'mmol P/mg C', &
         'minimum AG phosphorus quota', default=8.0e-5_rk, minimum=0.0_rk)
      call self%get_parameter(self%qp_max, 'qp_max', 'mmol P/mg C', &
         'maximum AG phosphorus quota', default=2.5e-4_rk, minimum=0.0_rk)
      call self%get_parameter(self%qn_bg, 'qn_bg', 'mmol N/mg C', &
         'fixed BG nitrogen quota', default=0.0035_rk, minimum=0.0_rk)
      call self%get_parameter(self%qp_bg, 'qp_bg', 'mmol P/mg C', &
         'fixed BG phosphorus quota', default=1.0e-4_rk, minimum=0.0_rk)

      call self%get_parameter(self%vmax_n, 'vmax_n', 'mmol N/mg C/d', &
         'maximum N uptake per unit AG carbon', default=6.0e-4_rk, minimum=0.0_rk)
      call self%get_parameter(self%vmax_p, 'vmax_p', 'mmol P/mg C/d', &
         'maximum P uptake per unit AG carbon', default=2.0e-5_rk, minimum=0.0_rk)
      call self%get_parameter(self%hN3, 'hN3', 'mmol N/m^3', &
         'half-saturation for leaf NO3 uptake', default=1.0_rk, minimum=1.0e-6_rk)
      call self%get_parameter(self%hN4, 'hN4', 'mmol N/m^3', &
         'half-saturation for leaf NH4 uptake', default=0.5_rk, minimum=1.0e-6_rk)
      call self%get_parameter(self%hP, 'hP', 'mmol P/m^3', &
         'half-saturation for leaf PO4 uptake', default=0.1_rk, minimum=1.0e-6_rk)
      call self%get_parameter(self%hKn, 'hKn', 'mmol N/m^2', &
         'half-saturation for root N uptake (layer amount)', default=20.0_rk, minimum=1.0e-6_rk)
      call self%get_parameter(self%hKp, 'hKp', 'mmol P/m^2', &
         'half-saturation for root P uptake (layer amount)', default=5.0_rk, minimum=1.0e-6_rk)
      call self%get_parameter(self%f_root, 'f_root', '-', &
         'root fraction of uptake capacity', default=0.5_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%isw_nupt, 'isw_nupt', '', &
         'nitrogen uptake rule (0: ammonium priority [legacy], 1: potentials capped by quota demand)', &
         default=0, minimum=0, maximum=1)
      call self%get_parameter(self%psiN4, 'psiN4', '-', &
         'ammonium affinity weight (isw_nupt=1 only)', default=1.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%psiN3, 'psiN3', '-', &
         'nitrate affinity weight (isw_nupt=1 only)', default=1.0_rk, minimum=0.0_rk)

      call self%get_parameter(self%srs_ag, 'srs_ag', '1/d', &
         'AG basal respiration at Topt', default=0.015_rk, minimum=0.0_rk)
      call self%get_parameter(self%srs_bg, 'srs_bg', '1/d', &
         'BG basal respiration at Topt', default=0.006_rk, minimum=0.0_rk)
      call self%get_parameter(self%r_nsc, 'r_nsc', '1/d', &
         'NSC maintenance respiration', default=0.002_rk, minimum=0.0_rk)
      call self%get_parameter(self%pu_ra, 'pu_ra', '-', &
         'respired fraction of gross production', default=0.25_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%hO2, 'hO2', 'mmol O_2/m^3', &
         'Monod half-saturation for AG respiration O2', default=15.625_rk, minimum=0.0_rk)
      call self%get_parameter(self%hG2o, 'hG2o', 'mmol O_2/m^2', &
         'Monod half-saturation for BG respiration on benthic O2', default=1.0_rk, minimum=0.0_rk)

      ! --- Oxygen stoichiometry (nippon-steel docs/51 s4, 2026-08-25) ------
      ! Until this build the module carried PQ = RQ = 1 as a COMMENT on the
      ! exchange line and no parameter: one mol O_2 per mol C fixed and per
      ! mol C respired. ERSEM's pelagic producers expose the same two
      ! constants as uB1c_O2 / urB1_O2 (mmol O_2/mg C); here they are mol/mol
      ! so that the defaults are exactly 1 and the pre-existing arithmetic is
      ! reproduced term by term (a product with 1.0_rk is exact and the
      ! order of the subtractions below is unchanged). Literature: PQ 1.0-1.4
      ! for a mixed community (nitrate-based growth at the top); RQ as
      ! CO_2:O_2 0.8-1.2, i.e. rq_o2c = 1/RQ in 0.83-1.25.
      call self%get_parameter(self%pq, 'pq', 'mol O_2/mol C', &
         'photosynthetic quotient: O_2 evolved per C fixed', default=1.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%rq_o2c, 'rq_o2c', 'mol O_2/mol C', &
         'O_2 consumed per C respired (reciprocal of the respiratory quotient CO_2:O_2)', &
         default=1.0_rk, minimum=0.0_rk)

      call self%get_parameter(self%tau_store, 'tau_store', '-', &
         'stored fraction of positive net AG production', default=0.3_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%tau_mob, 'tau_mob', '1/d', &
         'NSC remobilisation rate under light limitation', default=0.05_rk, minimum=0.0_rk)
      call self%get_parameter(self%k_bg, 'k_bg', '1/d', &
         'BG structural growth rate from NSC', default=0.02_rk, minimum=0.0_rk)
      ! v2 allocation control (2026-08-14): the v1 balance let the BG pool
      ! collapse over a full year (report_year.md). BG growth is now driven by
      ! the root:shoot deficit, fed both from positive net AG production
      ! (k_alloc) and from the reserve (k_bg), and storage stops once the
      ! reserve reaches its target fraction of BG structure.
      call self%get_parameter(self%rs_target, 'rs_target', '-', &
         'target BG:AG carbon ratio for allocation control', &
         default=0.75_rk, minimum=0.01_rk)
      call self%get_parameter(self%k_alloc, 'k_alloc', '-', &
         'maximum fraction of positive net AG production allocated to BG growth', &
         default=0.35_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%q_nsc, 'q_nsc', '-', &
         'target NSC reserve as a fraction of BG carbon', &
         default=0.15_rk, minimum=0.0_rk)

      call self%get_parameter(self%sd_ag, 'sd_ag', '1/d', &
         'AG background sloughing rate', default=0.004_rk, minimum=0.0_rk)
      call self%get_parameter(self%sd_bg, 'sd_bg', '1/d', &
         'BG background mortality rate', default=0.002_rk, minimum=0.0_rk)
      call self%get_parameter(self%sd_anx, 'sd_anx', '1/d', &
         'extra BG mortality under benthic anoxia', default=0.01_rk, minimum=0.0_rk)
      call self%get_parameter(self%nsc_starve, 'nsc_starve', '-', &
         'NSC:BG carbon ratio below which starvation mortality starts', &
         default=0.02_rk, minimum=0.0_rk)
      call self%get_parameter(self%sd_starve, 'sd_starve', '1/d', &
         'extra AG mortality under NSC starvation', default=0.02_rk, minimum=0.0_rk)
      call self%get_parameter(self%f_pel, 'f_pel', '-', &
         'fraction of AG sloughing routed to pelagic POM (R6)', &
         default=0.5_rk, minimum=0.0_rk, maximum=1.0_rk)

      ! --- Optional grazing closure on AG (docs/13, 2026-08-15) ----------
      ! Type-III (sigmoid) loss of AG carbon to external benthic grazers:
      !   F_gr = g_max * (grazer1c + grazer2c) * AGc^2 / (AGc^2 + h_ag^2)
      ! The sigmoid gives a low-biomass refuge that the density-independent
      ! sd_ag lacks (the year-run collapse of report_year_v6.md). A
      ! fraction pe_gr, limited by the plant's pelagic-O2 Monod factor, is
      ! respired by the grazers against pelagic O2/DIC with the module's
      ! standard nutrient and alkalinity return; the remainder is egested
      ! to the plant detritus pool at AG quota. Default g_max = 0 keeps the
      ! module bit-identical to the pre-grazing build (the code path is
      ! guarded; only the always-written 'graz' diagnostic is new).
      call self%get_parameter(self%g_max, 'g_max', '1/d', &
         'maximum grazing ration per unit grazer carbon', &
         default=0.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%h_ag, 'h_ag', 'mg C/m^2', &
         'AG carbon at the type-III grazing half-saturation', &
         default=500.0_rk, minimum=1.0e-6_rk)
      call self%get_parameter(self%pe_gr, 'pe_gr', '-', &
         'fraction of grazed carbon respired by the grazers', &
         default=0.4_rk, minimum=0.0_rk, maximum=1.0_rk)

      ! --- Own state variables (initial values must be set in fabm.yaml) --
      call self%register_state_variable(self%id_AGc, 'AGc', 'mg C/m^2', 'above-ground carbon', minimum=0.0_rk)
      call self%register_state_variable(self%id_AGn, 'AGn', 'mmol N/m^2', 'above-ground nitrogen', minimum=0.0_rk)
      call self%register_state_variable(self%id_AGp, 'AGp', 'mmol P/m^2', 'above-ground phosphorus', minimum=0.0_rk)
      call self%register_state_variable(self%id_BGc, 'BGc', 'mg C/m^2', 'below-ground carbon', minimum=0.0_rk)
      call self%register_state_variable(self%id_BGn, 'BGn', 'mmol N/m^2', 'below-ground nitrogen', minimum=0.0_rk)
      call self%register_state_variable(self%id_BGp, 'BGp', 'mmol P/m^2', 'below-ground phosphorus', minimum=0.0_rk)
      call self%register_state_variable(self%id_NSCc, 'NSCc', 'mg C/m^2', 'rhizome carbohydrate reserve', minimum=0.0_rk)

      ! Conservation bookkeeping
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_AGc, scale_factor=1._rk/CMass)
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_BGc, scale_factor=1._rk/CMass)
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_NSCc, scale_factor=1._rk/CMass)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_AGn)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_BGn)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_AGp)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_BGp)

      ! --- Dependencies ---------------------------------------------------
      call self%register_state_dependency(self%id_O2o, 'O2o', 'mmol O_2/m^3', 'pelagic oxygen')
      call self%register_state_dependency(self%id_O3c, 'O3c', 'mmol C/m^3', 'pelagic dissolved inorganic carbon')
      call self%register_state_dependency(self%id_TA, standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_state_dependency(self%id_N1p, 'N1p', 'mmol P/m^3', 'pelagic phosphate')
      call self%register_state_dependency(self%id_N3n, 'N3n', 'mmol N/m^3', 'pelagic nitrate')
      call self%register_state_dependency(self%id_N4n, 'N4n', 'mmol N/m^3', 'pelagic ammonium')
      call self%register_state_dependency(self%id_R6c, 'R6c', 'mg C/m^3', 'pelagic POM carbon')
      call self%register_state_dependency(self%id_R6n, 'R6n', 'mmol N/m^3', 'pelagic POM nitrogen')
      call self%register_state_dependency(self%id_R6p, 'R6p', 'mmol P/m^3', 'pelagic POM phosphorus')

      call self%register_state_dependency(self%id_K1p1, 'K1p1', 'mmol P/m^2', 'porewater phosphate, layer 1')
      call self%register_state_dependency(self%id_K1p2, 'K1p2', 'mmol P/m^2', 'porewater phosphate, layer 2')
      call self%register_state_dependency(self%id_K3n1, 'K3n1', 'mmol N/m^2', 'porewater nitrate, layer 1')
      call self%register_state_dependency(self%id_K3n2, 'K3n2', 'mmol N/m^2', 'porewater nitrate, layer 2')
      call self%register_state_dependency(self%id_K4n1, 'K4n1', 'mmol N/m^2', 'porewater ammonium, layer 1')
      call self%register_state_dependency(self%id_K4n2, 'K4n2', 'mmol N/m^2', 'porewater ammonium, layer 2')
      call self%register_state_dependency(self%id_G2o, 'G2o', 'mmol O_2/m^2', 'benthic oxygen, layer 1')
      call self%register_state_dependency(self%id_benTA, 'benTA', 'mEq/m^2', 'benthic alkalinity, layer 1')
      call self%register_state_dependency(self%id_Q6c, 'Q6c', 'mg C/m^2', 'plant detritus carbon')
      call self%register_state_dependency(self%id_Q6n, 'Q6n', 'mmol N/m^2', 'plant detritus nitrogen')
      call self%register_state_dependency(self%id_Q6p, 'Q6p', 'mmol P/m^2', 'plant detritus phosphorus')

      ! Grazer carbon (mg C/m^2); default zero_hz = no grazers. Couple to
      ! e.g. Y2/c and Y4/c in fabm.yaml to activate together with g_max.
      call self%register_dependency(self%id_gr1c, 'grazer1c', 'mg C/m^2', 'carbon of benthic grazer 1')
      call self%register_dependency(self%id_gr2c, 'grazer2c', 'mg C/m^2', 'carbon of benthic grazer 2')
      call self%request_coupling(self%id_gr1c, 'zero_hz')
      call self%request_coupling(self%id_gr2c, 'zero_hz')

      call self%register_dependency(self%id_ETW, standard_variables%temperature)
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

      ! --- Diagnostics ----------------------------------------------------
      call self%register_diagnostic_variable(self%id_gpp, 'GPP', 'mg C/m^2/d', &
         'gross primary production', domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_npp, 'NPP', 'mg C/m^2/d', &
         'net primary production', domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_resp, 'resp', 'mg C/m^2/d', &
         'total plant respiration', domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fT, 'fT', '-', &
         'CTMI temperature factor', domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_fI, 'fI', '-', &
         'canopy light factor', domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_graz, 'graz', 'mg C/m^2/d', &
         'grazing loss of AG carbon', domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_uN3l, 'upt_N3', 'mmol N/m^2/d', &
         'GROSS nitrate uptake from the water column by leaves', &
         domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_uN4l, 'upt_N4', 'mmol N/m^2/d', &
         'GROSS ammonium uptake from the water column by leaves', &
         domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_relN4, 'rel_N4', 'mmol N/m^2/d', &
         'respiratory ammonium return to the water column', &
         domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_uN3r, 'upt_N3_root', 'mmol N/m^2/d', &
         'nitrate uptake from porewater by roots', &
         domain=domain_bottom, source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_uN4r, 'upt_N4_root', 'mmol N/m^2/d', &
         'ammonium uptake from porewater by roots', &
         domain=domain_bottom, source=source_do_bottom)

   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_seagrass), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: AGc, AGn, AGp, BGc, BGn, BGp, NSCc
      real(rk) :: O2o, N1p, N3n, N4n
      real(rk) :: K1p1, K1p2, K3n1, K3n2, K4n1, K4n2, G2o
      real(rk) :: ETW, par
      real(rk) :: eT, lai, I_can, eI, qn, qp, eQ
      real(rk) :: tau
      real(rk) :: Pg, Ra_act, Ra_bas, Rn, Rb, fO2ag, fO2bg
      real(rk) :: Pnet_ag, T_st, T_mb, G_bg, G_bg_c
      real(rk) :: relE, nsc_gap, G_alloc, G_res, gscale
      real(rk) :: cap_n, cap_p, upt_leaf, upt_root
      real(rk) :: jN4_leaf, jN3_leaf, jP_leaf
      real(rk) :: jN4_root, jN3_root, jP_root, wsum, pot_n, nscale
      real(rk) :: M_ag, M_bg, starve
      real(rk) :: dAGc, dAGn, dAGp, dBGc, dBGn, dBGp, dNSC
      real(rk) :: gr1c, gr2c, F_gr, resp_gr, eges_gr

      _HORIZONTAL_LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_AGc, AGc)
         _GET_HORIZONTAL_(self%id_AGn, AGn)
         _GET_HORIZONTAL_(self%id_AGp, AGp)
         _GET_HORIZONTAL_(self%id_BGc, BGc)
         _GET_HORIZONTAL_(self%id_BGn, BGn)
         _GET_HORIZONTAL_(self%id_BGp, BGp)
         _GET_HORIZONTAL_(self%id_NSCc, NSCc)
         _GET_HORIZONTAL_(self%id_K1p1, K1p1)
         _GET_HORIZONTAL_(self%id_K1p2, K1p2)
         _GET_HORIZONTAL_(self%id_K3n1, K3n1)
         _GET_HORIZONTAL_(self%id_K3n2, K3n2)
         _GET_HORIZONTAL_(self%id_K4n1, K4n1)
         _GET_HORIZONTAL_(self%id_K4n2, K4n2)
         _GET_HORIZONTAL_(self%id_G2o, G2o)
         _GET_(self%id_O2o, O2o)
         _GET_(self%id_N1p, N1p)
         _GET_(self%id_N3n, N3n)
         _GET_(self%id_N4n, N4n)
         _GET_(self%id_ETW, ETW)
         _GET_(self%id_par, par)

         ! --- Environmental factors ---------------------------------------
         if (ETW <= self%Tmin .or. ETW >= self%Tmax) then
            eT = 0.0_rk
         else
            eT = (ETW - self%Tmin) * (ETW - self%Tmax) * (self%ctmi_a * ETW + self%ctmi_b)
            eT = max(0.0_rk, min(1.0_rk, eT))
         end if

         ! Optical depth of the canopy. k_can and a_lai enter ONLY as their
         ! product, so the guard must test the PRODUCT: testing `lai` alone
         ! divides by zero whenever k_can = 0 with a non-trivial a_lai, which
         ! the parameter reader allows (minimum=0.0). Found 2026-08-22.
         lai = self%a_lai * AGc
         tau = self%k_can * lai
         if (tau > 1.0e-8_rk) then
            I_can = max(0.0_rk, par) * (1.0_rk - exp(-tau)) / tau
         else
            I_can = max(0.0_rk, par)
         end if
         eI = tanh(self%alpha * I_can)

         qn = AGn / max(AGc, 1.0e-8_rk)
         qp = AGp / max(AGc, 1.0e-8_rk)
         eQ = min(max(0.0_rk, (qn - self%qn_min) / (self%qn_max - self%qn_min)), &
                  max(0.0_rk, (qp - self%qp_min) / (self%qp_max - self%qp_min)))
         eQ = min(1.0_rk, eQ)

         ! --- Production and respiration (mg C/m^2/d) ---------------------
         Pg = self%p_max * eT * eI * eQ * AGc
         fO2ag = max(0.0_rk, O2o) / (max(0.0_rk, O2o) + self%hO2)
         fO2bg = max(0.0_rk, G2o) / (max(0.0_rk, G2o) + self%hG2o)
         Ra_act = self%pu_ra * Pg
         Ra_bas = self%srs_ag * eT * AGc * fO2ag
         Rn = self%r_nsc * eT * NSCc * fO2ag
         Rb = self%srs_bg * eT * BGc * fO2bg

         ! --- Translocation and BG growth (v2 allocation control) ---------
         Pnet_ag = Pg - Ra_act - Ra_bas
         ! Root:shoot deficit: 1 when BG is absent, 0 at/above the target ratio
         relE = min(1.0_rk, max(0.0_rk, &
                1.0_rk - (BGc / max(AGc, 1.0e-8_rk)) / self%rs_target))
         ! Storage: bank surplus production until the reserve reaches its
         ! target fraction of BG structure
         nsc_gap = max(0.0_rk, 1.0_rk - NSCc / max(self%q_nsc * BGc, 1.0e-8_rk))
         T_st = self%tau_store * max(0.0_rk, Pnet_ag) * nsc_gap
         T_mb = self%tau_mob * NSCc * max(0.0_rk, 1.0_rk - eI)
         ! BG growth, driven by the deficit: direct allocation of net
         ! production (carbon from AG) plus reserve-fed growth (carbon from NSC)
         G_alloc = relE * self%k_alloc * max(0.0_rk, Pnet_ag)
         G_res = relE * self%k_bg * eT * NSCc
         G_bg = G_alloc + G_res
         ! BG structural growth needs N and P from the AG pools at fixed quota;
         ! scale it down when the AG pools cannot pay.
         G_bg_c = min(G_bg, &
                      0.5_rk * AGn / max(self%qn_bg, 1.0e-12_rk), &
                      0.5_rk * AGp / max(self%qp_bg, 1.0e-12_rk))
         gscale = G_bg_c / max(G_bg, 1.0e-12_rk)
         G_alloc = G_alloc * gscale
         G_res = G_res * gscale

         ! --- Nutrient uptake (mmol/m^2/d) --------------------------------
         cap_n = self%vmax_n * eT * AGc * max(0.0_rk, 1.0_rk - qn / self%qn_max)
         cap_p = self%vmax_p * eT * AGc * max(0.0_rk, 1.0_rk - qp / self%qp_max)

         upt_leaf = (1.0_rk - self%f_root) * cap_n
         upt_root = self%f_root * cap_n
         wsum = K4n1 + K4n2
         if (self%isw_nupt == 0) then
            ! LEGACY: ammonium first, nitrate takes the remainder.
            jN4_leaf = upt_leaf * N4n / (N4n + self%hN4)
            jN3_leaf = max(0.0_rk, upt_leaf - jN4_leaf) * N3n / (N3n + self%hN3)
            jN4_root = upt_root * wsum / (wsum + self%hKn)
            jN3_root = max(0.0_rk, upt_root - jN4_root) * (K3n1 + K3n2) / (K3n1 + K3n2 + self%hKn)
         else
            ! Substrate-limited potentials, then ONE quota cap. Each flux is
            ! zero when its own form is absent, the total tends to zero with
            ! total DIN, and no form receives a remainder it did not earn.
            jN4_leaf = upt_leaf * self%psiN4 * N4n / (N4n + self%hN4)
            jN3_leaf = upt_leaf * self%psiN3 * N3n / (N3n + self%hN3)
            jN4_root = upt_root * self%psiN4 * wsum / (wsum + self%hKn)
            jN3_root = upt_root * self%psiN3 * (K3n1 + K3n2) / (K3n1 + K3n2 + self%hKn)
            pot_n = jN4_leaf + jN3_leaf + jN4_root + jN3_root
            if (pot_n > 0.0_rk) then
               nscale = min(1.0_rk, cap_n / pot_n)
            else
               nscale = 0.0_rk
            end if
            jN4_leaf = jN4_leaf * nscale
            jN3_leaf = jN3_leaf * nscale
            jN4_root = jN4_root * nscale
            jN3_root = jN3_root * nscale
         end if
         jP_leaf = (1.0_rk - self%f_root) * cap_p * N1p / (N1p + self%hP)
         jP_root = self%f_root * cap_p * (K1p1 + K1p2) / (K1p1 + K1p2 + self%hKp)

         ! --- Mortality ----------------------------------------------------
         starve = 0.0_rk
         if (NSCc < self%nsc_starve * max(BGc, 1.0e-8_rk)) starve = self%sd_starve
         M_ag = (self%sd_ag + starve) * AGc
         M_bg = (self%sd_bg + self%sd_anx * (1.0_rk - fO2bg)) * BGc

         ! --- State ODEs (per day) ----------------------------------------
         dAGc = Pg - Ra_act - Ra_bas - T_st + T_mb - M_ag - G_alloc
         dNSC = T_st - T_mb - G_res - Rn
         dBGc = G_bg_c - Rb - M_bg
         dAGn = jN4_leaf + jN3_leaf + jN4_root + jN3_root &
                - self%qn_bg * G_bg_c - qn * M_ag
         dAGp = jP_leaf + jP_root - self%qp_bg * G_bg_c - qp * M_ag
         dBGn = self%qn_bg * G_bg_c - BGn / max(BGc, 1.0e-8_rk) * M_bg
         dBGp = self%qp_bg * G_bg_c - BGp / max(BGc, 1.0e-8_rk) * M_bg
         ! Respiration releases the associated N and P to the water (AG) as NH4
         ! and PO4 at the current quota, keeping the quota from drifting when
         ! carbon is respired.
         dAGn = dAGn - qn * (Ra_act + Ra_bas)
         dAGp = dAGp - qp * (Ra_act + Ra_bas)
         dBGn = dBGn - BGn / max(BGc, 1.0e-8_rk) * Rb
         dBGp = dBGp - BGp / max(BGc, 1.0e-8_rk) * Rb

         _SET_BOTTOM_ODE_(self%id_AGc, dAGc)
         _SET_BOTTOM_ODE_(self%id_AGn, dAGn)
         _SET_BOTTOM_ODE_(self%id_AGp, dAGp)
         _SET_BOTTOM_ODE_(self%id_BGc, dBGc)
         _SET_BOTTOM_ODE_(self%id_BGn, dBGn)
         _SET_BOTTOM_ODE_(self%id_BGp, dBGp)
         _SET_BOTTOM_ODE_(self%id_NSCc, dNSC)

         ! --- Exchanges with the water column -----------------------------
         ! Carbon and oxygen. Carbon is carbon; the oxygen carries the two
         ! quotients (pq, rq_o2c), both 1 by default, which reproduces the
         ! former "PQ = 1" line term by term.
         _SET_BOTTOM_EXCHANGE_(self%id_O3c, (-Pg + Ra_act + Ra_bas + Rn + Rb) / CMass)
         _SET_BOTTOM_EXCHANGE_(self%id_O2o, (self%pq * Pg - self%rq_o2c * Ra_act &
                                             - self%rq_o2c * Ra_bas - self%rq_o2c * Rn) / CMass)
         ! BG respiration draws benthic layer-1 oxygen instead
         _SET_BOTTOM_ODE_(self%id_G2o, -self%rq_o2c * Rb / CMass)

         ! Leaf nutrient uptake and respiratory return
         _SET_BOTTOM_EXCHANGE_(self%id_N4n, -jN4_leaf + qn * (Ra_act + Ra_bas) + BGn / max(BGc, 1.0e-8_rk) * Rb)
         _SET_BOTTOM_EXCHANGE_(self%id_N3n, -jN3_leaf)
         _SET_BOTTOM_EXCHANGE_(self%id_N1p, -jP_leaf + qp * (Ra_act + Ra_bas) + BGp / max(BGc, 1.0e-8_rk) * Rb)

         ! Root uptake from the porewater pools (split by availability)
         wsum = max(K4n1 + K4n2, 1.0e-8_rk)
         _SET_BOTTOM_ODE_(self%id_K4n1, -jN4_root * K4n1 / wsum)
         _SET_BOTTOM_ODE_(self%id_K4n2, -jN4_root * K4n2 / wsum)
         wsum = max(K3n1 + K3n2, 1.0e-8_rk)
         _SET_BOTTOM_ODE_(self%id_K3n1, -jN3_root * K3n1 / wsum)
         _SET_BOTTOM_ODE_(self%id_K3n2, -jN3_root * K3n2 / wsum)
         wsum = max(K1p1 + K1p2, 1.0e-8_rk)
         _SET_BOTTOM_ODE_(self%id_K1p1, -jP_root * K1p1 / wsum)
         _SET_BOTTOM_ODE_(self%id_K1p2, -jP_root * K1p2 / wsum)

         ! Alkalinity bookkeeping: +1 per NO3, -1 per NH4, +1 per PO4 taken up;
         ! respiratory NH4/PO4 return reverses the sign. Leaf terms on pelagic
         ! TA, root terms on the benthic alkalinity pool.
         _SET_BOTTOM_EXCHANGE_(self%id_TA, jN3_leaf - jN4_leaf + jP_leaf &
            + qn * (Ra_act + Ra_bas) + BGn / max(BGc, 1.0e-8_rk) * Rb &
            - qp * (Ra_act + Ra_bas) - BGp / max(BGc, 1.0e-8_rk) * Rb)
         _SET_BOTTOM_ODE_(self%id_benTA, jN3_root - jN4_root + jP_root)

         ! Mortality routing
         _SET_BOTTOM_EXCHANGE_(self%id_R6c, self%f_pel * M_ag)
         _SET_BOTTOM_EXCHANGE_(self%id_R6n, self%f_pel * qn * M_ag)
         _SET_BOTTOM_EXCHANGE_(self%id_R6p, self%f_pel * qp * M_ag)
         _SET_BOTTOM_ODE_(self%id_Q6c, (1.0_rk - self%f_pel) * M_ag + M_bg)
         _SET_BOTTOM_ODE_(self%id_Q6n, (1.0_rk - self%f_pel) * qn * M_ag &
            + BGn / max(BGc, 1.0e-8_rk) * M_bg)
         _SET_BOTTOM_ODE_(self%id_Q6p, (1.0_rk - self%f_pel) * qp * M_ag &
            + BGp / max(BGc, 1.0e-8_rk) * M_bg)

         ! --- Optional grazing closure (docs/13; inert when g_max = 0) ----
         F_gr = 0.0_rk
         if (self%g_max > 0.0_rk) then
            _GET_HORIZONTAL_(self%id_gr1c, gr1c)
            _GET_HORIZONTAL_(self%id_gr2c, gr2c)
            F_gr = self%g_max * (gr1c + gr2c) &
                   * AGc * AGc / (AGc * AGc + self%h_ag * self%h_ag)
            ! Grazer respiration against the near-bottom water, limited by
            ! the same pelagic-O2 Monod factor as the plant's AG respiration;
            ! the O2-suppressed remainder is egested with the faeces.
            resp_gr = self%pe_gr * fO2ag * F_gr
            eges_gr = F_gr - resp_gr
            _SET_BOTTOM_ODE_(self%id_AGc, -F_gr)
            _SET_BOTTOM_ODE_(self%id_AGn, -qn * F_gr)
            _SET_BOTTOM_ODE_(self%id_AGp, -qp * F_gr)
            _SET_BOTTOM_EXCHANGE_(self%id_O3c, resp_gr / CMass)
            _SET_BOTTOM_EXCHANGE_(self%id_O2o, -resp_gr / CMass)
            ! Nutrient and alkalinity return of the respired fraction
            ! (+1 eq per NH4 released, -1 eq per PO4 released)
            _SET_BOTTOM_EXCHANGE_(self%id_N4n, qn * resp_gr)
            _SET_BOTTOM_EXCHANGE_(self%id_N1p, qp * resp_gr)
            _SET_BOTTOM_EXCHANGE_(self%id_TA, (qn - qp) * resp_gr)
            ! Egestion to the plant detritus pool at AG quota
            _SET_BOTTOM_ODE_(self%id_Q6c, eges_gr)
            _SET_BOTTOM_ODE_(self%id_Q6n, qn * eges_gr)
            _SET_BOTTOM_ODE_(self%id_Q6p, qp * eges_gr)
         end if

         ! --- Diagnostics --------------------------------------------------
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_gpp, Pg)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_npp, Pg - Ra_act - Ra_bas - Rn - Rb)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resp, Ra_act + Ra_bas + Rn + Rb)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fT, eT)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fI, eI)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_graz, F_gr)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN3l, jN3_leaf)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN4l, jN4_leaf)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_relN4, &
            qn * (Ra_act + Ra_bas) + BGn / max(BGc, 1.0e-8_rk) * Rb)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN3r, jN3_root)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN4r, jN4_root)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module ersem_seagrass
