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
! 2026-09-29 (isw_fix = 1, default): the three accepted fixes of
! nippon-steel docs/21 s4.3 (Q10 maintenance with heat mortality; BG
! respiration products and root-uptake alkalinity in the sediment layers;
! no nightly light-coded mobilisation). isw_fix = 0: the equations below.
! v1 simplifications (documented in the design doc, section 6):
!   - all respired CO2 is returned to pelagic DIC (no porewater DIC)
!   - no epiphytes, no shoot demography, no root oxygen loss
!   - PAR received is the host-provided bottom PAR; in the 0-D box this
!     is the box PAR after background extinction
!
! 2026-10-07 (isw_uni = 1): the unified eelgrass formulation shared with MUSE (muse docs/EELGRASS_UNIFIED_SPEC_20261007.md):
! one law set and ONE parameter table in both models (tissue-specific uptake kinetics on pore-water concentration over the
! shared root profile, NSC-share maintenance on water O2, EMS allocation, resorption, leaf and root exudation, ...). The
! equivalence of the rate terms with MUSE is tested by muse/tests/eelgrass_equiv. isw_uni = 0 (default) is the module as it
! was, bit for bit (muse/tests/eelgrass_equiv/regress_isw0.py).
! Deactivation contract: if no instance of this model appears in
! fabm.yaml, nothing is registered and results are identical to a build
! without this module (same convention as benthic_cao with iswCaO=0).
!-----------------------------------------------------------------------

module ersem_seagrass

   use fabm_types
   use ersem_shared

   implicit none

   private

   ! jsasaki 2026-10-07: diagnostics of the unified formulation (isw_uni = 1): the rate terms the equivalence test compares
   integer, parameter :: nuni = 45
   character(len=12), parameter :: uni_name(nuni) = [character(len=12) :: 'Pg', 'Ract', 'RmA', 'RmB', 'phi', 'Tst', 'Alloc', 'Mob', &
      'ExuL', 'ExuR', 'MA', 'MB', 'MN', 'dAGc', 'dBGc', 'dNSC', 'L4', 'L3', 'LP', 'U4', 'U3', 'UP', 'w1', 'c4_1', 'c4_2', 'c3_1', &
      'c3_2', 'cP_1', 'cP_2', 'dic_w', 'o2_w', 'ta_w', 'nh4_w', 'po4_w', 'dic_pw', 'ta_pw', 'par_in', 'dAGn', 'dAGp', 'dBGn', 'dBGp', 'detN', 'detP', 'nh4_pw', 'po4_pw']

   ! units of the unified diagnostics (mg C or mmol m-2 d-1 as noted in docs/EELGRASS_UNIFIED_SPEC_20261007.md s7)
   character(len=16), parameter :: uni_unit(nuni) = [character(len=16) :: &
      'mg C/m^2/d', 'mg C/m^2/d', 'mg C/m^2/d', 'mg C/m^2/d', &
      '-', 'mg C/m^2/d', 'mg C/m^2/d', 'mg C/m^2/d', &
      'mg C/m^2/d', 'mg C/m^2/d', 'mg C/m^2/d', 'mg C/m^2/d', &
      'mg C/m^2/d', 'mg C/m^2/d', 'mg C/m^2/d', 'mg C/m^2/d', &
      'mmol N/m^2/d', 'mmol N/m^2/d', 'mmol P/m^2/d', 'mmol N/m^2/d', &
      'mmol N/m^2/d', 'mmol P/m^2/d', '-', 'mmol N/m^3', &
      'mmol N/m^3', 'mmol N/m^3', 'mmol N/m^3', 'mmol P/m^3', &
      'mmol P/m^3', 'mmol C/m^2/d', 'mmol O2/m^2/d', 'mmol eq/m^2/d', &
      'mmol N/m^2/d', 'mmol P/m^2/d', 'mmol C/m^2/d', 'mmol eq/m^2/d', &
      'W/m^2', 'mmol N/m^2/d', 'mmol P/m^2/d', 'mmol N/m^2/d', &
      'mmol P/m^2/d', 'mmol N/m^2/d', 'mmol P/m^2/d', 'mmol N/m^2/d', &
      'mmol P/m^2/d']

   type,extends(type_base_model),public :: type_ersem_seagrass
      ! Own bottom state variables
      type (type_bottom_state_variable_id) :: id_AGc, id_AGn, id_AGp
      type (type_bottom_state_variable_id) :: id_BGc, id_BGn, id_BGp
      type (type_bottom_state_variable_id) :: id_NSCc

      ! Pelagic dependencies (bottom cell)
      type (type_state_variable_id) :: id_O2o, id_O3c, id_TA
      type (type_state_variable_id) :: id_N1p, id_N3n, id_N4n
      type (type_state_variable_id) :: id_R6c, id_R6n, id_R6p
      type (type_state_variable_id) :: id_R2c   ! pelagic semi-labile DOC (f_exu > 0 only)

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

      ! ledger counters (isw_ledger = 1 only; nippon-steel docs/119 FC1)
      integer  :: isw_ledger
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_growth_O3c, id_ledger_growth_O2o
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_bgresp_G2o
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_leaf_N4n, id_ledger_leaf_N3n, id_ledger_leaf_N1p
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_root_K4n1, id_ledger_root_K4n2
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_root_K3n1, id_ledger_root_K3n2
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_root_K1p1, id_ledger_root_K1p2
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_leaf_TA, id_ledger_root_benTA
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_exu_R2c
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_slough_R6c, id_ledger_slough_R6n, id_ledger_slough_R6p
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_graze_O3c, id_ledger_graze_O2o
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_graze_N4n, id_ledger_graze_N1p, id_ledger_graze_TA

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
      real(rk) :: psi_inh   ! isw_nupt = 2: ammonium inhibition of nitrate uptake, (mmol N/m3)^-1
      real(rk) :: srs_ag, srs_bg, r_nsc, pu_ra, hO2, hG2o
      real(rk) :: f_rspn
      real(rk) :: pq, rq_o2c
      real(rk) :: tau_store, tau_mob, k_bg
      real(rk) :: f_exu
      real(rk) :: rs_target, k_alloc, q_nsc
      real(rk) :: sd_ag, sd_bg, sd_anx, nsc_starve, sd_starve, f_pel
      real(rk) :: g_max, h_ag, pe_gr
      ! the three accepted fixes of nippon-steel docs/21 s4.3 (isw_fix = 1, the default since 2026-09-29; 0 = the former
      ! equations, bit for bit): maintenance respiration on a Q10 law (q10_m, Tref_m) instead of the growth CTMI, with
      ! heat mortality sd_heat above T_heat; below-ground respiration products into the sediment layers (DIC to G3c1/2,
      ! NH4 to K4n1/2, PO4 to K1p1/2, TA to benTA/benTA2, split f_bg1 : 1 - f_bg1) and root-uptake alkalinity in the layer
      ! of each uptake; reserve mobilisation only when AG is below its share of BG (no nightly light switch), limited by
      ! the leaf quotas
      integer  :: isw_fix
      real(rk) :: q10_m, Tref_m, sd_heat, T_heat, f_bg1, K_dic
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_bgresp_G3c1, id_ledger_bgresp_G3c2, &
         id_ledger_bgresp_K4n1, id_ledger_bgresp_K4n2, id_ledger_bgresp_K1p1, id_ledger_bgresp_K1p2, &
         id_ledger_bgresp_benTA, id_ledger_bgresp_benTA2, id_ledger_root_benTA2
      type (type_bottom_state_variable_id) :: id_G3c1, id_G3c2, id_benTA2
      ! jsasaki 2026-10-07: unified eelgrass formulation shared with MUSE (muse docs/EELGRASS_UNIFIED_SPEC_20261007.md).
      ! isw_uni = 1 replaces the process laws by the unified set U1-U20 (needs isw_fix = 1); isw_uni = 0 (default) is the former module.
      integer  :: isw_uni
      real(rk) :: q10, Tref, k_tr, rs, k_mob, K_nsc, e_exu = 0.0_rk, e_leaf = 0.0_rk, f_recl, sd_hyp, z_p, z_max
      real(rk) :: V_l4, V_l3, V_lP, V_r4, V_r3, V_rP, K_l4, K_l3, K_lP, K_r4, K_r3, K_rP
      logical  :: no3_red
      ! jsasaki 2026-10-07: unified mat (microphytobenthos) switches, legacy path only (muse docs/UNIFY_MAT_SPEC_20261007.md), per instance.
      ! isw_no3red = 1: nitrate assimilation with its reductant, 2 AG C oxidised to DIC per NO3-N taken up, no O2 used (needs pq = rq_o2c = 1);
      ! f_dk < 1: dark factor on nitrate uptake, jN3 x (f_dk + (1 - f_dk) eI); isw_matdiag = 1: registers upt_P, mort_AG, red_C.
      ! Defaults (0, 1, 0) leave the module bit-identical to the former one.
      integer  :: isw_no3red, isw_matdiag
      real(rk) :: f_dk
      type (type_horizontal_diagnostic_variable_id) :: id_uPl, id_mAG, id_redC, id_xO3c, id_xO2, id_xTA, id_dAG, id_xG3c
      type (type_bottom_state_variable_id) :: id_Q1c
      type (type_horizontal_dependency_id) :: id_K1p1w, id_K1p2w, id_K3n1w, id_K3n2w, id_K4n1w, id_K4n2w
      type (type_horizontal_dependency_id) :: id_D1m, id_D2m, id_poro
      type (type_horizontal_diagnostic_variable_id) :: id_u(nuni)
   contains
      procedure :: initialize
      procedure :: do_bottom
      procedure :: do_bottom_uni
   end type

contains

   subroutine initialize(self, configunit)
      class (type_ersem_seagrass), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      real(rk) :: a, b, ab2
      logical  :: uni
      integer  :: iu

      ! Set time unit to d-1 (ERSEM convention)
      self%dt = 86400._rk

      ! jsasaki 2026-10-07: unified formulation switch (read first: it selects the defaults of the shared parameters)
      call self%get_parameter(self%isw_uni, 'isw_uni', '', 'unified eelgrass formulation shared with MUSE (0: this module as '// &
         'before [rfB13], 1: laws U1-U20 and the shared parameter set; requires isw_fix = 1)', default=0, minimum=0, maximum=1)
      uni = self%isw_uni == 1

      ! --- Parameters (defaults follow docs/10 section 5) -----------------
      call self%get_parameter(self%p_max, 'p_max', '1/d', &
         'maximum gross production at Topt', default=0.12_rk, minimum=0.0_rk)
      call self%get_parameter(self%alpha, 'alpha', '(W/m^2)^-1', &
         'initial slope of the P-I curve (tanh)', default=dflt(uni, 0.06_rk, 0.04_rk), minimum=0.0_rk)
      call self%get_parameter(self%a_lai, 'a_lai', 'm^2/mg C', &
         'leaf area per unit AG carbon', default=dflt(uni, 4.0e-5_rk, 4.8e-4_rk / CMass), minimum=0.0_rk)
      call self%get_parameter(self%k_can, 'k_can', '-', &
         'canopy attenuation coefficient per unit LAI', default=0.7_rk, minimum=0.0_rk)

      call self%get_parameter(self%Tmin, 'Tmin', 'degrees_Celsius', &
         'CTMI minimum temperature for growth', default=2.0_rk)
      call self%get_parameter(self%Topt, 'Topt', 'degrees_Celsius', &
         'CTMI optimal temperature for growth', default=dflt(uni, 18.0_rk, 25.0_rk))
      call self%get_parameter(self%Tmax, 'Tmax', 'degrees_Celsius', &
         'CTMI maximum temperature for growth', default=dflt(uni, 30.0_rk, 38.0_rk))
      call self%get_parameter(self%T_heat, 'T_heat', 'degrees_Celsius', 'temperature above which heat mortality '// &
         'acts (isw_fix = 1)', default=dflt(uni, self%Tmax, 30.0_rk))
      if (self%Tmin >= self%Topt) call self%fatal_error('initialize','CTMI requires Tmin < Topt')
      if (self%Topt >= self%Tmax) call self%fatal_error('initialize','CTMI requires Topt < Tmax')
      a = self%Topt - self%Tmin
      b = self%Topt - self%Tmax
      ab2 = (a * b)**2
      self%ctmi_a = -(a + b) / ab2
      self%ctmi_b = (a * b + (a + b) * self%Topt) / ab2

      call self%get_parameter(self%qn_min, 'qn_min', 'mmol N/mg C', &
         'minimum AG nitrogen quota', default=dflt(uni, 0.003_rk, 0.036_rk / CMass), minimum=0.0_rk)
      call self%get_parameter(self%qn_max, 'qn_max', 'mmol N/mg C', &
         'maximum AG nitrogen quota', default=dflt(uni, 0.008_rk, 0.096_rk / CMass), minimum=0.0_rk)
      call self%get_parameter(self%qp_min, 'qp_min', 'mmol P/mg C', &
         'minimum AG phosphorus quota', default=dflt(uni, 8.0e-5_rk, 0.00096_rk / CMass), minimum=0.0_rk)
      call self%get_parameter(self%qp_max, 'qp_max', 'mmol P/mg C', &
         'maximum AG phosphorus quota', default=dflt(uni, 2.5e-4_rk, 0.0030_rk / CMass), minimum=0.0_rk)
      call self%get_parameter(self%qn_bg, 'qn_bg', 'mmol N/mg C', &
         'fixed BG nitrogen quota', default=dflt(uni, 0.0035_rk, 0.042_rk / CMass), minimum=0.0_rk)
      call self%get_parameter(self%qp_bg, 'qp_bg', 'mmol P/mg C', &
         'fixed BG phosphorus quota', default=dflt(uni, 1.0e-4_rk, 0.0012_rk / CMass), minimum=0.0_rk)

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
         'nitrogen uptake rule (0: ammonium priority [legacy], 1: potentials capped by quota demand, 2: as 1 with ammonium inhibition of nitrate uptake)', &
         default=0, minimum=0, maximum=2)
      call self%get_parameter(self%psiN4, 'psiN4', '-', &
         'ammonium affinity weight (isw_nupt=1 only)', default=1.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%psiN3, 'psiN3', '-', &
         'nitrate affinity weight (isw_nupt=1 only)', default=1.0_rk, minimum=0.0_rk)
      ! isw_nupt = 2 (2026-09-08, nippon-steel docs/114 CX): the isw_nupt = 1 rule
      ! with the classical ammonium inhibition of nitrate uptake (Wroblewski 1977,
      ! exp(-psi_inh * NH4)). Motivation: a mat at its quota ceiling returns
      ! nitrogen by respiration and re-takes nitrate and ammonium in proportion to
      ! their concentrations, so the water's nitrate drains while ammonium stays;
      ! the tanks keep nitrate and hold ammonium low. isw_nupt = 1 is untouched.
      call self%get_parameter(self%psi_inh, 'psi_inh', '(mmol N/m3)^-1', &
         'ammonium inhibition of nitrate uptake (isw_nupt=2 only)', default=1.5_rk, minimum=0.0_rk)

      ! f_rspn (2026-09-09, nippon-steel docs/114 DH): the FRACTION of the
      ! above-ground respiratory nitrogen and phosphorus that the plant RETAINS
      ! instead of returning to the water. The code below was written to keep
      ! the quota from drifting when carbon is respired; physiologically basal
      ! respiration oxidises carbon skeletons while the nitrogen in protein and
      ! pigment is re-used, so the quota SHOULD rise. Default 0.0 reproduces the
      ! previous behaviour bit for bit; 1.0 retains everything. Carbon and
      ! oxygen fluxes are untouched, and so are the below-ground (Rb) terms.
      call self%get_parameter(self%f_rspn, 'f_rspn', '-', &
         'fraction of respiratory N and P retained by the plant (isw_uni = 1: structural respiration of AG and BG)', &
         default=0.0_rk, minimum=0.0_rk, maximum=1.0_rk)

      call self%get_parameter(self%srs_ag, 'srs_ag', '1/d', &
         'AG basal respiration (at Tref_m with isw_fix = 1, at Topt with isw_fix = 0)', default=0.015_rk, minimum=0.0_rk)
      call self%get_parameter(self%srs_bg, 'srs_bg', '1/d', &
         'BG basal respiration (at Tref_m with isw_fix = 1, at Topt with isw_fix = 0)', default=0.006_rk, minimum=0.0_rk)
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
      ! Exudation (docs/114 EV, 2026-09-12): a fraction of GROSS fixation
      ! released as dissolved organic carbon to the pelagic semi-labile pool
      ! (R2, carbon only) instead of entering the plant. DIC uptake and O2
      ! release are those of the fixation and do not change; only where the
      ! fixed carbon goes does. Default 0 = the former model, and the R2
      ! coupling is then not even requested (bit-identical).
      call self%get_parameter(self%f_exu, 'f_exu', '-', &
           'fraction of gross production exuded as DOC to pelagic R2 (isw_uni = 0 only; the unified formulation uses e_leaf)', default=0.0_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%tau_mob, 'tau_mob', '1/d', &
         'NSC remobilisation rate (isw_fix = 0: under light limitation; 1: when AG is below its share)', &
         default=0.05_rk, minimum=0.0_rk)
      call self%get_parameter(self%isw_fix, 'isw_fix', '', &
         'docs/21 s4.3 fixes (0: former equations, 1: Q10 maintenance, BG products in the sediment, no light-coded '// &
         'mobilisation)', default=1, minimum=0, maximum=1)
      call self%get_parameter(self%q10_m, 'q10_m', '-', 'Q10 of maintenance respiration (isw_fix = 1; Marsh et al. '// &
         '1986)', default=2.4_rk, minimum=1.0_rk)
      call self%get_parameter(self%Tref_m, 'Tref_m', 'degrees_Celsius', 'reference temperature of the maintenance '// &
         'rates (isw_fix = 1)', default=20.0_rk)
      call self%get_parameter(self%sd_heat, 'sd_heat', '1/d', 'heat mortality above T_heat (isw_fix = 1; assumption)', &
         default=1.0_rk / 30.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%f_bg1, 'f_bg1', '-', 'share of the below-ground actions in sediment layer 1 '// &
         '(isw_fix = 1; assumption)', default=0.2_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%K_dic, 'K_dic', 'mmol C/m^3', 'half-saturation of gross fixation in the water DIC '// &
         '(isw_fix = 1; a positivity guard: the fixation stops when the DIC is exhausted)', default=dflt(uni, 1.0_rk, 10.0_rk), &
         minimum=0.0_rk)
      if (self%isw_fix == 1 .and. .not. (self%K_dic > 0.0_rk)) &
         call self%fatal_error('initialize', 'isw_fix = 1 requires K_dic > 0 (review p12 #8)')
      ! positive O2 half-saturations (a zero one gives 0/0 in anoxia; review p13 #10), and finite values of the
      ! corrected formulation's parameters (a NaN T_heat would silently disable the heat mortality; review p13 #11)
      if (.not. (self%hO2 > 0.0_rk .and. self%hG2o > 0.0_rk)) &
         call self%fatal_error('initialize', 'hO2 and hG2o must be positive')
      if (self%isw_fix == 1 .and. .not. all(abs([self%q10_m, self%Tref_m, self%sd_heat, self%T_heat, self%f_bg1, &
          self%K_dic]) <= huge(1.0_rk))) call self%fatal_error('initialize', 'isw_fix = 1: q10_m, Tref_m, sd_heat, '// &
          'T_heat, f_bg1 and K_dic must be finite')
      if (self%isw_fix == 1 .and. (self%qn_min >= self%qn_max .or. self%qp_min >= self%qp_max)) &
         call self%fatal_error('initialize', 'isw_fix = 1 requires qn_min < qn_max and qp_min < qp_max')
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
         'AG background sloughing rate', default=dflt(uni, 0.004_rk, 0.015_rk), minimum=0.0_rk)
      call self%get_parameter(self%sd_bg, 'sd_bg', '1/d', &
         'BG background mortality rate', default=dflt(uni, 0.002_rk, 0.004_rk), minimum=0.0_rk)
      call self%get_parameter(self%sd_anx, 'sd_anx', '1/d', &
         'extra BG mortality under benthic anoxia', default=0.01_rk, minimum=0.0_rk)
      call self%get_parameter(self%nsc_starve, 'nsc_starve', '-', &
         'NSC:BG carbon ratio below which starvation mortality starts', &
         default=0.02_rk, minimum=0.0_rk)
      call self%get_parameter(self%sd_starve, 'sd_starve', '1/d', &
         'extra AG mortality under NSC starvation', default=dflt(uni, 0.02_rk, 1.0_rk / 30.0_rk), minimum=0.0_rk)
      call self%get_parameter(self%f_pel, 'f_pel', '-', &
         'fraction of AG sloughing routed to pelagic POM (R6)', &
         default=0.5_rk, minimum=0.0_rk, maximum=1.0_rk)

      ! jsasaki 2026-10-07: parameters of the unified formulation (isw_uni = 1 only; names as in MUSE &seagrass, per mg C here:
      ! the MUSE value per mmol C divided by CMass = 12.011 wherever carbon is in the unit denominator)
      if (uni) then
         if (self%isw_fix /= 1) call self%fatal_error('initialize', 'isw_uni = 1 requires isw_fix = 1')
         call self%get_parameter(self%q10, 'q10', '-', 'Q10 of maintenance respiration and nutrient uptake demand (not of exudation)', &
            default=2.4_rk, minimum=1.0_rk)
         call self%get_parameter(self%Tref, 'Tref', 'degrees_Celsius', 'reference temperature of q10', default=20.0_rk)
         call self%get_parameter(self%k_tr, 'k_tr', '1/d', 'relaxation rate of AG -> BG allocation to the BG:AG target', &
            default=0.033_rk, minimum=0.0_rk)
         call self%get_parameter(self%rs, 'rs', '-', 'target BG:AG carbon ratio', default=0.75_rk, minimum=0.01_rk)
         call self%get_parameter(self%k_mob, 'k_mob', '1/d', 'NSC remobilisation rate when AG is below its share', &
            default=0.05_rk, minimum=0.0_rk)
         call self%get_parameter(self%K_nsc, 'K_nsc', '-', 'NSC:BG ratio at which the reserve pays half of the maintenance', &
            default=0.02_rk, minimum=1.0e-9_rk)
         call self%get_parameter(self%e_exu, 'e_exu', '1/d', 'root exudation per unit BG carbon (to benthic Q1c, mg C m-2); light and '// &
            'temperature independent', default=0.0018_rk, minimum=0.0_rk)
         call self%get_parameter(self%e_leaf, 'e_leaf', '1/d', 'leaf exudation per unit AG carbon (to pelagic R2c); light and temperature '// &
            'independent', default=0.0014_rk, minimum=0.0_rk)
         call self%get_parameter(self%f_recl, 'f_recl', '-', 'fraction of leaf N, P resorbed at senescence', &
            default=0.2_rk, minimum=0.0_rk, maximum=1.0_rk)
         call self%get_parameter(self%sd_hyp, 'sd_hyp', '1/d', 'extra BG mortality per unit water-O2 deficit (1 - fO2)', &
            default=0.0_rk, minimum=0.0_rk)
         call self%get_parameter(self%z_p, 'z_p', 'm', 'depth of the root-profile maximum', default=0.03_rk, minimum=1.0e-4_rk, maximum=10.0_rk)
         call self%get_parameter(self%z_max, 'z_max', 'm', 'root depth', default=0.15_rk, minimum=1.0e-3_rk, maximum=100.0_rk)
         call self%get_parameter(self%V_l4, 'V_l4', 'mmol N/mg C/d', 'maximum leaf NH4 uptake', default=0.038_rk / CMass, minimum=0.0_rk)
         call self%get_parameter(self%V_l3, 'V_l3', 'mmol N/mg C/d', 'maximum leaf NO3 uptake', default=0.027_rk / CMass, minimum=0.0_rk)
         call self%get_parameter(self%V_lP, 'V_lP', 'mmol P/mg C/d', 'maximum leaf PO4 uptake', default=0.014_rk / CMass, minimum=0.0_rk)
         call self%get_parameter(self%V_r4, 'V_r4', 'mmol N/mg C/d', 'maximum root NH4 uptake', default=0.0161_rk / CMass, minimum=0.0_rk)
         call self%get_parameter(self%V_r3, 'V_r3', 'mmol N/mg C/d', 'maximum root NO3 uptake', default=0.02645_rk / CMass, minimum=0.0_rk)
         call self%get_parameter(self%V_rP, 'V_rP', 'mmol P/mg C/d', 'maximum root PO4 uptake', default=0.002645_rk / CMass, minimum=0.0_rk)
         call self%get_parameter(self%K_l4, 'K_l4', 'mmol N/m^3', 'leaf NH4 half-saturation', default=58.9_rk, minimum=1.0e-6_rk)
         call self%get_parameter(self%K_l3, 'K_l3', 'mmol N/m^3', 'leaf NO3 half-saturation', default=42.8_rk, minimum=1.0e-6_rk)
         call self%get_parameter(self%K_lP, 'K_lP', 'mmol P/m^3', 'leaf PO4 half-saturation', default=7.6_rk, minimum=1.0e-6_rk)
         call self%get_parameter(self%K_r4, 'K_r4', 'mmol N/m^3', 'root NH4 half-saturation (pore-water concentration)', &
            default=48.6_rk, minimum=1.0e-6_rk)
         call self%get_parameter(self%K_r3, 'K_r3', 'mmol N/m^3', 'root NO3 half-saturation (pore-water concentration)', &
            default=53.3_rk, minimum=1.0e-6_rk)
         call self%get_parameter(self%K_rP, 'K_rP', 'mmol P/m^3', 'root PO4 half-saturation (pore-water concentration)', &
            default=6.0_rk, minimum=1.0e-6_rk)
         call self%get_parameter(self%no3_red, 'no3_red', '', 'nitrate assimilation oxidises 2 C per N (leaves: water DIC; roots: pore DIC)', &
            default=.false.)
         if (.not. (self%Tmin < self%Topt .and. self%Topt < self%Tmax .and. 2.0_rk * self%Topt >= self%Tmin + self%Tmax)) &
            call self%fatal_error('initialize', 'isw_uni = 1: the cardinal-temperature model needs Tmin < Topt < Tmax and Topt >= (Tmin + Tmax)/2')
         if (.not. (self%K_nsc > 0.0_rk .and. self%rs > 0.0_rk .and. self%qn_bg > 0.0_rk .and. self%qp_bg > 0.0_rk .and. self%q_nsc > 0.0_rk)) &
            call self%fatal_error('initialize', 'isw_uni = 1: K_nsc, rs, qn_bg, qp_bg and q_nsc must be positive')
      end if

      ! jsasaki 2026-10-07: unified mat switches (legacy path; the unified eelgrass has its own no3_red and laws)
      call self%get_parameter(self%isw_no3red, 'isw_no3red', '', 'nitrate assimilation with its reductant (2 AG C oxidised to DIC '// &
         'per NO3-N, no O2; requires pq = rq_o2c = 1) [0: as before]', default=0, minimum=0, maximum=1)
      call self%get_parameter(self%f_dk, 'f_dk', '-', 'dark/light ratio of nitrate uptake: jN3 x (f_dk + (1 - f_dk) eI) [1: as before]', &
         default=1.0_rk, minimum=0.0_rk, maximum=1.0_rk)
      ! jsasaki 2026-10-07: review round 2 #7: FABM's bounds checks do not reject NaN
      if (.not. (self%f_dk >= 0.0_rk .and. self%f_dk <= 1.0_rk)) &
         call self%fatal_error('initialize', 'f_dk must be finite and in [0, 1]')
      call self%get_parameter(self%isw_matdiag, 'isw_matdiag', '', 'register the mat diagnostics upt_P, mort_AG, red_C [0: none]', &
         default=0, minimum=0, maximum=1)
      if (self%isw_uni == 1 .and. (self%isw_no3red /= 0 .or. self%f_dk < 1.0_rk .or. self%isw_matdiag /= 0)) &
         call self%fatal_error('initialize', 'isw_no3red, f_dk and isw_matdiag belong to the legacy path (isw_uni = 0)')
      if (self%f_dk < 1.0_rk .and. self%isw_nupt == 0) &
         call self%fatal_error('initialize', 'f_dk < 1 requires isw_nupt >= 1 (the legacy remainder rule has no nitrate potential)')
      if (self%isw_no3red == 1 .and. .not. (self%pq == 1.0_rk .and. self%rq_o2c == 1.0_rk)) &
         call self%fatal_error('initialize', 'isw_no3red = 1 requires pq = 1 and rq_o2c = 1 (CH2O electron balance)')

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
      if (uni .and. self%g_max > 0.0_rk) call self%fatal_error('initialize', 'isw_uni = 1 does not include the grazing closure (g_max must be 0)')
      call self%get_parameter(self%h_ag, 'h_ag', 'mg C/m^2', &
         'AG carbon at the type-III grazing half-saturation', &
         default=500.0_rk, minimum=1.0e-6_rk)
      call self%get_parameter(self%pe_gr, 'pe_gr', '-', &
         'fraction of grazed carbon respired by the grazers', &
         default=0.4_rk, minimum=0.0_rk, maximum=1.0_rk)

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
      ! The exudation (f_exu > 0) and grazing (g_max > 0) counters are always
      ! registered and written as zero when their branch is inactive.
      call self%get_parameter(self%isw_ledger, 'isw_ledger', '', &
           'ledger counters: diagnostics of the source terms applied (0: off, 1: on)', &
           default=0, minimum=0, maximum=1)
      if (uni .and. self%isw_ledger == 1) call self%fatal_error('initialize', 'isw_uni = 1 has its own diagnostics; use isw_ledger = 0')
      ! jsasaki 2026-10-07: the ledger counters do not carry the nitrate reductant
      if (self%isw_no3red == 1 .and. self%isw_ledger == 1) &
         call self%fatal_error('initialize', 'isw_no3red = 1 is not covered by the isw_ledger counters')
      if (self%isw_ledger == 1) then
         call self%register_diagnostic_variable(self%id_ledger_growth_O3c, 'ledger_growth_O3c', 'mmol C/m^2/d', &
              'ledger: gross fixation and AG, NSC and BG respiration -> pelagic DIC (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_growth_O2o, 'ledger_growth_O2o', 'mmol O_2/m^2/d', &
              'ledger: gross fixation and AG and NSC respiration -> pelagic oxygen (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_bgresp_G2o, 'ledger_bgresp_G2o', 'mmol O_2/m^2/d', &
              'ledger: BG respiration -> benthic oxygen layer 1', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_leaf_N4n, 'ledger_leaf_N4n', 'mmol N/m^2/d', &
              'ledger: leaf ammonium uptake and AG and BG respiratory return -> pelagic ammonium (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_leaf_N3n, 'ledger_leaf_N3n', 'mmol N/m^2/d', &
              'ledger: leaf nitrate uptake -> pelagic nitrate (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_leaf_N1p, 'ledger_leaf_N1p', 'mmol P/m^2/d', &
              'ledger: leaf phosphate uptake and AG and BG respiratory return -> pelagic phosphate (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_root_K4n1, 'ledger_root_K4n1', 'mmol N/m^2/d', &
              'ledger: root ammonium uptake -> porewater ammonium layer 1', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_root_K4n2, 'ledger_root_K4n2', 'mmol N/m^2/d', &
              'ledger: root ammonium uptake -> porewater ammonium layer 2', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_root_K3n1, 'ledger_root_K3n1', 'mmol N/m^2/d', &
              'ledger: root nitrate uptake -> porewater nitrate layer 1', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_root_K3n2, 'ledger_root_K3n2', 'mmol N/m^2/d', &
              'ledger: root nitrate uptake -> porewater nitrate layer 2', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_root_K1p1, 'ledger_root_K1p1', 'mmol P/m^2/d', &
              'ledger: root phosphate uptake -> porewater phosphate layer 1', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_root_K1p2, 'ledger_root_K1p2', 'mmol P/m^2/d', &
              'ledger: root phosphate uptake -> porewater phosphate layer 2', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_leaf_TA, 'ledger_leaf_TA', 'mmol eq/m^2/d', &
              'ledger: leaf N and P uptake and AG and BG respiratory N and P return -> pelagic alkalinity (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_root_benTA, 'ledger_root_benTA', 'mmol eq/m^2/d', &
              'ledger: root N and P uptake -> benthic alkalinity layer 1', &
              domain=domain_bottom, source=source_do_bottom)
         if (self%isw_fix == 1) then
            ! the sediment-side sources of fix 2, each equal to the source applied (the BG respiration products per
            ! layer, and the root-uptake alkalinity of layer 2)
            call self%register_diagnostic_variable(self%id_ledger_bgresp_G3c1, 'ledger_bgresp_G3c1', 'mmol C/m^2/d', &
                 'ledger: BG respiration -> benthic DIC layer 1', domain=domain_bottom, source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_bgresp_G3c2, 'ledger_bgresp_G3c2', 'mmol C/m^2/d', &
                 'ledger: BG respiration -> benthic DIC layer 2', domain=domain_bottom, source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_bgresp_K4n1, 'ledger_bgresp_K4n1', 'mmol N/m^2/d', &
                 'ledger: BG respiratory N -> porewater ammonium layer 1', domain=domain_bottom, source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_bgresp_K4n2, 'ledger_bgresp_K4n2', 'mmol N/m^2/d', &
                 'ledger: BG respiratory N -> porewater ammonium layer 2', domain=domain_bottom, source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_bgresp_K1p1, 'ledger_bgresp_K1p1', 'mmol P/m^2/d', &
                 'ledger: BG respiratory P -> porewater phosphate layer 1', domain=domain_bottom, &
                 source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_bgresp_K1p2, 'ledger_bgresp_K1p2', 'mmol P/m^2/d', &
                 'ledger: BG respiratory P -> porewater phosphate layer 2', domain=domain_bottom, &
                 source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_bgresp_benTA, 'ledger_bgresp_benTA', &
                 'mmol eq/m^2/d', 'ledger: BG respiratory N and P -> benthic alkalinity layer 1', &
                 domain=domain_bottom, source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_bgresp_benTA2, 'ledger_bgresp_benTA2', &
                 'mmol eq/m^2/d', 'ledger: BG respiratory N and P -> benthic alkalinity layer 2', &
                 domain=domain_bottom, source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_root_benTA2, 'ledger_root_benTA2', 'mmol eq/m^2/d', &
                 'ledger: root N and P uptake -> benthic alkalinity layer 2', domain=domain_bottom, &
                 source=source_do_bottom)
         end if
         call self%register_diagnostic_variable(self%id_ledger_exu_R2c, 'ledger_exu_R2c', 'mg C/m^2/d', &
              'ledger: exudation of gross fixation -> pelagic semi-labile DOC (exchange; 0 when f_exu = 0)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_slough_R6c, 'ledger_slough_R6c', 'mg C/m^2/d', &
              'ledger: AG sloughing routed to the water -> pelagic POM carbon (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_slough_R6n, 'ledger_slough_R6n', 'mmol N/m^2/d', &
              'ledger: AG sloughing routed to the water -> pelagic POM nitrogen (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_slough_R6p, 'ledger_slough_R6p', 'mmol P/m^2/d', &
              'ledger: AG sloughing routed to the water -> pelagic POM phosphorus (exchange)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_graze_O3c, 'ledger_graze_O3c', 'mmol C/m^2/d', &
              'ledger: grazer respiration -> pelagic DIC (exchange; 0 when g_max = 0)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_graze_O2o, 'ledger_graze_O2o', 'mmol O_2/m^2/d', &
              'ledger: grazer respiration -> pelagic oxygen (exchange; 0 when g_max = 0)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_graze_N4n, 'ledger_graze_N4n', 'mmol N/m^2/d', &
              'ledger: grazer respiratory return -> pelagic ammonium (exchange; 0 when g_max = 0)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_graze_N1p, 'ledger_graze_N1p', 'mmol P/m^2/d', &
              'ledger: grazer respiratory return -> pelagic phosphate (exchange; 0 when g_max = 0)', &
              domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_graze_TA, 'ledger_graze_TA', 'mmol eq/m^2/d', &
              'ledger: grazer respiratory N and P return -> pelagic alkalinity (exchange; 0 when g_max = 0)', &
              domain=domain_bottom, source=source_do_bottom)
      end if

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
      if (self%f_exu > 0.0_rk) &
         call self%register_state_dependency(self%id_R2c, 'R2c', 'mg C/m^3', 'pelagic semi-labile DOC carbon')
      if (uni) then
         if (self%e_leaf > 0.0_rk .and. .not. self%f_exu > 0.0_rk) &
            call self%register_state_dependency(self%id_R2c, 'R2c', 'mg C/m^3', 'pelagic semi-labile DOC carbon')
      end if
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
      if (self%isw_fix == 1) then
         ! the sediment side of the below-ground respiration and of the layer-2 root uptake (default couplings to the
         ! standard ERSEM benthic column; an ERSEM configuration without it fails explicitly)
         call self%register_state_dependency(self%id_G3c1, 'G3c1', 'mmol C/m^2', 'benthic DIC, layer 1')
         call self%register_state_dependency(self%id_G3c2, 'G3c2', 'mmol C/m^2', 'benthic DIC, layer 2')
         call self%register_state_dependency(self%id_benTA2, 'benTA2', 'mEq/m^2', 'benthic alkalinity, layer 2')
         call self%request_coupling(self%id_G3c1, 'G3/per_layer/c1')
         call self%request_coupling(self%id_G3c2, 'G3/per_layer/c2')
         call self%request_coupling(self%id_benTA2, 'G5/per_layer/a2')
      end if
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

      ! jsasaki 2026-10-07: mat diagnostics for the equivalence test (isw_matdiag = 1 only)
      if (self%isw_matdiag == 1) then
         call self%register_diagnostic_variable(self%id_uPl, 'upt_P', 'mmol P/m^2/d', 'phosphate uptake from the water column by leaves', &
            domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_mAG, 'mort_AG', 'mg C/m^2/d', 'AG mortality (before routing)', &
            domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_redC, 'red_C', 'mg C/m^2/d', 'AG carbon oxidised as nitrate reductant', &
            domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_xO3c, 'x_O3c', 'mmol C/m^2/d', 'DIC exchange with the water (leaf and AG terms)', &
            domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_xO2, 'x_O2', 'mmol O_2/m^2/d', 'oxygen exchange with the water', &
            domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_xTA, 'x_TA', 'mmol eq/m^2/d', 'alkalinity exchange with the water', &
            domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_xG3c, 'x_G3c', 'mmol C/m^2/d', 'reductant DIC into the pore-water layers (root nitrate)', &
            domain=domain_bottom, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_dAG, 'd_AGc', 'mg C/m^2/d', 'AG carbon tendency (without grazing)', &
            domain=domain_bottom, source=source_do_bottom)
      end if

      ! jsasaki 2026-10-07: couplings and diagnostics of the unified formulation (nothing is registered with isw_uni = 0)
      if (uni) then
         call self%register_dependency(self%id_K1p1w, 'K1p1_pw', 'mmol P/m^2', 'pore-water phosphate of layer 1')
         call self%register_dependency(self%id_K1p2w, 'K1p2_pw', 'mmol P/m^2', 'pore-water phosphate of layer 2')
         call self%register_dependency(self%id_K3n1w, 'K3n1_pw', 'mmol N/m^2', 'pore-water nitrate of layer 1')
         call self%register_dependency(self%id_K3n2w, 'K3n2_pw', 'mmol N/m^2', 'pore-water nitrate of layer 2')
         call self%register_dependency(self%id_K4n1w, 'K4n1_pw', 'mmol N/m^2', 'pore-water ammonium of layer 1')
         call self%register_dependency(self%id_K4n2w, 'K4n2_pw', 'mmol N/m^2', 'pore-water ammonium of layer 2')
         call self%register_dependency(self%id_D1m, depth_of_bottom_interface_of_layer_1)
         call self%register_dependency(self%id_D2m, depth_of_bottom_interface_of_layer_2)
         call self%register_dependency(self%id_poro, sediment_porosity)
         if (self%e_exu > 0.0_rk) call self%register_state_dependency(self%id_Q1c, 'Q1c', 'mg C/m^2', &
            'benthic dissolved organic carbon (root exudate)')
         do iu = 1, nuni
            call self%register_diagnostic_variable(self%id_u(iu), 'uni_'//trim(uni_name(iu)), trim(uni_unit(iu)), &
               'unified-formulation term '//trim(uni_name(iu)), domain=domain_bottom, source=source_do_bottom)
         end do
      end if

   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_seagrass), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: AGc, AGn, AGp, BGc, BGn, BGp, NSCc
      real(rk) :: O2o, N1p, N3n, N4n, O3c, fK1, fK3, fP1
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
      real(rk) :: M_ag, M_bg, starve, M_nsc
      real(rk) :: dAGc, dAGn, dAGp, dBGc, dBGn, dBGp, dNSC
      real(rk) :: gr1c, gr2c, F_gr, resp_gr, eges_gr
      real(rk) :: gTm, rbn, rbp, f1, rbw
      real(rk) :: Rred_l, Rred_r       ! jsasaki 2026-10-07: nitrate reductant carbon, leaf and root (isw_no3red = 1), mg C/m^2/d
      real(rk) :: xred_w, xred_s       ! jsasaki 2026-10-07: reductant DIC to the water and to the pore-water layers, mmol C/m^2/d (applied values)

      ! jsasaki 2026-10-07: the unified formulation (isw_uni = 1) has its own routine
      if (self%isw_uni == 1) then
         call self%do_bottom_uni(_ARGUMENTS_DO_BOTTOM_)
         return
      end if

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
         _GET_(self%id_O3c, O3c)
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
         ! no fixation from exhausted water DIC (isw_fix = 1; review p11 #17)
         if (self%isw_fix == 1) Pg = Pg * max(O3c, 0.0_rk) / (max(O3c, 0.0_rk) + self%K_dic)
         fO2ag = max(0.0_rk, O2o) / (max(0.0_rk, O2o) + self%hO2)
         fO2bg = max(0.0_rk, G2o) / (max(0.0_rk, G2o) + self%hG2o)
         Ra_act = self%pu_ra * Pg
         if (self%isw_fix == 1) then
            ! maintenance on a Q10 law, never switched off by the growth cardinal temperatures (fix 1)
            gTm = self%q10_m**((ETW - self%Tref_m) / 10.0_rk)
         else
            gTm = eT
         end if
         Ra_bas = self%srs_ag * gTm * AGc * fO2ag
         Rn = self%r_nsc * gTm * NSCc * fO2ag
         Rb = self%srs_bg * gTm * BGc * fO2bg

         ! --- Translocation and BG growth (v2 allocation control) ---------
         Pnet_ag = Pg - Ra_act - Ra_bas
         ! Root:shoot deficit: 1 when BG is absent, 0 at/above the target ratio
         relE = min(1.0_rk, max(0.0_rk, &
                1.0_rk - (BGc / max(AGc, 1.0e-8_rk)) / self%rs_target))
         ! Storage: bank surplus production until the reserve reaches its
         ! target fraction of BG structure
         nsc_gap = max(0.0_rk, 1.0_rk - NSCc / max(self%q_nsc * BGc, 1.0e-8_rk))
         T_st = self%tau_store * max(0.0_rk, Pnet_ag) * nsc_gap
         if (self%isw_fix == 1) then
            ! mobilisation only while AG is below its share of BG, as far as the leaf quotas support new structure (fix 3)
            T_mb = 0.0_rk
            ! the nutrient-surplus factor eQ is zero at the minimum quotas, so new structure never dilutes the leaf
            ! quotas below them (review p11 #2)
            if (BGc > 1.0e-8_rk) T_mb = self%tau_mob * NSCc * max(0.0_rk, 1.0_rk - AGc * self%rs_target / BGc) * eQ
         else
            T_mb = self%tau_mob * NSCc * max(0.0_rk, 1.0_rk - eI)
         end if
         ! BG growth, driven by the deficit: direct allocation of net
         ! production (carbon from AG) plus reserve-fed growth (carbon from NSC)
         G_alloc = relE * self%k_alloc * max(0.0_rk, Pnet_ag)
         G_res = relE * self%k_bg * eT * NSCc
         ! BG structure takes AG N, P at its own quotas: only from the leaves' surplus above their minimum quotas
         ! (isw_fix = 1; eQ is zero at the minimum quotas, so the leaves never cross them; review p12 #5)
         if (self%isw_fix == 1) then
            G_alloc = G_alloc * eQ; G_res = G_res * eQ
         end if
         G_bg = G_alloc + G_res
         ! BG structural growth needs N and P from the AG pools at fixed quota;
         ! scale it down when the AG pools cannot pay.
         G_bg_c = min(G_bg, &
                      0.5_rk * AGn / max(self%qn_bg, 1.0e-12_rk), &
                      0.5_rk * AGp / max(self%qp_bg, 1.0e-12_rk))
         if (self%isw_fix == 1) then
            ! the donors lose exactly what BG gains, also for trace growth (review p12 #7)
            gscale = 0.0_rk
            if (G_bg > 0.0_rk) gscale = G_bg_c / G_bg
         else
            gscale = G_bg_c / max(G_bg, 1.0e-12_rk)
         end if
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
            if (self%isw_nupt == 2) then
               ! Ammonium inhibition of nitrate uptake (Wroblewski 1977). The
               ! potentials above are unchanged for isw_nupt = 1 (bit-identical).
               jN3_leaf = jN3_leaf * exp(-self%psi_inh * max(N4n, 0.0_rk))
               jN3_root = jN3_root * exp(-self%psi_inh * max(wsum, 0.0_rk))
            end if
            ! jsasaki 2026-10-07: dark factor on nitrate uptake (f_dk < 1 only; nitrate reductase needs photosynthate), before the quota cap
            if (self%f_dk < 1.0_rk) then
               jN3_leaf = jN3_leaf * (self%f_dk + (1.0_rk - self%f_dk) * eI)
               jN3_root = jN3_root * (self%f_dk + (1.0_rk - self%f_dk) * eI)
            end if
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
         if (self%isw_fix == 1) then
            ! each root uptake drawn from the two layers in proportion to their non-negative inventories; nothing is
            ! taken from an empty pair, so the plant's gain always equals the porewater debit (review p11 #13)
            call split2(K4n1, K4n2, jN4_root, fK1)
            call split2(K3n1, K3n2, jN3_root, fK3)
            call split2(K1p1, K1p2, jP_root, fP1)
         end if

         ! --- Mortality ----------------------------------------------------
         starve = 0.0_rk
         if (NSCc < self%nsc_starve * max(BGc, 1.0e-8_rk)) starve = self%sd_starve
         M_ag = (self%sd_ag + starve) * AGc
         M_bg = (self%sd_bg + self%sd_anx * (1.0_rk - fO2bg)) * BGc
         if (self%isw_fix == 1 .and. ETW > self%T_heat) then            ! heat injury (fix 1, second half)
            M_ag = M_ag + self%sd_heat * AGc; M_bg = M_bg + self%sd_heat * BGc
         end if
         ! the reserve of the dying rhizomes dies with them, to plant detritus (isw_fix = 1; review p12 #6)
         M_nsc = 0.0_rk
         if (self%isw_fix == 1 .and. BGc > 0.0_rk) M_nsc = M_bg / BGc * NSCc

         ! --- State ODEs (per day) ----------------------------------------
         dAGc = Pg - Ra_act - Ra_bas - T_st + T_mb - M_ag - G_alloc - self%f_exu * Pg
         ! jsasaki 2026-10-07: nitrate assimilation oxidises 2 AG C per N to DIC, no O2 (MUSE no3_reductant); N and P stocks untouched
         Rred_l = 0.0_rk; Rred_r = 0.0_rk
         if (self%isw_no3red == 1) then
            Rred_l = 2.0_rk * CMass * jN3_leaf
            Rred_r = 2.0_rk * CMass * jN3_root
            dAGc = dAGc - Rred_l - Rred_r
         end if
         dNSC = T_st - T_mb - G_res - Rn - M_nsc
         dBGc = G_bg_c - Rb - M_bg
         dAGn = jN4_leaf + jN3_leaf + jN4_root + jN3_root &
                - self%qn_bg * G_bg_c - qn * M_ag
         dAGp = jP_leaf + jP_root - self%qp_bg * G_bg_c - qp * M_ag
         dBGn = self%qn_bg * G_bg_c - BGn / max(BGc, 1.0e-8_rk) * M_bg
         dBGp = self%qp_bg * G_bg_c - BGp / max(BGc, 1.0e-8_rk) * M_bg
         ! Respiration releases the associated N and P to the water (AG) as NH4
         ! and PO4 at the current quota, keeping the quota from drifting when
         ! carbon is respired.
         dAGn = dAGn - (1.0_rk - self%f_rspn) * qn * (Ra_act + Ra_bas)
         dAGp = dAGp - (1.0_rk - self%f_rspn) * qp * (Ra_act + Ra_bas)
         dBGn = dBGn - BGn / max(BGc, 1.0e-8_rk) * Rb
         dBGp = dBGp - BGp / max(BGc, 1.0e-8_rk) * Rb

         _SET_BOTTOM_ODE_(self%id_AGc, dAGc)
         _SET_BOTTOM_ODE_(self%id_AGn, dAGn)
         _SET_BOTTOM_ODE_(self%id_AGp, dAGp)
         _SET_BOTTOM_ODE_(self%id_BGc, dBGc)
         _SET_BOTTOM_ODE_(self%id_BGn, dBGn)
         _SET_BOTTOM_ODE_(self%id_BGp, dBGp)
         _SET_BOTTOM_ODE_(self%id_NSCc, dNSC)

         ! below-ground respiration products: to the water (isw_fix = 0) or into sediment layers 1 and 2 (fix 2)
         rbn = BGn / max(BGc, 1.0e-8_rk) * Rb; rbp = BGp / max(BGc, 1.0e-8_rk) * Rb
         f1 = self%f_bg1
         if (self%isw_fix == 1) then
            _SET_BOTTOM_ODE_(self%id_G3c1, f1 * Rb / CMass)
            _SET_BOTTOM_ODE_(self%id_G3c2, (1.0_rk - f1) * Rb / CMass)
            _SET_BOTTOM_ODE_(self%id_K4n1, f1 * rbn)
            _SET_BOTTOM_ODE_(self%id_K4n2, (1.0_rk - f1) * rbn)
            _SET_BOTTOM_ODE_(self%id_K1p1, f1 * rbp)
            _SET_BOTTOM_ODE_(self%id_K1p2, (1.0_rk - f1) * rbp)
            _SET_BOTTOM_ODE_(self%id_benTA, f1 * (rbn - rbp))
            _SET_BOTTOM_ODE_(self%id_benTA2, (1.0_rk - f1) * (rbn - rbp))
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_bgresp_G3c1, f1 * Rb / CMass)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_bgresp_G3c2, (1.0_rk - f1) * Rb / CMass)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_bgresp_K4n1, f1 * rbn)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_bgresp_K4n2, (1.0_rk - f1) * rbn)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_bgresp_K1p1, f1 * rbp)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_bgresp_K1p2, (1.0_rk - f1) * rbp)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_bgresp_benTA, f1 * (rbn - rbp))
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_bgresp_benTA2, (1.0_rk - f1) * (rbn - rbp))
            end if
            rbw = 0.0_rk                                  ! nothing of the BG respiration reaches the water
         else
            rbw = 1.0_rk
         end if

         ! --- Exchanges with the water column -----------------------------
         ! Carbon and oxygen. Carbon is carbon; the oxygen carries the two
         ! quotients (pq, rq_o2c), both 1 by default, which reproduces the
         ! former "PQ = 1" line term by term.
         _SET_BOTTOM_EXCHANGE_(self%id_O3c, (-Pg + Ra_act + Ra_bas + Rn + rbw * Rb) / CMass)
         _SET_BOTTOM_EXCHANGE_(self%id_O2o, (self%pq * Pg - self%rq_o2c * Ra_act &
                                             - self%rq_o2c * Ra_bas - self%rq_o2c * Rn) / CMass)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_growth_O3c, (-Pg + Ra_act + Ra_bas + Rn + rbw * Rb) / CMass)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_growth_O2o, (self%pq * Pg - self%rq_o2c * Ra_act &
                                             - self%rq_o2c * Ra_bas - self%rq_o2c * Rn) / CMass)
         end if
         ! jsasaki 2026-10-07: reductant DIC: leaf part to the water, root part to the pore-water layers of the uptake (isw_fix = 1) or the water
         xred_w = 0.0_rk; xred_s = 0.0_rk
         if (self%isw_no3red == 1) then
            if (self%isw_fix == 1) then
               xred_w = Rred_l / CMass
               xred_s = Rred_r / CMass
               _SET_BOTTOM_ODE_(self%id_G3c1, xred_s * fK3)
               _SET_BOTTOM_ODE_(self%id_G3c2, xred_s * (1.0_rk - fK3))
            else
               xred_w = (Rred_l + Rred_r) / CMass
            end if
            _SET_BOTTOM_EXCHANGE_(self%id_O3c, xred_w)
         end if
         ! BG respiration draws benthic layer-1 oxygen instead
         _SET_BOTTOM_ODE_(self%id_G2o, -self%rq_o2c * Rb / CMass)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_bgresp_G2o, -self%rq_o2c * Rb / CMass)
         end if

         ! Leaf nutrient uptake and respiratory return
         _SET_BOTTOM_EXCHANGE_(self%id_N4n, -jN4_leaf + (1.0_rk - self%f_rspn) * qn * (Ra_act + Ra_bas) + rbw * rbn)
         _SET_BOTTOM_EXCHANGE_(self%id_N3n, -jN3_leaf)
         _SET_BOTTOM_EXCHANGE_(self%id_N1p, -jP_leaf + (1.0_rk - self%f_rspn) * qp * (Ra_act + Ra_bas) + rbw * rbp)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_leaf_N4n, -jN4_leaf + (1.0_rk - self%f_rspn) * qn * (Ra_act + Ra_bas) + rbw * rbn)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_leaf_N3n, -jN3_leaf)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_leaf_N1p, -jP_leaf + (1.0_rk - self%f_rspn) * qp * (Ra_act + Ra_bas) + rbw * rbp)
         end if

         ! Root uptake from the porewater pools (split by availability)
         if (self%isw_fix == 1) then
            _SET_BOTTOM_ODE_(self%id_K4n1, -jN4_root * fK1)
            _SET_BOTTOM_ODE_(self%id_K4n2, -jN4_root * (1.0_rk - fK1))
            _SET_BOTTOM_ODE_(self%id_K3n1, -jN3_root * fK3)
            _SET_BOTTOM_ODE_(self%id_K3n2, -jN3_root * (1.0_rk - fK3))
            _SET_BOTTOM_ODE_(self%id_K1p1, -jP_root * fP1)
            _SET_BOTTOM_ODE_(self%id_K1p2, -jP_root * (1.0_rk - fP1))
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K4n1, -jN4_root * fK1)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K4n2, -jN4_root * (1.0_rk - fK1))
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K3n1, -jN3_root * fK3)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K3n2, -jN3_root * (1.0_rk - fK3))
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K1p1, -jP_root * fP1)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K1p2, -jP_root * (1.0_rk - fP1))
            end if
         else
            wsum = max(K4n1 + K4n2, 1.0e-8_rk)
            _SET_BOTTOM_ODE_(self%id_K4n1, -jN4_root * K4n1 / wsum)
            _SET_BOTTOM_ODE_(self%id_K4n2, -jN4_root * K4n2 / wsum)
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K4n1, -jN4_root * K4n1 / wsum)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K4n2, -jN4_root * K4n2 / wsum)
            end if
            wsum = max(K3n1 + K3n2, 1.0e-8_rk)
            _SET_BOTTOM_ODE_(self%id_K3n1, -jN3_root * K3n1 / wsum)
            _SET_BOTTOM_ODE_(self%id_K3n2, -jN3_root * K3n2 / wsum)
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K3n1, -jN3_root * K3n1 / wsum)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K3n2, -jN3_root * K3n2 / wsum)
            end if
            wsum = max(K1p1 + K1p2, 1.0e-8_rk)
            _SET_BOTTOM_ODE_(self%id_K1p1, -jP_root * K1p1 / wsum)
            _SET_BOTTOM_ODE_(self%id_K1p2, -jP_root * K1p2 / wsum)
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K1p1, -jP_root * K1p1 / wsum)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_K1p2, -jP_root * K1p2 / wsum)
            end if
         end if

         ! Alkalinity bookkeeping: +1 per NO3, -1 per NH4, +1 per PO4 taken up;
         ! respiratory NH4/PO4 return reverses the sign. Leaf terms on pelagic
         ! TA, root terms on the benthic alkalinity pool.
         _SET_BOTTOM_EXCHANGE_(self%id_TA, jN3_leaf - jN4_leaf + jP_leaf &
            + (1.0_rk - self%f_rspn) * qn * (Ra_act + Ra_bas) + rbw * rbn &
            - (1.0_rk - self%f_rspn) * qp * (Ra_act + Ra_bas) - rbw * rbp)
         if (self%isw_fix == 1) then
            ! the root-uptake alkalinity in the layer of each uptake (fix 2), each ledger equal to its source (p11 #6)
            _SET_BOTTOM_ODE_(self%id_benTA, jN3_root * fK3 - jN4_root * fK1 + jP_root * fP1)
            _SET_BOTTOM_ODE_(self%id_benTA2, jN3_root * (1.0_rk - fK3) - jN4_root * (1.0_rk - fK1) &
               + jP_root * (1.0_rk - fP1))
         else
            _SET_BOTTOM_ODE_(self%id_benTA, jN3_root - jN4_root + jP_root)
         end if
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_leaf_TA, jN3_leaf - jN4_leaf + jP_leaf &
            + (1.0_rk - self%f_rspn) * qn * (Ra_act + Ra_bas) + rbw * rbn &
            - (1.0_rk - self%f_rspn) * qp * (Ra_act + Ra_bas) - rbw * rbp)
            if (self%isw_fix == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_benTA, jN3_root * fK3 - jN4_root * fK1 + jP_root * fP1)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_benTA2, jN3_root * (1.0_rk - fK3) &
                  - jN4_root * (1.0_rk - fK1) + jP_root * (1.0_rk - fP1))
            else
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_root_benTA, jN3_root - jN4_root + jP_root)
            end if
         end if

         ! Exudate: the diverted share of gross fixation, carbon only
         if (self%f_exu > 0.0_rk) _SET_BOTTOM_EXCHANGE_(self%id_R2c, self%f_exu * Pg)
         if (self%isw_ledger == 1) then
            if (self%f_exu > 0.0_rk) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_exu_R2c, self%f_exu * Pg)
            else
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_exu_R2c, 0.0_rk)
            end if
         end if

         ! Mortality routing
         _SET_BOTTOM_EXCHANGE_(self%id_R6c, self%f_pel * M_ag)
         _SET_BOTTOM_EXCHANGE_(self%id_R6n, self%f_pel * qn * M_ag)
         _SET_BOTTOM_EXCHANGE_(self%id_R6p, self%f_pel * qp * M_ag)
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_slough_R6c, self%f_pel * M_ag)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_slough_R6n, self%f_pel * qn * M_ag)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_slough_R6p, self%f_pel * qp * M_ag)
         end if
         _SET_BOTTOM_ODE_(self%id_Q6c, (1.0_rk - self%f_pel) * M_ag + M_bg + M_nsc)
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
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_graze_O3c, resp_gr / CMass)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_graze_O2o, -resp_gr / CMass)
            end if
            ! Nutrient and alkalinity return of the respired fraction
            ! (+1 eq per NH4 released, -1 eq per PO4 released)
            _SET_BOTTOM_EXCHANGE_(self%id_N4n, qn * resp_gr)
            _SET_BOTTOM_EXCHANGE_(self%id_N1p, qp * resp_gr)
            _SET_BOTTOM_EXCHANGE_(self%id_TA, (qn - qp) * resp_gr)
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_graze_N4n, qn * resp_gr)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_graze_N1p, qp * resp_gr)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_graze_TA, (qn - qp) * resp_gr)
            end if
            ! Egestion to the plant detritus pool at AG quota
            _SET_BOTTOM_ODE_(self%id_Q6c, eges_gr)
            _SET_BOTTOM_ODE_(self%id_Q6n, qn * eges_gr)
            _SET_BOTTOM_ODE_(self%id_Q6p, qp * eges_gr)
         else if (self%isw_ledger == 1) then
            ! ledger counters of the inactive grazing branch (g_max = 0)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_graze_O3c, 0.0_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_graze_O2o, 0.0_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_graze_N4n, 0.0_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_graze_N1p, 0.0_rk)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_graze_TA, 0.0_rk)
         end if

         ! --- Diagnostics --------------------------------------------------
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_gpp, Pg)
         ! jsasaki 2026-10-07: review round 1 #2: the nitrate reductant carbon is respired carbon (Rred = 0 unless isw_no3red = 1)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_npp, Pg - Ra_act - Ra_bas - Rn - Rb - Rred_l - Rred_r)
         ! jsasaki 2026-10-07: review round 3 #4: total plant carbon respiration includes the carbon oxidised as nitrate reductant (Rred = 0 unless isw_no3red = 1)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resp, Ra_act + Ra_bas + Rn + Rb + Rred_l + Rred_r)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fT, eT)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fI, eI)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_graz, F_gr)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN3l, jN3_leaf)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN4l, jN4_leaf)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_relN4, &
            (1.0_rk - self%f_rspn) * qn * (Ra_act + Ra_bas) + rbw * rbn)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN3r, jN3_root)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN4r, jN4_root)
         ! jsasaki 2026-10-07: mat diagnostics of the equivalence test
         if (self%isw_matdiag == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uPl, jP_leaf)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_mAG, M_ag)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_redC, Rred_l + Rred_r)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_xO3c, (-Pg + Ra_act + Ra_bas + Rn + rbw * Rb) / CMass + xred_w)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_xO2, (self%pq * Pg - self%rq_o2c * (Ra_act + Ra_bas + Rn)) / CMass)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_xTA, jN3_leaf - jN4_leaf + jP_leaf &
               + (1.0_rk - self%f_rspn) * qn * (Ra_act + Ra_bas) + rbw * rbn &
               - (1.0_rk - self%f_rspn) * qp * (Ra_act + Ra_bas) - rbw * rbp)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dAG, dAGc)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_xG3c, xred_s)
         end if

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

   ! jsasaki 2026-10-07: the unified eelgrass formulation (isw_uni = 1), laws U1-U20 of muse docs/EELGRASS_UNIFIED_SPEC_20261007.md.
   ! Same laws as MUSE's seagrass_prepare/_rates/_water/_commit; the only structural difference is the root zone: two porewater
   ! layers (0 - D1m, D1m - D2m) whose root weights are the exact integrals of the shared profile r(z), renormalised over the two.
   subroutine do_bottom_uni(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_seagrass), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: AGc, AGn, AGp, BGc, BGn, BGp, NSCc, O2o, O3c, N1p, N3n, N4n, ETW, par
      real(rk) :: K1p(2), K3n(2), K4n(2), K1pw(2), K3nw(2), K4nw(2), D1m, D2m, poro
      real(rk) :: eT, gT, tau, I_can, eI, qn, qp, eQ, fC, fO2, Pg, Ract, RmA, RmB, phi, ExuL, ExuR, Pnet
      real(rk) :: Tst, Alloc, Mob, starve, heat, M_A, M_Ab, M_B, M_N, RmAr, RmBr, qnb, qpb
      real(rk) :: dn, dp, L4, L3, LP, scl, rest, su, sN, sP
      real(rk) :: dd(2), zt(2), zb(2), dF(2), wl(2), cn4(2), cn3(2), cp(2), u4(2), u3(2), up(2)
      real(rk) :: dAGc, dBGc, dNSC, dAGn, dAGp, dBGn, dBGp, qlost
      real(rk) :: red3l, red3r, rbn, rbp
      integer :: k
      real(rk), parameter :: epsC = 1.0e-8_rk * CMass      ! the MUSE carbon floor 1e-8 mmol C, in mg C

      _HORIZONTAL_LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_AGc, AGc);  _GET_HORIZONTAL_(self%id_AGn, AGn);  _GET_HORIZONTAL_(self%id_AGp, AGp)
         _GET_HORIZONTAL_(self%id_BGc, BGc);  _GET_HORIZONTAL_(self%id_BGn, BGn);  _GET_HORIZONTAL_(self%id_BGp, BGp)
         _GET_HORIZONTAL_(self%id_NSCc, NSCc)
         _GET_HORIZONTAL_(self%id_K1p1, K1p(1)); _GET_HORIZONTAL_(self%id_K1p2, K1p(2))
         _GET_HORIZONTAL_(self%id_K3n1, K3n(1)); _GET_HORIZONTAL_(self%id_K3n2, K3n(2))
         _GET_HORIZONTAL_(self%id_K4n1, K4n(1)); _GET_HORIZONTAL_(self%id_K4n2, K4n(2))
         _GET_HORIZONTAL_(self%id_K1p1w, K1pw(1)); _GET_HORIZONTAL_(self%id_K1p2w, K1pw(2))
         _GET_HORIZONTAL_(self%id_K3n1w, K3nw(1)); _GET_HORIZONTAL_(self%id_K3n2w, K3nw(2))
         _GET_HORIZONTAL_(self%id_K4n1w, K4nw(1)); _GET_HORIZONTAL_(self%id_K4n2w, K4nw(2))
         _GET_HORIZONTAL_(self%id_D1m, D1m); _GET_HORIZONTAL_(self%id_D2m, D2m); _GET_HORIZONTAL_(self%id_poro, poro)
         _GET_(self%id_O2o, O2o); _GET_(self%id_O3c, O3c); _GET_(self%id_N1p, N1p)
         _GET_(self%id_N3n, N3n); _GET_(self%id_N4n, N4n); _GET_(self%id_ETW, ETW); _GET_(self%id_par, par)

         AGc = max(AGc, 0.0_rk); BGc = max(BGc, 0.0_rk); NSCc = max(NSCc, 0.0_rk)

         ! U1, U2: temperature
         ! the cardinal-temperature model of Rosso et al. (1993, n = 2), as MUSE (review r2 #1); the cubic of the legacy path can have an interior cutoff
         if (ETW <= self%Tmin .or. ETW >= self%Tmax) then
            eT = 0.0_rk
         else
            eT = (ETW - self%Tmax) * (ETW - self%Tmin)**2 / ((self%Topt - self%Tmin) * ((self%Topt - self%Tmin) * (ETW - self%Topt) &
                 - (self%Topt - self%Tmax) * (self%Topt + self%Tmin - 2.0_rk * ETW)))
            eT = max(0.0_rk, min(1.0_rk, eT))
         end if
         gT = self%q10**((ETW - self%Tref) / 10.0_rk)
         ! U3: light with canopy self-shading
         tau = self%k_can * self%a_lai * AGc
         I_can = max(0.0_rk, par)
         if (tau > 1.0e-8_rk) I_can = I_can * (1.0_rk - exp(-tau)) / tau
         eI = tanh(self%alpha * I_can)
         ! U4: quota; U5: DIC and water O2
         qn = max(AGn, 0.0_rk) / max(AGc, epsC); qp = max(AGp, 0.0_rk) / max(AGc, epsC)
         eQ = min(max(0.0_rk, (qn - self%qn_min) / (self%qn_max - self%qn_min)), &
                  max(0.0_rk, (qp - self%qp_min) / (self%qp_max - self%qp_min)), 1.0_rk)
         fC = max(O3c, 0.0_rk) / (max(O3c, 0.0_rk) + self%K_dic)
         fO2 = max(O2o, 0.0_rk) / (max(O2o, 0.0_rk) + self%hO2)

         ! U6-U8: production, respiration, the reserve's share
         Pg = self%p_max * eT * eI * eQ * fC * AGc
         Ract = self%pu_ra * Pg
         RmA = self%srs_ag * gT * AGc * fO2
         RmB = self%srs_bg * gT * BGc * fO2
         phi = 0.0_rk
         if (NSCc > 0.0_rk) phi = NSCc / (NSCc + self%K_nsc * max(BGc, epsC))
         ExuL = self%e_leaf * AGc                                            ! U12 (light- and temperature-independent)
         Pnet = Pg - Ract - (1.0_rk - phi) * RmA - ExuL

         ! U9-U11: storage, allocation (relaxation to rs, quota-capped), mobilisation
         Tst = 0.0_rk
         if (BGc > 0.0_rk) Tst = self%tau_store * max(0.0_rk, Pnet) * max(0.0_rk, 1.0_rk - NSCc / (self%q_nsc * BGc))
         Alloc = self%k_tr * max(0.0_rk, (self%rs * AGc - BGc) / (1.0_rk + self%rs)) * eQ      ! U10: gated by the leaf nutrient surplus
         Alloc = min(Alloc, 0.5_rk * max(AGn, 0.0_rk) / self%qn_bg, 0.5_rk * max(AGp, 0.0_rk) / self%qp_bg)
         Mob = 0.0_rk
         if (BGc > 0.0_rk) Mob = self%k_mob * NSCc * max(0.0_rk, 1.0_rk - AGc * self%rs / BGc) * eQ

         ! U12, U13: root exudation, mortality
         ExuR = self%e_exu * BGc
         starve = 0.0_rk
         if (NSCc < self%nsc_starve * BGc) starve = self%sd_starve
         heat = 0.0_rk
         if (ETW > self%T_heat) heat = self%sd_heat
         M_Ab = self%sd_ag * AGc                                             ! the basal (senescence) loss: resorbed (U17)
         M_A = (self%sd_ag + starve + heat) * AGc
         M_B = (self%sd_bg + starve + heat + self%sd_hyp * (1.0_rk - fO2)) * BGc
         M_N = 0.0_rk
         if (BGc > 0.0_rk) M_N = M_B / BGc * NSCc

         ! U20: root weights of the two layers, exact integrals of the profile cut at z_max, renormalised over the two layers
         zt = [0.0_rk, min(D1m, self%z_max)]; zb = [min(D1m, self%z_max), min(D2m, self%z_max)]
         do k = 1, 2
            dF(k) = 0.0_rk
            if (zb(k) > zt(k)) dF(k) = root_int(zt(k), zb(k), self%z_p)
         end do
         if (dF(1) + dF(2) > 0.0_rk) then
            wl = dF / (dF(1) + dF(2))
         else
            wl = [1.0_rk, 0.0_rk]
         end if
         ! pore-water concentrations (mmol m-3 of pore water)
         dd = [D1m, D2m - D1m]
         do k = 1, 2
            if (poro * dd(k) > 0.0_rk) then
               cn4(k) = max(K4nw(k), 0.0_rk) / (poro * dd(k))
               cn3(k) = max(K3nw(k), 0.0_rk) / (poro * dd(k))
               cp(k)  = max(K1pw(k), 0.0_rk) / (poro * dd(k))
            else
               cn4(k) = 0.0_rk; cn3(k) = 0.0_rk; cp(k) = 0.0_rk
            end if
         end do

         ! U14: demand, leaves first, roots take the rest
         dn = self%vmax_n * gT * AGc * max(0.0_rk, 1.0_rk - qn / self%qn_max)
         dp = self%vmax_p * gT * AGc * max(0.0_rk, 1.0_rk - qp / self%qp_max)
         L4 = self%V_l4 * AGc * max(N4n, 0.0_rk) / (self%K_l4 + max(N4n, 0.0_rk))
         L3 = self%V_l3 * AGc * max(N3n, 0.0_rk) / (self%K_l3 + max(N3n, 0.0_rk))
         LP = self%V_lP * AGc * max(N1p, 0.0_rk) / (self%K_lP + max(N1p, 0.0_rk))
         if (L4 + L3 > dn .and. L4 + L3 > 0.0_rk) then
            scl = dn / (L4 + L3); L4 = L4 * scl; L3 = L3 * scl
         end if
         if (LP > dp .and. LP > 0.0_rk) LP = dp
         rest = max(0.0_rk, dn - L4 - L3)
         su = 0.0_rk
         do k = 1, 2
            u4(k) = self%V_r4 * BGc * wl(k) * cn4(k) / (self%K_r4 + cn4(k))
            u3(k) = self%V_r3 * BGc * wl(k) * cn3(k) / (self%K_r3 + cn3(k)) * exp(-self%psi_inh * cn4(k))
            su = su + u4(k) + u3(k)
         end do
         sN = 0.0_rk
         if (su > 0.0_rk) sN = min(1.0_rk, rest / su)
         u4 = sN * u4; u3 = sN * u3
         rest = max(0.0_rk, dp - LP)
         su = 0.0_rk
         do k = 1, 2
            up(k) = self%V_rP * BGc * wl(k) * cp(k) / (self%K_rP + cp(k))
            su = su + up(k)
         end do
         sP = 0.0_rk
         if (su > 0.0_rk) sP = min(1.0_rk, rest / su)
         up = sP * up

         ! U16, U17: respiratory return (structural share, times 1 - f_rspn); BG at the BG quotas
         RmAr = (1.0_rk - phi) * (1.0_rk - self%f_rspn) * RmA
         qnb = BGn / max(BGc, epsC); qpb = BGp / max(BGc, epsC)
         RmBr = (1.0_rk - phi) * (1.0_rk - self%f_rspn) * RmB
         rbn = RmBr * qnb; rbp = RmBr * qpb
         red3l = 0.0_rk; red3r = 0.0_rk
         if (self%no3_red) then
            red3l = 2.0_rk * L3                      ! mmol C per day, from the leaves
            red3r = 2.0_rk * (u3(1) + u3(2))         ! from the roots
         end if

         ! pools (per day)
         dAGc = Pg - Ract - (1.0_rk - phi) * RmA - ExuL - Tst - Alloc + Mob - M_A - red3l * CMass
         dNSC = Tst - Mob - phi * (RmA + RmB) - M_N
         dBGc = Alloc - (1.0_rk - phi) * RmB - ExuR - M_B - red3r * CMass
         dAGn = L4 + L3 + u4(1) + u4(2) + u3(1) + u3(2) - self%qn_bg * Alloc - RmAr * qn - (M_A - self%f_recl * M_Ab) * qn
         dAGp = LP + up(1) + up(2) - self%qp_bg * Alloc - RmAr * qp - (M_A - self%f_recl * M_Ab) * qp
         dBGn = self%qn_bg * Alloc - rbn - M_B * qnb
         dBGp = self%qp_bg * Alloc - rbp - M_B * qpb
         _SET_BOTTOM_ODE_(self%id_AGc, dAGc)
         _SET_BOTTOM_ODE_(self%id_AGn, dAGn)
         _SET_BOTTOM_ODE_(self%id_AGp, dAGp)
         _SET_BOTTOM_ODE_(self%id_BGc, dBGc)
         _SET_BOTTOM_ODE_(self%id_BGn, dBGn)
         _SET_BOTTOM_ODE_(self%id_BGp, dBGp)
         _SET_BOTTOM_ODE_(self%id_NSCc, dNSC)

         ! water column: leaf carbon and oxygen, leaf nutrients, structural-respiration return, alkalinity
         _SET_BOTTOM_EXCHANGE_(self%id_O3c, (-Pg + Ract + RmA + red3l * CMass) / CMass)
         _SET_BOTTOM_EXCHANGE_(self%id_O2o, (self%pq * Pg - self%rq_o2c * (Ract + RmA + RmB)) / CMass)
         _SET_BOTTOM_EXCHANGE_(self%id_N4n, -L4 + RmAr * qn)
         _SET_BOTTOM_EXCHANGE_(self%id_N3n, -L3)
         _SET_BOTTOM_EXCHANGE_(self%id_N1p, -LP + RmAr * qp)
         _SET_BOTTOM_EXCHANGE_(self%id_TA, L3 - L4 + LP + RmAr * (qn - qp))
         if (self%e_leaf > 0.0_rk) _SET_BOTTOM_EXCHANGE_(self%id_R2c, ExuL)

         ! sediment layers: root uptake, BG respiration products and exudation by the root weights
         _SET_BOTTOM_ODE_(self%id_K4n1, -u4(1) + wl(1) * rbn)
         _SET_BOTTOM_ODE_(self%id_K4n2, -u4(2) + wl(2) * rbn)
         _SET_BOTTOM_ODE_(self%id_K3n1, -u3(1))
         _SET_BOTTOM_ODE_(self%id_K3n2, -u3(2))
         _SET_BOTTOM_ODE_(self%id_K1p1, -up(1) + wl(1) * rbp)
         _SET_BOTTOM_ODE_(self%id_K1p2, -up(2) + wl(2) * rbp)
         _SET_BOTTOM_ODE_(self%id_G3c1, wl(1) * RmB / CMass + 2.0_rk * u3(1) * merge(1.0_rk, 0.0_rk, self%no3_red))
         _SET_BOTTOM_ODE_(self%id_G3c2, wl(2) * RmB / CMass + 2.0_rk * u3(2) * merge(1.0_rk, 0.0_rk, self%no3_red))
         _SET_BOTTOM_ODE_(self%id_benTA, u3(1) - u4(1) + up(1) + wl(1) * (rbn - rbp))
         _SET_BOTTOM_ODE_(self%id_benTA2, u3(2) - u4(2) + up(2) + wl(2) * (rbn - rbp))
         if (self%e_exu > 0.0_rk) _SET_BOTTOM_ODE_(self%id_Q1c, ExuR)

         ! mortality routing: AG split between pelagic POM and the plant detritus layer, BG and its reserve to the detritus
         _SET_BOTTOM_EXCHANGE_(self%id_R6c, self%f_pel * M_A)
         _SET_BOTTOM_EXCHANGE_(self%id_R6n, self%f_pel * qn * (M_A - self%f_recl * M_Ab))
         _SET_BOTTOM_EXCHANGE_(self%id_R6p, self%f_pel * qp * (M_A - self%f_recl * M_Ab))
         _SET_BOTTOM_ODE_(self%id_Q6c, (1.0_rk - self%f_pel) * M_A + M_B + M_N)
         _SET_BOTTOM_ODE_(self%id_Q6n, (1.0_rk - self%f_pel) * qn * (M_A - self%f_recl * M_Ab) + M_B * qnb)
         _SET_BOTTOM_ODE_(self%id_Q6p, (1.0_rk - self%f_pel) * qp * (M_A - self%f_recl * M_Ab) + M_B * qpb)

         ! standard diagnostics
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_gpp, Pg)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_npp, Pg - Ract - RmA - RmB - (red3l + red3r) * CMass)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_resp, Ract + RmA + RmB + (red3l + red3r) * CMass)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fT, eT)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fI, eI)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_graz, 0.0_rk)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN3l, L3)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN4l, L4)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_relN4, RmAr * qn)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN3r, u3(1) + u3(2))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uN4r, u4(1) + u4(2))
         ! the terms the equivalence test compares (mg C or mmol per m2 and day as in the pools; exchanges per mmol C)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(1), Pg);  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(2), Ract)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(3), RmA); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(4), RmB)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(5), phi); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(6), Tst)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(7), Alloc); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(8), Mob)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(9), ExuL); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(10), ExuR)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(11), M_A); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(12), M_B)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(13), M_N); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(14), dAGc)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(15), dBGc); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(16), dNSC)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(17), L4); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(18), L3)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(19), LP); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(20), u4(1) + u4(2))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(21), u3(1) + u3(2)); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(22), up(1) + up(2))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(23), wl(1)); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(24), cn4(1))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(25), cn4(2)); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(26), cn3(1))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(27), cn3(2)); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(28), cp(1))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(29), cp(2))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(30), (-Pg + Ract + RmA + red3l * CMass) / CMass)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(31), (self%pq * Pg - self%rq_o2c * (Ract + RmA + RmB)) / CMass)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(32), L3 - L4 + LP + RmAr * (qn - qp))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(33), -L4 + RmAr * qn)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(34), -LP + RmAr * qp)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(35), RmB / CMass + 2.0_rk * (u3(1) + u3(2)) * merge(1.0_rk, 0.0_rk, self%no3_red))
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(37), par)                      ! the bottom PAR the plant was given
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(38), dAGn); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(39), dAGp)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(40), dBGn); _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(41), dBGp)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(42), (M_A - self%f_recl * M_Ab) * qn + M_B * qnb)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(43), (M_A - self%f_recl * M_Ab) * qp + M_B * qpb)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(44), -(u4(1) + u4(2)) + rbn)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(45), -(up(1) + up(2)) + rbp)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_u(36), (u3(1) + u3(2)) - (u4(1) + u4(2)) + (up(1) + up(2)) + rbn - rbp)

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom_uni

   ! the integral over [za, zb] of the root profile (x/zp) e^(1 - x/zp), stably: e zp e^(-ua) [ua (1 - e^-d) + 1 - (1 + d) e^-d]
   pure real(rk) function root_int(za, zb, zp) result(F)
      real(rk), intent(in) :: za, zb, zp
      real(rk), parameter :: e1 = 2.718281828459045_rk
      real(rk) :: ua, d, g, om
      ua = za / zp; d = (zb - za) / zp
      if (d < 1.0e-3_rk) then
         g = d * d * (0.5_rk - d * (1.0_rk / 3.0_rk - d * 0.125_rk))
         om = d * (1.0_rk - d * (0.5_rk - d / 6.0_rk))      ! 1 - e^-d
      else
         g = 1.0_rk - (1.0_rk + d) * exp(-d)
         om = 1.0_rk - exp(-d)
      end if
      F = e1 * zp * exp(-ua) * (ua * om + g)
   end function

   ! jsasaki 2026-10-07: default of a parameter shared with the unified formulation: a0 (former module), a1 (isw_uni = 1)
   pure real(rk) function dflt(uni, a0, a1)
      logical,  intent(in) :: uni
      real(rk), intent(in) :: a0, a1
      if (uni) then
         dflt = a1
      else
         dflt = a0
      end if
   end function

   ! the layer-1 share f1 of an uptake j from two layers with inventories a, b (clipped at zero); j = 0 when both
   ! are empty
   pure subroutine split2(a, b, j, f1)
      real(rk), intent(in) :: a, b
      real(rk), intent(inout) :: j
      real(rk), intent(out) :: f1
      real(rk) :: t
      t = max(a, 0.0_rk) + max(b, 0.0_rk)
      if (t > 0.0_rk) then
         f1 = max(a, 0.0_rk) / t
      else
         f1 = 0.0_rk; j = 0.0_rk
      end if
   end subroutine

end module ersem_seagrass
