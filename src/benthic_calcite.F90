#include "fabm_driver.h"

! Benthic variable that supports resuspension and remineralization.
! Both processes return material to the pelagic.

module ersem_benthic_calcite

   use fabm_types

   use fabm_particle
   use ersem_shared
   use ersem_benthic_base
   ! jsasaki 2026-10-17: carbonate Step 2 CARB_OMEGA=pore: pore-water speciation and saturation routines
   use ersem_carbonate, only: CaCO3_Saturation
   use carbonate_engine, only: carbonate_engine_solve_porewater, calcium_of_salinity

   implicit none

!  default: all is private.
   private

   type,extends(type_ersem_benthic_base),public :: type_ersem_benthic_calcite
      type (type_horizontal_diagnostic_variable_id) :: id_dissolution
      type (type_horizontal_diagnostic_variable_id) :: id_precipitation
      type (type_dependency_id)                     :: id_Om_Cal
      type (type_dependency_id)                     :: id_par
      type (type_horizontal_dependency_id) :: id_par_bot   ! jsasaki 2026-10-11: Stage Z bed light
      integer :: isw_par_bot = 0   ! jsasaki 2026-10-11: Stage Z bed light

      ! Parameters
      real(rk) :: fdissmax, fdissmin, ndiss, KcalomX
      real(rk) :: K_prec, K_par_prec, n_prec
      integer  :: iswcal
      integer  :: iswdiss   ! jsasaki 2026-10-07: family E: 0 = the specific rate floor fdissmin (old); 1 = no floor (unified law with MUSE)
      ! E5 C1 (nippon-steel docs/125 s21): CO2*-promoted dissolution, default off
      real(rk) :: k_co2diss, co2_ref_diss
      type (type_dependency_id)                     :: id_CO2aq
      type (type_horizontal_diagnostic_variable_id) :: id_co2_dissolution
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_co2diss_c,id_ledger_co2diss_O3c,id_ledger_co2diss_TA
      ! nippon-steel docs/129 C2: reoxidation-driven dissolution in the oxic layer, default off
      real(rk) :: f_neut, c_half_neut
      type (type_horizontal_dependency_id)          :: id_acid_S, id_nrate
      type (type_bottom_state_variable_id)          :: id_benTA1, id_benDIC1
      type (type_horizontal_diagnostic_variable_id) :: id_acid_dissolution
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_acdiss_c, id_ledger_acdiss_benTA, id_ledger_acdiss_benDIC
      ! ledger counters (isw_ledger = 1 only; nippon-steel docs/119 FC1)
      integer  :: isw_ledger
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_diss_c,id_ledger_diss_O3c,id_ledger_diss_TA
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_prec_c,id_ledger_prec_O3c,id_ledger_prec_TA
      !---> jsasaki 2026-10-17: carbonate Step 2 CARB_OMEGA=pore (muse docs/CARBONATE_STEP2_CONTRACT_20261017.md): dissolution
      ! driven by the pore-water calcite saturation of each benthic layer, products into that layer's pore DIC and TA
      integer  :: isw_omega_pore = 0
      real(rk) :: z_cal
      type (type_dependency_id)                     :: id_ETW_p, id_X1X_p, id_dens_p
      type (type_horizontal_dependency_id)          :: id_poro_p, id_Dm_p(3)
      type (type_bottom_state_variable_id)          :: id_pDIC(3), id_pTA(3), id_pH2S(3), id_pNH4(3), id_pPO4(3), id_pSi(3)
      type (type_horizontal_diagnostic_variable_id) :: id_om_pore(3), id_pdiss(3), id_pfail
      type (type_horizontal_diagnostic_variable_id) :: id_ledger_pdiss_c, id_ledger_pdiss_benDIC(3), id_ledger_pdiss_benTA(3)
      !<--- jsasaki 2026-10-17

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
      ! jsasaki 2026-10-07: family E (docs/UNIFY_E_SPEC_20261007.md in muse). The floor fdissmin dissolves calcite also in
      ! supersaturated water (Omega > 1), which no mechanism supports; MUSE has none by default. iswdiss = 1 removes it, so the
      ! law is fdissmax (1 - Omega)^ndiss for Omega < 1 and zero above, as in MUSE calcite_coef. 0 (default) = old behaviour.
      call self%get_parameter(self%iswdiss,'iswdiss','','dissolution rate floor (0: fdissmin applies at every Omega, 1: no floor, unified law)',default=0,minimum=0,maximum=1)
      ! review round 2 #3: the unified power law needs a positive finite exponent (0**0 = 1 would dissolve a supersaturated bed)
      if (self%iswdiss == 1 .and. self%iswcal == 1) then
         if (.not. (self%ndiss > 0.0_rk .and. self%ndiss < huge(1.0_rk))) call self%fatal_error('initialize','iswdiss = 1 requires a positive finite ndiss')
      end if
      ! review round 1 #5: the hyperbolic law (iswcal = 2) is positive again above Omega = 1 + KcalomX, so it cannot be the unified law
      if (self%iswdiss == 1 .and. self%iswcal == 2) call self%fatal_error('initialize','iswdiss = 1 requires iswcal = 1 (power law) or 0; the hyperbolic law is not zero above saturation')
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
      ! jsasaki 2026-10-11: Stage Z bed light: optional bottom PAR (default 0 = off, bit-identical)
      call self%get_parameter(self%isw_par_bot, 'isw_par_bot', '', 'bed light: 0 = bottom-layer mean PAR, 1 = PAR at the bottom of the column (light/parEIR_bot)', default=0, minimum=0, maximum=1)
      if (self%isw_par_bot == 1) call self%register_dependency(self%id_par_bot, type_horizontal_standard_variable(name='ersem_bottom_downwelling_photosynthetic_radiative_flux'))
      call self%register_diagnostic_variable(self%id_precipitation,'precipitation','mg C/m^2/d','light-driven calcification',source=source_do_bottom)
      ! E5 C1 (2026-09-22, nippon-steel docs/125 s21, s21.9): acid-promoted
      ! bed-carbonate dissolution, a restricted comparator. CaCO3 + CO2 + H2O ->
      ! Ca2+ + 2HCO3-: per mol TA +2, DIC +1 (one carbon from the solid). Rate =
      ! k_co2diss * bL2c * [CO2*]/co2_ref_diss (mg C/m^2/d), added to the Omega
      ! law; the overlying water's CO2* is a declared proxy for pore-water
      ! acidity. Reported in its own diagnostic (co2_dissolution) so the
      ! Omega-law 'dissolution' is not double counted. k_co2diss = 0 (default)
      ! skips the term, its coupling and its diagnostic: bit-identical.
      call self%get_parameter(self%k_co2diss, 'k_co2diss', '1/d', &
         'E5 C1: CO2*-promoted specific dissolution rate at co2_ref_diss (0: off)', &
         default=0.0_rk, minimum=0.0_rk)
      call self%get_parameter(self%co2_ref_diss, 'co2_ref_diss', 'mmol/m^3', &
         'E5 C1: reference CO2* of the acid-promoted dissolution law', &
         default=15.0_rk, minimum=1.0e-6_rk)
      if (self%k_co2diss > 0.0_rk) then
         call self%register_dependency(self%id_CO2aq, 'CO2aq', 'mmol/m^3', &
            'carbonic acid concentration CO2* (E5 C1 driver)')
         call self%register_diagnostic_variable(self%id_co2_dissolution, 'co2_dissolution', 'mg C/m^2/d', &
            'E5 C1: CO2*-promoted dissolution', source=source_do_bottom)
      end if
      ! nippon-steel docs/129 C2 (2026-09-25): reoxidation-driven carbonate dissolution. The acid that sulfur
      ! reoxidation (2 H+ per S0 oxidised) and nitrification (2 H+ per NH4 oxidised) produce in the oxic layer is
      ! partly neutralised by the bed's CaCO3 instead of by the pore water's alkalinity (Green & Aller 2001:
      ! 40-100 % of coastal dissolution so driven). R = f_neut * (acid_S + 2 nrate) * c/(c + c_half_neut), in mmol
      ! CaCO3 per m^2 per day, f_neut = mol CaCO3 dissolved per mol H+ (0.5: the CO2 endpoint, CaCO3 + 2H+ -> Ca2+ +
      ! CO2, TA-neutral net; 1: the bicarbonate endpoint, CaCO3 + H+ -> Ca2+ + HCO3-, net +1 per H+). The solid
      ! loses R (mg C: R*CMass); the oxic layer's alkalinity gains 2R and its DIC R (the acid itself is already
      ! charged by the producing modules). f_neut = 0 (default) registers nothing: bit-identical.
      call self%get_parameter(self%f_neut, 'f_neut', 'mol/mol', &
         'C2: CaCO3 dissolved per mol of oxic-layer reoxidation acid (0: off; 0.5 CO2 endpoint, 1 HCO3 endpoint)', &
         default=0.0_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%c_half_neut, 'c_half_neut', 'mg C/m^2', &
         'C2: calcite stock at which the reoxidation-driven dissolution is halved', default=100.0_rk, minimum=1.0e-6_rk)
      if (self%f_neut > 0.0_rk) then
         call self%register_dependency(self%id_acid_S, 'acid_S', 'mmol eq/m^2/d', &
            'C2: oxic-layer acid production by sulfur reoxidation')
         call self%register_dependency(self%id_nrate, 'nrate', 'mmol N/m^2/d', 'C2: benthic nitrification rate')
         call self%register_state_dependency(self%id_benTA1, 'benTA1', 'mEq/m^2', 'C2: alkalinity of the oxic layer')
         call self%register_state_dependency(self%id_benDIC1, 'benDIC1', 'mmol C/m^2', 'C2: DIC of the oxic layer')
         call self%register_diagnostic_variable(self%id_acid_dissolution, 'acid_dissolution', 'mmol C/m^2/d', &
            'C2: reoxidation-driven CaCO3 dissolution', source=source_do_bottom)
      end if
      !---> jsasaki 2026-10-17: carbonate Step 2 CARB_OMEGA=pore. isw_omega_pore = 1: the power law fdissmax (1 - Omega)^ndiss
      ! (Omega < 1 only) reads the calcite saturation of the PORE WATER of each benthic layer instead of Om_Cal: speciation of the
      ! layer's DIC and TA (G3, G5) with its sulfide, ammonium, phosphate and silicate (carbonate_engine_solve_porewater: total
      ! scale, Lueker et al. 2000, borate of Lee et al. 2010, the option set of the water column) at the WATER's T and S, calcium
      ! from salinity. ERSEM keeps one inventory per dissolved column; its per_layer split is uniform (the same convention as the
      ! isw_fe laws), so the layer values differ only through what the layers hold. The stock is spread uniformly over the top
      ! z_cal metres; each layer dissolves its share at its own Omega, and the products (DIC +1, TA +2 per C) go into that
      ! layer's pore water instead of the overlying water. 0 (default) registers nothing: bit-identical.
      call self%get_parameter(self%isw_omega_pore,'isw_omega_pore','','dissolution driver: 0 = Om_Cal coupling (as before), 1 = pore-water calcite saturation per benthic layer',default=0,minimum=0,maximum=1)
      if (self%isw_omega_pore == 1) then
         if (self%iswcal /= 1) call self%fatal_error('initialize','isw_omega_pore = 1 requires iswcal = 1 (power law)')
         if (.not. (self%ndiss > 0.0_rk .and. self%ndiss < huge(1.0_rk))) call self%fatal_error('initialize','isw_omega_pore = 1 requires a positive finite ndiss')
         if (self%fdissmin /= 0.0_rk .and. self%iswdiss == 0) call self%fatal_error('initialize','isw_omega_pore = 1 requires no rate floor (fdissmin = 0 or iswdiss = 1)')
         call self%get_parameter(self%z_cal,'z_cal','m','depth over which the calcite stock is spread (pore-water driver)',default=0.01_rk,minimum=1.0e-6_rk)
         call self%register_dependency(self%id_ETW_p, standard_variables%temperature)
         call self%register_dependency(self%id_X1X_p, standard_variables%practical_salinity)
         call self%register_dependency(self%id_dens_p, standard_variables%density)
         call self%register_dependency(self%id_poro_p, sediment_porosity)
         call self%register_dependency(self%id_Dm_p(1), depth_of_bottom_interface_of_layer_1)
         call self%register_dependency(self%id_Dm_p(2), depth_of_bottom_interface_of_layer_2)
         call self%register_dependency(self%id_Dm_p(3), depth_of_sediment_column)
         call register_pore_layer(1, '1')
         call register_pore_layer(2, '2')
         call register_pore_layer(3, '3')
         call self%register_diagnostic_variable(self%id_pfail,'om_pore_fail','-','number of layers whose pore-water speciation failed (no dissolution there)',source=source_do_bottom)
      end if
      !<--- jsasaki 2026-10-17
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
         call self%register_diagnostic_variable(self%id_ledger_co2diss_c, 'ledger_co2diss_c', 'mg C/m^2/d', &
              'ledger: E5 C1 CO2*-promoted dissolution -> benthic calcite', source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_co2diss_O3c, 'ledger_co2diss_O3c', 'mmol C/m^2/d', &
              'ledger: E5 C1 CO2*-promoted dissolution -> pelagic dissolved inorganic carbon (bottom flux)', &
              source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_co2diss_TA, 'ledger_co2diss_TA', 'mmol eq/m^2/d', &
              'ledger: E5 C1 CO2*-promoted dissolution -> pelagic alkalinity (bottom flux)', source=source_do_bottom)
         if (self%f_neut > 0.0_rk) then
            call self%register_diagnostic_variable(self%id_ledger_acdiss_c, 'ledger_acdiss_c', 'mg C/m^2/d', &
                 'ledger: C2 reoxidation-driven dissolution -> benthic calcite', source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_acdiss_benTA, 'ledger_acdiss_benTA1', 'mmol eq/m^2/d', &
                 'ledger: C2 reoxidation-driven dissolution -> benthic alkalinity, layer 1', source=source_do_bottom)
            call self%register_diagnostic_variable(self%id_ledger_acdiss_benDIC, 'ledger_acdiss_benDIC1', 'mmol C/m^2/d', &
                 'ledger: C2 reoxidation-driven dissolution -> benthic DIC, layer 1', source=source_do_bottom)
         end if
      end if
      ! jsasaki 2026-10-17: carbonate Step 2 CARB_OMEGA=pore: its ledger counters (isw_ledger is read above)
      if (self%isw_ledger == 1 .and. self%isw_omega_pore == 1) then
         call self%register_diagnostic_variable(self%id_ledger_pdiss_c,'ledger_pdiss_c','mg C/m^2/d', &
              'ledger: pore-water-driven calcite dissolution -> benthic calcite', source=source_do_bottom)
         call register_pore_ledger(1, '1')
         call register_pore_ledger(2, '2')
         call register_pore_ledger(3, '3')
      end if

   contains

      ! jsasaki 2026-10-17: carbonate Step 2 CARB_OMEGA=pore: the pore-water inventories of one layer (couplings set in the yaml,
      ! e.g. benDIC1: G3/per_layer/c1, benTA1: G5/per_layer/a1, H2S1: G2_H2S/per_layer/h1_pw, NH41: K4/per_layer/n1_pw,
      ! PO41: K1/per_layer/p1_pw, Si1: K5/per_layer/s1_pw) and its diagnostics
      subroutine register_pore_layer(l, lab)
         integer,          intent(in) :: l
         character(len=*), intent(in) :: lab
         call self%register_state_dependency(self%id_pDIC(l),'benDIC'//lab,'mmol C/m^2','pore-water DIC of layer '//lab)
         call self%register_state_dependency(self%id_pTA(l),'benTA'//lab,'mEq/m^2','pore-water alkalinity of layer '//lab)
         call self%register_state_dependency(self%id_pH2S(l),'H2S'//lab,'mmol S/m^2','pore-water sulfide of layer '//lab)
         call self%register_state_dependency(self%id_pNH4(l),'NH4'//lab,'mmol N/m^2','pore-water ammonium of layer '//lab)
         call self%register_state_dependency(self%id_pPO4(l),'PO4'//lab,'mmol P/m^2','pore-water phosphate of layer '//lab)
         call self%register_state_dependency(self%id_pSi(l),'Si'//lab,'mmol Si/m^2','pore-water silicate of layer '//lab)
         call self%register_diagnostic_variable(self%id_om_pore(l),'om_cal_pore'//lab,'-','pore-water calcite saturation of layer '//lab,source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_pdiss(l),'pdissolution'//lab,'mg C/m^2/d','pore-water-driven dissolution in layer '//lab,source=source_do_bottom)
      end subroutine register_pore_layer

      subroutine register_pore_ledger(l, lab)
         integer,          intent(in) :: l
         character(len=*), intent(in) :: lab
         call self%register_diagnostic_variable(self%id_ledger_pdiss_benDIC(l),'ledger_pdiss_benDIC'//lab,'mmol C/m^2/d', &
              'ledger: pore-water-driven calcite dissolution -> benthic DIC, layer '//lab, source=source_do_bottom)
         call self%register_diagnostic_variable(self%id_ledger_pdiss_benTA(l),'ledger_pdiss_benTA'//lab,'mmol eq/m^2/d', &
              'ledger: pore-water-driven calcite dissolution -> benthic alkalinity, layer '//lab, source=source_do_bottom)
      end subroutine register_pore_ledger

   end subroutine

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)

      class (type_ersem_benthic_calcite),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: bL2c
      real(rk) :: Om_Cal
      real(rk) :: fdiss
      real(rk) :: par, F_prec
      real(rk) :: co2aq, F_co2
      real(rk) :: acid_S, nrate, R_ac
      ! jsasaki 2026-10-17: carbonate Step 2 CARB_OMEGA=pore
      real(rk) :: ETW_p, X1X_p, dens_p, poro_p, Dm_p(3), top_p, thk_p, w_p, vol_p, R_p(3), om_p(3), nfail_p
      real(rk) :: qDIC, qTA, qH2S, qNH4, qPO4, qSi, pH_p, fco2_p, h2co3_p, hco3_p, co3_p, k0_p, om_arg_p
      logical  :: ok_p
      integer  :: l_p

      _HORIZONTAL_LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_c, bL2c)
         if (self%iswcal>0 .or. self%K_prec>0._rk) then
            _GET_(self%id_om_cal, om_cal)
            om_cal=max(om_cal,0._rk)
         end if

         ! jsasaki 2026-10-17: carbonate Step 2 CARB_OMEGA=pore: the Om_Cal law unchanged when isw_omega_pore = 0
         if (self%isw_omega_pore == 0) then
            if (self%iswcal==0) then  ! NB select case would be cleaner but makes vectorization impossible for ifort 14
               fdiss = 0._rk
            elseif (self%iswcal==1) then
               fdiss = (max(1._rk-om_cal,0._rk))**self%ndiss
            else
               fdiss = max(0._rk,(1._rk-om_cal)/(1._rk-om_cal+self%KcalomX))
            end if

            if (self%iswdiss == 0 .or. self%iswcal == 0) then   ! iswcal = 0: fdissmin is the constant rate fdiss
               fdiss = max(fdiss * self%fdissmax, self%fdissmin)
            else
               fdiss = fdiss * self%fdissmax   ! jsasaki 2026-10-07: family E: unified law without the floor
            end if

            _SET_BOTTOM_ODE_(self%id_c, -fdiss*bL2c)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dissolution, -fdiss*bL2c)
            _SET_BOTTOM_EXCHANGE_(self%id_O3c, fdiss*bL2c/CMass)
            _SET_BOTTOM_EXCHANGE_(self%id_TA, 2*fdiss*bL2c/CMass)  ! Dissolution of CaCO3 increases alkalinity by 2 units
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_diss_c, -fdiss*bL2c)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_diss_O3c, fdiss*bL2c/CMass)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_diss_TA, 2*fdiss*bL2c/CMass)
            end if
         else
            !---> jsasaki 2026-10-17: carbonate Step 2 CARB_OMEGA=pore: each layer's share of the stock dissolves at its own
            ! pore-water calcite saturation; products into that layer's pore DIC and TA (nothing to the water directly)
            _GET_(self%id_ETW_p, ETW_p)
            _GET_(self%id_X1X_p, X1X_p)
            _GET_(self%id_dens_p, dens_p)
            _GET_HORIZONTAL_(self%id_poro_p, poro_p)
            do l_p = 1, 3
               _GET_HORIZONTAL_(self%id_Dm_p(l_p), Dm_p(l_p))
            end do
            R_p = 0.0_rk; om_p = 1.0_rk; nfail_p = 0.0_rk; top_p = 0.0_rk
            do l_p = 1, 3
               thk_p = Dm_p(l_p) - top_p
               w_p = max(0.0_rk, min(Dm_p(l_p), self%z_cal) - min(top_p, self%z_cal)) / self%z_cal
               top_p = Dm_p(l_p)
               vol_p = poro_p * thk_p
               if (.not. (vol_p > 1.0e-9_rk)) cycle
               _GET_HORIZONTAL_(self%id_pDIC(l_p), qDIC)
               _GET_HORIZONTAL_(self%id_pTA(l_p), qTA)
               _GET_HORIZONTAL_(self%id_pH2S(l_p), qH2S)
               _GET_HORIZONTAL_(self%id_pNH4(l_p), qNH4)
               _GET_HORIZONTAL_(self%id_pPO4(l_p), qPO4)
               _GET_HORIZONTAL_(self%id_pSi(l_p), qSi)
               ! layer inventories (mmol m-2) over the layer's pore volume -> mmol m-3 -> mol kg-1 (water density)
               call carbonate_engine_solve_porewater(ETW_p, X1X_p, 0.0_rk, &
                    max(qDIC, 0.0_rk)/vol_p/1.0e3_rk/dens_p, qTA/vol_p/1.0e3_rk/dens_p, &
                    max(qPO4, 0.0_rk)/vol_p/1.0e3_rk/dens_p, max(qSi, 0.0_rk)/vol_p/1.0e3_rk/dens_p, &
                    max(qH2S, 0.0_rk)/vol_p/1.0e3_rk/dens_p, max(qNH4, 0.0_rk)/vol_p/1.0e3_rk/dens_p, &
                    1, 10, 2, pH_p, fco2_p, h2co3_p, hco3_p, co3_p, k0_p, ok_p)
               if (.not. ok_p) then
                  nfail_p = nfail_p + 1.0_rk
                  cycle
               end if
               call CaCO3_Saturation(ETW_p, X1X_p, 0.0_rk, co3_p, om_p(l_p), om_arg_p, calcium_of_salinity(X1X_p))
               R_p(l_p) = self%fdissmax * (max(1.0_rk - om_p(l_p), 0.0_rk))**self%ndiss * w_p * bL2c
            end do
            fdiss = 0.0_rk
            if (bL2c > 0.0_rk) fdiss = sum(R_p) / bL2c
            _SET_BOTTOM_ODE_(self%id_c, -sum(R_p))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dissolution, -sum(R_p))
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pfail, nfail_p)
            do l_p = 1, 3
               _SET_BOTTOM_ODE_(self%id_pDIC(l_p), R_p(l_p)/CMass)
               _SET_BOTTOM_ODE_(self%id_pTA(l_p), 2.0_rk*R_p(l_p)/CMass)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_om_pore(l_p), om_p(l_p))
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pdiss(l_p), R_p(l_p))
            end do
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_diss_c, 0.0_rk)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_diss_O3c, 0.0_rk)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_diss_TA, 0.0_rk)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_pdiss_c, -sum(R_p))
               do l_p = 1, 3
                  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_pdiss_benDIC(l_p), R_p(l_p)/CMass)
                  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_pdiss_benTA(l_p), 2.0_rk*R_p(l_p)/CMass)
               end do
            end if
            !<--- jsasaki 2026-10-17
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
            ! jsasaki 2026-10-11: Stage Z bed light switch
            if (self%isw_par_bot == 1) then
               _GET_HORIZONTAL_(self%id_par_bot, par)
            else
               _GET_(self%id_par, par)
            end if
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
         ! E5 C1 acid-promoted dissolution (inert when k_co2diss = 0)
         F_co2 = 0.0_rk
         if (self%k_co2diss > 0.0_rk) then
            _GET_(self%id_CO2aq, co2aq)
            F_co2 = self%k_co2diss * bL2c * max(co2aq, 0.0_rk) / self%co2_ref_diss
            _SET_BOTTOM_ODE_(self%id_c, -F_co2)
            _SET_BOTTOM_EXCHANGE_(self%id_O3c, F_co2/CMass)
            _SET_BOTTOM_EXCHANGE_(self%id_TA, 2*F_co2/CMass)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_co2_dissolution, F_co2)
         end if
         if (self%isw_ledger == 1) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_co2diss_c, -F_co2)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_co2diss_O3c, F_co2/CMass)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_co2diss_TA, 2*F_co2/CMass)
         end if
         ! C2 reoxidation-driven dissolution (docs/129; inert when f_neut = 0)
         if (self%f_neut > 0.0_rk) then
            _GET_HORIZONTAL_(self%id_acid_S, acid_S)
            _GET_HORIZONTAL_(self%id_nrate, nrate)
            R_ac = self%f_neut * (max(acid_S, 0.0_rk) + 2.0_rk * max(nrate, 0.0_rk)) &
                 * max(bL2c, 0.0_rk) / (max(bL2c, 0.0_rk) + self%c_half_neut)
            _SET_BOTTOM_ODE_(self%id_c, -R_ac * CMass)
            _SET_BOTTOM_ODE_(self%id_benTA1, 2.0_rk * R_ac)
            _SET_BOTTOM_ODE_(self%id_benDIC1, R_ac)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_acid_dissolution, R_ac)
            if (self%isw_ledger == 1) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_acdiss_c, -R_ac * CMass)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_acdiss_benTA, 2.0_rk * R_ac)
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ledger_acdiss_benDIC, R_ac)
            end if
         end if

      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

end module
