! Shared iron / iron-sulfide laws of the ERSEM-MUSE unification, family D wave 2 part 2 (Fe2/Fe3/FeS/FeS2 retention and
! reoxidation).  jsasaki 2026-10-07 (unification plan annex D rows 16-19; muse/docs/UNIFY_D_SPEC_20261007.md section 9).
!
! THIS FILE EXISTS IN TWO REPOSITORIES WITH IDENTICAL CONTENT:
!   muse/src/benthos_iron_laws.F90   and   ersem/src/benthos_iron_laws.F90
! The checksums are compared by muse/tests/unify_d3/check_sync.sh; edit both or neither.  Every routine is pure and has no
! state.  The laws are those of muse/src/sed_network_c.F90 (reactions R9-R14, R16 and the Fe(III) and sulfate branches of the
! organic-matter cascade), written with the SAME operation order so that the equivalence test (muse/tests/unify_d3) compares
! like with like.  MUSE keeps its own copy of the expressions in sed_network_c (its default behaviour is not touched); the
! ERSEM benthic_sulfur_cycle (isw_fe = 1) calls this module.
!
! Units: per BULK volume (mmol m-3 d-1) for the extents; solutes enter as pore-water (phase) concentrations, solids and the
! sorbing Fe2 as bulk totals; temperatures in degC; the Q10 factor is th_q10 of benthos_shared_laws (Tref 20 degC).
module benthos_iron_laws
   use, intrinsic :: iso_fortran_env, only: real64
   use benthos_shared_laws, only: th_q10
   implicit none
   private
   integer, parameter :: rk = real64

   ! reactions of one cell (extent = one unit of the first-named species per extent, MUSE sed_network_c numbering in brackets)
   integer, parameter, public :: nfe = 8
   integer, parameter, public :: jFEOX = 1, jISP = 2, jFESOX = 3, jPYR = 4, jPYOX = 5, jSFE = 6, jISD = 7, jFE3RED = 8
   ! species of the sources matrix
   integer, parameter, public :: nsp_fe = 9
   integer, parameter, public :: sO2 = 1, sH2S = 2, sS0 = 3, sFe2 = 4, sFe3 = 5, sFeS = 6, sFeS2 = 7, sSO4 = 8, sTA = 9
   character(len=6), parameter, public :: fe_reac_name(nfe) = [character(len=6) :: 'R9', 'R10', 'R11', 'R12', 'R13', 'R14', 'R16', 'R3']
   character(len=4), parameter, public :: fe_spec_name(nsp_fe) = [character(len=4) :: 'O2', 'H2S', 'S0', 'Fe2', 'Fe3', 'FeS', 'FeS2', &
      'SO4', 'TA']

   type, public :: iron_params
      ! 20 degC rate constants, MUSE sed_network_c names and defaults (units as there: (mmol m-3)-1 d-1 on phase concentrations
      ! unless stated)
      real(rk) :: k_feox = 3.5_rk, k_fesox = 2.2e-3_rk, k_pyox = 3.0e-4_rk, k_sfe = 1.2e-4_rk
      real(rk) :: k_pyr0 = 1.73e-4_rk                      ! k_pyr(phi) = k_pyr0 (1 - phi)
      real(rk) :: k_isp = 2.74e4_rk, k_isd = 8.21e-3_rk    ! mmol m-3 solid d-1; 1/d
      real(rk) :: K_c = 130._rk                            ! (mmol m-3)^2 at the declared pH 7.5
      real(rk) :: Ks_Fe2 = 2.65_rk * 268._rk               ! sorption, rho_s K_D per solid volume
      ! Q10 in the order of the reactions jFEOX .. jISD (temperature table: Q10_fe2, Q10_isp, Q10_fes, Q10_pyr, Q10_pyox, Q10_sfe, Q10_isd)
      real(rk) :: q10(7) = [4.4_rk, 1.0_rk, 2.0_rk, 2.8_rk, 2.3_rk, 1.0_rk, 1.0_rk]
      ! oxidant cascade constants of the Fe(III) / sulfate branch (MUSE K_Fe3, K_SO4, Kin_*)
      real(rk) :: K_Fe3 = 5.e4_rk, K_SO4 = 1600._rk, Kin_O2 = 5._rk, Kin_NO3 = 5._rk, Kin_Fe3 = 5.e4_rk
   end type

   public :: iron_extents, iron_stoich, fe3_share, fe2_capacity, sulfur_ox_extents

contains

   ! the bulk capacity of Fe2 per unit dissolved concentration (dissolved + sorbed), phi + (1-phi) Ks
   pure real(rk) function fe2_capacity(p, phi) result(f)
      type(iron_params), intent(in) :: p
      real(rk), intent(in) :: phi
      f = phi + (1._rk - phi) * p%Ks_Fe2
   end function

   ! Reaction extents r(1:7) (mmol m-3 bulk d-1) of the seven iron-sulfur reactions at one state.
   !   cFe2 (total, dissolved + sorbed), cFe3, cFeS, cFeS2, cS0: bulk concentrations (mmol m-3 of bulk sediment)
   !   qH2S, qO2: pore-water (phase) concentrations (mmol m-3)
   ! r(jFE3RED) is not an extent of this routine (it is a fraction of the carbon flux, see fe3_share) and is returned 0.
   pure subroutine iron_extents(p, T, phi, cFe2, cFe3, cFeS, cFeS2, cS0, qH2S, qO2, r)
      type(iron_params), intent(in) :: p
      real(rk), intent(in) :: T, phi, cFe2, cFe3, cFeS, cFeS2, cS0, qH2S, qO2
      real(rk), intent(out) :: r(nfe)
      real(rk) :: th(7), sol, qFe2, qFe3, qFeS, qFeS2, qS0, cfe2b
      th = th_q10(T, p%q10, 20._rk)
      sol = 1._rk - phi
      cfe2b = max(cFe2, 0._rk)
      qFe2 = cfe2b / fe2_capacity(p, phi)
      qFe3 = max(cFe3, 0._rk) / sol
      qFeS = max(cFeS, 0._rk) / sol
      qFeS2 = max(cFeS2, 0._rk) / sol
      qS0 = max(cS0, 0._rk) / sol
      r = 0._rk
      r(jFEOX) = th(1) * p%k_feox * cfe2b * qO2                                                    ! R9
      r(jISP) = th(2) * p%k_isp * sol * max(0._rk, qFe2 * qH2S / p%K_c - 1._rk)                    ! R10
      r(jFESOX) = th(3) * p%k_fesox * qFeS * qO2 * sol                                             ! R11
      r(jPYR) = th(4) * p%k_pyr0 * sol * qFeS * qS0 * sol                                          ! R12
      r(jPYOX) = th(5) * p%k_pyox * qFeS2 * qO2 * sol                                              ! R13
      r(jSFE) = th(6) * p%k_sfe * qH2S * qFe3 * sol                                                ! R14
      r(jISD) = th(7) * p%k_isd * max(cFeS, 0._rk) * max(0._rk, 1._rk - qFe2 * qH2S / p%K_c)       ! R16
   end subroutine

   ! Aerobic sulfide and S0 oxidation of MUSE (R6, R7, R8; Q10 table: Q10_hs, Q10_s0) on pore-water concentrations, per bulk
   ! volume: r_aer = k_hs_tot qH2S qO2 phi;  r10 = th (1 - f_s0) r_aer (to SO4), r11 = th f_s0 r_aer (to S0),
   ! r12 = th k_s0ox qS0 qO2 (1 - phi) (cS0 bulk, qS0 = cS0/(1 - phi)).  The ERSEM benthic_sulfur_cycle (isw_fe = 1) uses
   ! these in place of the bottom-water-O2 / PAR gate.
   pure subroutine sulfur_ox_extents(k_hs_tot, f_s0, k_s0ox, q10_hs, q10_s0, T, phi, cS0, qH2S, qO2, r10, r11, r12)
      real(rk), intent(in) :: k_hs_tot, f_s0, k_s0ox, q10_hs, q10_s0, T, phi, cS0, qH2S, qO2
      real(rk), intent(out) :: r10, r11, r12
      real(rk) :: r_aer, sol, qS0
      sol = 1._rk - phi
      qS0 = max(cS0, 0._rk) / sol
      r_aer = k_hs_tot * qH2S * qO2 * phi
      r10 = th_q10(T, q10_hs, 20._rk) * (1._rk - f_s0) * r_aer
      r11 = th_q10(T, q10_hs, 20._rk) * f_s0 * r_aer
      r12 = th_q10(T, q10_s0, 20._rk) * k_s0ox * qS0 * qO2 * sol
   end subroutine

   ! The sources matrix S(species, reaction), hand-written as in sed_network_c (TA by the conservative-charge rule, Fe2 +2,
   ! SO4 -2).  Column jFE3RED is per mole of organic carbon oxidised by Fe(III) (CH2O + 4 Fe(III) -> CO2 + 4 Fe2, the
   ! DIC, N and P of the carbon are booked by the respiring module, as for sulfate reduction): Fe3 -4, Fe2 +4, TA +8.
   pure subroutine iron_stoich(s)
      real(rk), intent(out) :: s(nsp_fe, nfe)
      s = 0._rk
      ! R9 Fe2 + 0.25 O2 -> Fe3
      s(sFe2, jFEOX) = -1; s(sO2, jFEOX) = -0.25_rk; s(sFe3, jFEOX) = 1; s(sTA, jFEOX) = -2
      ! R10 Fe2 + H2S -> FeS
      s(sFe2, jISP) = -1; s(sH2S, jISP) = -1; s(sFeS, jISP) = 1; s(sTA, jISP) = -2
      ! R11 FeS + 2.25 O2 -> Fe3 + SO4
      s(sFeS, jFESOX) = -1; s(sO2, jFESOX) = -2.25_rk; s(sFe3, jFESOX) = 1; s(sSO4, jFESOX) = 1; s(sTA, jFESOX) = -2
      ! R12 FeS + S0 -> FeS2
      s(sFeS, jPYR) = -1; s(sS0, jPYR) = -1; s(sFeS2, jPYR) = 1
      ! R13 FeS2 + 3.75 O2 -> Fe3 + 2 SO4
      s(sFeS2, jPYOX) = -1; s(sO2, jPYOX) = -3.75_rk; s(sFe3, jPYOX) = 1; s(sSO4, jPYOX) = 2; s(sTA, jPYOX) = -4
      ! R14 H2S + 2 Fe3 -> S0 + 2 Fe2
      s(sH2S, jSFE) = -1; s(sFe3, jSFE) = -2; s(sS0, jSFE) = 1; s(sFe2, jSFE) = 2; s(sTA, jSFE) = 4
      ! R16 FeS -> Fe2 + H2S
      s(sFeS, jISD) = -1; s(sFe2, jISD) = 1; s(sH2S, jISD) = 1; s(sTA, jISD) = 2
      ! R3 (Fe(III) reduction by organic carbon, per mole C)
      s(sFe3, jFE3RED) = -4; s(sFe2, jFE3RED) = 4; s(sTA, jFE3RED) = 8
   end subroutine

   ! The share of the anaerobic carbon flux that runs on Fe(III) instead of sulfate when the two compete (MUSE cascade
   ! f(3), f(4) of sed_network_c before the normalisation; their ratio is what the normalisation leaves unchanged):
   !   f3 = qFe3/(K_Fe3 + qFe3) inO2 inNO3,   f4 = qSO4/(K_SO4 + qSO4) inO2 inNO3 inFe3,   share = f3/(f3 + f4).
   ! ERSEM's anaerobic respiration is controlled by the bacteria (a carbon flux that does not depend on the acceptor), so the
   ! cascade acts as a PARTITION of that flux; MUSE's additional slowing of the decay where sum(f) < 1 has no ERSEM
   ! counterpart.  qFe3 is the solid-phase Fe(III) concentration (bulk/(1-phi)), qSO4, qO2, qNO3 phase concentrations.
   pure real(rk) function fe3_share(p, qFe3, qSO4, qO2, qNO3) result(share)
      type(iron_params), intent(in) :: p
      real(rk), intent(in) :: qFe3, qSO4, qO2, qNO3
      real(rk) :: inO2, inNO3, inFe3, f3, f4, q3
      q3 = max(qFe3, 0._rk)
      inO2 = p%Kin_O2 / (p%Kin_O2 + max(qO2, 0._rk))
      inNO3 = p%Kin_NO3 / (p%Kin_NO3 + max(qNO3, 0._rk))
      inFe3 = p%Kin_Fe3 / (p%Kin_Fe3 + q3)
      f3 = q3 / (p%K_Fe3 + q3) * inO2 * inNO3
      f4 = max(qSO4, 0._rk) / (p%K_SO4 + max(qSO4, 0._rk)) * inO2 * inNO3 * inFe3
      if (f3 + f4 > 0._rk) then
         share = f3 / (f3 + f4)
      else
         share = 0._rk
      end if
   end function

end module benthos_iron_laws
