!-----------------------------------------------------------------------
! carbonate_engine.F90
!
! PyCO2SYS-compatible carbonate chemistry solver for ERSEM
! Solves DIC + TA -> pH + pCO2 using robust root-finding (Brent's method)
!
! This module provides an alternative carbonate chemistry engine that
! reproduces PyCO2SYS results without runtime Python dependencies.
!
! pH scales: opt_pH_scale = 1(Total), 2(SWS), 3(Free), 4(NBS)
! K1/K2:     opt_k_carbonic uses PyCO2SYS numbering:
!              10 = Lueker et al. (2000), Total scale
!              14 = Millero (2010), SWS scale
!            legacy_mode uses Mehrbach/DM87 on SWS scale
!
! Other constants:
!   K0: Weiss (1974)
!   KB: Dickson (1990) via Millero (1995), Total scale
!   KW: Millero (1995), SWS scale
!   KS: Dickson (1990), Free scale
!   KF: Perez & Fraga (1987), Free scale (native)
!   Total boron: Uppstrom (1974) or Lee et al. (2010)
!   fH: Takahashi-Williams-Benson (1982) for NBS scale
!
! Author: Claude Code
! Date: 2026
!-----------------------------------------------------------------------

module carbonate_engine

   use fabm_types, only: rk

   implicit none

   private

   ! Numeric tolerances for solver
   real(rk), parameter :: tol_ph = 1.0e-8_rk      ! pH convergence tolerance
   real(rk), parameter :: tol_f  = 1.0e-12_rk     ! Function value tolerance (mol/kg)
   integer,  parameter :: max_iter = 100          ! Maximum iterations

   ! pH bracket for root-finding
   real(rk), parameter :: pH_lo_default = 2.0_rk  ! Lower pH bound
   real(rk), parameter :: pH_hi_default = 12.0_rk ! Upper pH bound

   ! Universal gas constant (cm3 bar / mol K)
   real(rk), parameter :: Rgas = 83.14472_rk

   public :: carbonate_engine_solve

contains

   !-----------------------------------------------------------------------
   ! carbonate_engine_solve
   !
   ! Main solver interface: compute pH and pCO2 from DIC and TA
   !
   ! Inputs:
   !   T               - Temperature (degC)
   !   S               - Practical salinity (PSU)
   !   Pbar            - Pressure (bar), 0 for surface
   !   DIC_molkg       - Dissolved inorganic carbon (mol/kg)
   !   TA_molkg        - Total alkalinity (mol/kg)
   !   opt_pH_scale    - pH scale: 1=Total, 2=SWS, 3=Free, 4=NBS
   !   opt_k_carbonic  - K1/K2 (PyCO2SYS numbering): 10=Lueker, 14=Millero
   !   opt_total_borate - Total boron: 1=Uppstrom 1974, 2=Lee 2010
   !   legacy_mode     - .true. for backward-compat CO2DYN behavior
   !
   ! Outputs:
   !   pH        - pH on requested scale
   !   pCO2_atm  - Partial pressure of CO2 (atm)
   !   H2CO3     - Carbonic acid concentration (mol/kg)
   !   HCO3      - Bicarbonate concentration (mol/kg)
   !   CO3       - Carbonate concentration (mol/kg)
   !   K0        - Henry's law constant (mol/kg/atm)
   !   success   - .true. if solver converged
   !-----------------------------------------------------------------------
   subroutine carbonate_engine_solve(T, S, Pbar, DIC_molkg, TA_molkg, &
                                     opt_pH_scale, opt_k_carbonic, &
                                     opt_total_borate, legacy_mode, &
                                     pH, pCO2_atm, H2CO3, HCO3, CO3, &
                                     K0, success)
      real(rk), intent(in)  :: T, S, Pbar
      real(rk), intent(in)  :: DIC_molkg, TA_molkg
      integer,  intent(in)  :: opt_pH_scale
      integer,  intent(in)  :: opt_k_carbonic, opt_total_borate
      logical,  intent(in)  :: legacy_mode
      real(rk), intent(out) :: pH, pCO2_atm, H2CO3, HCO3, CO3, K0
      logical,  intent(out) :: success

      ! Equilibrium constants (on target pH scale after conversion)
      real(rk) :: K1, K2, KB, KW

      ! Concentrations
      real(rk) :: BT, ST, FT

      ! Solver variables
      real(rk) :: H_solution
      real(rk) :: Tclamp
      logical  :: include_OH_H

      ! Temperature floor: allow sub-zero for PyCO2SYS compatibility,
      ! but match legacy CO2DYN behavior for backward-compat mode
      if (legacy_mode) then
         Tclamp = max(T, 0.0_rk)
      else
         Tclamp = max(T, -2.0_rk)
      end if

      ! Legacy mode omits OH-H from alkalinity (matching CO2CLC)
      include_OH_H = (.not. legacy_mode)

      ! Calculate total concentrations from salinity
      call calc_total_concentrations(S, opt_total_borate, BT, ST, FT)

      ! Calculate equilibrium constants on target pH scale
      call calc_equilibrium_constants(Tclamp, S, Pbar, opt_pH_scale, &
                                      opt_k_carbonic, legacy_mode, &
                                      ST, FT, K0, K1, K2, KB, KW)

      ! Solve for [H+] using Brent's method
      call solve_pH_brent(DIC_molkg, TA_molkg, K1, K2, KB, KW, BT, &
                          H_solution, success, include_OH_H)

      if (.not. success) then
         pH = 0.0_rk
         pCO2_atm = 0.0_rk
         H2CO3 = 0.0_rk
         HCO3 = 0.0_rk
         CO3 = 0.0_rk
         return
      end if

      ! H_solution is on the target scale because all equilibrium
      ! constants were converted to target scale before solving
      pH = -log10(H_solution)

      ! Calculate carbonate species from [H+] and DIC
      call calc_carbonate_species(DIC_molkg, H_solution, K0, K1, K2, &
                                  H2CO3, HCO3, CO3, pCO2_atm)

   end subroutine carbonate_engine_solve

   !-----------------------------------------------------------------------
   ! calc_equilibrium_constants
   !
   ! Calculate all equilibrium constants on target pH scale.
   ! Strategy: compute each K on its native scale, apply pressure
   ! correction, then convert once to target scale via convert_K.
   !
   ! Native scales:
   !   KS  - Free (Dickson 1990)
   !   KF  - Free (Perez & Fraga 1987, native free scale)
   !   KB  - Total (Dickson 1990 via Millero 1995)
   !   KW  - SWS (Millero 1995)
   !   K1,K2 (opt=10) - Total (Lueker 2000)
   !   K1,K2 (opt=14) - SWS (Millero 2010)
   !   K1,K2 (legacy) - SWS (Mehrbach/DM87)
   !-----------------------------------------------------------------------
   subroutine calc_equilibrium_constants(T, S, Pbar, opt_pH_scale, &
                                         opt_k_carbonic, legacy_mode, &
                                         ST, FT, K0, K1, K2, KB, KW)
      real(rk), intent(in)  :: T, S, Pbar
      integer,  intent(in)  :: opt_pH_scale
      integer,  intent(in)  :: opt_k_carbonic
      logical,  intent(in)  :: legacy_mode
      real(rk), intent(in)  :: ST, FT
      real(rk), intent(out) :: K0, K1, K2, KB, KW

      real(rk) :: TK, TK100, lnTK, invTK
      real(rk) :: sqrtS, S15, S2
      real(rk) :: IS, sqrtIS
      real(rk) :: KS, KF_free
      real(rk) :: delta, kappa
      real(rk) :: free2total, free2sws, fH
      real(rk) :: pK1, pK2
      integer  :: k1k2_native

      ! Temperature in Kelvin
      TK = T + 273.15_rk
      TK100 = TK / 100.0_rk
      lnTK = log(TK)
      invTK = 1.0_rk / TK

      ! Salinity terms
      sqrtS = sqrt(S)
      S15 = S ** 1.5_rk
      S2 = S * S

      ! Ionic strength (Millero, 1982)
      IS = 19.924_rk * S / (1000.0_rk - 1.005_rk * S)
      sqrtIS = sqrt(IS)

      !-------------------------------------------------------------------
      ! KS = [H+][SO4--]/[HSO4-] on free scale
      ! Dickson (1990)
      !-------------------------------------------------------------------
      KS = exp(-4276.1_rk * invTK + 141.328_rk - 23.093_rk * lnTK &
               + (-13856.0_rk * invTK + 324.57_rk - 47.986_rk * lnTK) * sqrtIS &
               + (35474.0_rk * invTK - 771.54_rk + 114.723_rk * lnTK) * IS &
               - 2698.0_rk * invTK * IS**1.5_rk + 1776.0_rk * invTK * IS**2.0_rk &
               + log(1.0_rk - 0.001005_rk * S))

      ! Pressure correction for KS (on free scale)
      delta = -18.03_rk + 0.0466_rk * T + 0.000316_rk * T**2.0_rk
      kappa = -4.53_rk + 0.00009_rk * T
      KS = KS * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))

      !-------------------------------------------------------------------
      ! KF = [H+][F-]/[HF] on free scale
      ! Perez & Fraga (1987) - directly on free scale
      ! (PyCO2SYS kHF_FREE_PF87: no total->free conversion needed)
      !-------------------------------------------------------------------
      KF_free = exp(874.0_rk * invTK - 9.68_rk + 0.111_rk * sqrtS)

      ! Pressure correction for KF (on free scale)
      delta = -9.78_rk - 0.009_rk * T - 0.000942_rk * T**2.0_rk
      kappa = -3.91_rk + 0.000054_rk * T
      KF_free = KF_free * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar &
                               / (Rgas * TK))

      !-------------------------------------------------------------------
      ! Scale conversion factors at depth (using depth-corrected KS, KF)
      !-------------------------------------------------------------------
      free2total = 1.0_rk + ST / KS
      free2sws   = 1.0_rk + ST / KS + FT / KF_free

      ! Activity coefficient for NBS scale (Takahashi et al. 1982)
      fH = calc_fH(TK, S)

      !-------------------------------------------------------------------
      ! K0 = [CO2*]/pCO2 (mol/kg/atm) - scale-independent
      ! Weiss (1974)
      !-------------------------------------------------------------------
      K0 = exp(93.4517_rk / TK100 - 60.2409_rk + 23.3585_rk * log(TK100) &
               + S * (0.023517_rk - 0.023656_rk * TK100 &
                      + 0.0047036_rk * TK100**2.0_rk))

      if (Pbar > 0.0_rk) then
         K0 = K0 * exp((-Pbar * 32.3_rk) / (Rgas * TK))
      end if

      !-------------------------------------------------------------------
      ! KB = [H+][BO2-]/[HBO2] on total scale (native)
      ! Dickson (1990) via Millero (1995)
      !-------------------------------------------------------------------
      KB = exp((-8966.9_rk - 2890.53_rk * sqrtS - 77.942_rk * S &
                + 1.728_rk * S15 - 0.0996_rk * S2) / TK &
               + (148.0248_rk + 137.1942_rk * sqrtS + 1.62142_rk * S) &
               + (-24.4344_rk - 25.085_rk * sqrtS - 0.2474_rk * S) * lnTK &
               + 0.053105_rk * sqrtS * TK)

      ! Pressure correction (on native total scale)
      delta = -29.48_rk + 0.1622_rk * T - 0.002608_rk * T**2.0_rk
      kappa = -2.84_rk / 1000.0_rk
      KB = KB * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))

      ! Convert KB: total(1) -> target scale
      if (.not. legacy_mode) then
         call convert_K(KB, 1, opt_pH_scale, free2total, free2sws, fH)
      end if
      ! legacy_mode: keep KB on total (intentional: reproduces old bug)

      !-------------------------------------------------------------------
      ! K1, K2: Carbonic acid dissociation constants
      !-------------------------------------------------------------------
      if (legacy_mode) then
         ! Mehrbach (1973) refit by Dickson & Millero (1987) on SWS scale
         k1k2_native = 2  ! SWS
         pK1 = 3670.7_rk / TK - 62.008_rk + 9.7944_rk * lnTK &
               - 0.0118_rk * S + 0.000116_rk * S2
         K1 = 10.0_rk ** (-pK1)

         pK2 = 1394.7_rk / TK + 4.777_rk - 0.0184_rk * S + 0.000118_rk * S2
         K2 = 10.0_rk ** (-pK2)
      else
         select case (opt_k_carbonic)
         case (14)
            ! Millero (2010) on SWS scale (PyCO2SYS #14)
            k1k2_native = 2  ! SWS
            pK1 = -126.34048_rk + 6320.813_rk * invTK &
                  + 19.568224_rk * lnTK &
                  + (13.4038_rk * sqrtS + 0.03206_rk * S &
                     - 0.00005242_rk * S2) &
                  + (-530.659_rk * sqrtS - 5.8210_rk * S) * invTK &
                  + (-2.0664_rk * sqrtS) * lnTK
            K1 = 10.0_rk ** (-pK1)

            pK2 = -90.18333_rk + 5143.692_rk * invTK &
                  + 14.613358_rk * lnTK &
                  + (21.3728_rk * sqrtS + 0.1218_rk * S &
                     - 0.0003688_rk * S2) &
                  + (-788.289_rk * sqrtS - 19.189_rk * S) * invTK &
                  + (-3.374_rk * sqrtS) * lnTK
            K2 = 10.0_rk ** (-pK2)
         case default
            ! Lueker et al. (2000) on Total scale (PyCO2SYS #10, default)
            k1k2_native = 1  ! Total
            pK1 = 3633.86_rk * invTK - 61.2172_rk + 9.6777_rk * lnTK &
                  - 0.011555_rk * S + 0.0001152_rk * S2
            K1 = 10.0_rk ** (-pK1)

            pK2 = 471.78_rk * invTK + 25.929_rk - 3.16967_rk * lnTK &
                  - 0.01781_rk * S + 0.0001122_rk * S2
            K2 = 10.0_rk ** (-pK2)
         end select
      end if

      ! Pressure corrections for K1, K2 (on native scale)
      delta = -25.5_rk + 0.1271_rk * T
      kappa = (-3.08_rk + 0.0877_rk * T) / 1000.0_rk
      K1 = K1 * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))

      delta = -15.82_rk - 0.0219_rk * T
      kappa = (1.13_rk - 0.1475_rk * T) / 1000.0_rk
      K2 = K2 * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))

      ! Convert K1, K2: native -> target scale
      if (.not. legacy_mode) then
         call convert_K(K1, k1k2_native, opt_pH_scale, &
                        free2total, free2sws, fH)
         call convert_K(K2, k1k2_native, opt_pH_scale, &
                        free2total, free2sws, fH)
      end if
      ! legacy_mode: keep K1, K2 on SWS (no conversion)

      !-------------------------------------------------------------------
      ! KW = [H+][OH-] on SWS scale (native)
      ! Millero (1995)
      !-------------------------------------------------------------------
      KW = exp(148.9802_rk - 13847.26_rk * invTK - 23.6521_rk * lnTK &
               + (-5.977_rk + 118.67_rk * invTK + 1.0495_rk * lnTK) * sqrtS &
               - 0.01615_rk * S)

      ! Pressure correction (on native SWS scale)
      delta = -25.60_rk + 0.2324_rk * T - 0.0036246_rk * T**2.0_rk
      kappa = (-5.13_rk + 0.0794_rk * T) / 1000.0_rk
      KW = KW * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))

      ! Convert KW: SWS(2) -> target scale
      if (.not. legacy_mode) then
         call convert_K(KW, 2, opt_pH_scale, free2total, free2sws, fH)
      end if
      ! legacy_mode: keep KW on SWS (no conversion)

   end subroutine calc_equilibrium_constants

   !-----------------------------------------------------------------------
   ! calc_fH - Activity coefficient for NBS scale
   ! Takahashi-Williams-Benson (1982)
   !-----------------------------------------------------------------------
   function calc_fH(TK, S) result(fH)
      real(rk), intent(in) :: TK, S
      real(rk) :: fH

      fH = 1.2948_rk - 0.002036_rk * TK &
           + (0.0004607_rk - 0.000001475_rk * TK) * S * S
   end function calc_fH

   !-----------------------------------------------------------------------
   ! convert_K - Convert equilibrium constant between pH scales
   !
   ! Route: from_scale -> free -> to_scale
   !   from_scale: 1=Total, 2=SWS, 3=Free (native scale of constant)
   !   to_scale:   1=Total, 2=SWS, 3=Free, 4=NBS (target scale)
   !   free2total: 1 + ST/KS
   !   free2sws:   1 + ST/KS + FT/KF
   !   fH:         activity coefficient for NBS
   !-----------------------------------------------------------------------
   subroutine convert_K(K, from_scale, to_scale, free2total, free2sws, fH)
      real(rk), intent(inout) :: K
      integer,  intent(in)    :: from_scale, to_scale
      real(rk), intent(in)    :: free2total, free2sws, fH

      ! No conversion needed if scales match
      if (from_scale == to_scale) return

      ! Step 1: Convert from native scale to free scale
      select case (from_scale)
      case (1)  ! Total -> Free
         K = K / free2total
      case (2)  ! SWS -> Free
         K = K / free2sws
      case (3)  ! Already Free
         ! no-op
      end select

      ! Step 2: Convert from free scale to target scale
      select case (to_scale)
      case (1)  ! Free -> Total
         K = K * free2total
      case (2)  ! Free -> SWS
         K = K * free2sws
      case (3)  ! Stay Free
         ! no-op
      case (4)  ! Free -> NBS: K_NBS = K_free * free2sws * fH
         ! NBS via SWS: H_NBS = fH * H_SWS = fH * free2sws * H_free
         ! So K_NBS = K_free * free2sws * fH
         K = K * free2sws * fH
      end select

   end subroutine convert_K

   !-----------------------------------------------------------------------
   ! calc_total_concentrations
   !
   ! Calculate total boron, sulfate, and fluoride from salinity
   !
   ! opt_total_borate: 1 = Uppstrom (1974), 2 = Lee et al. (2010)
   !-----------------------------------------------------------------------
   subroutine calc_total_concentrations(S, opt_total_borate, BT, ST, FT)
      real(rk), intent(in)  :: S
      integer,  intent(in)  :: opt_total_borate
      real(rk), intent(out) :: BT, ST, FT

      real(rk) :: Cl

      ! Chloride from salinity
      Cl = S / 1.80655_rk

      ! Total boron (mol/kg)
      if (opt_total_borate == 1) then
         ! Uppstrom (1974)
         BT = 0.000416_rk * S / 35.0_rk
      else
         ! Lee et al. (2010) - default for opt_total_borate == 2
         BT = 0.0004326_rk * S / 35.0_rk
      end if

      ! Total sulfate (mol/kg) - Morris & Riley (1966)
      ST = 0.14_rk * Cl / 96.062_rk

      ! Total fluoride (mol/kg) - Riley (1965)
      FT = 0.000067_rk * Cl / 18.998_rk

   end subroutine calc_total_concentrations

   !-----------------------------------------------------------------------
   ! alkalinity_residual
   !
   ! Calculate the difference between modeled TA and input TA
   ! F(H) = TA_model - TA_input
   !
   ! TA_model = [HCO3-] + 2[CO3--] + [B(OH)4-] + [OH-] - [H+]
   !-----------------------------------------------------------------------
   function alkalinity_residual(H, DIC, TA_input, K1, K2, KB, KW, BT, &
                                include_OH_H) result(F)
      real(rk), intent(in) :: H, DIC, TA_input, K1, K2, KB, KW, BT
      logical,  intent(in) :: include_OH_H
      real(rk) :: F

      real(rk) :: denom, HCO3, CO3, BOH4
      real(rk) :: TA_model

      ! Carbonate species from DIC and [H+]
      denom = H * H + K1 * H + K1 * K2
      HCO3 = DIC * K1 * H / denom
      CO3 = DIC * K1 * K2 / denom

      ! Borate alkalinity
      BOH4 = BT * KB / (KB + H)

      ! Total alkalinity
      TA_model = HCO3 + 2.0_rk * CO3 + BOH4
      if (include_OH_H) then
         TA_model = TA_model + KW / H - H
      end if

      F = TA_model - TA_input

   end function alkalinity_residual

   !-----------------------------------------------------------------------
   ! solve_pH_brent
   !
   ! Solve for [H+] using Brent's method
   ! This is a robust root-finding algorithm that combines bisection,
   ! secant, and inverse quadratic interpolation
   !-----------------------------------------------------------------------
   subroutine solve_pH_brent(DIC, TA_input, K1, K2, KB, KW, BT, &
                             H_solution, success, include_OH_H)
      real(rk), intent(in)  :: DIC, TA_input, K1, K2, KB, KW, BT
      real(rk), intent(out) :: H_solution
      logical,  intent(out) :: success
      logical,  intent(in)  :: include_OH_H

      real(rk) :: a, b, c, d, e
      real(rk) :: fa, fb, fc
      real(rk) :: p, q, r, s
      real(rk) :: tol1, xm
      real(rk) :: H_lo, H_hi
      integer  :: iter

      ! Initialize with pH bounds converted to [H+]
      H_lo = 10.0_rk ** (-pH_hi_default)  ! Low [H+] = high pH
      H_hi = 10.0_rk ** (-pH_lo_default)  ! High [H+] = low pH

      a = H_lo
      b = H_hi
      fa = alkalinity_residual(a, DIC, TA_input, K1, K2, KB, KW, BT, &
                              include_OH_H)
      fb = alkalinity_residual(b, DIC, TA_input, K1, K2, KB, KW, BT, &
                              include_OH_H)

      ! Check if root is bracketed
      if (fa * fb > 0.0_rk) then
         ! Try to find a valid bracket by adjusting bounds
         ! First try wider pH range
         a = 10.0_rk ** (-14.0_rk)
         b = 10.0_rk ** (-1.0_rk)
         fa = alkalinity_residual(a, DIC, TA_input, K1, K2, KB, KW, BT, &
                              include_OH_H)
         fb = alkalinity_residual(b, DIC, TA_input, K1, K2, KB, KW, BT, &
                              include_OH_H)

         if (fa * fb > 0.0_rk) then
            success = .false.
            H_solution = 0.0_rk
            return
         end if
      end if

      ! Initialize Brent's method: c holds the previous iterate
      c = a
      fc = fa
      d = b - a
      e = d

      success = .false.

      do iter = 1, max_iter
         ! Ensure b is the best approximation (closest to root)
         if (abs(fc) < abs(fb)) then
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
         end if

         ! Convergence check
         tol1 = 2.0_rk * epsilon(1.0_rk) * abs(b) + 0.5_rk * tol_ph * b
         xm = 0.5_rk * (c - b)

         if (abs(xm) <= tol1 .or. abs(fb) < tol_f) then
            H_solution = b
            success = .true.
            return
         end if

         ! Try inverse quadratic interpolation or secant
         if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
            s = fb / fa
            if (abs(a - c) < epsilon(1.0_rk) * max(abs(a), abs(c))) then
               ! Linear interpolation (secant method)
               p = 2.0_rk * xm * s
               q = 1.0_rk - s
            else
               ! Inverse quadratic interpolation
               q = fa / fc
               r = fb / fc
               p = s * (2.0_rk * xm * q * (q - r) - (b - a) * (r - 1.0_rk))
               q = (q - 1.0_rk) * (r - 1.0_rk) * (s - 1.0_rk)
            end if

            ! Check if interpolation is acceptable
            if (p > 0.0_rk) then
               q = -q
            else
               p = -p
            end if

            s = e
            e = d

            if (2.0_rk * p < 3.0_rk * xm * q - abs(tol1 * q) .and. &
                p < abs(0.5_rk * s * q)) then
               d = p / q
            else
               d = xm
               e = d
            end if
         else
            ! Bisection
            d = xm
            e = d
         end if

         ! Update a to previous b
         a = b
         fa = fb

         ! Update b
         if (abs(d) > tol1) then
            b = b + d
         else
            if (xm >= 0.0_rk) then
               b = b + tol1
            else
               b = b - tol1
            end if
         end if

         fb = alkalinity_residual(b, DIC, TA_input, K1, K2, KB, KW, BT, &
                              include_OH_H)

         ! Ensure bracket is maintained
         if ((fb > 0.0_rk .and. fc > 0.0_rk) .or. &
             (fb < 0.0_rk .and. fc < 0.0_rk)) then
            c = a
            fc = fa
            e = b - a
            d = e
         end if
      end do

      ! Did not converge within max_iter
      success = .false.
      H_solution = 0.0_rk

   end subroutine solve_pH_brent

   !-----------------------------------------------------------------------
   ! calc_carbonate_species
   !
   ! Calculate carbonate species concentrations from [H+] and DIC
   !-----------------------------------------------------------------------
   subroutine calc_carbonate_species(DIC, H, K0, K1, K2, H2CO3, HCO3, CO3, pCO2_atm)
      real(rk), intent(in)  :: DIC, H, K0, K1, K2
      real(rk), intent(out) :: H2CO3, HCO3, CO3, pCO2_atm

      real(rk) :: denom

      ! Denominator for species fractions
      denom = H * H + K1 * H + K1 * K2

      ! Species concentrations (mol/kg)
      H2CO3 = DIC * H * H / denom        ! [CO2*] = [H2CO3] + [CO2(aq)]
      HCO3 = DIC * K1 * H / denom        ! [HCO3-]
      CO3 = DIC * K1 * K2 / denom        ! [CO3--]

      ! pCO2 in atm
      if (K0 > 0.0_rk) then
         pCO2_atm = H2CO3 / K0
      else
         pCO2_atm = 0.0_rk
      end if

   end subroutine calc_carbonate_species

end module carbonate_engine
