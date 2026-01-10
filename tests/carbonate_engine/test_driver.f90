!-----------------------------------------------------------------------
! test_driver.f90
!
! Offline verification harness for carbonate_engine module
!
! Reads reference cases from CSV file and compares against carbonate_engine
! solver results. Exits with non-zero status if any case exceeds tolerance.
!
! Tolerances:
!   pH:   1e-4 (absolute)
!   pCO2: 2 uatm = 2e-6 atm (absolute)
!
! Constants used (matching reference_cases.csv):
!   K0: Weiss (1974)
!   K1, K2: Mehrbach (1973) refit by Lueker et al. (2000) - Total scale
!   KB: Dickson (1990)
!   KW: Millero (1995) converted to Total scale
!   KS: Dickson (1990)
!   KF: Perez & Fraga (1987)
!   Total boron: Lee et al. (2010) - 432.6 * S / 35 umol/kg
!
! Usage:
!   ./test_driver [reference_cases.csv]
!
! If no argument given, defaults to "reference_cases.csv" in current directory.
!-----------------------------------------------------------------------

program test_driver
   use, intrinsic :: iso_fortran_env, only: real64, int32, error_unit, output_unit
   implicit none

   ! Use real64 to match typical double precision
   integer, parameter :: rk = real64

   ! Tolerances for comparison with reference data (very tight match)
   real(rk), parameter :: tol_pH = 1.0e-6_rk        ! 0.000001 pH units
   real(rk), parameter :: tol_pCO2_atm = 1.0e-8_rk  ! 0.01 uatm

   ! Variables for file I/O
   character(len=256) :: csv_file
   character(len=1024) :: line
   integer :: unit_num, ios
   integer :: n_cases, n_pass, n_fail
   logical :: file_exists

   ! Test case variables
   real(rk) :: T, S, Pbar, DIC_molkg, TA_molkg
   real(rk) :: expected_pH, expected_pCO2_atm
   real(rk) :: calc_pH, calc_pCO2_atm
   real(rk) :: calc_H2CO3, calc_HCO3, calc_CO3, calc_K0
   real(rk) :: err_pH, err_pCO2
   logical :: success
   integer :: phscale

   ! Get CSV file path from command line or use default
   if (command_argument_count() >= 1) then
      call get_command_argument(1, csv_file)
   else
      csv_file = 'reference_cases.csv'
   end if

   ! Check if file exists
   inquire(file=trim(csv_file), exist=file_exists)
   if (.not. file_exists) then
      write(error_unit, '(A)') 'ERROR: Reference file not found: ' // trim(csv_file)
      write(error_unit, '(A)') 'Generate it using: python generate_reference.py'
      stop 1
   end if

   ! Open CSV file
   open(newunit=unit_num, file=trim(csv_file), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      write(error_unit, '(A)') 'ERROR: Cannot open file: ' // trim(csv_file)
      stop 1
   end if

   ! Skip header line
   read(unit_num, '(A)', iostat=ios) line
   if (ios /= 0) then
      write(error_unit, '(A)') 'ERROR: Empty file or read error'
      close(unit_num)
      stop 1
   end if

   ! Initialize counters
   n_cases = 0
   n_pass = 0
   n_fail = 0
   phscale = 1  ! Total scale

   write(output_unit, '(A)') '=================================================='
   write(output_unit, '(A)') 'Carbonate Engine Test Driver'
   write(output_unit, '(A)') '=================================================='
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') 'Reading reference cases from: ' // trim(csv_file)
   write(output_unit, '(A)') ''

   ! Read and process each test case
   do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit  ! End of file or error

      ! Skip empty lines
      if (len_trim(line) == 0) cycle

      ! Parse CSV line: T,S,DIC_molkg,TA_molkg,pH,pCO2_atm
      read(line, *, iostat=ios) T, S, DIC_molkg, TA_molkg, expected_pH, expected_pCO2_atm
      if (ios /= 0) then
         write(error_unit, '(A)') 'WARNING: Failed to parse line: ' // trim(line)
         cycle
      end if

      n_cases = n_cases + 1

      ! Call carbonate_engine_solve
      call carbonate_engine_solve_wrapper(T, S, 0.0_rk, DIC_molkg, TA_molkg, phscale, &
                                          calc_pH, calc_pCO2_atm, calc_H2CO3, calc_HCO3, &
                                          calc_CO3, calc_K0, success)

      ! Calculate errors
      err_pH = abs(calc_pH - expected_pH)
      err_pCO2 = abs(calc_pCO2_atm - expected_pCO2_atm)

      ! Check if within tolerance
      if (.not. success) then
         n_fail = n_fail + 1
         write(output_unit, '(A,I4,A)') 'Case ', n_cases, ': FAIL - Solver did not converge'
         write(output_unit, '(A,F8.2,A,F8.3,A,ES12.5,A,ES12.5)') &
            '  Inputs: T=', T, ' S=', S, ' DIC=', DIC_molkg, ' TA=', TA_molkg
      else if (err_pH > tol_pH .or. err_pCO2 > tol_pCO2_atm) then
         n_fail = n_fail + 1
         write(output_unit, '(A,I4,A)') 'Case ', n_cases, ': FAIL - Exceeds tolerance'
         write(output_unit, '(A,F8.2,A,F8.3,A,ES12.5,A,ES12.5)') &
            '  Inputs: T=', T, ' S=', S, ' DIC=', DIC_molkg, ' TA=', TA_molkg
         write(output_unit, '(A,F10.6,A,F10.6,A,ES10.3)') &
            '  pH:   expected=', expected_pH, ' calculated=', calc_pH, ' error=', err_pH
         write(output_unit, '(A,ES12.5,A,ES12.5,A,ES10.3)') &
            '  pCO2: expected=', expected_pCO2_atm, ' calculated=', calc_pCO2_atm, ' error=', err_pCO2
      else
         n_pass = n_pass + 1
         write(output_unit, '(A,I4,A,F10.6,A,ES12.5,A)') &
            'Case ', n_cases, ': PASS (pH=', calc_pH, ', pCO2=', calc_pCO2_atm, ')'
      end if
   end do

   close(unit_num)

   ! Print summary
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '=================================================='
   write(output_unit, '(A)') 'Summary'
   write(output_unit, '(A)') '=================================================='
   write(output_unit, '(A,I4)') 'Total cases: ', n_cases
   write(output_unit, '(A,I4)') 'Passed:      ', n_pass
   write(output_unit, '(A,I4)') 'Failed:      ', n_fail
   write(output_unit, '(A)') ''

   if (n_fail > 0) then
      write(output_unit, '(A)') 'TEST FAILED'
      stop 1
   else if (n_cases == 0) then
      write(output_unit, '(A)') 'WARNING: No test cases found'
      stop 1
   else
      write(output_unit, '(A)') 'ALL TESTS PASSED'
      stop 0
   end if

contains

   !-----------------------------------------------------------------------
   ! carbonate_engine_solve_wrapper
   !
   ! Wrapper that implements the same interface as carbonate_engine_solve
   ! but is self-contained for standalone testing
   !-----------------------------------------------------------------------
   subroutine carbonate_engine_solve_wrapper(T, S, Pbar, DIC_molkg, TA_molkg, phscale, &
                                              pH, pCO2_atm, H2CO3, HCO3, CO3, K0, success)
      real(rk), intent(in)  :: T, S, Pbar
      real(rk), intent(in)  :: DIC_molkg, TA_molkg
      integer,  intent(in)  :: phscale
      real(rk), intent(out) :: pH, pCO2_atm, H2CO3, HCO3, CO3, K0
      logical,  intent(out) :: success

      ! Local constants
      real(rk), parameter :: Rgas = 83.14472_rk
      real(rk), parameter :: tol_ph_solver = 1.0e-8_rk
      real(rk), parameter :: tol_f_solver = 1.0e-12_rk
      integer,  parameter :: max_iter = 100

      ! Equilibrium constants
      real(rk) :: K1, K2, KB, KW, KS, KF
      real(rk) :: total2sws, sws2total

      ! Concentrations
      real(rk) :: BT, ST, FT

      ! Solver variables
      real(rk) :: H_solution, pH_total
      real(rk) :: Tmax

      ! Ensure temperature is non-negative
      Tmax = max(T, 0.0_rk)

      ! Calculate equilibrium constants
      call calc_constants(Tmax, S, Pbar, phscale, K0, K1, K2, KB, KW, KS, KF, total2sws)

      ! Calculate total concentrations
      call calc_totals(S, BT, ST, FT)

      ! Solve for [H+] using Brent's method
      call solve_brent(DIC_molkg, TA_molkg, K1, K2, KB, KW, BT, H_solution, success)

      if (.not. success) then
         pH = 0.0_rk
         pCO2_atm = 0.0_rk
         H2CO3 = 0.0_rk
         HCO3 = 0.0_rk
         CO3 = 0.0_rk
         return
      end if

      ! Calculate pH
      pH_total = -log10(H_solution)
      if (phscale == 1) then
         pH = pH_total
      else
         pH = -log10(H_solution * total2sws)
      end if

      ! Calculate species
      call calc_species(DIC_molkg, H_solution, K0, K1, K2, H2CO3, HCO3, CO3, pCO2_atm)

   end subroutine carbonate_engine_solve_wrapper

   subroutine calc_constants(T, S, Pbar, phscale, K0, K1, K2, KB, KW, KS, KF, total2sws)
      real(rk), intent(in)  :: T, S, Pbar
      integer,  intent(in)  :: phscale
      real(rk), intent(out) :: K0, K1, K2, KB, KW, KS, KF, total2sws

      real(rk), parameter :: Rgas = 83.14472_rk
      real(rk) :: TK, TK100, lnTK, invTK
      real(rk) :: sqrtS, S15, S2, IS, sqrtIS
      real(rk) :: Cl, ST, FT
      real(rk) :: delta, kappa
      real(rk) :: total2free, free2sws
      real(rk) :: pK1, pK2

      TK = T + 273.15_rk
      TK100 = TK / 100.0_rk
      lnTK = log(TK)
      invTK = 1.0_rk / TK

      ! Handle S=0 (freshwater) case
      if (S > 0.0_rk) then
         sqrtS = sqrt(S)
         S15 = S ** 1.5_rk
         IS = 19.924_rk * S / (1000.0_rk - 1.005_rk * S)
         sqrtIS = sqrt(IS)
      else
         sqrtS = 0.0_rk
         S15 = 0.0_rk
         IS = 0.0_rk
         sqrtIS = 0.0_rk
      end if
      S2 = S * S
      Cl = S / 1.80655_rk
      ST = 0.14_rk * Cl / 96.062_rk
      FT = 0.000067_rk * Cl / 18.998_rk

      ! KS - handle S=0 to avoid log(1) numerical issues
      if (S > 0.0_rk) then
         KS = exp(-4276.1_rk * invTK + 141.328_rk - 23.093_rk * lnTK &
                  + (-13856.0_rk * invTK + 324.57_rk - 47.986_rk * lnTK) * sqrtIS &
                  + (35474.0_rk * invTK - 771.54_rk + 114.723_rk * lnTK) * IS &
                  - 2698.0_rk * invTK * IS**1.5_rk + 1776.0_rk * invTK * IS**2.0_rk &
                  + log(1.0_rk - 0.001005_rk * S))
      else
         KS = 1.0_rk  ! Arbitrary; ST=0 so KS not used
      end if

      ! KF
      KF = exp(874.0_rk * invTK - 9.68_rk + 0.111_rk * sqrtS)

      total2free = 1.0_rk / (1.0_rk + ST / KS)
      free2sws = 1.0_rk + ST / KS + FT / KF
      total2sws = total2free * free2sws

      ! Pressure corrections for KS
      delta = -18.03_rk + 0.0466_rk * T + 0.000316_rk * T**2.0_rk
      kappa = -4.53_rk + 0.00009_rk * T
      KS = KS * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))

      ! Pressure corrections for KF
      delta = -9.78_rk - 0.009_rk * T - 0.000942_rk * T**2.0_rk
      kappa = -3.91_rk + 0.000054_rk * T
      KF = KF * total2free * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))
      total2free = 1.0_rk / (1.0_rk + ST / KS)
      KF = KF / total2free

      free2sws = 1.0_rk + ST / KS + FT / (KF * total2free)
      total2sws = total2free * free2sws

      ! K0 - Weiss (1974)
      K0 = exp(93.4517_rk / TK100 - 60.2409_rk + 23.3585_rk * log(TK100) &
               + S * (0.023517_rk - 0.023656_rk * TK100 + 0.0047036_rk * TK100**2.0_rk))
      ! Pressure correction only at depth (Pbar > 0)
      if (Pbar > 0.0_rk) then
         K0 = K0 * exp((-Pbar * 32.3_rk) / (Rgas * TK))
      end if

      ! KB
      KB = exp((-8966.9_rk - 2890.53_rk * sqrtS - 77.942_rk * S &
                + 1.728_rk * S15 - 0.0996_rk * S2) / TK &
               + (148.0248_rk + 137.1942_rk * sqrtS + 1.62142_rk * S) &
               + (-24.4344_rk - 25.085_rk * sqrtS - 0.2474_rk * S) * lnTK &
               + 0.053105_rk * sqrtS * TK)

      ! K1, K2 (Total scale - Lueker et al. 2000)
      ! Mehrbach (1973) refit by Lueker et al. (2000)
      pK1 = 3633.86_rk * invTK - 61.2172_rk + 9.6777_rk * lnTK &
            - 0.011555_rk * S + 0.0001152_rk * S2
      K1 = 10.0_rk ** (-pK1)

      pK2 = 471.78_rk * invTK + 25.929_rk - 3.16967_rk * lnTK &
            - 0.01781_rk * S + 0.0001122_rk * S2
      K2 = 10.0_rk ** (-pK2)

      ! Pressure corrections
      delta = -25.5_rk + 0.1271_rk * T
      kappa = (-3.08_rk + 0.0877_rk * T) / 1000.0_rk
      K1 = K1 * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))

      delta = -15.82_rk - 0.0219_rk * T
      kappa = (1.13_rk - 0.1475_rk * T) / 1000.0_rk
      K2 = K2 * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))

      delta = -29.48_rk + 0.1622_rk * T - 0.002608_rk * T**2.0_rk
      kappa = -2.84_rk / 1000.0_rk
      KB = KB * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))

      ! KW - Millero (1995) on SWS scale
      KW = exp(148.9802_rk - 13847.26_rk * invTK - 23.6521_rk * lnTK &
               + (-5.977_rk + 118.67_rk * invTK + 1.0495_rk * lnTK) * sqrtS &
               - 0.01615_rk * S)

      ! Convert KW from SWS to Total scale (phscale=1 assumed in test driver)
      ! KW_total = KW_sws * (1 + ST/KS) / (1 + ST/KS + FT/KF)
      ! At S=0, ST=FT=0 so conversion factor = 1
      if (S > 0.0_rk) then
         KW = KW * (1.0_rk + ST/KS) / (1.0_rk + ST/KS + FT/KF)
      end if

      delta = -25.60_rk + 0.2324_rk * T - 0.0036246_rk * T**2.0_rk
      kappa = (-5.13_rk + 0.0794_rk * T) / 1000.0_rk
      KW = KW * exp((-delta + 0.5_rk * kappa * Pbar) * Pbar / (Rgas * TK))

   end subroutine calc_constants

   subroutine calc_totals(S, BT, ST, FT)
      real(rk), intent(in)  :: S
      real(rk), intent(out) :: BT, ST, FT
      real(rk) :: Cl

      Cl = S / 1.80655_rk
      ! Total boron - Lee et al. (2010) equation for open ocean
      ! BT = 432.6 * S / 35 umol/kg
      BT = 0.0004326_rk * S / 35.0_rk
      ST = 0.14_rk * Cl / 96.062_rk
      FT = 0.000067_rk * Cl / 18.998_rk
   end subroutine calc_totals

   function alk_residual(H, DIC, TA_input, K1, K2, KB, KW, BT) result(F)
      real(rk), intent(in) :: H, DIC, TA_input, K1, K2, KB, KW, BT
      real(rk) :: F
      real(rk) :: denom, HCO3, CO3, BOH4, OH

      denom = H * H + K1 * H + K1 * K2
      HCO3 = DIC * K1 * H / denom
      CO3 = DIC * K1 * K2 / denom
      BOH4 = BT * KB / (KB + H)
      OH = KW / H
      F = HCO3 + 2.0_rk * CO3 + BOH4 + OH - H - TA_input
   end function alk_residual

   subroutine solve_brent(DIC, TA_input, K1, K2, KB, KW, BT, H_solution, success)
      real(rk), intent(in)  :: DIC, TA_input, K1, K2, KB, KW, BT
      real(rk), intent(out) :: H_solution
      logical,  intent(out) :: success

      real(rk), parameter :: tol_ph_s = 1.0e-8_rk
      real(rk), parameter :: tol_f_s = 1.0e-12_rk
      integer, parameter :: max_it = 100

      real(rk) :: a, b, c, d, e
      real(rk) :: fa, fb, fc
      real(rk) :: p, q, r, s
      real(rk) :: tol1, xm
      integer :: iter

      a = 10.0_rk ** (-12.0_rk)
      b = 10.0_rk ** (-2.0_rk)
      fa = alk_residual(a, DIC, TA_input, K1, K2, KB, KW, BT)
      fb = alk_residual(b, DIC, TA_input, K1, K2, KB, KW, BT)

      if (fa * fb > 0.0_rk) then
         a = 10.0_rk ** (-14.0_rk)
         b = 10.0_rk ** (-1.0_rk)
         fa = alk_residual(a, DIC, TA_input, K1, K2, KB, KW, BT)
         fb = alk_residual(b, DIC, TA_input, K1, K2, KB, KW, BT)
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

      do iter = 1, max_it
         ! Ensure b is the best approximation (closest to root)
         if (abs(fc) < abs(fb)) then
            a = b; b = c; c = a
            fa = fb; fb = fc; fc = fa
         end if

         tol1 = 2.0_rk * epsilon(1.0_rk) * abs(b) + 0.5_rk * tol_ph_s * b
         xm = 0.5_rk * (c - b)

         if (abs(xm) <= tol1 .or. abs(fb) < tol_f_s) then
            H_solution = b
            success = .true.
            return
         end if

         if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
            s = fb / fa
            if (abs(a - c) < epsilon(1.0_rk) * max(abs(a), abs(c))) then
               p = 2.0_rk * xm * s
               q = 1.0_rk - s
            else
               q = fa / fc
               r = fb / fc
               p = s * (2.0_rk * xm * q * (q - r) - (b - a) * (r - 1.0_rk))
               q = (q - 1.0_rk) * (r - 1.0_rk) * (s - 1.0_rk)
            end if

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
            d = xm
            e = d
         end if

         a = b
         fa = fb

         if (abs(d) > tol1) then
            b = b + d
         else
            if (xm >= 0.0_rk) then
               b = b + tol1
            else
               b = b - tol1
            end if
         end if

         fb = alk_residual(b, DIC, TA_input, K1, K2, KB, KW, BT)

         if ((fb > 0.0_rk .and. fc > 0.0_rk) .or. (fb < 0.0_rk .and. fc < 0.0_rk)) then
            c = a
            fc = fa
            e = b - a
            d = e
         end if
      end do

      success = .false.
      H_solution = 0.0_rk
   end subroutine solve_brent

   subroutine calc_species(DIC, H, K0, K1, K2, H2CO3, HCO3, CO3, pCO2_atm)
      real(rk), intent(in)  :: DIC, H, K0, K1, K2
      real(rk), intent(out) :: H2CO3, HCO3, CO3, pCO2_atm
      real(rk) :: denom

      denom = H * H + K1 * H + K1 * K2
      H2CO3 = DIC * H * H / denom
      HCO3 = DIC * K1 * H / denom
      CO3 = DIC * K1 * K2 / denom

      if (K0 > 0.0_rk) then
         pCO2_atm = H2CO3 / K0
      else
         pCO2_atm = 0.0_rk
      end if
   end subroutine calc_species

end program test_driver
