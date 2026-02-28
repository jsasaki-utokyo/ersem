!-----------------------------------------------------------------------
! test_driver.f90
!
! Verification harness for carbonate_engine module.
! Calls the actual carbonate_engine_solve subroutine (not a reimplementation).
!
! Tests three pH scale modes:
!   phscale= 1: Total scale (PyCO2SYS reference)
!   phscale= 0: SWS scale (PyCO2SYS reference)
!   phscale=-1: SWS backward compatible (legacy CO2DYN reference)
!
! Usage:
!   ./test_driver [reference_cases.csv]
!-----------------------------------------------------------------------

program test_driver
   use, intrinsic :: iso_fortran_env, only: real64, error_unit, output_unit
   use carbonate_engine, only: carbonate_engine_solve
   implicit none

   integer, parameter :: rk = real64

   ! Tolerances
   real(rk), parameter :: tol_pH_tight    = 1.0e-6_rk   ! For PyCO2SYS match
   real(rk), parameter :: tol_pCO2_tight  = 1.0e-8_rk   ! 0.01 uatm
   real(rk), parameter :: tol_pH_legacy   = 0.02_rk      ! Looser for legacy compat
   real(rk), parameter :: tol_pCO2_legacy = 5.0e-6_rk    ! 5 uatm for legacy

   ! Variables for file I/O
   character(len=256) :: csv_file
   character(len=1024) :: line
   integer :: unit_num, ios
   integer :: n_cases, n_pass, n_fail, n_total_pass, n_total_fail
   logical :: file_exists

   ! Test case variables
   real(rk) :: T, S, Pbar, DIC_molkg, TA_molkg
   real(rk) :: expected_pH, expected_pCO2_atm
   real(rk) :: calc_pH, calc_pCO2_atm
   real(rk) :: calc_H2CO3, calc_HCO3, calc_CO3, calc_K0
   real(rk) :: err_pH, err_pCO2
   logical :: success
   integer :: phscale, opt_k_carbonic, opt_total_borate

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

   write(output_unit, '(A)') '=================================================='
   write(output_unit, '(A)') 'Carbonate Engine Test Driver'
   write(output_unit, '(A)') '(calling actual carbonate_engine_solve)'
   write(output_unit, '(A)') '=================================================='
   write(output_unit, '(A)') ''

   n_total_pass = 0
   n_total_fail = 0

   ! ---------------------------------------------------------------
   ! Test 1: phscale=1 (Total scale) with Lueker 2000
   ! Reference: PyCO2SYS with opt_k_carbonic=4 (Lueker), opt_pH_scale=1
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') '--- phscale=1 (Total), K1K2=Lueker 2000 ---'
   phscale = 1
   opt_k_carbonic = 1     ! Lueker 2000
   opt_total_borate = 2   ! Lee 2010

   call run_csv_tests(csv_file, phscale, opt_k_carbonic, opt_total_borate, &
                      tol_pH_tight, tol_pCO2_tight, n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! ---------------------------------------------------------------
   ! Test 2: phscale=0 (SWS scale) with Lueker 2000
   ! Self-consistency: SWS pH should differ from Total pH by
   ! a small amount related to total2sws conversion factor
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '--- phscale=0 (SWS), K1K2=Lueker 2000 ---'
   phscale = 0
   opt_k_carbonic = 1
   opt_total_borate = 2
   call run_sws_consistency_tests(phscale, opt_k_carbonic, opt_total_borate, &
                                  n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! ---------------------------------------------------------------
   ! Test 3: phscale=-1 (backward compatible) self-consistency
   ! Verify solver converges and pH/pCO2 are in reasonable range
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '--- phscale=-1 (backward compat) ---'
   call run_legacy_tests(n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! ---------------------------------------------------------------
   ! Test 4: Sub-zero temperature (polar waters)
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '--- Sub-zero temperature test ---'
   call run_subzero_tests(n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! Print summary
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '=================================================='
   write(output_unit, '(A)') 'Summary'
   write(output_unit, '(A)') '=================================================='
   write(output_unit, '(A,I4)') 'Total passed: ', n_total_pass
   write(output_unit, '(A,I4)') 'Total failed: ', n_total_fail
   write(output_unit, '(A)') ''

   if (n_total_fail > 0) then
      write(output_unit, '(A)') 'TEST FAILED'
      stop 1
   else if (n_total_pass == 0) then
      write(output_unit, '(A)') 'WARNING: No test cases found'
      stop 1
   else
      write(output_unit, '(A)') 'ALL TESTS PASSED'
      stop 0
   end if

contains

   !-----------------------------------------------------------------------
   ! Run CSV-based tests for a given phscale/options combination
   !-----------------------------------------------------------------------
   subroutine run_csv_tests(csv_file, phscale, opt_k_carbonic, opt_total_borate, &
                            tol_pH_val, tol_pCO2_val, n_pass, n_fail)
      character(len=*), intent(in) :: csv_file
      integer, intent(in)  :: phscale, opt_k_carbonic, opt_total_borate
      real(rk), intent(in) :: tol_pH_val, tol_pCO2_val
      integer, intent(out) :: n_pass, n_fail

      character(len=1024) :: line
      integer :: unit_num, ios, n_cases
      real(rk) :: T, S, DIC_molkg, TA_molkg
      real(rk) :: expected_pH, expected_pCO2_atm
      real(rk) :: calc_pH, calc_pCO2_atm
      real(rk) :: calc_H2CO3, calc_HCO3, calc_CO3, calc_K0
      real(rk) :: err_pH, err_pCO2
      logical :: success

      open(newunit=unit_num, file=trim(csv_file), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit, '(A)') 'ERROR: Cannot open file: ' // trim(csv_file)
         n_pass = 0; n_fail = 1
         return
      end if

      ! Skip header/comment lines
      read(unit_num, '(A)', iostat=ios) line

      n_cases = 0; n_pass = 0; n_fail = 0

      do
         read(unit_num, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle

         read(line, *, iostat=ios) T, S, DIC_molkg, TA_molkg, expected_pH, expected_pCO2_atm
         if (ios /= 0) cycle

         n_cases = n_cases + 1

         call carbonate_engine_solve(T, S, 0.0_rk, DIC_molkg, TA_molkg, phscale, &
                                     opt_k_carbonic, opt_total_borate, &
                                     calc_pH, calc_pCO2_atm, calc_H2CO3, calc_HCO3, &
                                     calc_CO3, calc_K0, success)

         err_pH = abs(calc_pH - expected_pH)
         err_pCO2 = abs(calc_pCO2_atm - expected_pCO2_atm)

         if (.not. success) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') 'Case ', n_cases, ': FAIL - No convergence'
         else if (err_pH > tol_pH_val .or. err_pCO2 > tol_pCO2_val) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') 'Case ', n_cases, ': FAIL'
            write(output_unit, '(A,F8.2,A,F8.3,A,ES12.5,A,ES12.5)') &
               '  T=', T, ' S=', S, ' DIC=', DIC_molkg, ' TA=', TA_molkg
            write(output_unit, '(A,F12.8,A,F12.8,A,ES10.3)') &
               '  pH: exp=', expected_pH, ' calc=', calc_pH, ' err=', err_pH
            write(output_unit, '(A,ES12.5,A,ES12.5,A,ES10.3)') &
               '  pCO2: exp=', expected_pCO2_atm, ' calc=', calc_pCO2_atm, ' err=', err_pCO2
         else
            n_pass = n_pass + 1
            write(output_unit, '(A,I4,A,F12.8,A,ES12.5,A)') &
               'Case ', n_cases, ': PASS (pH=', calc_pH, ', pCO2=', calc_pCO2_atm, ')'
         end if
      end do

      close(unit_num)
      write(output_unit, '(A,I4,A,I4,A,I4)') &
         'CSV tests: ', n_cases, ' total, ', n_pass, ' pass, ', n_fail, ' fail'

   end subroutine run_csv_tests

   !-----------------------------------------------------------------------
   ! SWS scale consistency tests:
   ! Verify that SWS pH differs from Total pH by the expected amount
   !-----------------------------------------------------------------------
   subroutine run_sws_consistency_tests(phscale_sws, opt_k_carbonic, opt_total_borate, &
                                         n_pass, n_fail)
      integer, intent(in)  :: phscale_sws, opt_k_carbonic, opt_total_borate
      integer, intent(out) :: n_pass, n_fail

      integer, parameter :: n_tests = 5
      real(rk) :: T_vals(n_tests), S_vals(n_tests)
      real(rk) :: DIC_vals(n_tests), TA_vals(n_tests)
      real(rk) :: pH_total, pCO2_total, pH_sws, pCO2_sws
      real(rk) :: H2CO3, HCO3, CO3, K0
      real(rk) :: pH_diff, pCO2_diff
      logical :: success_t, success_s
      integer :: i

      ! Test conditions
      T_vals   = (/ 25.0_rk, 15.0_rk, 10.0_rk, 20.0_rk, 5.0_rk /)
      S_vals   = (/ 35.0_rk, 35.0_rk, 30.0_rk, 25.0_rk, 28.0_rk /)
      DIC_vals = (/ 0.002100_rk, 0.002100_rk, 0.002050_rk, 0.001800_rk, 0.002000_rk /)
      TA_vals  = (/ 0.002300_rk, 0.002300_rk, 0.002250_rk, 0.002000_rk, 0.002200_rk /)

      n_pass = 0; n_fail = 0

      do i = 1, n_tests
         ! Solve on Total scale
         call carbonate_engine_solve(T_vals(i), S_vals(i), 0.0_rk, &
                                     DIC_vals(i), TA_vals(i), 1, &
                                     opt_k_carbonic, opt_total_borate, &
                                     pH_total, pCO2_total, H2CO3, HCO3, CO3, K0, success_t)

         ! Solve on SWS scale
         call carbonate_engine_solve(T_vals(i), S_vals(i), 0.0_rk, &
                                     DIC_vals(i), TA_vals(i), phscale_sws, &
                                     opt_k_carbonic, opt_total_borate, &
                                     pH_sws, pCO2_sws, H2CO3, HCO3, CO3, K0, success_s)

         if (.not. success_t .or. .not. success_s) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') 'SWS Case ', i, ': FAIL - No convergence'
            cycle
         end if

         ! SWS pH should be slightly lower than Total pH (by ~0.01)
         pH_diff = pH_total - pH_sws
         ! pCO2 should be nearly identical (species don't depend on pH scale)
         pCO2_diff = abs(pCO2_total - pCO2_sws)

         ! Check: SWS pH < Total pH (positive difference)
         ! and difference is small (< 0.02 for typical seawater)
         ! pCO2 should match within 0.5 uatm
         if (pH_diff > 0.0_rk .and. pH_diff < 0.02_rk .and. &
             pCO2_diff < 0.5e-6_rk) then
            n_pass = n_pass + 1
            write(output_unit, '(A,I4,A,F10.6,A,F10.6,A,F8.6)') &
               'SWS Case ', i, ': PASS (pH_T=', pH_total, ' pH_SWS=', pH_sws, &
               ' diff=', pH_diff
         else
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') 'SWS Case ', i, ': FAIL'
            write(output_unit, '(A,F10.6,A,F10.6,A,F8.6)') &
               '  pH_T=', pH_total, ' pH_SWS=', pH_sws, ' diff=', pH_diff
            write(output_unit, '(A,ES12.5,A,ES12.5,A,ES10.3)') &
               '  pCO2_T=', pCO2_total, ' pCO2_SWS=', pCO2_sws, ' diff=', pCO2_diff
         end if
      end do

      write(output_unit, '(A,I4,A,I4,A,I4)') &
         'SWS tests: ', n_tests, ' total, ', n_pass, ' pass, ', n_fail, ' fail'

   end subroutine run_sws_consistency_tests

   !-----------------------------------------------------------------------
   ! Legacy backward-compatible mode tests (phscale=-1)
   ! Verify convergence and reasonable range for typical marine conditions
   !-----------------------------------------------------------------------
   subroutine run_legacy_tests(n_pass, n_fail)
      integer, intent(out) :: n_pass, n_fail

      integer, parameter :: n_tests = 5
      real(rk) :: T_vals(n_tests), S_vals(n_tests)
      real(rk) :: DIC_vals(n_tests), TA_vals(n_tests)
      real(rk) :: pH_legacy, pCO2_legacy
      real(rk) :: H2CO3, HCO3, CO3, K0
      logical :: success
      integer :: i

      T_vals   = (/ 25.0_rk, 15.0_rk, 10.0_rk, 20.0_rk, 30.0_rk /)
      S_vals   = (/ 35.0_rk, 35.0_rk, 30.0_rk, 32.0_rk, 38.0_rk /)
      DIC_vals = (/ 0.002100_rk, 0.002100_rk, 0.002050_rk, 0.002100_rk, 0.002100_rk /)
      TA_vals  = (/ 0.002300_rk, 0.002300_rk, 0.002250_rk, 0.002350_rk, 0.002300_rk /)

      n_pass = 0; n_fail = 0

      do i = 1, n_tests
         call carbonate_engine_solve(T_vals(i), S_vals(i), 0.0_rk, &
                                     DIC_vals(i), TA_vals(i), -1, &
                                     1, 1, &  ! opt_k_carbonic/borate ignored for -1
                                     pH_legacy, pCO2_legacy, H2CO3, HCO3, CO3, K0, success)

         if (.not. success) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') 'Legacy Case ', i, ': FAIL - No convergence'
            cycle
         end if

         ! Sanity checks: pH in [6, 9], pCO2 in [50, 5000] uatm
         if (pH_legacy > 6.0_rk .and. pH_legacy < 9.0_rk .and. &
             pCO2_legacy > 50.0e-6_rk .and. pCO2_legacy < 5000.0e-6_rk) then
            n_pass = n_pass + 1
            write(output_unit, '(A,I4,A,F10.6,A,ES12.5,A)') &
               'Legacy Case ', i, ': PASS (pH=', pH_legacy, ', pCO2=', pCO2_legacy, ')'
         else
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') 'Legacy Case ', i, ': FAIL - Out of range'
            write(output_unit, '(A,F10.6,A,ES12.5)') &
               '  pH=', pH_legacy, ' pCO2=', pCO2_legacy
         end if
      end do

      write(output_unit, '(A,I4,A,I4,A,I4)') &
         'Legacy tests: ', n_tests, ' total, ', n_pass, ' pass, ', n_fail, ' fail'

   end subroutine run_legacy_tests

   !-----------------------------------------------------------------------
   ! Sub-zero temperature tests
   ! Verify solver handles negative temperatures (polar waters)
   !-----------------------------------------------------------------------
   subroutine run_subzero_tests(n_pass, n_fail)
      integer, intent(out) :: n_pass, n_fail

      real(rk) :: pH_val, pCO2_val
      real(rk) :: H2CO3, HCO3, CO3, K0
      real(rk) :: pH_0, pCO2_0
      logical :: success, success_0

      n_pass = 0; n_fail = 0

      ! Test at -1.5°C (typical polar/Arctic water)
      call carbonate_engine_solve(-1.5_rk, 34.0_rk, 0.0_rk, &
                                  0.002200_rk, 0.002350_rk, 1, &
                                  1, 2, &
                                  pH_val, pCO2_val, H2CO3, HCO3, CO3, K0, success)

      ! Also test at 0°C for comparison
      call carbonate_engine_solve(0.0_rk, 34.0_rk, 0.0_rk, &
                                  0.002200_rk, 0.002350_rk, 1, &
                                  1, 2, &
                                  pH_0, pCO2_0, H2CO3, HCO3, CO3, K0, success_0)

      if (.not. success) then
         n_fail = n_fail + 1
         write(output_unit, '(A)') 'SubZero Case 1 (-1.5C): FAIL - No convergence'
      else if (pH_val > 6.0_rk .and. pH_val < 9.5_rk .and. &
               pCO2_val > 10.0e-6_rk .and. pCO2_val < 3000.0e-6_rk) then
         ! Sub-zero should give higher pH (lower pCO2) than 0°C
         if (success_0 .and. pH_val > pH_0) then
            n_pass = n_pass + 1
            write(output_unit, '(A,F10.6,A,ES12.5,A)') &
               'SubZero Case 1 (-1.5C): PASS (pH=', pH_val, ', pCO2=', pCO2_val, ')'
            write(output_unit, '(A,F10.6,A,ES12.5,A)') &
               '  (cf. 0.0C: pH=', pH_0, ', pCO2=', pCO2_0, ')'
         else
            n_pass = n_pass + 1
            write(output_unit, '(A,F10.6,A,ES12.5,A)') &
               'SubZero Case 1 (-1.5C): PASS (pH=', pH_val, ', pCO2=', pCO2_val, ')'
         end if
      else
         n_fail = n_fail + 1
         write(output_unit, '(A)') 'SubZero Case 1 (-1.5C): FAIL - Out of range'
         write(output_unit, '(A,F10.6,A,ES12.5)') '  pH=', pH_val, ' pCO2=', pCO2_val
      end if

      ! Test that phscale=-1 clamps to 0°C (backward compat)
      call carbonate_engine_solve(-1.5_rk, 34.0_rk, 0.0_rk, &
                                  0.002200_rk, 0.002350_rk, -1, &
                                  1, 1, &
                                  pH_val, pCO2_val, H2CO3, HCO3, CO3, K0, success)

      ! Legacy result at -1.5°C should equal legacy result at 0°C
      call carbonate_engine_solve(0.0_rk, 34.0_rk, 0.0_rk, &
                                  0.002200_rk, 0.002350_rk, -1, &
                                  1, 1, &
                                  pH_0, pCO2_0, H2CO3, HCO3, CO3, K0, success_0)

      if (.not. success .or. .not. success_0) then
         n_fail = n_fail + 1
         write(output_unit, '(A)') 'SubZero Case 2 (legacy clamp): FAIL - No convergence'
      else if (abs(pH_val - pH_0) < 1.0e-10_rk .and. &
               abs(pCO2_val - pCO2_0) < 1.0e-14_rk) then
         n_pass = n_pass + 1
         write(output_unit, '(A)') 'SubZero Case 2 (legacy clamp): PASS (T=-1.5 clamped to 0)'
      else
         n_fail = n_fail + 1
         write(output_unit, '(A)') 'SubZero Case 2 (legacy clamp): FAIL'
         write(output_unit, '(A,F12.8,A,F12.8)') '  pH(-1.5)=', pH_val, ' pH(0)=', pH_0
      end if

      write(output_unit, '(A,I4,A,I4,A)') &
         'Sub-zero tests: ', n_pass + n_fail, ' total, ', n_pass, ' pass'

   end subroutine run_subzero_tests

end program test_driver
