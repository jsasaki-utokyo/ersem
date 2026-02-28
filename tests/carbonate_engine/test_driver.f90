!-----------------------------------------------------------------------
! test_driver.f90
!
! Verification harness for carbonate_engine module.
! Calls the actual carbonate_engine_solve subroutine.
!
! Tests:
!   1. opt_pH_scale=1 (Total): PyCO2SYS reference CSV
!   2. Cross-scale consistency: Total vs SWS vs Free vs NBS
!   3. legacy_mode=.true.: backward-compat (convergence + range)
!   4. Sub-zero temperature (polar waters)
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
   real(rk), parameter :: tol_pH_tight    = 1.0e-4_rk   ! pH match
   real(rk), parameter :: tol_pCO2_tight  = 1.0e-6_rk   ! pCO2 in atm

   ! Variables
   character(len=256) :: csv_file
   integer :: n_pass, n_fail, n_total_pass, n_total_fail
   logical :: file_exists

   ! Get CSV file path
   if (command_argument_count() >= 1) then
      call get_command_argument(1, csv_file)
   else
      csv_file = 'reference_cases.csv'
   end if

   inquire(file=trim(csv_file), exist=file_exists)
   if (.not. file_exists) then
      write(error_unit, '(A)') 'ERROR: Reference file not found: ' &
                                // trim(csv_file)
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
   ! Test 1: opt_pH_scale=1 (Total), Lueker 2000, Lee 2010 boron
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') '--- Total scale (opt_pH_scale=1), Lueker 2000 ---'
   call run_csv_tests(csv_file, 1, 10, 2, .false., &
                      tol_pH_tight, tol_pCO2_tight, n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! ---------------------------------------------------------------
   ! Test 2: Cross-scale consistency (Total vs SWS vs Free vs NBS)
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '--- Cross-scale consistency (4 scales) ---'
   call run_cross_scale_tests(n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! ---------------------------------------------------------------
   ! Test 3: legacy_mode=.true. (backward compatible)
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '--- legacy_mode=.true. (backward compat) ---'
   call run_legacy_tests(n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! ---------------------------------------------------------------
   ! Test 4: Sub-zero temperature
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '--- Sub-zero temperature test ---'
   call run_subzero_tests(n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! Summary
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
   ! Run CSV-based tests against PyCO2SYS reference values
   !-----------------------------------------------------------------------
   subroutine run_csv_tests(csv_file, opt_pH_scale, opt_k_carbonic, &
                            opt_total_borate, legacy_mode, &
                            tol_pH_val, tol_pCO2_val, n_pass, n_fail)
      character(len=*), intent(in) :: csv_file
      integer, intent(in)  :: opt_pH_scale, opt_k_carbonic, opt_total_borate
      logical, intent(in)  :: legacy_mode
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

      open(newunit=unit_num, file=trim(csv_file), status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit, '(A)') 'ERROR: Cannot open: ' // trim(csv_file)
         n_pass = 0; n_fail = 1
         return
      end if

      n_cases = 0; n_pass = 0; n_fail = 0

      do
         read(unit_num, '(A)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle

         read(line, *, iostat=ios) T, S, DIC_molkg, TA_molkg, &
                                   expected_pH, expected_pCO2_atm
         if (ios /= 0) cycle

         n_cases = n_cases + 1

         call carbonate_engine_solve(T, S, 0.0_rk, DIC_molkg, TA_molkg, &
                                     opt_pH_scale, opt_k_carbonic, &
                                     opt_total_borate, legacy_mode, &
                                     calc_pH, calc_pCO2_atm, &
                                     calc_H2CO3, calc_HCO3, &
                                     calc_CO3, calc_K0, success)

         err_pH = abs(calc_pH - expected_pH)
         err_pCO2 = abs(calc_pCO2_atm - expected_pCO2_atm)

         if (.not. success) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') &
               'Case ', n_cases, ': FAIL - No convergence'
         else if (err_pH > tol_pH_val .or. err_pCO2 > tol_pCO2_val) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') 'Case ', n_cases, ': FAIL'
            write(output_unit, '(A,F8.2,A,F8.3,A,ES12.5,A,ES12.5)') &
               '  T=', T, ' S=', S, ' DIC=', DIC_molkg, ' TA=', TA_molkg
            write(output_unit, '(A,F12.8,A,F12.8,A,ES10.3)') &
               '  pH: exp=', expected_pH, ' calc=', calc_pH, &
               ' err=', err_pH
            write(output_unit, '(A,ES12.5,A,ES12.5,A,ES10.3)') &
               '  pCO2: exp=', expected_pCO2_atm, &
               ' calc=', calc_pCO2_atm, ' err=', err_pCO2
         else
            n_pass = n_pass + 1
            write(output_unit, '(A,I4,A,F12.8,A,ES12.5,A)') &
               'Case ', n_cases, ': PASS (pH=', calc_pH, &
               ', pCO2=', calc_pCO2_atm, ')'
         end if
      end do

      close(unit_num)
      write(output_unit, '(A,I4,A,I4,A,I4)') &
         'CSV tests: ', n_cases, ' total, ', n_pass, ' pass, ', &
         n_fail, ' fail'

   end subroutine run_csv_tests

   !-----------------------------------------------------------------------
   ! Cross-scale consistency tests:
   ! Same input -> 4 pH scales should satisfy:
   !   pH_NBS > pH_free > pH_total > pH_sws
   !   pCO2 identical across all scales
   !-----------------------------------------------------------------------
   subroutine run_cross_scale_tests(n_pass, n_fail)
      integer, intent(out) :: n_pass, n_fail

      integer, parameter :: n_tests = 5
      real(rk) :: T_vals(n_tests), S_vals(n_tests)
      real(rk) :: DIC_vals(n_tests), TA_vals(n_tests)
      real(rk) :: pH_t, pH_s, pH_f, pH_n
      real(rk) :: pCO2_t, pCO2_s, pCO2_f, pCO2_n
      real(rk) :: H2CO3, HCO3, CO3, K0
      logical :: suc_t, suc_s, suc_f, suc_n
      integer :: i
      logical :: ok

      T_vals   = (/ 25.0_rk, 15.0_rk, 10.0_rk, 20.0_rk,  5.0_rk /)
      S_vals   = (/ 35.0_rk, 35.0_rk, 30.0_rk, 25.0_rk, 28.0_rk /)
      DIC_vals = (/ 0.002100_rk, 0.002100_rk, 0.002050_rk, &
                    0.001800_rk, 0.002000_rk /)
      TA_vals  = (/ 0.002300_rk, 0.002300_rk, 0.002250_rk, &
                    0.002000_rk, 0.002200_rk /)

      n_pass = 0; n_fail = 0

      do i = 1, n_tests
         ! Total (1)
         call carbonate_engine_solve(T_vals(i), S_vals(i), 0.0_rk, &
              DIC_vals(i), TA_vals(i), 1, 10, 2, .false., &
              pH_t, pCO2_t, H2CO3, HCO3, CO3, K0, suc_t)
         ! SWS (2)
         call carbonate_engine_solve(T_vals(i), S_vals(i), 0.0_rk, &
              DIC_vals(i), TA_vals(i), 2, 10, 2, .false., &
              pH_s, pCO2_s, H2CO3, HCO3, CO3, K0, suc_s)
         ! Free (3)
         call carbonate_engine_solve(T_vals(i), S_vals(i), 0.0_rk, &
              DIC_vals(i), TA_vals(i), 3, 10, 2, .false., &
              pH_f, pCO2_f, H2CO3, HCO3, CO3, K0, suc_f)
         ! NBS (4)
         call carbonate_engine_solve(T_vals(i), S_vals(i), 0.0_rk, &
              DIC_vals(i), TA_vals(i), 4, 10, 2, .false., &
              pH_n, pCO2_n, H2CO3, HCO3, CO3, K0, suc_n)

         if (.not.(suc_t .and. suc_s .and. suc_f .and. suc_n)) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') &
               'XScale Case ', i, ': FAIL - No convergence'
            cycle
         end if

         ok = .true.
         ! pH ordering: NBS > Free > Total > SWS
         ! H_total = H_free*(1+ST/KS), H_sws = H_free*(1+ST/KS+FT/KF)
         ! H_nbs = fH*H_free where fH < 1, so H_nbs < H_free
         if (.not. (pH_f > pH_t .and. pH_t > pH_s)) then
            ok = .false.
            write(output_unit, '(A,I4,A)') &
               'XScale Case ', i, ': FAIL - pH ordering violated'
            write(output_unit, '(A,F10.6,A,F10.6,A,F10.6)') &
               '  pH_F=', pH_f, ' pH_T=', pH_t, ' pH_SWS=', pH_s
         end if

         ! NBS pH should be highest (fH < 1 -> H_nbs < H_free)
         if (.not. (pH_n > pH_f)) then
            ok = .false.
            write(output_unit, '(A,I4,A)') &
               'XScale Case ', i, ': FAIL - NBS should be > Free'
            write(output_unit, '(A,F10.6,A,F10.6)') &
               '  pH_NBS=', pH_n, ' pH_F=', pH_f
         end if

         ! pCO2 should be nearly identical across all scales
         if (abs(pCO2_t - pCO2_s) > 0.5e-6_rk .or. &
             abs(pCO2_t - pCO2_f) > 0.5e-6_rk .or. &
             abs(pCO2_t - pCO2_n) > 0.5e-6_rk) then
            ok = .false.
            write(output_unit, '(A,I4,A)') &
               'XScale Case ', i, ': FAIL - pCO2 mismatch'
            write(output_unit, '(A,ES12.5,A,ES12.5)') &
               '  pCO2_T=', pCO2_t, ' pCO2_SWS=', pCO2_s
            write(output_unit, '(A,ES12.5,A,ES12.5)') &
               '  pCO2_F=', pCO2_f, ' pCO2_NBS=', pCO2_n
         end if

         ! Free-Total difference should be reasonable (< 0.15)
         ! Total-SWS difference should be small (< 0.02)
         if (pH_f - pH_t > 0.15_rk .or. pH_t - pH_s > 0.02_rk) then
            ok = .false.
            write(output_unit, '(A,I4,A,F8.6,A,F8.6)') &
               'XScale Case ', i, &
               ': FAIL - scale diffs: F-T=', pH_f - pH_t, &
               ' T-SWS=', pH_t - pH_s
         end if

         if (ok) then
            n_pass = n_pass + 1
            write(output_unit, '(A,I4,A)') 'XScale Case ', i, ': PASS'
            write(output_unit, '(A,F10.6,A,F10.6,A,F10.6,A,F10.6)') &
               '  pH: T=', pH_t, ' SWS=', pH_s, &
               ' Free=', pH_f, ' NBS=', pH_n
         else
            n_fail = n_fail + 1
         end if
      end do

      write(output_unit, '(A,I4,A,I4,A,I4)') &
         'Cross-scale tests: ', n_tests, ' total, ', &
         n_pass, ' pass, ', n_fail, ' fail'

   end subroutine run_cross_scale_tests

   !-----------------------------------------------------------------------
   ! Legacy mode tests (legacy_mode=.true.)
   ! Verify convergence and reasonable range
   !-----------------------------------------------------------------------
   subroutine run_legacy_tests(n_pass, n_fail)
      integer, intent(out) :: n_pass, n_fail

      integer, parameter :: n_tests = 5
      real(rk) :: T_vals(n_tests), S_vals(n_tests)
      real(rk) :: DIC_vals(n_tests), TA_vals(n_tests)
      real(rk) :: pH_val, pCO2_val
      real(rk) :: H2CO3, HCO3, CO3, K0
      logical :: success
      integer :: i

      T_vals   = (/ 25.0_rk, 15.0_rk, 10.0_rk, 20.0_rk, 30.0_rk /)
      S_vals   = (/ 35.0_rk, 35.0_rk, 30.0_rk, 32.0_rk, 38.0_rk /)
      DIC_vals = (/ 0.002100_rk, 0.002100_rk, 0.002050_rk, &
                    0.002100_rk, 0.002100_rk /)
      TA_vals  = (/ 0.002300_rk, 0.002300_rk, 0.002250_rk, &
                    0.002350_rk, 0.002300_rk /)

      n_pass = 0; n_fail = 0

      do i = 1, n_tests
         ! legacy_mode=.true., opt_pH_scale=2 (SWS), opt_k_carbonic
         ! is ignored (Mehrbach/DM87 forced)
         call carbonate_engine_solve(T_vals(i), S_vals(i), 0.0_rk, &
              DIC_vals(i), TA_vals(i), 2, 10, 1, .true., &
              pH_val, pCO2_val, H2CO3, HCO3, CO3, K0, success)

         if (.not. success) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') &
               'Legacy Case ', i, ': FAIL - No convergence'
            cycle
         end if

         ! Sanity: pH in [6, 9], pCO2 in [50, 5000] uatm
         if (pH_val > 6.0_rk .and. pH_val < 9.0_rk .and. &
             pCO2_val > 50.0e-6_rk .and. pCO2_val < 5000.0e-6_rk) then
            n_pass = n_pass + 1
            write(output_unit, '(A,I4,A,F10.6,A,ES12.5,A)') &
               'Legacy Case ', i, ': PASS (pH=', pH_val, &
               ', pCO2=', pCO2_val, ')'
         else
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') &
               'Legacy Case ', i, ': FAIL - Out of range'
            write(output_unit, '(A,F10.6,A,ES12.5)') &
               '  pH=', pH_val, ' pCO2=', pCO2_val
         end if
      end do

      write(output_unit, '(A,I4,A,I4,A,I4)') &
         'Legacy tests: ', n_tests, ' total, ', &
         n_pass, ' pass, ', n_fail, ' fail'

   end subroutine run_legacy_tests

   !-----------------------------------------------------------------------
   ! Sub-zero temperature tests
   !-----------------------------------------------------------------------
   subroutine run_subzero_tests(n_pass, n_fail)
      integer, intent(out) :: n_pass, n_fail

      real(rk) :: pH_val, pCO2_val, pH_0, pCO2_0
      real(rk) :: H2CO3, HCO3, CO3, K0
      logical :: success, success_0

      n_pass = 0; n_fail = 0

      ! Test 1: -1.5C on Total scale (should work with T >= -2)
      call carbonate_engine_solve(-1.5_rk, 34.0_rk, 0.0_rk, &
           0.002200_rk, 0.002350_rk, 1, 10, 2, .false., &
           pH_val, pCO2_val, H2CO3, HCO3, CO3, K0, success)

      call carbonate_engine_solve(0.0_rk, 34.0_rk, 0.0_rk, &
           0.002200_rk, 0.002350_rk, 1, 10, 2, .false., &
           pH_0, pCO2_0, H2CO3, HCO3, CO3, K0, success_0)

      if (.not. success) then
         n_fail = n_fail + 1
         write(output_unit, '(A)') &
            'SubZero Case 1 (-1.5C): FAIL - No convergence'
      else if (pH_val > 6.0_rk .and. pH_val < 9.5_rk .and. &
               pCO2_val > 10.0e-6_rk .and. pCO2_val < 3000.0e-6_rk) then
         n_pass = n_pass + 1
         write(output_unit, '(A,F10.6,A,ES12.5,A)') &
            'SubZero Case 1 (-1.5C): PASS (pH=', pH_val, &
            ', pCO2=', pCO2_val, ')'
      else
         n_fail = n_fail + 1
         write(output_unit, '(A)') &
            'SubZero Case 1 (-1.5C): FAIL - Out of range'
      end if

      ! Test 2: legacy_mode clamps T to 0C
      call carbonate_engine_solve(-1.5_rk, 34.0_rk, 0.0_rk, &
           0.002200_rk, 0.002350_rk, 2, 10, 1, .true., &
           pH_val, pCO2_val, H2CO3, HCO3, CO3, K0, success)

      call carbonate_engine_solve(0.0_rk, 34.0_rk, 0.0_rk, &
           0.002200_rk, 0.002350_rk, 2, 10, 1, .true., &
           pH_0, pCO2_0, H2CO3, HCO3, CO3, K0, success_0)

      if (.not. success .or. .not. success_0) then
         n_fail = n_fail + 1
         write(output_unit, '(A)') &
            'SubZero Case 2 (legacy clamp): FAIL - No convergence'
      else if (abs(pH_val - pH_0) < 1.0e-10_rk .and. &
               abs(pCO2_val - pCO2_0) < 1.0e-14_rk) then
         n_pass = n_pass + 1
         write(output_unit, '(A)') &
            'SubZero Case 2 (legacy clamp): PASS (T=-1.5 clamped to 0)'
      else
         n_fail = n_fail + 1
         write(output_unit, '(A)') &
            'SubZero Case 2 (legacy clamp): FAIL'
         write(output_unit, '(A,F12.8,A,F12.8)') &
            '  pH(-1.5)=', pH_val, ' pH(0)=', pH_0
      end if

      write(output_unit, '(A,I4,A,I4,A)') &
         'Sub-zero tests: ', n_pass + n_fail, ' total, ', &
         n_pass, ' pass'

   end subroutine run_subzero_tests

end program test_driver
