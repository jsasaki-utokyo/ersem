!-----------------------------------------------------------------------
! test_driver.f90
!
! Verification harness for carbonate_engine module.
! Calls the actual carbonate_engine_solve subroutine.
!
! Tests:
!   1. All 4 pH scales vs PyCO2SYS reference CSV (9-column format)
!   2. Cross-scale consistency: Total vs SWS vs Free vs NBS
!   3. legacy_mode=.true.: backward-compat (convergence + range)
!   4. Sub-zero temperature (polar waters)
!   5. opt_k_carbonic=14 (Millero 2010) vs PyCO2SYS
!   6. Pressure correction (Pbar > 0) vs PyCO2SYS
!   7. convert_pH_scale round-trip consistency
!   8. engine=0 phscale x opt_pH_scale regression (carbonate.F90 logic)
!
! CSV format (9 columns):
!   T, S, DIC, TA, pH_total, pH_sws, pH_free, pH_nbs, pCO2_atm
!
! Usage:
!   ./test_driver [reference_cases.csv]
!-----------------------------------------------------------------------

program test_driver
   use, intrinsic :: iso_fortran_env, only: real64, error_unit, output_unit
   use carbonate_engine, only: carbonate_engine_solve, convert_pH_scale
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
   ! Test 1: All 4 pH scales vs PyCO2SYS reference (Lueker 2000)
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') '--- 4-scale PyCO2SYS reference comparison ---'
   call run_csv_tests_4scale(csv_file, 10, 2, &
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

   ! ---------------------------------------------------------------
   ! Test 5: opt_k_carbonic=14 (Millero 2010)
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '--- Millero 2010 (opt_k_carbonic=14) ---'
   call run_millero2010_tests(n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! ---------------------------------------------------------------
   ! Test 6: Pressure correction (Pbar > 0)
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '--- Pressure correction (Pbar > 0) ---'
   call run_pressure_tests(n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! ---------------------------------------------------------------
   ! Test 7: convert_pH_scale round-trip
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '--- convert_pH_scale round-trip ---'
   call run_convert_pH_tests(n_pass, n_fail)
   n_total_pass = n_total_pass + n_pass
   n_total_fail = n_total_fail + n_fail

   ! ---------------------------------------------------------------
   ! Test 8: engine=0 pH_total/pH_selected logic (carbonate.F90)
   ! ---------------------------------------------------------------
   write(output_unit, '(A)') ''
   write(output_unit, '(A)') '--- engine=0 phscale x opt_pH_scale ---'
   call run_engine0_scale_tests(n_pass, n_fail)
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
   ! Run CSV-based tests against PyCO2SYS reference values (all 4 scales)
   ! CSV format: T, S, DIC, TA, pH_total, pH_sws, pH_free, pH_nbs, pCO2
   !-----------------------------------------------------------------------
   subroutine run_csv_tests_4scale(csv_file, opt_k_carbonic, &
                                   opt_total_borate, &
                                   tol_pH_val, tol_pCO2_val, n_pass, n_fail)
      character(len=*), intent(in) :: csv_file
      integer, intent(in)  :: opt_k_carbonic, opt_total_borate
      real(rk), intent(in) :: tol_pH_val, tol_pCO2_val
      integer, intent(out) :: n_pass, n_fail

      character(len=1024) :: line
      character(len=8) :: scale_name(4)
      integer :: unit_num, ios, n_cases, iscale
      real(rk) :: T, S, DIC_molkg, TA_molkg
      real(rk) :: exp_pH(4), exp_pCO2_atm
      real(rk) :: calc_pH, calc_pCO2_atm
      real(rk) :: calc_H2CO3, calc_HCO3, calc_CO3, calc_K0
      real(rk) :: err_pH, err_pCO2
      logical :: success, case_ok

      scale_name(1) = 'Total'
      scale_name(2) = 'SWS'
      scale_name(3) = 'Free'
      scale_name(4) = 'NBS'

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

         ! Read 9-column CSV
         read(line, *, iostat=ios) T, S, DIC_molkg, TA_molkg, &
                                   exp_pH(1), exp_pH(2), exp_pH(3), &
                                   exp_pH(4), exp_pCO2_atm
         if (ios /= 0) cycle

         n_cases = n_cases + 1
         case_ok = .true.

         ! Test each of the 4 pH scales
         do iscale = 1, 4
            call carbonate_engine_solve(T, S, 0.0_rk, &
                 DIC_molkg, TA_molkg, iscale, opt_k_carbonic, &
                 opt_total_borate, .false., &
                 calc_pH, calc_pCO2_atm, &
                 calc_H2CO3, calc_HCO3, calc_CO3, calc_K0, success)

            if (.not. success) then
               case_ok = .false.
               write(output_unit, '(A,I3,A,A,A)') &
                  'Case ', n_cases, ' ', trim(scale_name(iscale)), &
                  ': FAIL - No convergence'
               cycle
            end if

            err_pH = abs(calc_pH - exp_pH(iscale))
            err_pCO2 = abs(calc_pCO2_atm - exp_pCO2_atm)

            if (err_pH > tol_pH_val .or. err_pCO2 > tol_pCO2_val) then
               case_ok = .false.
               write(output_unit, '(A,I3,A,A,A)') &
                  'Case ', n_cases, ' ', &
                  trim(scale_name(iscale)), ': FAIL'
               write(output_unit, '(A,F8.2,A,F8.3)') &
                  '  T=', T, ' S=', S
               write(output_unit, '(A,F12.8,A,F12.8,A,ES10.3)') &
                  '  pH: exp=', exp_pH(iscale), &
                  ' calc=', calc_pH, ' err=', err_pH
               write(output_unit, '(A,ES12.5,A,ES12.5,A,ES10.3)') &
                  '  pCO2: exp=', exp_pCO2_atm, &
                  ' calc=', calc_pCO2_atm, ' err=', err_pCO2
            end if
         end do

         if (case_ok) then
            n_pass = n_pass + 1
            write(output_unit, '(A,I3,A)') &
               'Case ', n_cases, ': PASS (all 4 scales)'
         else
            n_fail = n_fail + 1
         end if
      end do

      close(unit_num)
      write(output_unit, '(A,I4,A,I4,A,I4)') &
         '4-scale CSV tests: ', n_cases, ' total, ', n_pass, &
         ' pass, ', n_fail, ' fail'

   end subroutine run_csv_tests_4scale

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

   !-----------------------------------------------------------------------
   ! Millero 2010 (opt_k_carbonic=14) tests vs PyCO2SYS
   !-----------------------------------------------------------------------
   subroutine run_millero2010_tests(n_pass, n_fail)
      integer, intent(out) :: n_pass, n_fail

      integer, parameter :: n_tests = 5
      real(rk) :: T_v(n_tests), S_v(n_tests)
      real(rk) :: DIC_v(n_tests), TA_v(n_tests)
      real(rk) :: exp_pH(n_tests), exp_fCO2(n_tests)
      real(rk) :: calc_pH, calc_pCO2, H2CO3, HCO3, CO3, K0
      real(rk) :: err_pH, err_pCO2
      logical :: success
      integer :: i

      ! PyCO2SYS reference: opt_k_carbonic=14, opt_k_fluoride=2,
      ! opt_total_borate=2, opt_pH_scale=1(Total), Pbar=0
      T_v   = (/ 25.0_rk, 15.0_rk, 20.0_rk, 10.0_rk, 25.0_rk /)
      S_v   = (/ 35.0_rk, 35.0_rk, 32.0_rk, 30.0_rk, 10.0_rk /)
      DIC_v = (/ 0.002100_rk, 0.002100_rk, 0.002100_rk, &
                 0.002050_rk, 0.002000_rk /)
      TA_v  = (/ 0.002300_rk, 0.002300_rk, 0.002350_rk, &
                 0.002250_rk, 0.002200_rk /)
      exp_pH  = (/ 7.8562326840_rk, 8.0033184745_rk, 8.0599744209_rk, &
                   8.1494226599_rk, 8.2920775215_rk /)
      exp_fCO2 = (/ 6.7311234012e-04_rk, 4.5472896597e-04_rk, &
                    4.0791169220e-04_rk, 3.1055875674e-04_rk, &
                    2.8860492098e-04_rk /)

      n_pass = 0; n_fail = 0

      do i = 1, n_tests
         call carbonate_engine_solve(T_v(i), S_v(i), 0.0_rk, &
              DIC_v(i), TA_v(i), 1, 14, 2, .false., &
              calc_pH, calc_pCO2, H2CO3, HCO3, CO3, K0, success)

         err_pH = abs(calc_pH - exp_pH(i))
         err_pCO2 = abs(calc_pCO2 - exp_fCO2(i))

         if (.not. success) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') &
               'M2010 Case ', i, ': FAIL - No convergence'
         else if (err_pH > tol_pH_tight .or. &
                  err_pCO2 > tol_pCO2_tight) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') 'M2010 Case ', i, ': FAIL'
            write(output_unit, '(A,F12.8,A,F12.8,A,ES10.3)') &
               '  pH: exp=', exp_pH(i), ' calc=', calc_pH, &
               ' err=', err_pH
            write(output_unit, '(A,ES12.5,A,ES12.5,A,ES10.3)') &
               '  fCO2: exp=', exp_fCO2(i), ' calc=', calc_pCO2, &
               ' err=', err_pCO2
         else
            n_pass = n_pass + 1
            write(output_unit, '(A,I4,A,F10.6,A,ES12.5,A)') &
               'M2010 Case ', i, ': PASS (pH=', calc_pH, &
               ', fCO2=', calc_pCO2, ')'
         end if
      end do

      write(output_unit, '(A,I4,A,I4,A,I4)') &
         'Millero 2010 tests: ', n_tests, ' total, ', &
         n_pass, ' pass, ', n_fail, ' fail'

   end subroutine run_millero2010_tests

   !-----------------------------------------------------------------------
   ! Pressure correction tests (Pbar > 0) vs PyCO2SYS
   !-----------------------------------------------------------------------
   subroutine run_pressure_tests(n_pass, n_fail)
      integer, intent(out) :: n_pass, n_fail

      integer, parameter :: n_tests = 5
      real(rk) :: T_v(n_tests), S_v(n_tests)
      real(rk) :: DIC_v(n_tests), TA_v(n_tests), Pbar_v(n_tests)
      real(rk) :: exp_pH(n_tests), exp_fCO2(n_tests)
      real(rk) :: calc_pH, calc_pCO2, H2CO3, HCO3, CO3, K0
      real(rk) :: err_pH, err_pCO2
      logical :: success
      integer :: i

      ! PyCO2SYS reference: opt_k_carbonic=10, Pbar > 0
      T_v    = (/ 25.0_rk, 15.0_rk, 10.0_rk, 5.0_rk, 25.0_rk /)
      S_v    = (/ 35.0_rk, 35.0_rk, 30.0_rk, 34.0_rk, 35.0_rk /)
      DIC_v  = (/ 0.002100_rk, 0.002100_rk, 0.002050_rk, &
                  0.002200_rk, 0.002100_rk /)
      TA_v   = (/ 0.002300_rk, 0.002300_rk, 0.002250_rk, &
                  0.002350_rk, 0.002300_rk /)
      Pbar_v = (/ 10.0_rk, 10.0_rk, 50.0_rk, 100.0_rk, 500.0_rk /)
      exp_pH   = (/ 7.8491574775_rk, 7.9978894874_rk, &
                    8.1321585842_rk, 8.0084727413_rk, &
                    7.6832256078_rk /)
      exp_fCO2 = (/ 6.7132856591e-04_rk, 4.4849834020e-04_rk, &
                    3.0286963838e-04_rk, 3.9439259266e-04_rk, &
                    6.3706467418e-04_rk /)

      n_pass = 0; n_fail = 0

      do i = 1, n_tests
         call carbonate_engine_solve(T_v(i), S_v(i), Pbar_v(i), &
              DIC_v(i), TA_v(i), 1, 10, 2, .false., &
              calc_pH, calc_pCO2, H2CO3, HCO3, CO3, K0, success)

         err_pH = abs(calc_pH - exp_pH(i))
         err_pCO2 = abs(calc_pCO2 - exp_fCO2(i))

         if (.not. success) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A)') &
               'Depth Case ', i, ': FAIL - No convergence'
         else if (err_pH > tol_pH_tight .or. &
                  err_pCO2 > tol_pCO2_tight) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I4,A,F6.0,A)') &
               'Depth Case ', i, ' (Pbar=', Pbar_v(i), '): FAIL'
            write(output_unit, '(A,F12.8,A,F12.8,A,ES10.3)') &
               '  pH: exp=', exp_pH(i), ' calc=', calc_pH, &
               ' err=', err_pH
            write(output_unit, '(A,ES12.5,A,ES12.5,A,ES10.3)') &
               '  fCO2: exp=', exp_fCO2(i), ' calc=', calc_pCO2, &
               ' err=', err_pCO2
         else
            n_pass = n_pass + 1
            write(output_unit, '(A,I4,A,F6.0,A,F10.6,A)') &
               'Depth Case ', i, ' (Pbar=', Pbar_v(i), &
               '): PASS (pH=', calc_pH, ')'
         end if
      end do

      write(output_unit, '(A,I4,A,I4,A,I4)') &
         'Pressure tests: ', n_tests, ' total, ', &
         n_pass, ' pass, ', n_fail, ' fail'

   end subroutine run_pressure_tests

   !-----------------------------------------------------------------------
   ! convert_pH_scale round-trip tests
   ! Solve on Total, convert to each scale, convert back -> should match
   !-----------------------------------------------------------------------
   subroutine run_convert_pH_tests(n_pass, n_fail)
      integer, intent(out) :: n_pass, n_fail

      real(rk) :: pH_total, pH_conv, pH_back
      real(rk) :: pCO2, H2CO3, HCO3, CO3, K0
      logical :: success
      integer :: iscale
      character(len=8) :: sname(4)

      sname(1) = 'Total'
      sname(2) = 'SWS'
      sname(3) = 'Free'
      sname(4) = 'NBS'

      n_pass = 0; n_fail = 0

      ! Get reference total pH
      call carbonate_engine_solve(25.0_rk, 35.0_rk, 0.0_rk, &
           0.002100_rk, 0.002300_rk, 1, 10, 2, .false., &
           pH_total, pCO2, H2CO3, HCO3, CO3, K0, success)

      if (.not. success) then
         n_fail = n_fail + 4
         write(output_unit, '(A)') 'Convert: FAIL - No convergence'
         return
      end if

      do iscale = 1, 4
         ! Total -> iscale -> Total
         call convert_pH_scale(25.0_rk, 35.0_rk, 0.0_rk, &
              pH_total, 1, iscale, 2, pH_conv)
         call convert_pH_scale(25.0_rk, 35.0_rk, 0.0_rk, &
              pH_conv, iscale, 1, 2, pH_back)

         if (abs(pH_back - pH_total) < 1.0e-10_rk) then
            n_pass = n_pass + 1
            write(output_unit, '(A,A,A,F10.6,A,F10.6,A)') &
               'Convert T->', trim(sname(iscale)), '->T: PASS (', &
               pH_conv, ' ->', pH_back, ')'
         else
            n_fail = n_fail + 1
            write(output_unit, '(A,A,A,ES10.3)') &
               'Convert T->', trim(sname(iscale)), &
               '->T: FAIL err=', abs(pH_back - pH_total)
         end if
      end do

      write(output_unit, '(A,I4,A,I4,A)') &
         'Convert tests: ', n_pass + n_fail, ' total, ', &
         n_pass, ' pass'

   end subroutine run_convert_pH_tests

   !-----------------------------------------------------------------------
   ! engine=0 phscale x opt_pH_scale regression test
   !
   ! Simulates carbonate.F90 do-routine logic for engine=0:
   !   pH_total is derived from CO2DYN output (scale = phscale)
   !   pH_selected is derived from pH_total via convert_pH_scale
   !
   ! Verifies that pH_total is always correct total-scale pH regardless
   ! of the phscale/opt_pH_scale combination.
   !-----------------------------------------------------------------------
   subroutine run_engine0_scale_tests(n_pass, n_fail)
      integer, intent(out) :: n_pass, n_fail

      real(rk), parameter :: T = 25.0_rk, S = 35.0_rk, Pbar = 0.0_rk
      real(rk), parameter :: DIC = 0.002100_rk, TA = 0.002300_rk
      integer, parameter :: opt_tb = 2  ! Lee 2010

      real(rk) :: pH_ref_total, pH_ref_sws
      real(rk) :: pCO2, H2CO3, HCO3, CO3, K0
      logical  :: success

      ! Variables mirroring carbonate.F90 logic
      real(rk) :: pH_co2dyn     ! simulated CO2DYN output
      real(rk) :: pH_total      ! standard variable (must be total)
      real(rk) :: pH_selected   ! selected-scale diagnostic
      real(rk) :: err

      integer :: phscale, opt_pH_scale
      character(len=60) :: label

      n_pass = 0; n_fail = 0

      ! Get reference values: true total and SWS pH
      call carbonate_engine_solve(T, S, Pbar, DIC, TA, &
           1, 10, opt_tb, .false., &
           pH_ref_total, pCO2, H2CO3, HCO3, CO3, K0, success)
      call carbonate_engine_solve(T, S, Pbar, DIC, TA, &
           2, 10, opt_tb, .false., &
           pH_ref_sws, pCO2, H2CO3, HCO3, CO3, K0, success)

      ! ---- Case 1: Consistent, phscale=1 opt_pH_scale=1 ----
      ! CO2DYN returns Total. pH_total = pH. No pH_selected.
      phscale = 1; opt_pH_scale = 1
      label = 'phscale=1, opt_pH_scale=1 (consistent)'
      pH_co2dyn = pH_ref_total  ! CO2DYN with phscale=1 returns Total

      ! carbonate.F90 logic:
      if (phscale == 1) then
         pH_total = pH_co2dyn
      else
         call convert_pH_scale(T, S, Pbar, pH_co2dyn, &
              2, 1, opt_tb, pH_total)
      end if

      err = abs(pH_total - pH_ref_total)
      if (err < 1.0e-10_rk) then
         n_pass = n_pass + 1
         write(output_unit, '(A,A,A)') 'E0 ', trim(label), ': PASS'
      else
         n_fail = n_fail + 1
         write(output_unit, '(A,A,A,ES10.3)') &
            'E0 ', trim(label), ': FAIL err=', err
      end if

      ! ---- Case 2: Consistent, phscale=0 opt_pH_scale=2 ----
      ! CO2DYN returns SWS. Convert SWS->Total for pH_total.
      phscale = 0; opt_pH_scale = 2
      label = 'phscale=0, opt_pH_scale=2 (consistent)'
      pH_co2dyn = pH_ref_sws  ! CO2DYN with phscale=0 returns SWS

      if (phscale == 1) then
         pH_total = pH_co2dyn
      else
         call convert_pH_scale(T, S, Pbar, pH_co2dyn, &
              2, 1, opt_tb, pH_total)
      end if

      err = abs(pH_total - pH_ref_total)
      if (err < 1.0e-4_rk) then
         n_pass = n_pass + 1
         write(output_unit, '(A,A,A)') 'E0 ', trim(label), ': PASS'
      else
         n_fail = n_fail + 1
         write(output_unit, '(A,A,A,ES10.3)') &
            'E0 ', trim(label), ': FAIL err=', err
      end if

      ! ---- Case 3: Inconsistent, phscale=1 opt_pH_scale=3 (Free) ----
      ! CO2DYN returns Total (phscale=1). opt_pH_scale=3 only affects
      ! pH_selected diagnostic. pH_total must still be correct Total.
      phscale = 1; opt_pH_scale = 3
      label = 'phscale=1, opt_pH_scale=3 (inconsistent)'
      pH_co2dyn = pH_ref_total  ! CO2DYN with phscale=1 returns Total

      ! carbonate.F90 logic: check phscale (not opt_pH_scale!)
      if (phscale == 1) then
         pH_total = pH_co2dyn
      else
         call convert_pH_scale(T, S, Pbar, pH_co2dyn, &
              2, 1, opt_tb, pH_total)
      end if

      ! pH_selected: convert Total -> Free
      call convert_pH_scale(T, S, Pbar, pH_total, &
           1, opt_pH_scale, opt_tb, pH_selected)

      err = abs(pH_total - pH_ref_total)
      if (err < 1.0e-10_rk) then
         n_pass = n_pass + 1
         write(output_unit, '(A,A,A)') 'E0 ', trim(label), &
            ': PASS (pH_total correct)'
      else
         n_fail = n_fail + 1
         write(output_unit, '(A,A,A,ES10.3)') &
            'E0 ', trim(label), ': FAIL pH_total err=', err
      end if

      ! Verify pH_selected is actually on Free scale
      ! Free pH > Total pH for seawater (less H+ counted on free)
      if (pH_selected > pH_total .and. &
          abs(pH_selected - pH_total) > 0.05_rk .and. &
          abs(pH_selected - pH_total) < 0.15_rk) then
         n_pass = n_pass + 1
         write(output_unit, '(A,F8.4,A,F8.4,A)') &
            'E0 pH_selected(Free)=', pH_selected, &
            ' vs pH_total=', pH_total, ': PASS'
      else
         n_fail = n_fail + 1
         write(output_unit, '(A,F8.4,A,F8.4,A)') &
            'E0 pH_selected(Free)=', pH_selected, &
            ' vs pH_total=', pH_total, ': FAIL'
      end if

      write(output_unit, '(A,I4,A,I4,A,I4)') &
         'engine=0 scale tests: ', n_pass + n_fail, ' total, ', &
         n_pass, ' pass, ', n_fail, ' fail'

   end subroutine run_engine0_scale_tests

end program test_driver
