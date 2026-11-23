!-----------------------------------------------------------------------
! Shared sulfur laws of the water (ERSEM sulfur_cycle, isw_law = 1) and the bed (MUSE sed_network_c rates C2, C6):
! mass action in both reactants with a per-reaction Q10 (TB-FVCOM box-ersem/docs/sulfur_unification.md sections 3
! and 4.2, 2026-11-23). Pure procedures without FABM, so a test program can compare them with MUSE (gate G2).
!
!   r_HS  = theta_hs  k_hs_tot H2S O2              total aerobic sulfide oxidation (f_s0 to S0, the rest to SO4)
!   r_S0  = theta_s0  k_s0ox   S0  O2              S0 oxidation
!   r_dis = theta_dis k_disp   S0  max(0, 1 - H2S/h2s_stop)   S0 disproportionation (per S0)
!
! Integration over a host step h: ns symmetric (Strang) sub-passes of the exact extents of the isolated reactions,
! disproportionation tau/2, S0 oxidation tau/2, sulfide oxidation tau, S0 oxidation tau/2, disproportionation tau/2,
! each on the stocks left by the previous stage; the summed extents form every stock change with the stoichiometry,
! so S, O2 - (2 H2S + 1.5 S0) and TA + 2 SO4 are conserved exactly and no reactant is overdrawn.
!-----------------------------------------------------------------------
module ersem_sulfur_laws

   use, intrinsic :: iso_fortran_env, only: real64
   use, intrinsic :: iso_c_binding, only: c_double

   implicit none
   private

   integer, parameter, public :: srk = real64

   ! expm1 of the C library (libm; no Fortran intrinsic): accurate 1 - exp(-z) for small z
   interface
      pure real(c_double) function c_expm1(x) bind(C, name='expm1')
         import :: c_double
         real(c_double), value :: x
      end function
   end interface

   public :: theta_q10, extent_bimolecular, sulfur_extents, sulfur_nu

contains

   ! Q10^((T - 20)/10), exactly 1 at 20 degC (as MUSE theta)
   pure elemental real(srk) function theta_q10(q10, temp)
      real(srk), intent(in) :: q10, temp
      if (temp == 20._srk) then
         theta_q10 = 1._srk
      else
         theta_q10 = q10**((temp - 20._srk) / 10._srk)
      end if
   end function

   ! O2 per S of total aerobic sulfide oxidation: 0.5 to S0, 2 to SO4
   pure elemental real(srk) function sulfur_nu(f_s0)
      real(srk), intent(in) :: f_s0
      sulfur_nu = 2._srk * (1._srk - f_s0) + 0.5_srk * f_s0
   end function

   ! extent x over tau of dx/dt = k (A - nu x)(B - x), x(0) = 0 (closed form; 0 <= x <= min(A/nu, B))
   pure elemental real(srk) function extent_bimolecular(k, a_in, nu, b_in, tau) result(x)
      real(srk), intent(in) :: k, a_in, nu, b_in, tau
      real(srk) :: a, m, mm, z, g
      x = 0._srk
      if (k <= 0._srk .or. tau <= 0._srk .or. a_in <= 0._srk .or. b_in <= 0._srk) return
      a = a_in / nu
      m = min(a, b_in)
      mm = max(a, b_in)
      z = k * nu * (mm - m) * tau
      if (z >= 1.e-6_srk) then
         g = -c_expm1(-z) / (mm - m)
      else
         g = k * nu * tau * (1._srk - z / 2._srk + z * z / 6._srk)
      end if
      x = min(m, m * mm * g / (1._srk + m * g))
   end function

   ! summed extents (mmol S m-3) of the three reactions over h (d) by ns Strang sub-passes
   pure subroutine sulfur_extents(h2s0, s00, o20, k_hs, f_s0, k_s0, k_dis, h2s_stop, h, ns, e_hs, e_s0, e_dis)
      real(srk), intent(in)  :: h2s0, s00, o20   ! start stocks (mmol m-3), >= 0
      real(srk), intent(in)  :: k_hs, k_s0       ! theta x k (x water factor), (mmol m-3)-1 d-1
      real(srk), intent(in)  :: f_s0             ! share of sulfide oxidation to S0
      real(srk), intent(in)  :: k_dis, h2s_stop  ! theta x k_disp (d-1); sulfide stop (mmol m-3)
      real(srk), intent(in)  :: h                ! host step (d), > 0
      integer,   intent(in)  :: ns               ! number of sub-passes
      real(srk), intent(out) :: e_hs, e_s0, e_dis
      real(srk) :: hs, s0, o2, tau, nu, x
      integer :: i
      hs = max(h2s0, 0._srk); s0 = max(s00, 0._srk); o2 = max(o20, 0._srk)
      e_hs = 0._srk; e_s0 = 0._srk; e_dis = 0._srk
      tau = h / ns
      nu = sulfur_nu(f_s0)
      do i = 1, ns
         call stage_disp(k_dis, h2s_stop, 0.5_srk * tau, hs, s0, e_dis)
         call stage_s0ox(k_s0, 0.5_srk * tau, s0, o2, e_s0)
         x = extent_bimolecular(k_hs, o2, nu, hs, tau)
         hs = max(hs - x, 0._srk); s0 = s0 + f_s0 * x; o2 = max(o2 - nu * x, 0._srk)
         e_hs = e_hs + x
         call stage_s0ox(k_s0, 0.5_srk * tau, s0, o2, e_s0)
         call stage_disp(k_dis, h2s_stop, 0.5_srk * tau, hs, s0, e_dis)
      end do
   end subroutine

   ! one S0-oxidation stage over t on the running stocks
   pure subroutine stage_s0ox(k_s0, t, s0, o2, e_s0)
      real(srk), intent(in)    :: k_s0, t
      real(srk), intent(inout) :: s0, o2, e_s0
      real(srk) :: y
      y = extent_bimolecular(k_s0, o2, 1.5_srk, s0, t)
      s0 = max(s0 - y, 0._srk); o2 = max(o2 - 1.5_srk * y, 0._srk)
      e_s0 = e_s0 + y
   end subroutine

   ! one disproportionation stage over t on the running stocks (evolving sulfide stop)
   pure subroutine stage_disp(k_dis, h2s_stop, t, hs, s0, e_dis)
      real(srk), intent(in)    :: k_dis, h2s_stop, t
      real(srk), intent(inout) :: hs, s0, e_dis
      real(srk) :: y
      y = 0._srk
      if (hs < h2s_stop .and. h2s_stop > 0._srk) &
         y = extent_bimolecular(k_dis / h2s_stop, h2s_stop - hs, 0.75_srk, s0, t)
      s0 = max(s0 - y, 0._srk); hs = hs + 0.75_srk * y
      e_dis = e_dis + y
   end subroutine

end module ersem_sulfur_laws
