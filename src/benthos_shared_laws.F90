! Shared benthic process laws of the ERSEM-MUSE unification, family C (sediment OM decay, temperature laws,
! interface resistance, bioturbation / bioirrigation enhancement).  jsasaki 2026-10-07 (unification plan,
! muse/docs/UNIFY_C_SPEC_20261007.md).
!
! THIS FILE EXISTS IN TWO REPOSITORIES WITH IDENTICAL CONTENT:
!   muse/src/benthos_shared_laws.F90   and   ersem/src/benthos_shared_laws.F90
! The checksums are compared by muse/tests/unify_C/check_sync.sh; edit both or neither.  Every function is pure,
! elemental and has no state, so that both models and the equivalence test evaluate literally the same text.
! Units: temperatures degC, rates d-1, thicknesses m, diffusivities m2 d-1, resistance d m-1, activities in the
! units of the caller (the ratio Y/(Y+h) only needs Y and h in the same unit).
module benthos_shared_laws
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   integer, parameter :: rk = real64

   public :: th_q10, eT_mod, tur_enh, irr_enh, r_dbl_from, dbl_from_r, r_dbl_species, om_decay_rate, k_from_affinity

   public :: o2_supply_extent   ! jsasaki 2026-10-07: family B wave 2 (amendment 3): realised respiration extent under an O2 supply constraint

contains

   ! Pure Q10 temperature factor, th = Q10**((T-Tref)/10); exactly 1 at T = Tref.  Used for every OM decay class,
   ! acceptor path and chemical oxidation (per-reaction Q10).  Tref = 20 degC in both models.
   pure elemental function th_q10(T, q10, Tref) result(th)
      real(rk), intent(in) :: T, q10, Tref
      real(rk) :: th
      if (T == Tref) then
         th = 1._rk
      else
         th = q10**((T - Tref) / 10._rk)
      end if
   end function

   ! Modified Q10 of biological uptake and growth (fauna, bacteria; ERSEM benthic_fauna / benthic_bacteria):
   ! eT = max(0, Q10**((T-Tref)/10) - Q10**((T-Tcut)/wcut)), Tcut 32 degC, wcut 3 degC.  Not the same as th_q10.
   pure elemental function eT_mod(T, q10, Tref, Tcut, wcut) result(eT)
      real(rk), intent(in) :: T, q10, Tref, Tcut, wcut
      real(rk) :: eT
      eT = max(0._rk, q10**((T - Tref) / 10._rk) - q10**((T - Tcut) / wcut))
   end function

   ! Bioturbation enhancement of the particulate diffusivity: D_b = Etur * tur_enh(Y, mtur, htur).
   pure elemental function tur_enh(Y, mtur, htur) result(f)
      real(rk), intent(in) :: Y, mtur, htur
      real(rk) :: f
      f = 1.0_rk + mtur * Y / (Y + htur)
   end function

   ! Bioirrigation enhancement: irr_min + mirr * Y/(Y+hirr).  ERSEM multiplies the layer diffusivities EDZ_i by it
   ! (irr_min = 2 in the 711-14a configuration); MUSE multiplies the surface exchange rate irr_surf by it with
   ! irr_min = 1 (the floor is already contained in irr_surf).  mirr and hirr are the same quantities.
   pure elemental function irr_enh(Y, irr_min, mirr, hirr) result(f)
      real(rk), intent(in) :: Y, irr_min, mirr, hirr
      real(rk) :: f
      f = irr_min + mirr * Y / (Y + hirr)
   end function

   ! jsasaki 2026-10-07: family B wave 2 (amendment 3), the explicit O2 SUPPLY constraint of fauna respiration.
   ! D = demanded O2 consumption, S = O2 that can be supplied to the same place in the same time (same unit, per day).
   ! The REALISED extent R = D S/(D + S) (harmonic combination; R -> D for S >> D, R -> S for S << D, R <= min(D, S)) is the
   ! single number the caller applies to the carbon loss, the DIC production and the O2 consumption (one extent, no negative-O2
   ! clipping), and it is smooth in S so that Newton iterations stay well conditioned as S -> 0.  D <= 0 gives 0; S <= 0 gives 0.
   pure elemental function o2_supply_extent(D, S) result(R)
      real(rk), intent(in) :: D, S
      real(rk) :: R
      if (D > 0._rk .and. S > 0._rk) then
         R = D * S / (D + S)
      else
         R = 0._rk
      end if
   end function

   ! Interface (diffusive boundary layer) resistance R = thickness / diffusivity, d m-1.
   ! ERSEM EDZ_mix = R_dbl; MUSE dbl/d0w = R_dbl.
   pure elemental function r_dbl_from(dbl, d0w) result(r)
      real(rk), intent(in) :: dbl, d0w
      real(rk) :: r
      r = dbl / d0w
   end function

   ! The thickness that gives the resistance R_dbl for the diffusivity d0w (MUSE keeps d0w for the pore diffusion).
   pure elemental function dbl_from_r(r_dbl, d0w) result(dbl)
      real(rk), intent(in) :: r_dbl, d0w
      real(rk) :: dbl
      dbl = r_dbl * d0w
   end function

   ! Species-dependent interface resistance (plan 2a item 6): the diffusive boundary layer of solute i is thinner for a
   ! smaller diffusivity, delta_i = delta_O2 (D_i/D_O2)^(1/3) (Schmidt-number scaling), so R_i = delta_i / D_i
   ! = delta_O2 (D_i/D_O2)^(1/3) / D_i.  Equals r_dbl_from(delta_O2, D_O2) for i = O2.  Not yet used by either model
   ! (both keep one R for all solutes in wave 1); tested and specified for the per-solute step.
   pure elemental function r_dbl_species(delta_o2, d_o2, d_i) result(r)
      real(rk), intent(in) :: delta_o2, d_o2, d_i
      real(rk) :: r
      r = delta_o2 * (d_i / d_o2)**(1._rk / 3._rk) / d_i
   end function

   ! Unified OM decay kernel of one class: R = (k * OM) * f_ox * th * acc  (OM in any amount unit; R in the same unit
   ! per day; th from th_q10).  f_ox is the acceptor limitation (0..1), acc the acceptor-specific multiplier.  The
   ! multiplication order is that of MUSE sed_network_c rates (om * f * th * acc_fac), so the two agree bit for bit.
   pure elemental function om_decay_rate(k, OM, f_ox, th, acc) result(R)
      real(rk), intent(in) :: k, OM, f_ox, th, acc
      real(rk) :: R
      R = k * max(OM, 0._rk)
      R = R * f_ox * th * acc
   end function

   ! First-order constant at the reference temperature that a bacteria-mediated law  (su + suf eN) eT eOx H  (ERSEM
   ! benthic_bacteria, per unit substrate) reduces to for a quasi-steady biomass H_ref (mg C m-2) and unit eT, eOx:
   ! k = (su + suf eN) * H_ref.  su is per mg bacterial C per day.
   pure elemental function k_from_affinity(su, suf, eN, H_ref) result(k)
      real(rk), intent(in) :: su, suf, eN, H_ref
      real(rk) :: k
      k = (su + suf * eN) * H_ref
   end function

end module benthos_shared_laws
