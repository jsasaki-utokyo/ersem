#include "fabm_driver.h"

! Diel vertical migration of a phytoplankton functional type by swimming (design v3, after Codex review r2,
! TB-FVCOM box-ersem docs/phytoplankton_types.md section 11.3; jsasaki 2026-11-22). Intended for a raphidophyte such as
! Heterosigma akashiwo: up toward the light in the day phase, down toward the nutrients in the night phase when the
! cells are nitrogen- or phosphorus-limited, at a bounded swimming speed.
!
!   phyto_migration_target (do_column)
!     phase: an entrained clock in local time, day phase t_up <= local time < t_down (ascent starts before sunrise,
!            descent before sunset); local time in whole seconds of the day, from the host yearday (UTC) + utc_offset.
!     day:   z_target = the shallowest cell-centre depth whose PAR <= I_opt (the surface if the surface PAR is
!            already <= I_opt; the deepest cell if PAR > I_opt everywhere); strength = 1.
!     night: column-integral quota factors iNn, iNp (ERSEM primary_producer formulas), iN = min(iNn, iNp);
!            strength s = min(1, max(0, (iN_trig - iN)/iN_width)) (replete cells, s = 0, do not move);
!            the limiting nutrient is N if iNn <= iNp, else P; its attainable uptake per unit carbon with the PFT's
!            configured linear affinities is b(z) = qun3*N3n + qun4*N4n (N) or qurp*N1p (P); z_target = the depth of the
!            maximum of b among the cells BELOW the biomass centroid z_c and within the reach of the strength-scaled tanh
!            law, z_c <= z <= z_c + L asinh(sinh 1 exp(s v_max t_left / L)) (v3.1, loop C r1 #1; the reach holds the
!            current strength and candidate fixed; both are recomputed every step); the score at z_c is that of the
!            layer containing z_c; if b(z_target) does not exceed b(z_c) the strength is 0 (no
!            benefit, no descent). One strength per column: the column's cells move as one population, so the trigger
!            uses the column-integral quotas (a per-cell strength would change quotas by differential transport).
!   phyto_migration_move (get_vertical_movement)
!     w = strength * (v_max/3600) * tanh((z - z_target)/L)  (m s-1, positive up; z depth, positive down; v_max in m h-1),
!     at night w = min(w, 0): no ascent in the night phase (v3.1, loop C r1 #2; coupling day_phase),
!     the same velocity for every constituent coupled as target1..ntarget (c, n, p, Chl). The host transports it with
!     its vertical-movement scheme (box driver: implicit upwind, zero flux through the surface and the bed): each
!     constituent is conserved per column and all are carried by the same velocity, |w| <= v_max/3600 at every
!     interface; local quotas still change where the moved population mixes with the resident one, and the upwind
!     scheme spreads the profile numerically.
!
! Not used by any existing configuration.

module ersem_phyto_migration

   use fabm_types

   implicit none

   private

   integer, parameter :: max_target = 5
   real(rk), parameter :: qnRPIcX_pm = 1.26E-02_rk   ! = ersem_shared qnRPIcX (Redfield N:C, mmol N/mg C)
   real(rk), parameter :: qpRPIcX_pm = 7.86E-04_rk   ! = ersem_shared qpRPIcX (Redfield P:C, mmol P/mg C)

   type, extends(type_base_model), public :: type_ersem_phyto_migration_target
      type (type_state_variable_id)                 :: id_c, id_n, id_p, id_N3n, id_N4n, id_N1p
      type (type_dependency_id)                     :: id_par, id_thickness, id_depth
      type (type_global_dependency_id)              :: id_yday
      type (type_bottom_diagnostic_variable_id)     :: id_ztarget, id_strength, id_iN, id_dayphase
      real(rk) :: utc_offset, t_up, t_down, I_opt, qnlc, qplc, xqcn, xqcp, iN_trig, iN_width, qun3, qun4, qurp, v_max, L
   contains
      procedure :: initialize => target_initialize
      procedure :: do_column  => target_do_column
   end type

   type, extends(type_base_model), public :: type_ersem_phyto_migration_move
      type (type_state_variable_id)        :: id_target(max_target)
      type (type_dependency_id)            :: id_depth
      type (type_bottom_dependency_id)     :: id_ztarget, id_strength, id_dayphase   ! day phase: loop C r1 #2
      integer  :: ntarget
      real(rk) :: v_max, L
   contains
      procedure :: initialize            => move_initialize
      procedure :: get_vertical_movement => move_get_vertical_movement
   end type

contains

   subroutine target_initialize(self, configunit)
      class (type_ersem_phyto_migration_target), intent(inout), target :: self
      integer,                                   intent(in)            :: configunit

      call self%get_parameter(self%utc_offset, 'utc_offset', 'h', 'local time minus UTC', default=9.0_rk)
      call self%get_parameter(self%t_up, 't_up', 'h', 'local hour at which the day (ascent) phase starts', &
                              default=4.0_rk, minimum=0.0_rk, maximum=24.0_rk)
      call self%get_parameter(self%t_down, 't_down', 'h', 'local hour at which the night (descent) phase starts', &
                              default=16.0_rk, minimum=0.0_rk, maximum=24.0_rk)
      if (self%t_down <= self%t_up) call self%fatal_error('initialize', 't_down must be later than t_up')
      call self%get_parameter(self%I_opt, 'I_opt', 'W/m^2', 'PAR sought in the day phase', default=150.0_rk, &
                              minimum=0.0_rk)
      call self%get_parameter(self%qnlc, 'qnlc', 'mmol N/mg C', 'minimum N:C quota (as in the PFT)', default=0.0042_rk)
      call self%get_parameter(self%qplc, 'qplc', 'mmol P/mg C', 'minimum P:C quota (as in the PFT)', default=0.00015_rk)
      call self%get_parameter(self%xqcn, 'xqcn', '-', 'N limitation threshold relative to Redfield (as in the PFT)', &
                              default=1.0_rk)
      call self%get_parameter(self%xqcp, 'xqcp', '-', 'P limitation threshold relative to Redfield (as in the PFT)', &
                              default=1.0_rk)
      call self%get_parameter(self%iN_trig, 'iN_trig', '-', 'nutrient limitation factor below which cells descend', &
                              default=0.8_rk, minimum=0.0_rk, maximum=1.0_rk)
      call self%get_parameter(self%iN_width, 'iN_width', '-', 'width of the descent trigger', default=0.3_rk, &
                              minimum=1.0e-6_rk)
      call self%get_parameter(self%qun3, 'qun3', 'm^3/mg C/d', 'nitrate affinity (as in the PFT)', default=0.002_rk)
      call self%get_parameter(self%qun4, 'qun4', 'm^3/mg C/d', 'ammonium affinity (as in the PFT)', default=0.002_rk)
      call self%get_parameter(self%qurp, 'qurp', 'm^3/mg C/d', 'phosphate affinity (as in the PFT)', default=0.002_rk)
      call self%get_parameter(self%v_max, 'v_max', 'm/h', 'swimming speed (reach of the night search)', &
                              default=0.5_rk, minimum=0.0_rk)
      call self%get_parameter(self%L, 'L', 'm', 'slow-down depth scale of the movement model (and of the reach, v3.1)', &
                              default=1.0_rk, minimum=0.0_rk)

      call self%register_state_dependency(self%id_c, 'c', 'mg C/m^3', 'carbon of the migrating PFT')
      call self%register_state_dependency(self%id_n, 'n', 'mmol N/m^3', 'nitrogen of the migrating PFT')
      call self%register_state_dependency(self%id_p, 'p', 'mmol P/m^3', 'phosphorus of the migrating PFT')
      call self%register_state_dependency(self%id_N3n, 'N3n', 'mmol N/m^3', 'nitrate')
      call self%register_state_dependency(self%id_N4n, 'N4n', 'mmol N/m^3', 'ammonium')
      call self%register_state_dependency(self%id_N1p, 'N1p', 'mmol P/m^3', 'phosphate')
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
      call self%register_dependency(self%id_depth, standard_variables%depth)
      call self%register_dependency(self%id_yday, standard_variables%number_of_days_since_start_of_the_year)
      call self%register_diagnostic_variable(self%id_ztarget, 'z_target', 'm', 'target depth of the migration', &
         missing_value=0.0_rk, act_as_state_variable=.true., source=source_do_column)
      call self%register_diagnostic_variable(self%id_strength, 'strength', '-', 'fraction of v_max swum', &
         missing_value=0.0_rk, act_as_state_variable=.true., source=source_do_column)
      call self%register_diagnostic_variable(self%id_iN, 'iN', '-', 'column nutrient limitation factor of the PFT', &
         missing_value=1.0_rk, source=source_do_column)
      call self%register_diagnostic_variable(self%id_dayphase, 'day_phase', '-', '1 in the day phase, 0 at night', &
         missing_value=0.0_rk, source=source_do_column)
   end subroutine target_initialize

   subroutine target_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_ersem_phyto_migration_target), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: yday, par, dz, z, c, n, p, x, b, t_left
      integer  :: isec, isec_up, isec_down
      logical  :: nlim
      real(rk) :: cint, nsum, pint, czint, zc, zreach, iNn, iNp, iN, s, zt, zdeep, b_best, b_zc, dzc, zsurf
      logical  :: day, found

      _GET_GLOBAL_(self%id_yday, yday)
      ! local time as whole seconds of the day (review r2 #15: exact switch hours)
      isec = modulo(nint(86400.0_rk*(yday - floor(yday))) + nint(3600.0_rk*self%utc_offset), 86400)
      isec_up = nint(3600.0_rk*self%t_up)
      isec_down = nint(3600.0_rk*self%t_down)
      day = isec >= isec_up .and. isec < isec_down

      ! pass 1 (surface to bed): integrals, centroid, the day target
      cint = 0.0_rk; nsum = 0.0_rk; pint = 0.0_rk; czint = 0.0_rk
      found = .false.; zt = 0.0_rk; zdeep = 0.0_rk; zsurf = -1.0_rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_thickness, dz)
         _GET_(self%id_depth, z)
         _GET_(self%id_c, c)
         _GET_(self%id_n, n)
         _GET_(self%id_p, p)
         _GET_(self%id_par, par)
         if (zsurf < 0.0_rk) zsurf = z
         cint = cint + max(0.0_rk, c)*dz
         nsum = nsum + max(0.0_rk, n)*dz
         pint = pint + max(0.0_rk, p)*dz
         czint = czint + max(0.0_rk, c)*dz*z
         zdeep = z
         if (.not. found .and. par <= self%I_opt) then
            zt = z
            found = .true.
         end if
      _DOWNWARD_LOOP_END_
      if (.not. found) zt = zdeep
      if (cint > 0.0_rk) then
         zc = czint/cint
         iNn = min(1.0_rk, max(0.0_rk, (nsum/cint - self%qnlc)/(self%xqcn*qnRPIcX_pm - self%qnlc)))
         iNp = min(1.0_rk, max(0.0_rk, (pint/cint - self%qplc)/(self%xqcp*qpRPIcX_pm - self%qplc)))
      else
         zc = zsurf
         iNn = 1.0_rk
         iNp = 1.0_rk
      end if
      iN = min(iNn, iNp)
      nlim = iNn <= iNp

      if (day) then
         s = 1.0_rk
      else
         s = min(1.0_rk, max(0.0_rk, (self%iN_trig - iN)/self%iN_width))
         ! pass 2: the best attainable-uptake depth below the centroid and within reach of the time left until t_up
         t_left = real(modulo(isec_up - isec, 86400), rk)/3600.0_rk
         ! loop C r1 #1: reach of the strength-scaled tanh law: travel time to within L of a target at distance D is
         ! (L/(s v)) ln(sinh(D/L)/sinh(1)); reachable if <= t_left, i.e. D <= L asinh(sinh(1) exp(s v t_left / L))
         zreach = zc + self%L*asinh(sinh(1.0_rk)*exp(min(50.0_rk, s*self%v_max*t_left/self%L)))
         b_best = -1.0_rk; b_zc = -1.0_rk; dzc = huge(1.0_rk); zt = zc
         _DOWNWARD_LOOP_BEGIN_
            _GET_(self%id_depth, z)
            if (nlim) then
               _GET_(self%id_N3n, x)
               b = self%qun3*max(0.0_rk, x)
               _GET_(self%id_N4n, x)
               b = b + self%qun4*max(0.0_rk, x)
            else
               _GET_(self%id_N1p, x)
               b = self%qurp*max(0.0_rk, x)
            end if
            if (z >= zc .and. z <= zreach .and. b > b_best) then
               b_best = b
               zt = z
            end if
            _GET_(self%id_thickness, dz)
            ! loop C r1 #3: the centroid score is that of the layer containing the centroid
            if (abs(z - zc) <= 0.5_rk*dz .and. dzc > 0.0_rk) then
               dzc = 0.0_rk
               b_zc = b
            else if (dzc > 0.0_rk .and. abs(z - zc) < dzc) then
               dzc = abs(z - zc)
               b_zc = b
            end if
         _DOWNWARD_LOOP_END_
         if (b_best <= b_zc) then
            s = 0.0_rk
            zt = zc
         end if
      end if

      _MOVE_TO_BOTTOM_
      _SET_BOTTOM_DIAGNOSTIC_(self%id_ztarget, zt)
      _SET_BOTTOM_DIAGNOSTIC_(self%id_strength, s)
      _SET_BOTTOM_DIAGNOSTIC_(self%id_iN, iN)
      _SET_BOTTOM_DIAGNOSTIC_(self%id_dayphase, merge(1.0_rk, 0.0_rk, day))
   end subroutine target_do_column

   subroutine move_initialize(self, configunit)
      class (type_ersem_phyto_migration_move), intent(inout), target :: self
      integer,                                 intent(in)            :: configunit
      integer :: i
      character(len=8) :: idx

      call self%get_parameter(self%ntarget, 'ntarget', '', 'number of constituents moved (target1..)', &
                              default=4, minimum=1, maximum=max_target)
      call self%get_parameter(self%v_max, 'v_max', 'm/h', 'maximum swimming speed', default=0.5_rk, minimum=0.0_rk)
      call self%get_parameter(self%L, 'L', 'm', 'depth scale over which the speed drops to zero at the target', &
                              default=1.0_rk, minimum=1.0e-3_rk)
      call self%register_dependency(self%id_depth, standard_variables%depth)
      call self%register_dependency(self%id_ztarget, 'z_target', 'm', 'target depth of the migration')
      call self%register_dependency(self%id_strength, 'strength', '-', 'fraction of v_max swum')
      call self%register_dependency(self%id_dayphase, 'day_phase', '-', '1 in the day phase, 0 at night')   ! loop C r1 #2
      do i = 1, self%ntarget
         write (idx,'(i0)') i
         call self%register_state_dependency(self%id_target(i), 'target'//trim(idx), '', 'constituent '//trim(idx))
      end do
   end subroutine move_initialize

   subroutine move_get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_ersem_phyto_migration_move), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      real(rk) :: z, zt, s, w, ph
      integer  :: i

      _LOOP_BEGIN_
         _GET_(self%id_depth, z)
         _GET_BOTTOM_(self%id_ztarget, zt)
         _GET_BOTTOM_(self%id_strength, s)
         _GET_BOTTOM_(self%id_dayphase, ph)
         w = s*self%v_max/3600.0_rk*tanh((z - zt)/self%L)
         ! loop C r1 #2: no ascent at night (the convergence law would lift the layers below the target)
         if (ph < 0.5_rk) w = min(w, 0.0_rk)
         do i = 1, self%ntarget
            _ADD_VERTICAL_VELOCITY_(self%id_target(i), w)
         end do
      _LOOP_END_
   end subroutine move_get_vertical_movement

end module ersem_phyto_migration
