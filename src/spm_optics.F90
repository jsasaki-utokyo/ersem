#include "fabm_driver.h"
!-----------------------------------------------------------------------
! Simple suspended particulate matter (SPM) model for light attenuation.
!
! State variables:
!   c   (mg/m^3, interior)  — SPM concentration in water column
!   bed (g/m^2,  bottom)    — sediment bed mass per unit area
!
! Processes:
!   Settling:     vertical_movement = -w_s (FABM handles advection)
!   Deposition:   bottom flux  = w_s * c_bottom (water → bed)
!   Erosion:      Partheniades = M * max(0, tau/tau_cr - 1) (bed → water)
!
! Provides:
!   absorption_of_silt (1/m) = k_abs * c  (consumed by light_iop)
!
! River/OBC input:
!   FVCOM reads SPM_c from river/OBC NetCDF files automatically.
!   Variables not present default to 0 (no dilution).
!
! References:
!   Partheniades (1965) erosion formula
!   Lee et al. (2005) Eq. 11 for IOP-based Kd (in light_iop.F90)
!
! jsasaki 2026-03-09
!-----------------------------------------------------------------------
module ersem_spm_optics

   use fabm_types
   use fabm_standard_variables

   implicit none

   private

   type, extends(type_base_model), public :: type_ersem_spm_optics
      ! State variable IDs
      type(type_state_variable_id)        :: id_c
      type(type_bottom_state_variable_id) :: id_bed

      ! Dependency IDs
      type(type_horizontal_dependency_id) :: id_taub

      ! Diagnostic IDs
      type(type_diagnostic_variable_id)   :: id_abs_silt

      ! Parameters
      real(rk) :: w_s      ! settling velocity (m/d)
      real(rk) :: M_ero    ! erosion rate constant (g/m^2/s)
      real(rk) :: tau_cr   ! critical shear stress for erosion (Pa)
      real(rk) :: k_abs    ! specific absorption coefficient (m^2/mg)
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
   end type type_ersem_spm_optics

contains

   subroutine initialize(self, configunit)
      class(type_ersem_spm_optics), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      ! Parameters
      call self%get_parameter(self%w_s,    'w_s',    'm/d',     'settling velocity',               default=1.0_rk)
      call self%get_parameter(self%M_ero,  'M_ero',  'g/m^2/s', 'erosion rate constant',           default=5.0e-4_rk)
      call self%get_parameter(self%tau_cr, 'tau_cr', 'Pa',      'critical shear stress for erosion',default=0.1_rk)
      call self%get_parameter(self%k_abs,  'k_abs',  'm^2/mg',  'specific absorption coefficient', default=1.5e-5_rk)

      ! Interior state variable: SPM concentration
      ! vertical_movement < 0 means sinking (FABM convention)
      call self%register_state_variable(self%id_c, 'c', 'mg/m^3', &
         'suspended particulate matter concentration', &
         minimum=0.0_rk, &
         vertical_movement=-self%w_s / 86400.0_rk)

      ! Bottom state variable: bed sediment mass
      call self%register_bottom_state_variable(self%id_bed, 'bed', 'g/m^2', &
         'bed sediment mass per unit area', &
         minimum=0.0_rk)

      ! Environmental dependency: bottom stress (Pa)
      ! Already provided by FVCOM via TAUBM_N in PREPARE_ENVIRONMENT
      call self%register_horizontal_dependency(self%id_taub, standard_variables%bottom_stress)

      ! Diagnostic: absorption of silt (1/m)
      ! This satisfies the dependency in light_iop.F90 (line 65)
      ! replacing the previous bulk_constant silt_absorption
      call self%register_diagnostic_variable(self%id_abs_silt, 'abs_silt', '1/m', &
         'absorption coefficient of suspended silt', &
         standard_variable=type_bulk_standard_variable(name='absorption_of_silt'))
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_ersem_spm_optics), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: c

      _LOOP_BEGIN_
         _GET_(self%id_c, c)
         ! Dynamic light absorption: replaces the former constant silt_absorption
         _SET_DIAGNOSTIC_(self%id_abs_silt, self%k_abs * c)
      _LOOP_END_
   end subroutine do

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class(type_ersem_spm_optics), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: c_bot, bed, taub
      real(rk) :: dep, ero

      _BOTTOM_LOOP_BEGIN_
         _GET_(self%id_c, c_bot)
         _GET_HORIZONTAL_(self%id_bed, bed)
         _GET_HORIZONTAL_(self%id_taub, taub)

         ! Deposition: settling flux at water-bed interface
         ! vertical_movement handles sinking WITHIN the water column;
         ! this handles the EXIT flux at the bottom boundary.
         dep = (self%w_s / 86400.0_rk) * c_bot   ! [mg/m^2/s]

         ! Erosion: Partheniades formula E = M * max(0, tau/tau_cr - 1)
         ero = self%M_ero * max(0.0_rk, taub / self%tau_cr - 1.0_rk) &
               * 1.0e3_rk                        ! [g/m^2/s -> mg/m^2/s]
         ! Limit erosion by available bed mass
         ero = min(ero, bed * 1.0e3_rk / 86400.0_rk)

         ! Water column bottom flux: positive = into water
         _ADD_BOTTOM_FLUX_(self%id_c, ero - dep)

         ! Bed source: positive = into bed
         _ADD_BOTTOM_SOURCE_(self%id_bed, (dep - ero) * 1.0e-3_rk)  ! [mg -> g]
      _BOTTOM_LOOP_END_
   end subroutine do_bottom

end module ersem_spm_optics
