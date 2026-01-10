!-----------------------------------------------------------------------
! carbonate_engine.F90
!
! PyCO2SYS-compatible carbonate chemistry solver for ERSEM
! Solves DIC + TA -> pH + pCO2 using robust root-finding (Brent's method)
!
! This module provides an alternative carbonate chemistry engine that
! reproduces PyCO2SYS results without runtime Python dependencies.
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
   real(rk), parameter :: pH_lo_default = 4.0_rk  ! Lower pH bound
   real(rk), parameter :: pH_hi_default = 10.0_rk ! Upper pH bound

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
   !   T         - Temperature (degC)
   !   S         - Practical salinity (PSU)
   !   Pbar      - Pressure (bar), 0 for surface
   !   DIC_molkg - Dissolved inorganic carbon (mol/kg)
   !   TA_molkg  - Total alkalinity (mol/kg)
   !   phscale   - pH scale (1: total, 0: SWS, -1: SWS backward compat)
   !
   ! Outputs:
   !   pH        - pH on requested scale
   !   pCO2_atm  - Partial pressure of CO2 (atm)
   !   H2CO3     - Carbonic acid concentration (mol/kg)
   !   HCO3      - Bicarbonate concentration (mol/kg)
   !   CO3       - Carbonate concentration (mol/kg)
   !   K0        - Henry's law constant
   !   success   - .true. if solver converged
   !-----------------------------------------------------------------------
   subroutine carbonate_engine_solve(T, S, Pbar, DIC_molkg, TA_molkg, phscale, &
                                     pH, pCO2_atm, H2CO3, HCO3, CO3, K0, success)
      real(rk), intent(in)  :: T, S, Pbar
      real(rk), intent(in)  :: DIC_molkg, TA_molkg
      integer,  intent(in)  :: phscale
      real(rk), intent(out) :: pH, pCO2_atm, H2CO3, HCO3, CO3, K0
      logical,  intent(out) :: success

      ! STUB IMPLEMENTATION - always fails, triggers fallback to legacy solver
      ! Full implementation to be added in next commit

      success = .false.
      pH = 0.0_rk
      pCO2_atm = 0.0_rk
      H2CO3 = 0.0_rk
      HCO3 = 0.0_rk
      CO3 = 0.0_rk
      K0 = 0.0_rk

   end subroutine carbonate_engine_solve

end module carbonate_engine
