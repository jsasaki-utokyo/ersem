!-----------------------------------------------------------------------
! fabm_types_mock.F90
!
! Minimal mock of fabm_types for standalone testing of carbonate_engine.
! Provides only the rk kind parameter needed by carbonate_engine.
!-----------------------------------------------------------------------
module fabm_types
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   integer, parameter :: rk = real64
end module fabm_types
