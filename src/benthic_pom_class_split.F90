#include "fabm_driver.h"

! Splitter of benthic POM sources over the OM classes (ERSEM-MUSE unification, family C wave 2; jsasaki 2026-10-07).
! muse/docs/UNIFY_C2_SPEC_20261007.md section 1.4.
!
! An organism that writes faeces, dead matter or detritus to ONE POM sink (coupling Q, Q6 or Q6c/n/p/s) can be coupled to this model
! instead of a benthic_column_particulate_matter_layer. Like that layer, it registers one "state-like" diagnostic per element (c, n, p, s),
! so that its sources are collected by FABM (variables c_sms_tot ...). A child processor then hands every source on, element by element,
! to up to three class layers (target1: fast class Q6f, target2: slow class Q6s, target3: refractory class Q7), with the shares
!    qx{c,n,p,s}{i}  (i = 2..ntarget; target 1 takes 1 - sum, the rule and names of pelagic_base deposition).
! The class layers are ordinary benthic_column_particulate_matter_layer instances (one per class, each with its own Q and depth bounds), so
! the change of the class stock AND of the class penetration depth is computed by the existing layer machinery. Mass is conserved
! exactly: the shares of every element sum to 1.
!
! Nothing else in ERSEM is edited; a configuration without this model is bit-identical to the version before it existed.

module ersem_benthic_pom_class_split

   use fabm_types
   use fabm_particle

   use ersem_shared

   implicit none

   private

   integer, parameter :: max_target = 3

   ! Submodel that converts the sources collected by the parent into changes of the class layers.
   type, extends(type_base_model) :: type_pom_split_processor
      type (type_horizontal_dependency_id) :: id_sms(4)                       ! collected source of c, n, p, s (per second)
      type (type_bottom_state_variable_id) :: id_target(max_target,4)         ! target(class, element)
      logical  :: has(4) = .false.                                            ! element present in the composition
      logical  :: tgt(max_target,4) = .false.                                 ! target registered for this element
      real(rk) :: share(max_target,4) = 0.0_rk
      integer  :: ntarget = 2
   contains
      procedure :: do_bottom => processor_do_bottom
   end type

   type, extends(type_particle_model), public :: type_ersem_benthic_pom_class_split
      type (type_horizontal_diagnostic_variable_id) :: id_local(4)
      type (type_bottom_state_variable_id) :: id_target(max_target,4)
      logical :: has(4) = .false.
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_ersem_benthic_pom_class_split), intent(inout), target :: self
      integer,                                    intent(in)            :: configunit

      character(len=10) :: composition, composition3
      integer           :: ntarget, iel, itarget
      character(len=1), parameter :: elname(4) = ['c','n','p','s']
      character(len=8), parameter :: elunits(4) = ['mg C   ','mmol N ','mmol P ','mmol Si']
      character(len=12),parameter :: ellong(4) = ['carbon      ','nitrogen    ','phosphorus  ','silicate    ']
      character(len=16) :: num
      real(rk)          :: qx
      class (type_pom_split_processor), pointer :: proc

      call self%get_parameter(composition,'composition','','elemental composition (subset of cnps)',default='cnps')
      call self%get_parameter(ntarget,'ntarget','','number of class targets (2: fast, slow; 3: plus refractory)',default=2,minimum=2,maximum=max_target)
      call self%get_parameter(composition3,'composition3','','elements accepted by the third target (refractory class)',default='cnp')

      allocate(proc)
      proc%dt = 86400._rk
      proc%ntarget = ntarget
      call self%add_child(proc,'processor',configunit=-1)

      do iel=1,4
         self%has(iel) = index(composition,elname(iel))/=0
         proc%has(iel) = self%has(iel)
         if (.not.self%has(iel)) cycle

         ! Collector: a diagnostic that acts as a state variable (as in benthic_column_particulate_matter_layer).
         call self%register_diagnostic_variable(self%id_local(iel),elname(iel),trim(elunits(iel))//'/m^2','source collector for '//trim(ellong(iel)), &
            act_as_state_variable=.true.,domain=domain_bottom,output=output_none,source=source_do_bottom)
         select case (elname(iel))
         case ('c')
            call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_local(iel),1.0_rk/CMass)
         case ('n')
            call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_local(iel))
         case ('p')
            call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_local(iel))
         case ('s')
            call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_local(iel))
         end select

         call proc%register_dependency(proc%id_sms(iel),elname(iel)//'_sms',trim(elunits(iel))//'/m^2/s','sources collected for '//trim(ellong(iel)))
         call proc%request_coupling(proc%id_sms(iel),'../'//elname(iel)//'_sms_tot')

         ! Shares: target 1 takes the remainder (rule of pelagic_base).
         proc%share(1,iel) = 1.0_rk
         do itarget=2,ntarget
            write (num,'(i0)') itarget
            if (itarget==3 .and. index(composition3,elname(iel))==0) then
               qx = 0.0_rk
            else
               call self%get_parameter(qx,'qx'//elname(iel)//trim(num),'-','share of '//trim(ellong(iel))//' to target '//trim(num),default=0.0_rk,minimum=0.0_rk,maximum=1.0_rk)
            end if
            proc%share(itarget,iel) = qx
            proc%share(1,iel) = proc%share(1,iel) - qx
         end do
         if (proc%share(1,iel) < -1.0e-12_rk) call self%fatal_error('initialize','the shares of '//trim(ellong(iel))//' to targets 2.. exceed 1')
         proc%share(1,iel) = max(0.0_rk,proc%share(1,iel))

         ! Class targets (state dependencies, coupled as a whole model: target<i>: <layer instance>).
         do itarget=1,ntarget
            write (num,'(i0)') itarget
            if (itarget==3 .and. index(composition3,elname(iel))==0) cycle
            proc%tgt(itarget,iel) = .true.
            call self%register_state_dependency(self%id_target(itarget,iel),'target'//trim(num)//'_'//elname(iel),trim(elunits(iel))//'/m^2','class target '//trim(num)//' for '//trim(ellong(iel)))
            call self%request_coupling_to_model(self%id_target(itarget,iel),'target'//trim(num),elname(iel))
            call proc%register_state_dependency(proc%id_target(itarget,iel),'target'//trim(num)//'_'//elname(iel),trim(elunits(iel))//'/m^2','class target '//trim(num)//' for '//trim(ellong(iel)))
            call proc%request_coupling(proc%id_target(itarget,iel),'../target'//trim(num)//'_'//elname(iel))
         end do
      end do
   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_ersem_benthic_pom_class_split), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer :: iel

      ! The collectors carry no content of their own (they only receive sources).
      _HORIZONTAL_LOOP_BEGIN_
         do iel=1,4
            if (self%has(iel)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_local(iel),0.0_rk)
         end do
      _HORIZONTAL_LOOP_END_
   end subroutine do_bottom

   subroutine processor_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_pom_split_processor), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer  :: iel, itarget
      real(rk) :: sms

      _HORIZONTAL_LOOP_BEGIN_
         do iel=1,4
            if (.not.self%has(iel)) cycle
            _GET_HORIZONTAL_(self%id_sms(iel),sms)
            ! per second (FABM) -> per day (the time unit of ERSEM's ODE macros, dt = 86400)
            sms = sms*self%dt
            do itarget=1,self%ntarget
               if (self%tgt(itarget,iel)) _SET_BOTTOM_ODE_(self%id_target(itarget,iel),self%share(itarget,iel)*sms)
            end do
         end do
      _HORIZONTAL_LOOP_END_
   end subroutine processor_do_bottom

end module ersem_benthic_pom_class_split
