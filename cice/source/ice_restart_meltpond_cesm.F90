!  SVN:$Id$
!=======================================================================
!
! CESM meltpond parameterization
!
! This meltpond parameterization was developed for use with the delta-
! Eddington radiation scheme, and only affects the radiation budget in
! the model.  That is, although the pond volume is tracked, that liquid
! water is not used elsewhere in the model for mass budgets or other
! physical processes.
!
! authors David A. Bailey (NCAR)
!         Marika M. Holland (NCAR)
!         Elizabeth C. Hunke (LANL)

      module ice_restart_meltpond_cesm

      use ice_kinds_mod

      implicit none
      private
      public :: write_restart_pond_cesm, read_restart_pond_cesm

      logical (kind=log_kind), public :: & 
         restart_pond_cesm ! if .true., read meltponds restart file

      real (kind=dbl_kind), public :: &
         hs0               ! snow depth for transition to bare sea ice (m)

!=======================================================================

      contains

!=======================================================================
!
! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine write_restart_pond_cesm()

      use ice_domain_size, only: ncat
      use ice_fileunits, only: nu_dump_pond
      use ice_state, only: trcrn, nt_apnd, nt_hpnd
      use ice_restart, only: write_restart_field

      ! local variables

      logical (kind=log_kind) :: diag

      diag = .true.

      call write_restart_field(nu_dump_pond,0,trcrn(:,:,nt_apnd,:,:),'ruf8', &
                               'apnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0,trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                               'hpnd',ncat,diag)

      end subroutine write_restart_pond_cesm

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine read_restart_pond_cesm()

      use ice_domain_size, only: ncat
      use ice_communicate, only: my_task, master_task
      use ice_fileunits, only: nu_diag, nu_restart_pond 
      use ice_state, only: trcrn, nt_apnd, nt_hpnd
      use ice_restart, only: read_restart_field

      ! local variables

      logical (kind=log_kind) :: &
         diag

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) 'min/max cesm ponds'

      call read_restart_field(nu_restart_pond,0,trcrn(:,:,nt_apnd,:,:),'ruf8', &
                              'apnd',ncat,diag)
      call read_restart_field(nu_restart_pond,0,trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                              'hpnd',ncat,diag)

      end subroutine read_restart_pond_cesm

!=======================================================================

      end module ice_restart_meltpond_cesm

!=======================================================================
