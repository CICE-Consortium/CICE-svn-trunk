!  SVN:$Id$
!=======================================================================
! authors Elizabeth Hunke

      module ice_restart_age

      use ice_kinds_mod

      implicit none
      private
      public :: write_restart_age, read_restart_age

      logical (kind=log_kind), public :: & 
         restart_age      ! if .true., read age tracer restart file

!=======================================================================

      contains

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_age()

      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: ncat
      use ice_fileunits, only: nu_diag, nu_dump_age
      use ice_state, only: trcrn, nt_iage
      use ice_restart,only: write_restart_field

      ! local variables

      logical (kind=log_kind) :: diag

      diag = .true.

      !-----------------------------------------------------------------

      call write_restart_field(nu_dump_age,0,trcrn(:,:,nt_iage,:,:),'ruf8', &
                               'iage',ncat,diag)

      end subroutine write_restart_age

!=======================================================================

! Reads all values needed for an ice age restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_age()

      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: ncat
      use ice_fileunits, only: nu_diag, nu_restart_age
      use ice_state, only: trcrn, nt_iage
      use ice_restart,only: read_restart_field

      ! local variables

      logical (kind=log_kind) :: &
         diag

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) 'min/max age (s)'

      call read_restart_field(nu_restart_age,0,trcrn(:,:,nt_iage,:,:),'ruf8', &
                              'iage',ncat,diag)

      end subroutine read_restart_age

!=======================================================================

      end module ice_restart_age

!=======================================================================
