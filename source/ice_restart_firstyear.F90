!  SVN:$Id$
!=======================================================================

! First year concentration tracer restart files
! author: E. C. Hunke, LANL

      module ice_restart_firstyear

      use ice_kinds_mod

      implicit none
      private
      public :: write_restart_FY, read_restart_FY

      logical (kind=log_kind), public :: & 
         restart_FY      ! if .true., read FY tracer restart file

!=======================================================================

      contains

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_FY()

      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: ncat
      use ice_fileunits, only: nu_diag, nu_dump_FY
      use ice_flux, only: frz_onset
      use ice_state, only: trcrn, nt_FY
      use ice_restart, only: write_restart_field

      ! local variables

      logical (kind=log_kind) :: diag

      diag = .true.

      !-----------------------------------------------------------------

      call write_restart_field(nu_dump_FY,0,trcrn(:,:,nt_FY,:,:),'ruf8', &
                               'FY',ncat,diag)
      call write_restart_field(nu_dump_FY,0,frz_onset,'ruf8', &
                               'frz_onset',1,diag)

      end subroutine write_restart_FY

!=======================================================================

! Reads all values needed for an ice FY restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_FY()

      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: ncat
      use ice_fileunits, only: nu_diag, nu_restart_FY
      use ice_flux, only: frz_onset
      use ice_state, only: trcrn, nt_FY
      use ice_restart, only: read_restart_field

      ! local variables

      logical (kind=log_kind) :: &
         diag

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) 'min/max first-year ice area'

      call read_restart_field(nu_restart_FY,0,trcrn(:,:,nt_FY,:,:),'ruf8', &
                              'FY',ncat,diag)

      if (my_task == master_task) write(nu_diag,*) 'min/max frz_onset'

      call read_restart_field(nu_restart_FY,0,frz_onset,'ruf8', &
                              'frz_onset',1,diag)

      end subroutine read_restart_FY

!=======================================================================

      end module ice_restart_firstyear

!=======================================================================
