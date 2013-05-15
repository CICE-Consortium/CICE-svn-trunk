!=======================================================================
!
!BOP
!
! !MODULE: ice_restart_firstyear - First year concentration tracer for sea ice
!
! !DESCRIPTION:
!
! see 
! Armour, K. C., C. M. Bitz, L. Thompson and E. C. Hunke (2011). Controls
! on Arctic sea ice from first-year and multi-year ice survivability.
! J. Climate, 24, 23782390. doi: 10.1175/2010JCLI3823.1.
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors C. Bitz, University of Washington, modified from ice_age module
!
! 2012: E. Hunke adopted from CESM into CICE, changed name from ice_FY.F90
!
! !INTERFACE:
!
      module ice_restart_firstyear
!
! !USES:
!
      use ice_kinds_mod
!
!EOP
!
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
!
!BOP
!
! !IROUTINE: write_restart_FY - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_FY(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for restarting
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, year_init
      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: ncat
      use ice_fileunits, only: nu_diag, nu_rst_pointer, nu_dump_FY
      use ice_flux, only: frz_onset
      use ice_read_write, only: ice_open, ice_write
      use ice_restart, only: lenstr, restart_dir, restart_file
      use ice_state, only: trcrn, nt_FY
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         iyear = nyr + year_init - 1
         imonth = month
         iday = mday
         
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.FY.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_FY,filename,0)

      if (my_task == master_task) then
        write(nu_dump_FY) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do n = 1, ncat
         call ice_write(nu_dump_FY,0,trcrn(:,:,nt_FY,n,:),'ruf8',diag)
         call ice_write(nu_dump_FY,0,frz_onset,'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_dump_FY)

      end subroutine write_restart_FY

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_FY - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_FY(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for an ice FY restart
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: ncat
      use ice_calendar, only: istep1, time, time_forc
      use ice_fileunits, only: nu_diag, nu_rst_pointer, nu_restart_FY
      use ice_flux, only: frz_onset
      use ice_read_write, only: ice_open, ice_read
      use ice_restart, only: lenstr, restart_file, pointer_file
      use ice_state, only: trcrn, nt_FY
      use ice_exit, only: abort_ice
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
         i, j, k, n, it, iblk, & ! counting indices
         iyear, imonth, iday , & ! year, month, day
         iignore                 ! dummy variable

      real (kind=real_kind) :: &
         rignore                 ! dummy variable

      character(len=char_len_long) :: &
         filename, filename0, string1, string2

      logical (kind=log_kind) :: &
         diag

      if (my_task == master_task) then
         ! reconstruct path/file
         if (present(filename_spec)) then
            filename = filename_spec
         else
            open(nu_rst_pointer,file=pointer_file)
            read(nu_rst_pointer,'(a)') filename0
            filename = trim(filename0)
            close(nu_rst_pointer)

            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice('FY restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.FY', &
               string2(1:lenstr(string2))
         endif
      endif ! master_task

      call ice_open(nu_restart_FY,filename,0)

      if (my_task == master_task) then
!        read (nu_restart_FY) istep1,time,time_forc
        read (nu_restart_FY) iignore,rignore,rignore
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do n = 1, ncat
         call ice_read(nu_restart_FY,0,trcrn(:,:,nt_FY,n,:),'ruf8',diag)
         call ice_read(nu_restart_FY,0,frz_onset,'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_restart_FY)

      end subroutine read_restart_FY

!=======================================================================

      end module ice_restart_firstyear

!=======================================================================
