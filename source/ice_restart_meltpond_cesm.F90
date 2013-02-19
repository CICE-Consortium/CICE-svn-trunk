!=======================================================================
!
!BOP
!
! !MODULE: ice_restart_meltpond_cesm - CESM meltpond parameterization
!
! !DESCRIPTION:
!
! This meltpond parameterization was developed for use with the delta-
! Eddington radiation scheme, and only affects the radiation budget in
! the model.  That is, although the pond volume is tracked, that liquid
! water is not used elsewhere in the model for mass budgets or other
! physical processes.
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors David A. Bailey (NCAR)
!         Marika M. Holland (NCAR)
!         Elizabeth C. Hunke (LANL)
!
! !INTERFACE:
!
      module ice_restart_meltpond_cesm
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_fileunits
      use ice_read_write
      use ice_restart, only: lenstr, restart_dir, restart_file, &
                             pointer_file, runtype
      use ice_communicate, only: my_task, master_task
!
!EOP
!
      implicit none

      logical (kind=log_kind) :: & 
         restart_pond_cesm ! if .true., read meltponds restart file

      real (kind=dbl_kind) :: &
         hs0               ! snow depth for transition to bare sea ice (m)

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_pond - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_pond_cesm(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for restarting
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
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
              restart_file(1:lenstr(restart_file)),'.cesmpond.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_pond,filename,0)

      if (my_task == master_task) then
        write(nu_dump_pond) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      diag = .true.

      do n = 1, ncat
         call ice_write(nu_dump_pond,0,trcrn(:,:,nt_apnd,n,:),'ruf8',diag)
         call ice_write(nu_dump_pond,0,trcrn(:,:,nt_hpnd,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_dump_pond)

      end subroutine write_restart_pond_cesm

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_pond - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_pond_cesm(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for a meltpond volume restart
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
      use ice_exit, only: abort_ice
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: &
         filename, filename0, string1, string2

      logical (kind=log_kind) :: &
         diag

      if (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         read(nu_rst_pointer,'(a)') filename0
         filename = trim(filename0)
         close(nu_rst_pointer)

         ! reconstruct path/file
         n = index(filename0,trim(restart_file))
         if (n == 0) call abort_ice('cesm pond restart: filename discrepancy')
         string1 = trim(filename0(1:n-1))
         string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
         write(filename,'(a,a,a,a)') &
            string1(1:lenstr(string1)), &
            restart_file(1:lenstr(restart_file)),'.cesmpond', &
            string2(1:lenstr(string2))
      endif ! master_task

      call ice_open(nu_restart_pond,filename,0)

      if (my_task == master_task) then
        read(nu_restart_pond) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      do n = 1, ncat
         call ice_read(nu_restart_pond,0,trcrn(:,:,nt_apnd,n,:),'ruf8',diag)
         call ice_read(nu_restart_pond,0,trcrn(:,:,nt_hpnd,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_restart_pond)

      end subroutine read_restart_pond_cesm

!=======================================================================

      end module ice_restart_meltpond_cesm

!=======================================================================
