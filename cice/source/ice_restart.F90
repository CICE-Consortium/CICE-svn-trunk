!=======================================================================
!
!BOP
!
! !MODULE: ice_restart - ice model restart files
!
! !DESCRIPTION:
!
! Read and write ice model restart files
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb LANL
!
! 2004-05: Block structure added by William Lipscomb
!          Restart module separated from history module
! 2006 ECH: Accepted some CCSM code into mainstream CICE
!           Converted to free source form (F90) 
! 
! !INTERFACE:
!
      module ice_restart
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_blocks
      use ice_read_write
      use ice_fileunits
      use ice_timers
!
!EOP
!
      implicit none
      save

      logical (kind=log_kind) :: &
         restart ! if true, initialize using restart file instead of defaults

      character (len=char_len) :: &
         restart_file  , & ! output file for restart dump
         runtype           ! initial, continue, hybrid or branch

      character (len=char_len_long) :: &
         restart_dir   , & ! directory name for restart dump
         runid             ! identifier for CCSM coupled run

      character (len=char_len_long) :: &
         pointer_file      ! input pointer file for restarts

!=======================================================================

      contains

!=======================================================================

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: dumpfile - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine dumpfile(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for a restart
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_flux
      use ice_grid
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
      use ice_dyn_evp
      use ice_work, only: work1
      use ice_ocean, only: oceanmixed_ice
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
              restart_file(1:lenstr(restart_file)),'.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! write pointer (path/file)
      if (my_task == master_task) then
        open(nu_rst_pointer,file=pointer_file)
        write(nu_rst_pointer,'(a)') filename
        close(nu_rst_pointer)
      endif

      ! begin writing restart data
      call ice_open(nu_dump,filename,0)

      if (my_task == master_task) then
        write(nu_dump) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
        write(nu_diag,*) 'Restart written ',istep1,time,time_forc
      endif

      diag = .true.

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------

      do n=1,ncat
         call ice_write(nu_dump,0,aicen(:,:,n,:),'ruf8',diag)
         call ice_write(nu_dump,0,vicen(:,:,n,:),'ruf8',diag)
         call ice_write(nu_dump,0,vsnon(:,:,n,:),'ruf8',diag)
         do it = 1, ntrcr
            call ice_write(nu_dump,0,trcrn(:,:,it,n,:),'ruf8',diag)
         enddo
      enddo

      do k=1,ntilyr
         call ice_write(nu_dump,0,eicen(:,:,k,:),'ruf8',diag)
      enddo

      do k=1,ntslyr
         call ice_write(nu_dump,0,esnon(:,:,k,:),'ruf8',diag)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,uvel,'ruf8',diag)
      call ice_write(nu_dump,0,vvel,'ruf8',diag)

      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,fresh,'ruf8',diag)
      call ice_write(nu_dump,0,fsalt,'ruf8',diag)
      call ice_write(nu_dump,0,fhocn,'ruf8',diag)

      !-----------------------------------------------------------------
      ! ocean stress (for bottom heat flux in thermo)
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,strocnxT,'ruf8',diag)
      call ice_write(nu_dump,0,strocnyT,'ruf8',diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,stressp_1,'ruf8',diag)
      call ice_write(nu_dump,0,stressm_1,'ruf8',diag)
      call ice_write(nu_dump,0,stress12_1,'ruf8',diag)

      call ice_write(nu_dump,0,stressp_2,'ruf8',diag)
      call ice_write(nu_dump,0,stressm_2,'ruf8',diag)
      call ice_write(nu_dump,0,stress12_2,'ruf8',diag)

      call ice_write(nu_dump,0,stressp_3,'ruf8',diag)
      call ice_write(nu_dump,0,stressm_3,'ruf8',diag)
      call ice_write(nu_dump,0,stress12_3,'ruf8',diag)

      call ice_write(nu_dump,0,stressp_4,'ruf8',diag)
      call ice_write(nu_dump,0,stressm_4,'ruf8',diag)
      call ice_write(nu_dump,0,stress12_4,'ruf8',diag)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work1(i,j,iblk) = c0
            if (iceumask(i,j,iblk)) work1(i,j,iblk) = c1
         enddo
         enddo
      enddo
      call ice_write(nu_dump,0,work1,'ruf8',diag)

      ! for mixed layer model
      if (oceanmixed_ice) then
         call ice_write(nu_dump,0,sst,'ruf8',diag)
         call ice_write(nu_dump,0,frzmlt,'ruf8',diag)
      endif

      if (my_task == master_task) close(nu_dump)

      end subroutine dumpfile

!=======================================================================
!BOP
!
! !IROUTINE: restartfile  - restarts from a dumpfile
!
! !INTERFACE:
!
      subroutine restartfile
!
! !DESCRIPTION:
!
! Restarts from a dump
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_broadcast
      use ice_boundary
      use ice_domain_size
      use ice_domain
      use ice_calendar, only: istep0, istep1, time, time_forc
      use ice_flux
      use ice_state
      use ice_grid, only: tmask
      use ice_itd
      use ice_ocean, only: oceanmixed_ice
      use ice_work, only: work1
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, k, n, it, iblk ! counting indices

      character(len=char_len_long) :: &
         filename, filename0

      logical (kind=log_kind) :: &
         diag, hit_eof

      if (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         read(nu_rst_pointer,'(a)') filename0
         filename = trim(filename0)
         close(nu_rst_pointer)
      endif

      call ice_open(nu_restart,filename,0)

      if (my_task == master_task) then
         read (nu_restart) istep0,time,time_forc
         write(nu_diag,*) 'read ',pointer_file(1:lenstr(pointer_file))
         write(nu_diag,*) 'restart read at istep=',istep0,time,time_forc
      endif

      call broadcast_scalar(istep0,master_task)

      istep1 = istep0

      call broadcast_scalar(time,master_task)
      call broadcast_scalar(time_forc,master_task)

      diag = .true.     ! write min/max diagnostics for field

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------
      do n=1,ncat
         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' min/max area, vol ice, vol snow, Tsfc'

         call ice_read(nu_restart,0,aicen(:,:,n,:),'ruf8',diag)
         call ice_read(nu_restart,0,vicen(:,:,n,:),'ruf8',diag)
         call ice_read(nu_restart,0,vsnon(:,:,n,:),'ruf8',diag)
         do it = 1, ntrcr
            call ice_read(nu_restart,0,trcrn(:,:,it,n,:),'ruf8',diag)
         enddo
      enddo

      if (my_task == master_task) &
           write(nu_diag,*) 'min/max eicen for each layer'
      do k=1,ntilyr
         call ice_read(nu_restart,0,eicen(:,:,k,:),'ruf8',diag)
      enddo

      if (my_task == master_task) &
           write(nu_diag,*) 'min/max esnon for each layer'
      do k=1,ntslyr
         call ice_read(nu_restart,0,esnon(:,:,k,:),'ruf8',diag)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max velocity components'

      call ice_read(nu_restart,0,uvel,'ruf8',diag)
      call ice_read(nu_restart,0,vvel,'ruf8',diag)

      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux
      !-----------------------------------------------------------------
      if (my_task == master_task) &
         write(nu_diag,*) 'min/max fresh water and heat flux components'

      call ice_read(nu_restart,0,fresh,'ruf8',diag)
      call ice_read(nu_restart,0,fsalt,'ruf8',diag)
      call ice_read(nu_restart,0,fhocn,'ruf8',diag)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max ocean stress components'

      call ice_read(nu_restart,0,strocnxT,'ruf8',diag)
      call ice_read(nu_restart,0,strocnyT,'ruf8',diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      if (my_task == master_task) write(nu_diag,*) &
           'internal stress components'
      
      call ice_read(nu_restart,0,stressp_1,'ruf8',diag)
      call ice_read(nu_restart,0,stressm_1,'ruf8',diag)
      call ice_read(nu_restart,0,stress12_1,'ruf8',diag)

      call ice_read(nu_restart,0,stressp_2,'ruf8',diag)
      call ice_read(nu_restart,0,stressm_2,'ruf8',diag)
      call ice_read(nu_restart,0,stress12_2,'ruf8',diag)

      call ice_read(nu_restart,0,stressp_3,'ruf8',diag)
      call ice_read(nu_restart,0,stressm_3,'ruf8',diag)
      call ice_read(nu_restart,0,stress12_3,'ruf8',diag)

      call ice_read(nu_restart,0,stressp_4,'ruf8',diag)
      call ice_read(nu_restart,0,stressm_4,'ruf8',diag)
      call ice_read(nu_restart,0,stress12_4,'ruf8',diag)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'ice mask for dynamics'

      call ice_read(nu_restart,0,work1,'ruf8',diag)

      iceumask(:,:,:) = .false.
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (work1(i,j,iblk) > p5) iceumask(i,j,iblk) = .true.
         enddo
         enddo
      enddo

      ! for mixed layer model
      if (oceanmixed_ice) then

         if (my_task == master_task) &
              write(nu_diag,*) 'min/max sst, frzmlt'

         call ice_read(nu_restart,0,sst,'ruf8',diag)
         call ice_read(nu_restart,0,frzmlt,'ruf8',diag)
      endif

      if (my_task == master_task) close(nu_restart)

      !-----------------------------------------------------------------
      ! update boundary conditions
      !-----------------------------------------------------------------

      call ice_timer_start(timer_bound)

      call bound_state (aicen, trcrn, &
                        vicen, vsnon, &
                        eicen, esnon)

      call update_ghost_cells (uvel,               bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (vvel,               bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stressp_1,               bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stressm_1,               bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stress12_1,              bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stressp_2,               bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stressm_2,               bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stress12_2,              bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stressp_3,               bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stressm_3,               bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stress12_3,              bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stressp_4,               bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stressm_4,               bndy_info, &
                               field_loc_NEcorner, field_type_vector)
      call update_ghost_cells (stress12_4,              bndy_info, &
                               field_loc_NEcorner, field_type_vector)

      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! Ensure ice is binned in correct categories
      ! (should not be necessary unless restarting from a run with
      !  different category boundaries).
      !
      ! If called, this subroutine does not give exact restart.
      !-----------------------------------------------------------------
!!!      call cleanup_itd

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

      do iblk = 1, nblocks

         call aggregate (nx_block, ny_block, &
                         aicen(:,:,:,iblk),  &
                         trcrn(:,:,:,:,iblk),&
                         vicen(:,:,:,iblk),  &
                         vsnon(:,:,:,iblk),  &
                         eicen(:,:,:,iblk),  &
                         esnon(:,:,:,iblk),  &
                         aice (:,:,  iblk),  &
                         trcr (:,:,:,iblk),  &
                         vice (:,:,  iblk),  &
                         vsno (:,:,  iblk),  &
                         eice (:,:,  iblk),  &
                         esno (:,:,  iblk),  &
                         aice0(:,:,  iblk),  &
                         tmask(:,:,  iblk),  &
                         trcr_depend)

         aice_init(:,:,iblk) = aice(:,:,iblk)

      enddo

      end subroutine restartfile

!=======================================================================
!BOP
!
! !IROUTINE: integer function lenstr(label) - compute length string
!
! !INTERFACE:
!
      integer function lenstr(label)
!
! !DESCRIPTION:
!
! Compute length of string by finding first non-blank
! character from the right.
!
! !REVISION HISTORY:
!
! author:   ?
!
! !INPUT/OUTPUT PARAMETERS:
!
      character*(*) label
!
!EOP
!
      integer (kind=int_kind) :: &
         length, & ! length of character string
         n         ! loop index

      length = len(label)
      do n=length,1,-1
        if( label(n:n) /= ' ' ) exit
      enddo
      lenstr = n

      end function lenstr

!=======================================================================

      end module ice_restart

!=======================================================================
