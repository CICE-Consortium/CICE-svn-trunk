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
! 2008 ECH: Rearranged order in which internal stresses are written and read
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
      ! Tsfc is the only tracer written to this file.  All other
      ! tracers are written to their own dump/restart files.
      !-----------------------------------------------------------------

      do n=1,ncat
         call ice_write(nu_dump,0,aicen(:,:,n,:),'ruf8',diag)
         call ice_write(nu_dump,0,vicen(:,:,n,:),'ruf8',diag)
         call ice_write(nu_dump,0,vsnon(:,:,n,:),'ruf8',diag)
         call ice_write(nu_dump,0,trcrn(:,:,nt_Tsfc,n,:),'ruf8',diag)
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
      ! radiation fields
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,scale_factor,'ruf8',diag)

      call ice_write(nu_dump,0,swvdr,'ruf8',diag)
      call ice_write(nu_dump,0,swvdf,'ruf8',diag)
      call ice_write(nu_dump,0,swidr,'ruf8',diag)
      call ice_write(nu_dump,0,swidf,'ruf8',diag)

      !-----------------------------------------------------------------
      ! ocean stress (for bottom heat flux in thermo)
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,strocnxT,'ruf8',diag)
      call ice_write(nu_dump,0,strocnyT,'ruf8',diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,stressp_1,'ruf8',diag)
      call ice_write(nu_dump,0,stressp_3,'ruf8',diag)
      call ice_write(nu_dump,0,stressp_2,'ruf8',diag)
      call ice_write(nu_dump,0,stressp_4,'ruf8',diag)

      call ice_write(nu_dump,0,stressm_1,'ruf8',diag)
      call ice_write(nu_dump,0,stressm_3,'ruf8',diag)
      call ice_write(nu_dump,0,stressm_2,'ruf8',diag)
      call ice_write(nu_dump,0,stressm_4,'ruf8',diag)

      call ice_write(nu_dump,0,stress12_1,'ruf8',diag)
      call ice_write(nu_dump,0,stress12_3,'ruf8',diag)
      call ice_write(nu_dump,0,stress12_2,'ruf8',diag)
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
      subroutine restartfile (ice_ic)
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
      use ice_calendar, only: istep0, istep1, time, time_forc, calendar
      use ice_flux
      use ice_state
      use ice_grid, only: tmask
      use ice_itd
      use ice_ocean, only: oceanmixed_ice
      use ice_work, only: work1, work_g1, work_g2
      use ice_gather_scatter, only: scatter_global_stress
!
! !INPUT/OUTPUT PARAMETERS:

      character (*), optional :: ice_ic
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, k, n, it, iblk ! counting indices

      character(len=char_len_long) :: &
         filename, filename0

      logical (kind=log_kind) :: &
         diag

      if (present(ice_ic)) then
         filename = ice_ic
      elseif (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         read(nu_rst_pointer,'(a)') filename0
         filename = trim(filename0)
         close(nu_rst_pointer)
         write(nu_diag,*) 'Read ',pointer_file(1:lenstr(pointer_file))
      endif

      call ice_open(nu_restart,filename,0)

      if (my_task == master_task) then
         write(nu_diag,*) 'Using restart dump=', trim(filename)
         read (nu_restart) istep0,time,time_forc
         write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc
      endif
      call calendar(time)

      call broadcast_scalar(istep0,master_task)

      istep1 = istep0

      call broadcast_scalar(time,master_task)
      call broadcast_scalar(time_forc,master_task)

      diag = .true.     ! write min/max diagnostics for field

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      do n=1,ncat
         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' min/max area, vol ice, vol snow, Tsfc'

         call ice_read(nu_restart,0,aicen(:,:,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,vicen(:,:,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,vsnon(:,:,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,trcrn(:,:,nt_Tsfc,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      enddo

      if (my_task == master_task) &
           write(nu_diag,*) 'min/max eicen for each layer'
      do k=1,ntilyr
         call ice_read(nu_restart,0,eicen(:,:,k,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      enddo

      if (my_task == master_task) &
           write(nu_diag,*) 'min/max esnon for each layer'
      do k=1,ntslyr
         call ice_read(nu_restart,0,esnon(:,:,k,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max velocity components'

      call ice_read(nu_restart,0,uvel,'ruf8',diag, &
                       field_loc_NEcorner, field_type_vector)
      call ice_read(nu_restart,0,vvel,'ruf8',diag, &
                       field_loc_NEcorner, field_type_vector)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      if (my_task == master_task) &
         write(nu_diag,*) 'radiation fields'

      call ice_read(nu_restart,0,scale_factor,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swvdr,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swvdf,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swidr,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swidf,'ruf8',diag, &
                    field_loc_center, field_type_scalar)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max ocean stress components'

      call ice_read(nu_restart,0,strocnxT,'ruf8',diag, &
                       field_loc_center, field_type_vector)
      call ice_read(nu_restart,0,strocnyT,'ruf8',diag, &
                       field_loc_center, field_type_vector)

      !-----------------------------------------------------------------
      ! internal stress
      ! The stress tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
      if (my_task == master_task) write(nu_diag,*) &
           'internal stress components'
      
      allocate (work_g1(nx_global,ny_global), &
                work_g2(nx_global,ny_global))

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressp_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressp_3
      call scatter_global_stress(stressp_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressp_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressp_4
      call scatter_global_stress(stressp_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressm_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressm_3
      call scatter_global_stress(stressm_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressm_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressm_4
      call scatter_global_stress(stressm_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stress12_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stress12_3
      call scatter_global_stress(stress12_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stress12_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stress12_4
      call scatter_global_stress(stress12_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      deallocate (work_g1, work_g2)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'ice mask for dynamics'

      call ice_read(nu_restart,0,work1,'ruf8',diag, &
                    field_loc_center, field_type_scalar)

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

         call ice_read(nu_restart,0,sst,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,frzmlt,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      endif

      if (my_task == master_task) close(nu_restart)

      !-----------------------------------------------------------------
      ! Ensure unused stress values in west and south ghost cells are 0
      !-----------------------------------------------------------------
      do iblk = 1, nblocks
         do j = 1, nghost
         do i = 1, nx_block
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
         do j = 1, ny_block
         do i = 1, nghost
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
      enddo

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
