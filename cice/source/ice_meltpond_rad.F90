!=======================================================================
!
!BOP
!
! !MODULE: ice_meltpond_rad - Meltpond parameterization
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
!
! !INTERFACE:
!
      module ice_meltpond_rad
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
      use ice_exit, only: abort_ice
!
!EOP
!
      implicit none

      logical (kind=log_kind) :: & 
         restart_pond     ! if .true., read meltponds restart file

      real (kind=dbl_kind) :: &
         hs0              ! snow depth for transition to bare sea ice (m)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_meltponds
!
! !DESCRIPTION:
!
!  Initialize melt ponds.
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_meltponds_rad
!
! !USES:
!
      use ice_domain_size
      use ice_blocks
      use ice_domain
      use ice_flux
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      if (trim(runtype) == 'continue') restart_pond = .true.

      if (restart_pond) then
         call read_restart_pond_rad
      else
         trcrn(:,:,nt_apnd,:,:) = c0
         trcrn(:,:,nt_hpnd,:,:) = c0
!         trcrn(:,:,nt_volp,:,:) = c0
      endif

      end subroutine init_meltponds_rad

!=======================================================================
!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
!
      subroutine compute_ponds_rad(nx_block,ny_block,          &
                                   ilo, ihi, jlo, jhi,         &
                                   rfrac, meltt, melts,  frain,&
                                   aicen, vicen,  vsnon,       &
                                   trcrn)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
!      use ice_state, only: nt_Tsfc, nt_volp, nt_apnd, nt_hpnd
      use ice_state, only: nt_Tsfc, nt_apnd, nt_hpnd
      use ice_calendar, only: dt
      use ice_domain_size, only: max_ntrcr
      use ice_itd, only: hi_min

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         rfrac, &    ! water fraction retained for melt ponds
         meltt, &
         melts, &
         frain, &
         aicen, &
         vicen, &
         vsnon

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_ntrcr), &
         intent(inout) :: &
         trcrn

!     local temporary variables

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         volpn, &
         Tsfcn

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
         indxi, indxj     ! compressed indices for cells with ice melting

      integer (kind=int_kind) :: i,j,ij,icells

      real (kind=dbl_kind) :: &
         hi                     , & ! ice thickness (m)
         hs                     , & ! snow depth (m)
         dTs                    , & ! surface temperature diff for freeze-up (C)
         Tp                     , & ! pond freezing temperature (C)
         asnow                  , & ! area fraction of snow on ice
         apondn, &
         hpondn   

      real (kind=dbl_kind), parameter :: &
         Td       = c2          , & ! temperature difference for freeze-up (C)
         rexp     = p01         , & ! pond contraction scaling
         dpthhi   = 0.9_dbl_kind, & ! ratio of pond depth to ice thickness
         dpthfrac = 0.8_dbl_kind    ! ratio of pond depth to pond fraction

      !-----------------------------------------------------------------
      ! Initialize 
      !-----------------------------------------------------------------
      Tsfcn(:,:) = trcrn(:,:,nt_Tsfc)
      volpn(:,:) = trcrn(:,:,nt_hpnd) * trcrn(:,:,nt_apnd) * aicen(:,:)
!      volpn(:,:) = trcrn(:,:,nt_volp) * aicen(:,:)

      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (aicen(i,j) > puny) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo                     ! i
      enddo                     ! j

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hi = vicen(i,j)/aicen(i,j)
         hs = vsnon(i,j)/aicen(i,j)

         if (hi < hi_min) then

         !--------------------------------------------------------------
         ! Remove ponds on thin ice
         !--------------------------------------------------------------
            apondn = c0
            hpondn = c0
            volpn (i,j) = c0

         else

            !-----------------------------------------------------------
            ! Update pond volume
            !-----------------------------------------------------------
            volpn(i,j) = volpn(i,j) &
                       + rfrac(i,j)/rhofresh*(meltt(i,j)*rhoi &
                       +                      melts(i,j)*rhos &
                       +                      frain(i,j)*  dt)&
                       * aicen(i,j)

            !-----------------------------------------------------------
            ! Shrink pond volume under freezing conditions
            !-----------------------------------------------------------
            Tp = Timelt - Td
            dTs = max(Tp - Tsfcn(i,j),c0)
            volpn(i,j) = volpn(i,j) * exp(rexp*dTs/Tp)
            volpn(i,j) = max(volpn(i,j), c0)

            ! fraction of ice covered by ponds
            apondn = min (sqrt(volpn(i,j)/(dpthfrac*aicen(i,j))), c1)
            hpondn = dpthfrac * apondn
            ! fraction of grid cell covered by ponds
            apondn = apondn * aicen(i,j)

            !-----------------------------------------------------------
            ! Limit pond depth
            !-----------------------------------------------------------
             hpondn = min(hpondn, dpthhi*hi)
             volpn(i,j) = hpondn*apondn

            !-----------------------------------------------------------
            ! If surface is freezing or has snow, do not change albedo
            !-----------------------------------------------------------
!            if (Tsfcn(i,j) < Tp) apondn = c0
!            if (hs > puny) apondn = c0

            !-----------------------------------------------------------
            ! reduce pond area if there is snow, preserving volume
            !-----------------------------------------------------------
!echmod - do this in ice_shortwave.F90
!            if (hs >= hsmin) then
!               asnow = min(hs/hs0, c1) ! delta-Eddington formulation
!               apondn = (c1 - asnow) * apondn
!            endif

            !-----------------------------------------------------------
            ! Reload tracer array
            !-----------------------------------------------------------
            trcrn(i,j,nt_apnd) = apondn / aicen(i,j)
            trcrn(i,j,nt_hpnd) = hpondn
!            trcrn(i,j,nt_volp) = volpn(i,j) / aicen(i,j)

         endif

      enddo

      end subroutine compute_ponds_rad

!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_pond - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_pond_rad(filename_spec)
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
!              restart_file(1:lenstr(restart_file)),'.volpn.', &
              restart_file(1:lenstr(restart_file)),'.pond.', &
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
!         call ice_write(nu_dump_pond,0,trcrn(:,:,nt_volp,n,:),'ruf8',diag)
         call ice_write(nu_dump_pond,0,trcrn(:,:,nt_hpnd,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_dump_pond)

      end subroutine write_restart_pond_rad

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_pond - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_pond_rad(filename_spec)
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
!         if (n == 0) call abort_ice('volpn restart: filename discrepancy')
         if (n == 0) call abort_ice('pond restart: filename discrepancy')
         string1 = trim(filename0(1:n-1))
         string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
         write(filename,'(a,a,a,a)') &
            string1(1:lenstr(string1)), &
!            restart_file(1:lenstr(restart_file)),'.volpn', &
            restart_file(1:lenstr(restart_file)),'.pond', &
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
!         call ice_read(nu_restart_pond,0,trcrn(:,:,nt_volp,n,:),'ruf8',diag)
         call ice_read(nu_restart_pond,0,trcrn(:,:,nt_hpnd,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_restart_pond)

      end subroutine read_restart_pond_rad

!=======================================================================

      end module ice_meltpond_rad

!=======================================================================
