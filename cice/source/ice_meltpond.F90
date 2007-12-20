!=======================================================================
!
!BOP
!
! !MODULE: ice_meltpond - Meltpond parameterization
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors David A. Bailey (NCAR)
!         Marika M. Holland (NCAR)
!
! !INTERFACE:
!
      module ice_meltpond
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
         tr_pond,       & ! if .true., use explicit meltponds
         restart_pond    ! if .true., read meltponds restart file

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
      subroutine init_meltponds
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
!     local temporary variables

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      integer (kind=int_kind) :: i, j, ij, n, iblk, ilo, ihi, jlo, jhi

      ! Need to compute albedos before init_cpl in CCSM

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      restart_pond = .false.
      if (trim(runtype) == 'continue') restart_pond = .true.

      if (restart_pond) then
         call read_restart_pond
      else
         trcrn(:,:,nt_volpn,:,:) = c0
         apondn(:,:,:,:) = c0
         hpondn(:,:,:,:) = c0
      endif

      end subroutine init_meltponds

!=======================================================================
!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
!
      subroutine compute_ponds(nx_block,ny_block,nghost,   &
                               meltt, melts,  frain,       &
                               aicen, vicen,  vsnon,       &
                               trcrn, apondn, hpondn)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_state, only: nt_Tsfc, nt_volpn
      use ice_calendar, only: dt
      use ice_domain_size, only: ntrcr

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         nghost

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         meltt, &
         melts, &
         frain, &
         aicen, &
         vicen, &
         vsnon

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         apondn, &
         hpondn

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         volpn, &
         Tsfcn

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
         indxi, indxj     ! compressed indices for cells with ice melting

      integer (kind=int_kind) :: i,j,ij,icells,ilo,ihi,jlo,jhi

      real (kind=dbl_kind) :: hi,hs,dTs

      real (kind=dbl_kind), parameter :: &
         hi_min = p1

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      Tsfcn(:,:) = trcrn(:,:,nt_Tsfc)
      volpn(:,:) = trcrn(:,:,nt_volpn)

      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt.
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
         dTs = Timelt - Tsfcn(i,j)

         volpn(i,j) = volpn(i,j) + (c1-rfrac) &
            * (meltt(i,j)*(rhoi/rhofresh) + melts(i,j)*(rhos/rhofresh) &
            + frain(i,j)*dt/rhofresh)

!        Use exponential decrease in pond volume DAB
         if (Tsfcn(i,j) .lt. Timelt-c2) then
            volpn(i,j) = volpn(i,j) &
               *exp(p01*(dTs-c2)/(Timelt-c2))
         endif

         volpn(i,j) = max(volpn(i,j), c0)

         apondn(i,j) = min ( sqrt ( volpn(i,j) * 1.25 ), c1)
         hpondn(i,j) = 0.8 * apondn(i,j)

!        Limit pond depth to 90% of ice thickness
         if (hpondn(i,j) .gt. 0.9*hi) then
            hpondn(i,j) = 0.9*hi
            volpn(i,j) = hpondn(i,j)*apondn(i,j)
!           apondn(i,j) = min(volpn(i,j)/hpondn(i,j), 1.)
         endif

!        If we are freezing, do not change albedo
!        if (Tsfcn(i,j) .lt. Timelt-0.15) apondn(i,j) = c0

!        If we have snow, do not change albedo
         if ( hs .gt. puny ) apondn(i,j) = c0

!        remove ponds if ice becomes very thin
         if (hi .lt. hi_min) then
            apondn(i,j) = c0
            hpondn(i,j) = c0
            volpn(i,j)  = c0
         endif

         trcrn(i,j,nt_volpn) = volpn(i,j)

      enddo

      end subroutine compute_ponds

!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_pond - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_pond(filename_spec)
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
              restart_file(1:lenstr(restart_file)),'.volpn.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_pond,filename,0)

      if (my_task == master_task) then
        write(nu_dump_pond) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do n = 1, ncat
         call ice_write(nu_dump_pond,0,trcrn(:,:,nt_volpn,n,:),'ruf8',diag)
         call ice_write(nu_dump_pond,0,apondn(:,:,n,:),'ruf8',diag)
         call ice_write(nu_dump_pond,0,hpondn(:,:,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_dump_pond)

      end subroutine write_restart_pond

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_pond - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_pond(filename_spec)
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
         if (n == 0) call abort_ice('volpn restart: filename discrepancy')
         string1 = trim(filename0(1:n-1))
         string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
         write(filename,'(a,a,a,a)') &
            string1(1:lenstr(string1)), &
            restart_file(1:lenstr(restart_file)),'.volpn', &
            string2(1:lenstr(string2))
      endif ! master_task

      call ice_open(nu_restart_pond,filename,0)

      if (my_task == master_task) then
        read(nu_restart_pond) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do n = 1, ncat
         call ice_read(nu_restart_pond,0,trcrn(:,:,nt_volpn,n,:),'ruf8',diag)
         call ice_read(nu_restart_pond,0,apondn(:,:,n,:),'ruf8',diag)
         call ice_read(nu_restart_pond,0,hpondn(:,:,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_restart_pond)

      end subroutine read_restart_pond

!=======================================================================

      end module ice_meltpond

!=======================================================================
