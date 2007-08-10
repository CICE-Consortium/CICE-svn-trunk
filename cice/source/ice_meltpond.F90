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
!
!EOP
!
      implicit none

      integer (kind=int_kind) :: & ! defined in namelist 
         kpond          ! 1 = explicit meltponds

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

      ! Local variables to keep track of melt for ponds
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         melts_old, &
         meltt_old, &
         melts_tmp, &
         meltt_tmp

      integer (kind=int_kind) :: i, j, ij, n, iblk, ilo, ihi, jlo, jhi

      ! Need to compute albedos before init_cpl in CCSM

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do iblk=1,nblocks
      do n=1,ncat

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n,iblk) > puny) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

      !-----------------------------------------------------------------
      ! Melt pond initialization
      !-----------------------------------------------------------------

         apondn(:,:,n,iblk) = c0
         hpondn(:,:,n,iblk) = c0

         if (kpond == 1) then

            melts_tmp = c0
            meltt_tmp = c0

            call compute_ponds(nx_block, ny_block, nghost,              &
                               meltt_tmp, melts_tmp, frain(:,:,iblk),   &
                               aicen (:,:,n,iblk), vicen (:,:,n,iblk),  &
                               vsnon (:,:,n,iblk), trcrn (:,:,:,n,iblk),&
                               apondn(:,:,n,iblk), hpondn(:,:,n,iblk))

         endif

      enddo
      enddo

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

      end module ice_meltpond

!=======================================================================
