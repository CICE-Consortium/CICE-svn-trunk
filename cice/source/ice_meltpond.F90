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
! !ROUTINE: 
!
! !INTERFACE:
!
      subroutine compute_ponds(nx_block,ny_block,nghost,   &
                               meltt, melts,               &
                               aicen, vicen, vsnon, Tsfcn, &
                               volpn, apondn, hpondn)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         nghost

      real (kind=dbl_kind),intent(in) :: &
         meltt(nx_block,ny_block), &
         melts(nx_block,ny_block), &
         aicen(nx_block,ny_block), &
         vicen(nx_block,ny_block), &
         vsnon(nx_block,ny_block), &
         Tsfcn(nx_block,ny_block)

      real (kind=dbl_kind),intent(inout) :: &
         volpn(nx_block,ny_block),  &
         apondn(nx_block,ny_block), &
         hpondn(nx_block,ny_block)

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
            * (meltt(i,j)*(rhoi/rhofresh) + melts(i,j)*(rhos/rhofresh))

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

      enddo

      end subroutine compute_ponds

!=======================================================================

      end module ice_meltpond

!=======================================================================
