!=======================================================================
!
!BOP
!
! !MODULE: ice_meltpond_cesm - CESM meltpond parameterization
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
      module ice_meltpond_cesm
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
!
!EOP
!
      implicit none

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
      subroutine init_meltponds_cesm
!
! !USES:
!
      use ice_state, only: trcrn, nt_apnd, nt_hpnd
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      trcrn(:,:,nt_apnd,:,:) = c0
      trcrn(:,:,nt_hpnd,:,:) = c0

      end subroutine init_meltponds_cesm

!=======================================================================
!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
!
      subroutine compute_ponds_cesm(nx_block,ny_block,          &
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
      use ice_state, only: nt_Tsfc, nt_apnd, nt_hpnd
      use ice_calendar, only: dt
      use ice_domain_size, only: max_ntrcr
      use ice_itd, only: hi_min
      use ice_restart_meltpond_lvl, only: pndaspect

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
         dpthhi   = 0.9_dbl_kind    ! ratio of pond depth to ice thickness

      !-----------------------------------------------------------------
      ! Initialize 
      !-----------------------------------------------------------------
      Tsfcn(:,:) = trcrn(:,:,nt_Tsfc)
      volpn(:,:) = trcrn(:,:,nt_hpnd) * trcrn(:,:,nt_apnd) * aicen(:,:)

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
            apondn = min (sqrt(volpn(i,j)/(pndaspect*aicen(i,j))), c1)
            hpondn = pndaspect * apondn
            ! fraction of grid cell covered by ponds
            apondn = apondn * aicen(i,j)

            !-----------------------------------------------------------
            ! Limit pond depth
            !-----------------------------------------------------------
             hpondn = min(hpondn, dpthhi*hi)

         endif

         !-----------------------------------------------------------
         ! Reload tracer array
         !-----------------------------------------------------------
         trcrn(i,j,nt_apnd) = apondn / aicen(i,j)
         trcrn(i,j,nt_hpnd) = hpondn

      enddo

      end subroutine compute_ponds_cesm

!=======================================================================
!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
!
      subroutine compute_ponds_simple(nx_block,ny_block,          &
                                      ilo, ihi, jlo, jhi,         &
                                      rfrac, meltt, melts, frain, &
                                      aicen, vicen, vsnon,        &
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
      use ice_calendar, only: dt, istep
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
         hicemin  = c0!hi_min          ! minimum ice thickness with ponds (m)

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

         if (hi < hicemin) then

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

            apondn = c1
            hpondn = volpn(i,j)

            !-----------------------------------------------------------
            ! Reload tracer array
            !-----------------------------------------------------------
            trcrn(i,j,nt_apnd) = apondn / aicen(i,j)
            trcrn(i,j,nt_hpnd) = hpondn
!            trcrn(i,j,nt_volp) = volpn(i,j) / aicen(i,j)

         endif

      enddo

      end subroutine compute_ponds_simple

!=======================================================================

      end module ice_meltpond_cesm

!=======================================================================
