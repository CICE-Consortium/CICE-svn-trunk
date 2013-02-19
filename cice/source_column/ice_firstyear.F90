!=======================================================================
!
!BOP
!
! !MODULE: ice_firstyear - First year concentration tracer for sea ice
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
      module ice_firstyear
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
! !ROUTINE: init_FY
!
! !DESCRIPTION:
!
!  Initialize ice FY tracer (call prior to reading restart data)
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_FY 
!
! !USES:
!
      use ice_state, only: nt_FY, trcrn
!
!EOP
!
      trcrn(:,:,nt_FY,:,:) = c0

      end subroutine init_FY

!=======================================================================

!BOP
!
! !ROUTINE: update_FYarea 
!
! !DESCRIPTION:
!
!  Zero ice FY tracer on fixed day of year. Zeroing FY ice tracer promotes
!  ice to MY ice. Unfortunately some frazil ice may grow before the 
!  zeroing date and thus get promoted to MY ice too soon.
!  Bummer.
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine update_FYarea (nx_block, ny_block, &
                                dt,       icells,   &
                                indxi,    indxj,    &
                                nhmask,   shmask,   &
                                FYarea)
!
! !USES:
!
      use ice_calendar, only: secday, yday
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      logical (kind=log_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         nhmask, shmask

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         FYarea
!
!  local variables
!
      integer (kind=int_kind) :: i, j, ij
!
!EOP
!
     if ((yday >= 259._dbl_kind) .and. &
           (yday <  259._dbl_kind+dt/secday)) then
        do ij = 1, icells
           i = indxi(ij)
           j = indxj(ij)
           if (nhmask(i,j)) FYarea(i,j) = c0;
        enddo
      endif

     if ((yday >= 75._dbl_kind) .and. &
           (yday <  75._dbl_kind+dt/secday)) then
        do ij = 1, icells
           i = indxi(ij)
           j = indxj(ij)
           if (shmask(i,j)) FYarea(i,j) = c0;
        enddo
      endif

      end subroutine update_FYarea

!=======================================================================

      end module ice_firstyear

!=======================================================================
