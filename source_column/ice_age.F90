!=======================================================================
!
!BOP
!
! !MODULE: ice_age - Age tracer for sea ice
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_age
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
! !ROUTINE: init_age
!
! !DESCRIPTION:
!
!  Initialize ice age tracer (call prior to reading restart data)
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_age 
!
! !USES:
!
      use ice_state, only: nt_iage, trcrn
!
!EOP
!
         trcrn(:,:,nt_iage,:,:) = c0

      end subroutine init_age

!=======================================================================

!BOP
!
! !ROUTINE: increment_age 
!
! !DESCRIPTION:
!
!  Increase ice age tracer by timestep length.
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine increment_age (nx_block, ny_block, &
                                dt,       icells,   &
                                indxi,    indxj,    &
                                iage)
!
! !USES:
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

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         iage
!
!  local variables
!
      integer (kind=int_kind) :: i, j, ij
!
!EOP
!
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         iage(i,j) = iage(i,j) + dt 
      enddo

      end subroutine increment_age

!=======================================================================

      end module ice_age

!=======================================================================
