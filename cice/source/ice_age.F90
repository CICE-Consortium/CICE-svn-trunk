!  SVN:$Id$
!=======================================================================
!
! authors Elizabeth Hunke

      module ice_age

      use ice_kinds_mod
      use ice_constants

      implicit none

      private
      public :: init_age, increment_age

!=======================================================================

      contains

!=======================================================================

!  Initialize ice age tracer (call prior to reading restart data)

      subroutine init_age(nx_block, ny_block, ncat, iage)

      integer(kind=int_kind), intent(in) :: &
             nx_block , &
             ny_block , &
             ncat

      real(kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
             intent(out) :: iage

      iage(:,:,:) = c0

      end subroutine init_age

!=======================================================================

!  Increase ice age tracer by timestep length.

      subroutine increment_age (nx_block, ny_block, &
                                dt,       icells,   &
                                indxi,    indxj,    &
                                iage)

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

      !  local variables

      integer (kind=int_kind) :: i, j, ij

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         iage(i,j) = iage(i,j) + dt 
      enddo

      end subroutine increment_age

!=======================================================================

      end module ice_age

!=======================================================================
