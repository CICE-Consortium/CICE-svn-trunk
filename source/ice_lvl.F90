!  SVN:$Id$
!=======================================================================

! Ridged ice tracers for sea ice
!
! authors Elizabeth Hunke

      module ice_lvl

      use ice_kinds_mod

      implicit none

      private
      public :: init_lvl

!=======================================================================

      contains

!=======================================================================

!  Initialize ice lvl tracers (call prior to reading restart data)

      subroutine init_lvl(nx_block, ny_block, ncat, alvl, vlvl) 

      use ice_constants, only: c1

      integer(kind=int_kind), intent(in) :: &
             nx_block , &
             ny_block , &
             ncat

      real(kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
             intent(out) :: &
             alvl , & ! level ice area fraction
             vlvl     ! level ice volume

      alvl(:,:,:) = c1 ! level ice area fraction
      vlvl(:,:,:) = c1 ! level ice volume

      end subroutine init_lvl

!=======================================================================

      end module ice_lvl

!=======================================================================
