!=======================================================================
!
!BOP
!
! !MODULE: ice_lvl - Ridged ice tracers for sea ice
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
      module ice_lvl
!
! !USES:
!
      use ice_kinds_mod
!
!EOP
!
      implicit none

      private
      public :: init_lvl

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_lvl
!
! !DESCRIPTION:
!
!  Initialize ice lvl tracers (call prior to reading restart data)
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_lvl(nx_block, ny_block, ncat, alvl, vlvl) 
!
! !USES:
!
      use ice_constants, only: c1
!
        integer(kind=int_kind), intent(in) :: &
             nx_block , &
             ny_block , &
             ncat

        real(kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
             intent(out) :: &
             alvl , & ! level ice area fraction
             vlvl     ! level ice volume
!
!EOP
!
      alvl(:,:,:) = c1 ! level ice area fraction
      vlvl(:,:,:) = c1 ! level ice volume

      end subroutine init_lvl

!=======================================================================

      end module ice_lvl

!=======================================================================
