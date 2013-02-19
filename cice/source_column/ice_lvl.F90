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
      subroutine init_lvl 
!
! !USES:
!
      use ice_state, only: nt_alvl, nt_vlvl, trcrn, aicen, vicen
!
!EOP
!
 
      ! assume all ice is level
      trcrn(:,:,nt_alvl,:,:) = c1 ! level ice area fraction
      trcrn(:,:,nt_vlvl,:,:) = c1 ! level ice volume

      end subroutine init_lvl

!=======================================================================

      end module ice_lvl

!=======================================================================
