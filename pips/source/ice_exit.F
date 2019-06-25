!=======================================================================
!
!BOP
!
! !MODULE: ice_exit - exit the model
!
! !DESCRIPTION:
!
! Exit the model.  Logically, this routine should be used for "normal"
! exit of the ice model (called from CICE.F), but there would be only 
! one such call, and creating a subroutine here to accomplish that
! causes circular dependencies due to the coupler exit strategy.  Hence
! this module is used only to isolate CCSM "shr-" code calls during an
! abort.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb (LANL)
!         Elizabeth C. Hunke (LANL)
!
! !INTERFACE:
!
      module ice_exit
!
! !USES:
!
      use ice_kinds_mod
!
!EOP
!
      implicit none

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: abort_ice - abort the model
!
! !INTERFACE:
!
      subroutine abort_ice(error_message)
!
! !DESCRIPTION:
!
!  This routine aborts the ice model and prints an error message.
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_domain
      use ice_fileunits
      use ice_mpi_internal
#ifdef coupled
      use shr_sys_mod
#endif
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (len=*), intent(in) :: error_message
!
!EOP
!

#ifdef coupled
      call shr_sys_abort(error_message)
#else
      write (nu_diag,*) error_message
      call end_run
      stop
#endif

      end subroutine abort_ice

!=======================================================================

      end module ice_exit

!=======================================================================
