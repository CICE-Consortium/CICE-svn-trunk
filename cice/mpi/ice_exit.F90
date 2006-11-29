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
!  SVN:$Id$
!
! authors William H. Lipscomb (LANL)
!         Elizabeth C. Hunke (LANL)
! 2006 ECH: separated serial and mpi functionality
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
      use ice_fileunits
      use ice_communicate
#ifdef CCSM
      use shr_sys_mod
#endif

      include 'mpif.h'   ! MPI Fortran include file
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (len=*), intent(in) :: error_message
!
!EOP
!
      integer (int_kind) :: ierr ! MPI error flag

#ifdef CCSM
      call shr_sys_abort(error_message)
#else
      write (nu_diag,*) error_message
      call MPI_ABORT(MPI_COMM_WORLD, ierr)
      stop
#endif

      end subroutine abort_ice

!=======================================================================
!BOP
!
! !IROUTINE: end_run - ends run
!
! !INTERFACE:
!
      subroutine end_run
!
! !DESCRIPTION:
!
! Ends run by calling MPI_FINALIZE.
!
! !REVISION HISTORY:
!
! author: ?
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!

      integer (int_kind) :: ierr ! MPI error flag

      call MPI_FINALIZE(ierr)
!
!EOP
!
      end subroutine end_run

!=======================================================================

      end module ice_exit

!=======================================================================
