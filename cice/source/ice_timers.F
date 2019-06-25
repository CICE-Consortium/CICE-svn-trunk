! $Id: ice_timers.F,v 1.5 2004/02/25 17:35:14 eclare Exp $
!=======================================================================
!
!BOP
!
! !MODULE: ice_timers
!
! !DESCRIPTION:
!
! Timing routines
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !INTERFACE:
!
      module ice_timers
!
! !DESCRIPTION:
!
! Timing routines
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
!
!EOP
!
      implicit none

      integer (kind=int_kind), parameter, private ::
     &   nb_chronos = 10        ! total number of timers

      real (kind=dbl_kind), private ::
     &   tchrono(0:nb_chronos)  ! elapsed time
     &,  uchrono(0:nb_chronos)  ! initial clock time

      character (len=8), private :: timern(0:nb_chronos)  ! timer labels

      data timern /'Total'   ,
     &  'TimeLoop','Dynamics','Advectn' ,'Column'  ,'Thermo',  
     &  'Ridging' ,'Cat Conv','Coupling','ReadWrit','Bound' /

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: ice_timer_clear(n) - initialize timer n to 0
!
! !INTERFACE:
!
      subroutine ice_timer_clear(n)
!
! !DESCRIPTION:
!
! Initialize timer n to 0 \\
! if n = -1 initialize all timers
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: 
     &   n                      ! timer number
!
!EOP
!
      integer (kind=int_kind) :: l

      if (n.eq.-1) then
        do l=0,nb_chronos
          tchrono(l)=c0
          uchrono(l)=c0
        enddo
      else
        tchrono(n)=c0
        uchrono(n)=c0
      endif

      end subroutine ice_timer_clear

!=======================================================================
!BOP
!
! !IROUTINE: ice_timer_start(n) - begin timing with timer n
!
! !INTERFACE:
!
      subroutine ice_timer_start(n)
!
! !DESCRIPTION:
!
! Begin timing with timer n
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: 
     &    n                     ! timer number
!
!EOP
!
      call timers(uchrono(n))

      end subroutine ice_timer_start

!=======================================================================
!BOP
!
! !IROUTINE: ice_timer_stop(n) - end (or pause) timing with timer n
!
! !INTERFACE:
!
      subroutine ice_timer_stop(n)
!
! !DESCRIPTION:
!
! End (or pause) timing with timer n
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: 
     &    n                     ! timer number
!
!EOP
!
      real (kind=dbl_kind) :: end_time

      call timers(end_time)
      tchrono(n) = tchrono(n) + (end_time - uchrono(n))

      end subroutine ice_timer_stop

!=======================================================================
!BOP
!
! !IROUTINE: ice_timer_print - print timing results of timer n 
!
! !INTERFACE:
!
      subroutine ice_timer_print(n)
!
! !DESCRIPTION:
!
! Print timing results of timer n \\
! if n = -1 print timing results of all timers
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
      use ice_domain 
      use ice_mpi_internal
      use ice_fileunits
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: 
     &    n                     ! timer number
!
!EOP
!
      integer (kind=int_kind) :: l
      real (kind=dbl_kind) :: amin,amax
      real (kind=dbl_kind) :: timerw(1)

      if (n.eq.-1) then   ! print all timers
        do l=0,nb_chronos
        timerw(1)=tchrono(l)
        amin = ice_global_real_minval(1,timerw)
        amax = ice_global_real_maxval(1,timerw)
        if (my_task.eq.master_task) then
          write (nu_diag,100) l,timern(l),tchrono(l),amin,amax
        endif
        enddo
      else                ! print timer n
        timerw(1)=tchrono(n)
        amin = ice_global_real_minval(1,timerw)
        amax = ice_global_real_maxval(1,timerw)
        if (my_task.eq.master_task) then
          write (nu_diag,100) n,timern(n),tchrono(n),amin,amax
        endif
      endif
  100 format ('Timer number',i3,2x,a8,' =',f9.2,' seconds',
     &  '  min/max = ',2f9.2)

      end subroutine ice_timer_print

!=======================================================================
!BOP
!
! !IROUTINE: timers - do the work
!
! !INTERFACE:
!
      subroutine timers(t1)
!
! !DESCRIPTION:
!
! Do the work
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(out) :: t1
!
!EOP
!

#ifdef _MPI
      include "mpif.h"         ! MPI library definitions
#else
#ifdef Linux
      integer (kind=int_kind) ::
     &     count        !  current count of the system clock
     &,    count_rate   !  number of clock ticks per second
#else
      real (kind=dbl_kind) :: rtc
      integer (kind=int_kind) :: rate,irtc_rate
#endif
#endif

#ifdef _MPI
      t1=MPI_WTIME()
#else
#ifdef IRIX64
      rate=irtc_rate()
      t1=rtc()/rate
#else
# ifdef OSF1
      t1=c0
# else
#  ifdef Linux
      t1=c0
#  else
      t1=rtc()
#  endif
# endif
#endif
#endif

      end subroutine timers

!=======================================================================

      end module ice_timers

!=======================================================================

