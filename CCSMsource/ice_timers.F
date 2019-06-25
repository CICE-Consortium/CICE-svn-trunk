! $Id: ice_timers.F,v 1.10 2006/06/08 22:16:07 eclare Exp $
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
      save

#ifdef CCSMcoupled
      integer (kind=int_kind), parameter, private ::
     &   nb_chronos = 20        ! total number of timers
#else
      integer (kind=int_kind), parameter, private ::
     &   nb_chronos = 10        ! total number of timers
#endif

#ifndef _MPI
      integer (kind=int_kind) ::
     &   cycles_max   !  number of clock ticks per second

      real (kind=dbl_kind) ::
     &   clock_rate   !  number of clock ticks per second
#endif

      real (kind=dbl_kind), private ::
     &   tchrono(0:nb_chronos)  ! elapsed time
     &,  uchrono(0:nb_chronos)  ! initial clock time

      character (len=8), private :: timern(0:nb_chronos)  ! timer labels

#ifdef CCSMcoupled
      data timern /'Total'   ,
     &  'TimeLoop','Dynamics','Advectn' ,'Column'  ,'Thermo',  
     &  'Ridging' ,'Cat Conv','Coupling','ReadWrit','Bound',
     &  'Pre-cpl' ,'MPI-send','MPI-recv','Snd->Rcv','Rcv->Snd',
     &  'Cpl-recv','CR-unpck','CS-pack' ,'Cpl-send',' ' / 
#else
      data timern /'Total'   ,
     &  'TimeLoop','Dynamics','Advectn' ,'Column'  ,'Thermo',  
     &  'Ridging' ,'Cat Conv','Coupling','ReadWrit','Bound' /
#endif

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: ice_timer_init
!
! !INTERFACE:
!
      subroutine ice_timer_init
!
! !DESCRIPTION:
!
! initialize F90 intrinsic system_clock if MPI timer not available
!
! !REVISION HISTORY:
!
! author: Elizabeth Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!

#ifndef _MPI
      integer (kind=int_kind) :: cycles

      cycles = 0
      clock_rate = c0

      ! F90 intrinsic
      call system_clock(count_rate=cycles, count_max=cycles_max)

      if (cycles /= 0) then
        clock_rate = c1/real(cycles,kind=dbl_kind)
      else
        clock_rate = c0
        write(6,*) ' '
        write(6,'(a33)') '--- No system clock available ---'
        write(6,*) ' '
      endif
#endif

      call ice_timer_clear(-1) ! initialize all timers

      end subroutine ice_timer_init

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

      if (n == -1) then
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
!         Elizabeth Hunke, LANL, modified for non-MPI runs 5/2006
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

#ifdef _MPI
      tchrono(n) = tchrono(n) + (end_time - uchrono(n))
#else
      if (end_time >= uchrono(n)) then
        tchrono(n) = tchrono(n) 
     &             + clock_rate*(end_time - uchrono(n))
      else  ! correct wrap-around
            ! NOTE this only works for one wrap-around
        tchrono(n) = tchrono(n) 
     &             + clock_rate*(cycles_max + end_time - uchrono(n))
      endif
#endif

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

      if (n == -1) then   ! print all timers
        if (my_task == master_task) write (nu_diag,*) ' '
        do l=0,nb_chronos
        timerw(1)=tchrono(l)
        amin = ice_global_real_minval(1,timerw)
        amax = ice_global_real_maxval(1,timerw)
        if (my_task == master_task) then
          write (nu_diag,100) l,timern(l),tchrono(l),amin,amax
        endif
        enddo
      else                ! print timer n
        timerw(1)=tchrono(n)
        amin = ice_global_real_minval(1,timerw)
        amax = ice_global_real_maxval(1,timerw)
        if (my_task == master_task) then
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

      t1=MPI_WTIME()

#else
      integer (int_kind) :: cycles

      ! F90 intrinsic
      call system_clock(count=cycles)
      t1 = real(cycles,kind=dbl_kind)

#endif

      end subroutine timers

!=======================================================================

      end module ice_timers

!=======================================================================

