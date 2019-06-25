c $Id: timers.F,v 2.0 1999/12/09 22:25:36 eclare Exp $

c-----------------------------------------------------------------------
c.. timing routines
c-----------------------------------------------------------------------

      subroutine timer_clear(n)

      implicit none
      include "timers.h"       ! ice common blocks for timers

      integer n,l

      if (n.eq.-1) then
        do l=0,nb_chronos
          tchrono(l)=0.
          uchrono(l)=0.
        enddo
      else
        tchrono(n)=0.
        uchrono(n)=0.
      endif

      end
 
************************************************************************
                    
      subroutine timer_start(n)

      implicit none
      include "timers.h"       ! ice common blocks for timers

      integer n
      real elapse

      call timers(uchrono(n), elapse)

      end
 
************************************************************************

      subroutine timer_stop(n)

      implicit none
      include "timers.h"       ! ice common blocks for timers

      integer n
      real end_time, elapse

      call timers(end_time, elapse)
      tchrono(n) = tchrono(n) + (end_time - uchrono(n))

      end
 
************************************************************************

      subroutine timer_print(n)

      implicit none
      include "ice.h"          ! ice code common blocks
      include "timers.h"       ! ice common blocks for timers

      integer n,l
      real timer,amin,amax,rate

      if (n.eq.-1) then
        do l=0,nb_chronos
        timer=tchrono(l)
        amin = ice_global_real_minval(1,timer)
        amax = ice_global_real_maxval(1,timer)
        if (my_task.eq.master_task) then
          write (6,100) l,timern(l),tchrono(l),amin,amax
        endif
        enddo
      else
        timer=tchrono(n)
        amin = ice_global_real_minval(1,timer)
        amax = ice_global_real_maxval(1,timer)
        if (my_task.eq.master_task) then
          write (6,100) n,timern(n),tchrono(n),amin,amax
        endif
      endif
  100 format ('Timer number',i3,2x,a8,' =',f9.2,' seconds',
     &  '  min/max = ',2f9.2)

      end

************************************************************************
      subroutine timers(t1,e1)
      implicit none
#ifdef _MPI
      include "mpif.h"         ! MPI library definitions
#else
      real rtc
      integer(kind=8) rate,irtc_rate
#endif
      real t1,e1
#ifdef _MPI
      t1=MPI_WTIME()
#else
#ifdef IRIX64
      rate=irtc_rate()
      t1=rtc()/rate
#else
      t1=rtc()
#endif
#endif
      return
      end


