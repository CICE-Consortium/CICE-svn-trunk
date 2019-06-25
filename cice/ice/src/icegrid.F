c $Id: icegrid.F,v 1.4 1998/04/14 14:59:53 eclare Exp $

c     program icegrid.F
c   spatial grids, masks, boundary conditions and calendar
c-----------------------------------------------------------------------
      subroutine grid
c.. Horizontal grid initialization
c..   output or read-in variables:
c..     HT{N,E} = cell widths on {N,E} sides of T cell
c..     U{LAT,LONG} = true {latitude,longitude} of U points
c..     D{X,Y}{T,U} = {x,y} spacing centered at {T,U} points
      implicit none
#include "ice.h"
      integer i, j
      real xymnp

c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        hm(i,j) = 0.             ! initialize T-cell mask
       enddo
      enddo

      call popgrid               ! read POP grid lengths directly
c**      call rectgrid           ! regular rectangular grid

      call bound(HTN)
      call bound(HTE)
      call makemask              ! velocity mask

c.. construct T-grid cell and U-grid cell widths

c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
        dxt(i,j) = 0.5*(HTN(i,j) + HTN(i,j-1))
        dyt(i,j) = 0.5*(HTE(i,j) + HTE(i-1,j))

        dxtr(i,j) = 1./dxt(i,j)
        dytr(i,j) = 1./dyt(i,j)

        tarea(i,j) = dxt(i,j)*dyt(i,j)

        dxu(i,j) = 0.5*(HTN(i,j) + HTN(i+1,j))
        dyu(i,j) = 0.5*(HTE(i,j) + HTE(i,j+1))

       enddo
      enddo

      call bound(dxt)
      call bound(dyt)
      call bound(dxu)
      call bound(dyu)
      call bound(tarea)

      xymnp = 1.e+10       ! initialize xymin to big value
c$DOACROSS LOCAL(i,j), REDUCTION(xymnp)
      do j=1,jmt-1
       do i=1,imt-1
         xymnp = min(xymnp,dxt(i,j),dyt(i,j))
       enddo
      enddo
      xymin = xymnp**2
      print*,'min(dx,dy) = ',xymnp

      return
      end
c-----------------------------------------------------------------------
      subroutine popgrid
c.. POP grid and land mask
c..      rec no.         field         units
c..      -------         -----         -----
c..         1             KMT
c..         2            ULAT         radians
c..         3            ULON         radians
c..         4             HTN           cm
c..         5             HTE           cm
c..         6             HUS           cm
c..         7             HUW           cm
c..         8            ANGLE        radians
      implicit none
#include "ice.h"
      integer i, j
      integer KMT  (imt,jmt)             ! boundary mask for thickness (T-cell)
      real    HTEp (imt,jmt)             ! length of eastern edge of T-cell, cm
      real    HTNp (imt,jmt)             ! length of northern edge of T-cell

      open (11,file=grid_file,access='direct',form='unformatted',
     &         recl=imt*jmt*8, status='unknown') ! double precision data
      read(11, rec=1) KMT
      read(11, rec=2) ULAT
      read(11, rec=3) ULONG
      read(11, rec=4) HTNp
      read(11, rec=5) HTEp
      read(11, rec=8) ANGLE
      close (11) 

c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
        hm(i,j)  = real(KMT(i,j))
        if (hm(i,j).ge.1.) hm(i,j) = 1.
        HTE(i,j) = HTEp(i,j)
        HTN(i,j) = HTNp(i,j)
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rectgrid
c.. Regular rectangular grid and mask
c.. Default polar grid with constant long-lat spacing based on 
c..   POP code:  xy_grid
      implicit none
#include "ice.h"
      integer i, j
      real dxc, dyc, xdeg(0:imt), ydeg(0:jmt), radian, radius

      radian = 180./acos(-1.)            ! pi = acos(-1.)
      radius = 6.37e8                    ! radius of Earth, cm

      dxc = 360./real(imt)
      dyc = 180./real(jmt)
      xdeg(0) = -180.
      ydeg(0) = -90.

      do i=1,imt
         xdeg(i) = dxc + xdeg(i-1)       ! constant longitude spacing
      enddo
      do j=1,jmt
         ydeg(j) = dyc + ydeg(j-1)       ! constant latitude  spacing
      enddo

c.. calculate various geometric 2d arrays
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
         ULAT(i,j)  = ydeg(j)/radian     ! radians
         ULONG(i,j) = xdeg(i)/radian     ! radians
         HTN  (i,j) = dxc*radius/radian  ! constant longitude spacing, cm
         HTE (i,j)  = dyc*radius/radian  ! constant latitude  spacing, cm
         ANGLE(i,j) = 0.                 ! "square with the world"
       enddo
      enddo

c.. construct T-cell land mask
c        do j=0,jmt+1                    ! open
c         do i=0,imt+1                   ! open
c$DOACROSS LOCAL(i,j)
        do j=2,jmt-1                     ! closed
         do i=2,imt-1                    ! closed
          hm(i,j) = 1.
         enddo
        enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine makemask
c.. Sets the boundary values for the T cell land mask (hm) and
c.. makes the land mask for U cells (uvm) from hm
      implicit none
#include "ice.h"
      integer i, j

c$DOACROSS LOCAL(i,j)
      do j=0,jmt
       do i=0,imt
        uvm(i,j) = 0.
       enddo
      enddo

      call bound(hm)

c.. construct velocity mask uvm
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
        uvm(i,j) = min(hm(i,j),hm(i+1,j),hm(i,j+1),hm(i+1,j+1))
       enddo
      enddo
      call boundu(uvm)

      return
      end
c-----------------------------------------------------------------------
      subroutine calendar
c.. determine the date and time
      implicit none
#include "ice.h"

      integer k
      real cal(12)                              ! 365-day calendar
      data cal/0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334./
      integer nyr                               ! year number
      integer mday                              ! day of the month
      real fmday                                ! day of the month
      real tday                                 ! absolute day number
      real dayyr                                ! number of days per year
      dayyr = secyr/86400. 

      if (mod(istep,iplot).eq.0) print*, istep

      istep1 = istep + istep0
      time = real(istep1)*dt
      tday = time/86400.                        ! absolute day number (decimal)
      yday = mod(tday,dayyr)                    ! day of the year (decimal)
      monthl = month                            ! (save value - monthly data)
      do k = 1, 12
        if (yday .GE. cal(k)) month = k         ! month
      enddo
      fmday = yday - cal(month)                 ! day of the month (decimal)
      mday = int(fmday) + 1                     ! day of the month (integer)
      nyr = int(tday/dayyr)                     ! year number - 1

      idate = (year+nyr)*10000 + month*100 + mday  ! date (yyyymmdd)
      sec = (fmday - (mday - 1))*86400.         ! elapsed seconds into date

      if (mod(istep,iplot).eq.0) print*,'date: ',idate

      return
      end
c-----------------------------------------------------------------------
      subroutine bound(work1)
c.. Periodic/Neumann boundary conditions for T cells
c.. The only difference between subroutine bound and subroutine boundu
c.. is the size of the array work1.
      implicit none
#include "ice.h"
      integer i, j
      real work1(0:imt+1,0:jmt+1)

c.. Neumann conditions (POP grid land points)
c$DOACROSS LOCAL(i)
      do i=1,imt
        work1(i,0) = work1(i,1)
        work1(i,jmt+1) = work1(i,jmt)
      enddo

c.. Periodic conditions
c$DOACROSS LOCAL(j)
      do j=1,jmt
       work1(0,j) = work1(imt,j)
       work1(imt+1,j) = work1(1,j)
      enddo

c   nearest neighbor corners
      work1(0,0) = work1(1,1)
      work1(imt+1,0) = work1(imt,1)
      work1(0,jmt+1) = work1(1,jmt)
      work1(imt+1,jmt+1) = work1(imt,jmt)

      return
      end
c-----------------------------------------------------------------------
      subroutine boundu(work1)
c.. Periodic/Neumann boundary conditions for U cells
c.. The only difference between subroutine bound and subroutine boundu
c.. is the size of the array work1.
      implicit none
#include "ice.h"
      integer i, j
      real work1(0:imt,0:jmt)

c.. Neumann conditions (POP grid land points)
c$DOACROSS LOCAL(i)
      do i=1,imt-1
        work1(i,0) = work1(i,1)
        work1(i,jmt) = work1(i,jmt-1)
      enddo

c.. Periodic conditions
c$DOACROSS LOCAL(j)
      do j=1,jmt-1
       work1(0,j) = work1(imt-1,j)
       work1(imt,j) = work1(1,j)
      enddo

c   nearest neighbor corners
      work1(0,0) = work1(1,1)
      work1(imt,0) = work1(imt-1,1)
      work1(0,jmt) = work1(1,jmt-1)
      work1(imt,jmt) = work1(imt-1,jmt-1)

      return
      end
