c $Id: flux.F,v 2.0 1999/12/09 22:25:34 eclare Exp $
c-----------------------------------------------------------------------
c.. reads data from files and interpolates if needed
c-----------------------------------------------------------------------
      subroutine getflux
c.. data for stand-alone runs
      implicit none
      include "ice.h"          ! ice code common blocks
      integer nsec

      if (kdata.eq.1) then 

       if (istep.eq.1.and.my_task.eq.master_task) then  

         open(41,file='forcing.dat',status='unknown')
         read(41,*) height_file       ! NOTE that you also need
         read(41,*) uwind_file        ! ocean currents
         read(41,*) vwind_file        !
         read(41,*) potT_file         ! This is only meant to
         read(41,*) tair_file         ! be an illustrative
         read(41,*) humid_file        ! example
         read(41,*) rhoa_file         !
         read(41,*) fsw_file
         read(41,*) flw_file
         read(41,*) rain_file
         read(41,*) sst_file
         read(41,*) sss_file
         close (41)
             write (6,*) height_file
             write (6,*) uwind_file
             write (6,*) vwind_file
             write (6,*) potT_file
             write (6,*) tair_file
             write (6,*) humid_file
             write (6,*) rhoa_file
             write (6,*) fsw_file
             write (6,*) flw_file
             write (6,*) rain_file
             write (6,*) sst_file
             write (6,*) sss_file
       endif
       call monthly_dat  ! interpolate monthly to 1x daily

      endif

      time_forc=time

      return
      end
c-----------------------------------------------------------------------
      subroutine monthly_dat
c.. monthly data
      implicit none
      include "ice.h"          ! ice code common blocks
      include "thermw.h"       ! thermodynamics common blocks
      integer i, j, im,ip,ipp

c   Compute index in table
      im  = mod(month+10,12) + 1
      ip  = mod(month,12) + 1
      ipp = mod(month+1,12) + 1

      call interpolate         ! interpolation coefficients

c....... Read and interpolate
c.. Atmo data
      call intp_data(0,im,ip,ipp,height_file,zlvl)
      call intp_data(0,im,ip,ipp, uwind_file,uatm)
      call intp_data(0,im,ip,ipp, vwind_file,vatm)
      call intp_data(0,im,ip,ipp,  potT_file,potT)
      call intp_data(0,im,ip,ipp,  tair_file,Tair)
      call intp_data(0,im,ip,ipp, humid_file,Qa)
      call intp_data(0,im,ip,ipp,  rhoa_file,rhoa)
      call intp_data(0,im,ip,ipp,   fsw_file,fsw)
      call intp_data(0,im,ip,ipp,   flw_file,flw)
      call intp_data(0,im,ip,ipp,  rain_file,snow)
c.. Ocean data
      call intp_data(0,im,ip,ipp,sst_file,sst)
      call intp_data(0,im,ip,ipp,sss_file,sss)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        snow (i,j) = snow   (i,j)/rhos       ! snowfall rate, m/s 
        wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s
       enddo
      enddo

      call complete_getflux_ocn

      return
      end
c-----------------------------------------------------------------------
      subroutine interpolate
c   NOTE this assumes 12 months of equal length
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real tt, t0, t1, t2, t3, secmo, s3

c  make time cyclic
      tt = mod(time,secyr)
 
c   Find neighboring times
      secmo = secyr/12.  
      t1 = (real(month)-0.5)*secmo        ! midpoint, current month
      t0 = t1 - secmo                     !  - 1 month
      t2 = t1 + secmo                     !  + 1 month
      t3 = t1 + 2.*secmo                  !  + 2 months
 
c   Compute coefficients
      s3 = 2.*secmo**3
c     c0intp = - (tt - t1)*(tt - t2)*(tt - t3)/(3.*s3)
      c1intp =   (tt - t0)*(tt - t2)*(tt - t3)/s3
      c2intp = - (tt - t0)*(tt - t1)*(tt - t3)/s3
      c3intp =   (tt - t0)*(tt - t1)*(tt - t2)/(3.*s3)
      c0intp =   1.0 - c1intp - c2intp - c3intp
  
      return
      end
c-----------------------------------------------------------------------
      subroutine intp_data(recd,im,ip,ipp,data_file,field)
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j, nbits, recd,im,ip,ipp, nrec
      character (72) ::  data_file
      real fieldm(ilo:ihi,jlo:jhi)
      real field (ilo:ihi,jlo:jhi)
      real fieldp(ilo:ihi,jlo:jhi)
      real fieldpp(ilo:ihi,jlo:jhi)

      if (my_task.eq.master_task) write (6,*) '  ',data_file

        nbits = 64

      !  monthly data
      call ice_open (49,data_file,nbits)
      ! read data
      nrec = recd+im
      call ice_read(49,nrec,fieldm,'rda8')
      nrec = recd+month
      call ice_read(49,nrec,field,'rda8')
      nrec = recd+ip
      call ice_read(49,nrec,fieldp,'rda8')
      nrec = recd+ipp
      call ice_read(49,nrec,fieldpp,'rda8')
      ! close file
      if (my_task.eq.master_task) close (49)

      ! interpolate
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
         field(i,j) = c0intp*fieldm(i,j) + c1intp*field(i,j) 
     &              + c2intp*fieldp(i,j) + c3intp*fieldpp(i,j)
       enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine complete_getflux_ocn
c.. compute remaining ocean forcing fields
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real cht, tscale

      cht = 4.19e4            ! heat capacity of water (J kg-1 deg-1)
      tscale = 1.e7           ! time scale for fw, s
      tscale = max(tscale,dt)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
c....... Initialize freezing temperature
         Tf     (i,j) = -mu_Tf*sss(i,j)   ! freezing temp, deg C

c....... Compute oceanic heat flux
         fw     (i,j) = -(sst(i,j)-Tf(i,j))*cht*30./tscale
         if (fw(i,j).gt.0) then
           frzmlt(i,j) = -(sst(i,j)-Tf(i,j))*cht*30./dt
         else
           frzmlt(i,j) = 0.
         endif
         fw  (i,j) = min(fw(i,j),0.)      ! heat flux, W/m^2
       enddo
      enddo

      call t2ugrid(gwatx)
      call t2ugrid(gwaty)
c      call t2ugrid(tiltx)
c      call t2ugrid(tilty)

      return
      end
c-----------------------------------------------------------------------
