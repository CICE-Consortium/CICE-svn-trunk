c $Id: transp.F,v 2.1 1999/12/10 22:26:53 eclare Exp $

c-----------------------------------------------------------------------
c   ice advection, thermodynamic growth and type conversions
!!! This is currently hardwired for ncat=2 !!!
c-----------------------------------------------------------------------
      subroutine transport
c.. computes the transport equations for one timestep, 
c.. including the thermodynamics
      implicit none
      include "ice.h"                   ! ice code common blocks
      include "thermw.h"                ! thermodynamics common blocks
      integer i, j, k, n
      real fice(ilo:ihi,jlo:jhi)        ! d(hi)/dt from thermodynamics (m/s)
      real fsnow(ilo:ihi,jlo:jhi)       ! d(hs)/dt from thermodynamics (m/s)
      real dqice(ilo:ihi,jlo:jhi,nlyr)  ! dq/dt from thermodynamics 
      real dTsfc(ilo:ihi,jlo:jhi)   ! d(Tsfc)/dt from thermodynamics (deg/s)
      real ctarea, omctarea, ttmp, qtmp(nlyr)
      real hi_eff(imt_local,jmt_local,ncat), 
     1 hs_eff(imt_local,jmt_local,ncat), 
     1 work(imt_local,jmt_local),
     1 cold(imt_local,jmt_local,ncat), trhs(ilo:ihi,jlo:jhi), 
     1 vice(imt_local,jmt_local,ncat),
     1 vsnow(imt_local,jmt_local,ncat),
     1 work0(ilo:ihi,jlo:jhi), hepsilon

      call timer_start(4)           ! times advection and category conversions

c..  predict compactness  ............................

c$OMP PARALLEL DO PRIVATE(i,j,n)
      do n=1,ncat
      do j=1,jmt_local
      do i=1,imt_local
        cold(i,j,n) = compact(i,j,n)
      enddo
      enddo  
      enddo  

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
      do i=1,imt_local
        work(i,j) = compact(i,j,ncat)
      enddo
      enddo  

      call mpdata(work)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
      do i=1,imt_local
        compact(i,j,ncat) = work(i,j)
      enddo
      enddo  

      call timer_stop(4)

c..  predict ice and snow thicknesses .....................

         do n=1,ncat

      ! make thermodynamic prediction for growth rates
      call timer_start(3)
      call thermo_winton(n,fice,fsnow,dqice,dTsfc)
      call timer_stop(3)

      call timer_start(4)

      ! transport ice thickness
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        work(i,j) = cold(i,j,n)*hice(i,j,n)
       enddo
      enddo

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        trhs(i,j) = cold(i,j,n)*fice(i,j)
c**        trhs(i,j) = 0.      ! for testing transport w/o thermo
       enddo
      enddo

      call euler(work,trhs)
      call mpdata(work)
      call euler(work,trhs)

      ! ice volume 
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        vice(i,j,n) = work(i,j)*tarea(i,j)
        hi_eff(i,j,n) = work(i,j)
       enddo
      enddo

      ! transport snow depth
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        work(i,j) = cold(i,j,n)*hsnow(i,j,n)
       enddo
      enddo

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        trhs(i,j) = cold(i,j,n)*fsnow(i,j)
c**        trhs(i,j) = 0.      ! for testing transport w/o thermo
       enddo
      enddo

      call euler(work,trhs)
      call mpdata(work)
      call euler(work,trhs)

      ! snow volume
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        vsnow(i,j,n) = work(i,j)*tarea(i,j)
        hs_eff(i,j,n) = work(i,j)
       enddo
      enddo

      ! transport surface temperature
c  NOTE Tsfc < 0 but we must transport positive quantities...
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        work(i,j) = -cold(i,j,n)*hice(i,j,n)*Tsfc(i,j,n)
       enddo
      enddo

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        trhs(i,j) = -cold(i,j,n)*hice(i,j,n)*dTsfc(i,j)
c**        trhs(i,j) = 0.      ! for testing transport w/o thermo
       enddo
      enddo

      call euler(work,trhs)
      call mpdata(work)
      call euler(work,trhs)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        if (hi_eff(i,j,n).gt.tiny) then
         Tsfc(i,j,n) = -work(i,j)/hi_eff(i,j,n)
        else
         Tsfc(i,j,n) = 0.
        endif
       enddo
      enddo

      ! transport enthalpy
      do k=1,nlyr

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        work(i,j) = cold(i,j,n)*hice(i,j,n)*qice(i,j,n,k)
       enddo
      enddo

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        trhs(i,j) = cold(i,j,n)*hice(i,j,n)*dqice(i,j,k)
c**        trhs(i,j) = 0.      ! for testing transport w/o thermo
       enddo
      enddo

      call euler(work,trhs)
      call mpdata(work)
      call euler(work,trhs)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        if (hi_eff(i,j,n).gt.tiny) then
         qice(i,j,n,k) = work(i,j)/hi_eff(i,j,n)
        else
         qice(i,j,n,k) = 0.
        endif
       enddo
      enddo
      enddo ! enthalpy layers

      call timer_stop(4)

          enddo  ! ice categories

      call timer_start(5)

c.. conversions among categories ...........................
c.. NOTE this is hardwired for 2 ice categories ..

c.. frazil ice growth
      work0 = 0.
      call frazil(work0)  
c.. accumulate/add contents of category 0 to category 1 or 2  
c$OMP PARALLEL DO PRIVATE(i,j,k)
      do j=jlo,jhi
       do i=ilo,ihi

        ! accumulate with frazil ice left from previous timesteps
        work0(i,j) = work0(i,j) + c_max*hice(i,j,0)*tarea(i,j)

        if (work0(i,j).gt.tiny) then

        if (vice(i,j,1).gt.himin*c_min*tarea(i,j)) then
        ! add to category 1 if category 1 already contains ice
c           write (6,*) 'frazil to cat 1',i,j
         do k=1,nlyr
          qice(i,j,1,k) = qice(i,j,1,k)*vice(i,j,1) 
     1                  + qice(i,j,0,k)*work0(i,j)
         enddo
         vice(i,j,1) = vice(i,j,1) + work0(i,j)
         do k=1,nlyr
          qice(i,j,1,k) = qice(i,j,1,k)/vice(i,j,1) 
         enddo

        elseif (vice(i,j,2).gt.himin*c_min*tarea(i,j)) then
        ! add to category 2 if category 2 already contains ice
        ! (increase concentration)
c           write (6,*) 'frazil to cat 2',i,j
         ttmp = vice(i,j,2)
         do k=1,nlyr
          qice(i,j,2,k) = qice(i,j,2,k)*ttmp
     1                  + qice(i,j,0,k)*work0(i,j)
         enddo
         vice(i,j,2) = ttmp + work0(i,j)
         do k=1,nlyr
          qice(i,j,2,k) = qice(i,j,2,k)/vice(i,j,2)
         enddo
         compact(i,j,2) = min(c_max,
     &     compact(i,j,2)*(1.+lateral*work0(i,j)/ttmp))
         compact(i,j,1) = 1. - compact(i,j,2)

        else
        ! if the total accumulated volume is large enough,
        ! create category 1 ice,
        ! otherwise store it in category 0 
         ttmp = work0(i,j)/(c_max*tarea(i,j))   ! m
         if (ttmp.gt.himin) then
c           write (6,*) 'frazil to cat 1, new',i,j
          vice(i,j,1) = work0(i,j) 
          compact(i,j,1) = c_max
          compact(i,j,2) = c_min
          do k=1,nlyr                     ! some error here if frazil
           qice(i,j,1,k) = qice(i,j,0,k)  ! is accumulated and Tf 
          enddo                           ! changes, but small
          Tsfc(i,j,1) = Tsfc(i,j,0)
          compact(i,j,0) = 0.
          hice(i,j,0) = 0.
          call zero(i,j,0)
         else
          compact(i,j,0) = c_max
          hice(i,j,0) = hice(i,j,0) + ttmp
         endif
        endif

        endif
       enddo
      enddo

c$OMP PARALLEL DO PRIVATE(i,j,k,ctarea,omctarea,ttmp)
      do 100 j=jlo,jhi
       do 100 i=ilo,ihi

          if (hm(i,j).gt.tiny) then

c**       if (hi_eff(i,j,2).le.tiny.or.compact(i,j,2).lt.c_min) then
c**           compact(i,j,2) = c_min
       if (hi_eff(i,j,2).le.tiny) then
           compact(i,j,2) = 0.
           hi_eff(i,j,2) = 0.
           hs_eff(i,j,2) = 0.
           vice(i,j,2) = 0.
           vsnow(i,j,2) = 0.
           call zero(i,j,2)
       endif

       compact(i,j,1) = 1.-compact(i,j,2)   ! for now !

       if (compact(i,j,2).ge.c_min.and.compact(i,j,2).lt.c_max) then
        hice(i,j,1) = hi_eff(i,j,1)/compact(i,j,1)
        hice(i,j,2) = hi_eff(i,j,2)/compact(i,j,2)
        hsnow(i,j,1) = hs_eff(i,j,1)/compact(i,j,1)
        hsnow(i,j,2) = hs_eff(i,j,2)/compact(i,j,2)

        if (hice(i,j,2).gt.hzero(1)) then
         if (hice(i,j,1).gt.hzero(1)) then
          compact(i,j,2) = (vice(i,j,2) + vice(i,j,1) 
     1      - hzero(1)*tarea(i,j))
     1      /((hice(i,j,2) - hzero(1))*tarea(i,j))
          ctarea = compact(i,j,2)*tarea(i,j)
          omctarea = tarea(i,j) - ctarea
          ttmp = vice(i,j,1) - omctarea*hzero(1)
          do k=1,nlyr
           qice(i,j,2,k) = qice(i,j,2,k)*vice(i,j,2)  ! the factor of 2
     1                   + qice(i,j,1,k)*ttmp         ! cancels below
          enddo
          vice(i,j,2) = hice(i,j,2)*ctarea
          hice(i,j,1) = hzero(1)
          vice(i,j,1) = omctarea*hzero(1)
          vsnow(i,j,2) = vsnow(i,j,2) + vsnow(i,j,1) 
     1        - omctarea*hsnow(i,j,1)
          vsnow(i,j,1) = hsnow(i,j,1)*omctarea
          hsnow(i,j,2) = vsnow(i,j,2)/ctarea
          hsnow(i,j,1) = vsnow(i,j,1)/omctarea
          do k=1,nlyr
           qice(i,j,2,k) = qice(i,j,2,k)/vice(i,j,2)
          enddo
         endif

        else
c            write (6,*)  'H < h0'
         ttmp = vice(i,j,2) + vice(i,j,1)
         do k=1,nlyr
          if (ttmp.gt.tiny) then
           qtmp(k) = (qice(i,j,1,k)*vice(i,j,1)
     1       + qice(i,j,2,k)*vice(i,j,2))/ttmp
          else
           qtmp(k) = 0.
          endif
         enddo
         vice(i,j,1) = ttmp
         vsnow(i,j,1) = vsnow(i,j,2) + vsnow(i,j,1)
         if (vice(i,j,1).gt.hzero(1)*tarea(i,j)) then
          compact(i,j,2) = compact(i,j,1)
          ctarea = compact(i,j,1)*tarea(i,j)
          omctarea = tarea(i,j) - ctarea
          hice(i,j,2) = hice(i,j,1)
          hice(i,j,1) = vice(i,j,2)/tarea(i,j)
          vice(i,j,2) = hice(i,j,2)*ctarea
          vice(i,j,1) = hice(i,j,1)*omctarea
          hsnow(i,j,2) = hsnow(i,j,1)
          hsnow(i,j,1) = vsnow(i,j,2)/tarea(i,j)
          vsnow(i,j,2) = hsnow(i,j,2)*ctarea
          vsnow(i,j,1) = hsnow(i,j,1)*omctarea
          do k=1,nlyr   
           ttmp = qice(i,j,2,k) 
           qice(i,j,2,k) = qice(i,j,1,k)
           qice(i,j,1,k) = ttmp
          enddo
         else
          hice(i,j,1) = vice(i,j,1)/tarea(i,j)
          hice(i,j,2) = 0.
          hsnow(i,j,1) = vsnow(i,j,1)/tarea(i,j)
          hsnow(i,j,2) = 0.
          compact(i,j,2) = c_min
          vice(i,j,2) = 0.
          vsnow(i,j,2) = 0.
          do k=1,nlyr
           qice(i,j,1,k) = qtmp(k)
          enddo
          call zero(i,j,2)
         endif
        endif
       endif

       if (compact(i,j,2).ge.c_max) then
         ttmp = vice(i,j,2) + vice(i,j,1)
         compact(i,j,2) = c_max
         hice(i,j,2) = ttmp/(c_max*tarea(i,j))
         hice(i,j,1) = 0.
         hsnow(i,j,2) = (vsnow(i,j,2) + vsnow(i,j,1))
     1                       /(c_max*tarea(i,j))
         hsnow(i,j,1) = 0.
         do k=1,nlyr
          qice(i,j,2,k) = (qice(i,j,1,k)*vice(i,j,1)
     1       + qice(i,j,2,k)*vice(i,j,2))/ttmp
         enddo
         if (hice(i,j,2).lt.hzero(1)) then
c            write (6,*)  'c = 1, H < h0'
           hice(i,j,1) = hice(i,j,2)
           hice(i,j,2) = 0.
           compact(i,j,2) = c_min
           hsnow(i,j,1) = hsnow(i,j,2)
           hsnow(i,j,2) = 0.
           do k=1,nlyr
            qice(i,j,1,k) = qice(i,j,2,k)
           enddo
           call zero(i,j,2)
         else
           call zero(i,j,1)
         endif

       elseif (compact(i,j,2).lt.c_min) then
         compact(i,j,2) = c_min
         hice(i,j,2) = 0.
         hsnow(i,j,2) = 0.
         ttmp = vice(i,j,2) + vice(i,j,1)
         if (ttmp.gt.tiny) then
         vice(i,j,1) = ttmp
         vsnow(i,j,1) = vsnow(i,j,2) + vsnow(i,j,1)
         hice(i,j,1) = vice(i,j,1)/tarea(i,j)
         hsnow(i,j,1) = vsnow(i,j,1)/tarea(i,j)
         do k=1,nlyr
          qice(i,j,1,k) = (qice(i,j,1,k)*vice(i,j,1)
     1       + qice(i,j,2,k)*vice(i,j,2))/ttmp
         enddo
         else
         vice(i,j,1) = 0.
         vsnow(i,j,1) = 0.
         hice(i,j,1) = 0.
         hsnow(i,j,1) = 0.
         do k=1,nlyr
          qice(i,j,1,k) = 0.
         enddo
         endif                 

         if (hice(i,j,1).gt.hzero(1)) then
c            write (6,*)  'c = 0, h > h0'
           compact(i,j,2) = c_max
           hice(i,j,2) = vice(i,j,1)/(c_max*tarea(i,j))
           hice(i,j,1) = 0.
           hsnow(i,j,2) = vsnow(i,j,1)/(c_max*tarea(i,j))
           hsnow(i,j,1) = 0.
           do k=1,nlyr
           qice(i,j,2,k) = qice(i,j,1,k)
           enddo
           call zero(i,j,1)
         else
           call zero(i,j,2)
         endif
       endif

       compact(i,j,1) = 1.-compact(i,j,2)     ! for now !

       endif
 100   continue

c..  limit snow on thin ice                   ! This is often a first
c$OMP PARALLEL DO PRIVATE(i,j)                ! indication that 
      do j=jlo,jhi                            ! something is terribly
       do i=ilo,ihi                           ! wrong and the solution
        if (hsnow(i,j,1).gt.10.) then         ! is about to blow up
         hsnow(i,j,1) = 10.
         write (6,*)  'Limiting snow on thin ice',istep,i,j
        endif
       enddo
      enddo

 555   continue

c..  mask
c$OMP PARALLEL DO PRIVATE(i,j,n)
      do n=1,ncat
       do j=1,jmt_local
        do i=1,imt_local
        hice(i,j,n) = hice(i,j,n)*hm(i,j) 
        hsnow(i,j,n) = hsnow(i,j,n)*hm(i,j) 
        compact(i,j,n) = compact(i,j,n)*hm(i,j) 
        Tsfc(i,j,n) = Tsfc(i,j,n)*hm(i,j)
        qice(i,j,n,1) = qice(i,j,n,1)*hm(i,j)
        qice(i,j,n,2) = qice(i,j,n,2)*hm(i,j)
        enddo
       enddo
      enddo

      do n = 0,ncat
        call bound(hice(:,:,n))
        call bound(hsnow(:,:,n))
        call bound(compact(:,:,n))
        call bound(Tsfc(:,:,n))
        call bound(qice(:,:,n,1))
        call bound(qice(:,:,n,2))
      end do

      call ice_fraction        ! total ice area

      call timer_stop(5)

      return
      end
c-----------------------------------------------------------------------
      subroutine euler(phi,trhs)
c.. incorporates the inhomogeneous terms in the transport equations,
c.. using a time-splitting method
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real phi(imt_local,jmt_local), trhs(ilo:ihi,jlo:jhi), dtime

      dtime = 0.5*dt
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        phi(i,j) = phi(i,j) + dtime*trhs(i,j)
        if (phi(i,j).lt.tiny) phi(i,j) = 0.
       enddo
      enddo  
      call bound(phi)

      return
      end
c-----------------------------------------------------------------------
      subroutine mpdata(phi)
c.. Smolarkiewicz, P. K., 1984:  A fully multidimensional positive 
c.. definite advection transport algorithm with small implicit diffusion,
c.. J. Comput. Phys., 54, 325-362.
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j, k, ix, iy, iord
      real dive, divn, ph, xnu, ynu, betax, betay
      real donor, y1, y2, a, h, eps
      real phi(imt_local,jmt_local), phiavg(imt_local,jmt_local),
     1 uee(imt_local,jmt_local), vnn(imt_local,jmt_local),
     1 ua(imt_local,jmt_local), va(imt_local,jmt_local),
     1 fx(imt_local,jmt_local), fy(imt_local,jmt_local)

      data eps/1.e-15/
      donor(y1,y2,a,h)=0.5*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
       uee(i,j)=0.5*(u(i,j)+u(i,j-1))
       vnn(i,j)=0.5*(v(i,j)+v(i-1,j))
      enddo
      enddo
      call bound(uee)
      call bound(vnn)
 
c.. upwind 
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jhi
      do i=1,ihi
       fx(i,j)=donor(phi(i,j),phi(i+1,j),uee(i,j),HTE(i,j))
       fy(i,j)=donor(phi(i,j),phi(i,j+1),vnn(i,j),HTN(i,j))
      enddo
      enddo

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
       phi(i,j)=phi(i,j)-(fx(i,j)-fx(i-1,j)+fy(i,j)-fy(i,j-1))
     1                          /tarea(i,j)
      enddo
      enddo
      call bound(phi)

      iord=3                             ! hard-wired MPDATA parameter 
c      if (iord.eq.0) goto 555           ! Donor Cell

      do 300 k=1,iord

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jhi
      do i=1,ihi
       phiavg(i,j) = 0.25*(phi(i,j)+phi(i+1,j)+phi(i+1,j+1)+phi(i,j+1))
      enddo
      enddo

c.. antidiffusive velocities
c$OMP PARALLEL DO PRIVATE(i,j,dive,divn)
      do 100 j=jlo,jhi
      do 100 i=ilo,ihi
       dive = ((dyt(i+1,j)*(uee(i+1,j)+uee(i,j))*phi(i+1,j)
     1  - dyt(i,j)*(uee(i-1,j)+uee(i,j))*phi(i,j))
     1  /(phi(i+1,j)+phi(i,j)+eps)
     2  + (dxu(i,j)*(vnn(i+1,j)+vnn(i,j))*phiavg(i,j)
     2  - dxu(i,j-1)*(vnn(i+1,j-1)+vnn(i,j-1))*phiavg(i,j-1))
     2  /(phiavg(i,j)+phiavg(i,j-1)+eps))
     3         /(HTE(i,j)*(dxu(i,j)+dxu(i,j-1)))

       divn = ((dxt(i,j+1)*(vnn(i,j+1)+vnn(i,j))*phi(i,j+1)
     1  - dxt(i,j)*(vnn(i,j-1)+vnn(i,j))*phi(i,j))
     1  /(phi(i,j+1)+phi(i,j)+eps)
     2  + (dyu(i,j)*(uee(i,j+1)+uee(i,j))*phiavg(i,j)
     2  - dyu(i-1,j)*(uee(i-1,j)+uee(i-1,j+1))*phiavg(i-1,j))
     2  /(phiavg(i,j)+phiavg(i-1,j)+eps))
     3         /(HTN(i,j)*(dyu(i,j)+dyu(i-1,j)))

       ua(i,j) = abs(uee(i,j))*(phi(i+1,j)-phi(i,j))
     1   /(phi(i+1,j)+phi(i,j)+eps) - dt*uee(i,j)*dive

       va(i,j) = abs(vnn(i,j))*(phi(i,j+1)-phi(i,j))
     1   /(phi(i,j+1)+phi(i,j)+eps) - dt*vnn(i,j)*divn


 100  continue
      call bound_sw(ua)
      call bound_sw(va)

c.. upwind with antidiffusive velocities
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jhi
      do i=1,ihi
       fx(i,j)=donor(phi(i,j),phi(i+1,j),ua(i,j),HTE(i,j))
       fy(i,j)=donor(phi(i,j),phi(i,j+1),va(i,j),HTN(i,j))
      enddo
      enddo

      ix=-1
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
       phi(i,j)=phi(i,j)-(fx(i,j)-fx(i-1,j)+fy(i,j)-fy(i,j-1))
     1                          /tarea(i,j)
        if (phi(i,j).lt.-1.e-12) then 
c          write (6,*)  i,j,phi(i,j)
           ix=i
           iy=j
c           ph=phi(i,j)
        elseif (phi(i,j).lt.0.) then
          phi(i,j) = 0. 
        endif
c      xnu=dt*ua(i,j)/htn(i,j)      ! error checking
c      ynu=dt*va(i,j)/hte(i,j)
cc       if (xnu.ge.0.5) write (6,*)  i,j,'  xnu=',xnu
cc       if (ynu.ge.0.5) write (6,*)  i,j,'  ynu=',ynu
c      betax=(phiavg(i,j)-phiavg(i,j-1))/(phiavg(i,j)+phiavg(i,j-1)+eps)
c      betay=(phiavg(i,j)-phiavg(i-1,j))/(phiavg(i,j)+phiavg(i-1,j)+eps)
cc       if (betax.gt.1.00001) write (6,*)  i,j,'  betax=',betax
cc       if (betay.gt.1.00001) write (6,*)  i,j,'  betay=',betay
      enddo
      enddo
      if (ix.ge.0) then
        write (6,*)  ix,iy,' transport unstable ',istep,k
      endif

      call bound(phi)

c.. update advective velocities
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
       uee(i,j)=ua(i,j)
       vnn(i,j)=va(i,j)
      enddo
      enddo
      call bound(uee)
      call bound(vnn)

 300  continue
 555  continue

      return
      end   
c-----------------------------------------------------------------------
      subroutine merge_state(work1,work2)
c.. merges ice state properties based on category concentrations,
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j, n
      real work1(ilo:ihi,jlo:jhi,ncat),work2(ilo:ihi,jlo:jhi)

c$OMP PARALLEL DO PRIVATE(i,j,n)
      do j=jlo,jhi
       do i=ilo,ihi
        work2(i,j) = 0.
         do n = 1,ncat
          work2(i,j) = work2(i,j) + compact(i,j,n)*work1(i,j,n)
         enddo
       enddo
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine zero(i,j,n)
c.. reinitialize if ice disappears
      implicit none
      include "ice.h"          ! ice code common blocks

      integer i,j              ! grid indices
      integer n                ! category index
      integer k                ! layer index

      Tsfc(i,j,n) = 0.
      do k = 1, nlyr
         qice(i,j,n,k) = 0.
      enddo

      return
      end 
c-----------------------------------------------------------------------
      subroutine ice_fraction
c.. determine total ice fraction (excludes category 0)
      implicit none
      include "ice.h"          ! ice code common blocks
      include "thermw.h"       ! thermodynamics common blocks
      integer i,j              ! grid indices
      integer n                ! category index

c$OMP PARALLEL DO PRIVATE(i,j,n)
      do j = jlo,jhi
      do i = ilo,ihi
        ifrc(i,j) = 0.                              ! ice fraction
        do n=1,ncat
         if (hice(i,j,n).gt.himin)
     &     ifrc(i,j) = ifrc(i,j) + compact(i,j,n)
        enddo
      enddo
      enddo

      return
      end 
c-----------------------------------------------------------------------

