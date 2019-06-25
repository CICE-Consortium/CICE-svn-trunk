c $Id: transp.F,v 1.3.2.1 1998/04/14 16:52:26 eclare Exp $

c     program transp.F
c   ice advection, thermodynamic growth and type conversions
c-----------------------------------------------------------------------
      subroutine transport
c.. computes the transport equations for one timestep, 
c.. including the thermodynamics
      implicit none
#include "ice.h"
      integer i, j, k, nk, nkm1
      real oneminusc, ctarea, omctarea, htup, htdown, 
     1 enk, con1, con2
      real work1(0:imt+1,0:jmt+1), work2(0:imt+1,0:jmt+1), 
     1 work1s(0:imt+1,0:jmt+1), work2s(0:imt+1,0:jmt+1), 
     1 cold(0:imt+1,0:jmt+1), trhs(imt,jmt), fsi(imt,jmt),
     1 vithick(0:imt+1,0:jmt+1), vithin(0:imt+1,0:jmt+1),
     1 vsthick(0:imt+1,0:jmt+1), vsthin(0:imt+1,0:jmt+1),
     1 qtmp(0:imt+1,0:jmt+1),
     1 t(0:nkmax), tp(nkmax), z(0:nkmax), zp(0:nkmax)

c..  predict compactness  ............................

c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        cold(i,j) = compact(i,j)
       enddo
      enddo  

      call mpdata(compact)

c..  predict thick ice and snow thicknesses .....................
      kth = 1
c$DOACROSS LOCAL(i,j)
      do j = 0,jmt+1
       do i = 0,imt+1
        hiold(i,j) = hithick(i,j)
        hsold(i,j) = hsthick(i,j)
       enddo
      enddo

      call thermo_semtner

c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        work1(i,j) = cold(i,j)*hithick(i,j)
       enddo
      enddo

c$DOACROSS LOCAL(i,j,qdrown)
      do j=1,jmt
       do i=1,imt
c..        snow to ice conversion rate
        qdrown = max(0.,rhosiw*(rhos*hsthick(i,j)
     1                            - rhoiw*hithick(i,j)))
        fsi(i,j) = gamma*hsthick(i,j) + qdrown

        trhs(i,j) = cold(i,j)*(fice(i,j)+fsi(i,j))
       enddo
      enddo

      call euler(work1,trhs)
      call mpdata(work1)
      call euler(work1,trhs)

c..  volume of thick ice
c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        vithick(i,j) = work1(i,j)*tarea(i,j)
       enddo
      enddo

c..   snow
c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        work1s(i,j) = cold(i,j)*hsthick(i,j)
       enddo
      enddo

c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
        trhs(i,j) = cold(i,j)*(fsnow(i,j)-fsi(i,j))
       enddo
      enddo

      call euler(work1s,trhs)
      call mpdata(work1s)
      call euler(work1s,trhs)

c..  volume of snow on thick ice
c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        vsthick(i,j) = work1s(i,j)*tarea(i,j)
       enddo
      enddo

c..  predict thin ice thickness  ............................
      kth = 0
c$DOACROSS LOCAL(i,j)
      do j = 0,jmt+1
       do i = 0,imt+1
        hiold(i,j) = hithin(i,j)
        hsold(i,j) = hsthin(i,j)
       enddo
      enddo

      call thermo_semtner

c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        work2(i,j) = (1.-cold(i,j))*hithin(i,j)
       enddo
      enddo

c$DOACROSS LOCAL(i,j,qdrown)
      do j=1,jmt
       do i=1,imt
c..        snow to ice conversion rate
        qdrown = max(0.,rhosiw*(rhos*hsthin(i,j)
     1                            - rhoiw*hithin(i,j)))
        fsi(i,j) = gamma*hsthin(i,j) + qdrown

        trhs(i,j) = (1.-cold(i,j))*(fice(i,j)+fsi(i,j))
       enddo
      enddo

      call euler(work2,trhs)
      call mpdata(work2)
      call euler(work2,trhs)

c..  volume of thin ice
c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        vithin(i,j) = work2(i,j)*tarea(i,j)
       enddo
      enddo

c..  snow
c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        work2s(i,j) = (1.-cold(i,j))*hsthin(i,j)
       enddo
      enddo

c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
        trhs(i,j) = (1.-cold(i,j))*(fsnow(i,j)-fsi(i,j))
       enddo
      enddo

      call euler(work2s,trhs)
      call mpdata(work2s)
      call euler(work2s,trhs)

c..  volume of snow on thin ice
c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        vsthin(i,j) = work2s(i,j)*tarea(i,j)
       enddo
      enddo

c.. conversions between thick/thin ice/snow ...........................
c$DOACROSS LOCAL(i,j,oneminusc,ctarea,omctarea)
      do 100 j=1,jmt
       do 100 i=1,imt

          if (hm(i,j).gt.tiny) then

       oneminusc = 1.-compact(i,j)

       if (compact(i,j).ge.tiny.and.compact(i,j).lt.max_comp) then
        hithin(i,j) = work2(i,j)/oneminusc
        hithick(i,j) = work1(i,j)/compact(i,j)
        hsthin(i,j) = work2s(i,j)/oneminusc
        hsthick(i,j) = work1s(i,j)/compact(i,j)

        if (hithick(i,j).gt.hzero) then
         if (hithin(i,j).gt.hzero) then
          compact(i,j) = (vithick(i,j) + vithin(i,j)
     1      - hzero*tarea(i,j))/((hithick(i,j) - hzero)*tarea(i,j))
          ctarea = compact(i,j)*tarea(i,j)
          omctarea  = tarea(i,j) - ctarea
          vithick(i,j) = hithick(i,j)*ctarea
          hithin(i,j) = hzero
          vithin(i,j) = omctarea*hzero
          vsthick(i,j) = vsthick(i,j) + vsthin(i,j) 
     1        - omctarea*hsthin(i,j)
          vsthin(i,j) = hsthin(i,j)*omctarea
          hsthick(i,j) = vsthick(i,j)/ctarea
          hsthin(i,j) = vsthin(i,j)/omctarea
         endif

        else
c            print*, 'H < h0'
         vithin(i,j) = vithick(i,j) + vithin(i,j)
         vsthin(i,j) = vsthick(i,j) + vsthin(i,j)
         if (vithin(i,j).gt.hzero*tarea(i,j)) then
          compact(i,j) = oneminusc
          ctarea = oneminusc*tarea(i,j)
          omctarea = tarea(i,j) - ctarea
          hithick(i,j) = hithin(i,j)
          hithin(i,j) = vithick(i,j)/tarea(i,j)
          vithick(i,j) = hithick(i,j)*ctarea
          vithin(i,j) = hithin(i,j)*omctarea
          hsthick(i,j) = hsthin(i,j)
          hsthin(i,j) = vsthick(i,j)/tarea(i,j)
          vsthick(i,j) = hsthick(i,j)*ctarea
          vsthin(i,j) = hsthin(i,j)*omctarea
         else
          hithin(i,j) = vithin(i,j)/tarea(i,j)
          hithick(i,j) = 0.
          hsthin(i,j) = vsthin(i,j)/tarea(i,j)
          hsthick(i,j) = 0.
          compact(i,j) = 0.
          vithick(i,j) = 0.
          vsthick(i,j) = 0.
          nkij(i,j,1) = -1
         endif
        endif
       endif

       if (compact(i,j).ge.max_comp) then
         compact(i,j) = max_comp
         hithick(i,j) = (vithick(i,j) + vithin(i,j))
     1                      /(max_comp*tarea(i,j))
         hithin(i,j) = 0.
         hsthick(i,j) = (vsthick(i,j) + vsthin(i,j))
     1                       /(max_comp*tarea(i,j))
         hsthin(i,j) = 0.
         if (hithick(i,j).lt.hzero) then
c            print*, 'c = 1, H < h0'
           hithin(i,j) = hithick(i,j)
           hithick(i,j) = 0.
           compact(i,j) = 0.
           hsthin(i,j) = hsthick(i,j)
           hsthick(i,j) = 0.
           nkij(i,j,1) = -1
         else
           nkij(i,j,0) = -1
         endif

       elseif (compact(i,j).lt.tiny) then
         compact(i,j) = 0.
         vithin(i,j) = vithick(i,j) + vithin(i,j)
         vsthin(i,j) = vsthick(i,j) + vsthin(i,j)
         hithin(i,j) = vithin(i,j)/tarea(i,j)
         hithick(i,j) = 0.
         hsthin(i,j) = vsthin(i,j)/tarea(i,j)
         hsthick(i,j) = 0.
         if (hithin(i,j).gt.hzero) then
c            print*, 'c = 0, h > h0'
           compact(i,j) = max_comp
           hithick(i,j) = vithin(i,j)/(max_comp*tarea(i,j))
           hithin(i,j) = 0.
           hsthick(i,j) = vsthin(i,j)/(max_comp*tarea(i,j))
           hsthin(i,j) = 0.
           nkij(i,j,0) = -1
         else
           nkij(i,j,1) = -1
         endif
       endif

       if (qstorij(i,j,1).lt.tiny) qstorij(i,j,1) = 0.
       if (qstorij(i,j,0).lt.tiny) qstorij(i,j,0) = 0.
       endif
 100   continue

c..  limit snow on thin ice                   ! This is often a first
c$DOACROSS LOCAL(i,j)                         ! indication that 
      do j=1,jmt                              ! something is terribly
       do i=1,imt                             ! wrong and the solution
        if (hsthin(i,j).gt.1.e3) then         ! is about to blow up
         hsthin(i,j) = 1.e3
         print*, 'Limiting snow on thin ice'
        endif
       enddo
      enddo

 555   continue

c..  mask
c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        compact(i,j) = compact(i,j)*hm(i,j) 
        hithick(i,j) = hithick(i,j)*hm(i,j) 
        hithin(i,j) = hithin(i,j)*hm(i,j) 
        hsthick(i,j) = hsthick(i,j)*hm(i,j) 
        hsthin(i,j) = hsthin(i,j)*hm(i,j) 
       enddo
      enddo

c..  update information for thermodynamics ............
c$DOACROSS LOCAL(i,j,kth,hice,nk,k,htup,htdown,enk,con1,con2,t,tp)
         do 200 kth=0,1
        do 200 j=1,jmt
       do 200 i=1,imt
      if (kth.eq.0) hice = hithin(i,j)
      if (kth.eq.1) hice = hithick(i,j)
      nk = nkij(i,j,kth)
      do k=0,nk
       t(k) = tij(i,j,kth,k)
      enddo
      if (hice.eq.0.) then
       nk = -1
       tsfcij(i,j,kth) = tf(i,j)
      endif
      z(0) = 0.
      zp(0) = 0.

      if (nk.ge.0) then
c..   test whether to change number of layers  ! mostly Semtner's code 
      htup = real(nk+1)*hmin
      htdown = real(nk)*hmin
c.................................................
      if (hice.ge.htup .and. nk.lt.nkmax) then
c..     increase the number of ice layers
      nk = nk+1
      nkm1 = nk-1
      enk = real(nk)
      tp(1) = t(1)
      tp(nk) = t(nkm1)
      if (nk.gt.2) then
        con1 = hice/real(nkm1)
        con2 = hice/enk
        do k = 1,nk
          z(k) = z(k-1) - con1
          zp(k) = zp(k-1) - con2
        enddo
        do k=2,nkm1
         tp(k) = ((zp(k-1)-z(k-1))*t(k-1) + (z(k-1)-zp(k))*t(k))/con2
        enddo
      endif
      do k=1,nk
        t(k) = tp(k)
      enddo
c.................................................
      elseif (hice.le.htdown) then
c..     decrease the number of ice layers
      nk = nk-1
      nkm1 = nk-1
      enk = real(nk)
      if (nk.gt.0) then
        con1 = hice/real(nk+1)
        con2 = hice/enk
        do k=1,nk
          z(k) = z(k-1) - con1
          zp(k) = zp(k-1) - con2
          tp(k) = ((zp(k-1)-z(k))*t(k) + (z(k)-zp(k))*t(k+1))/con2
          t(k) = tp(k)
        enddo
      endif
      endif        ! changing number of layers  ! end Semtner's code
c.................................................
      elseif (hice.gt.0.) then     ! nk=-1  ==>  nk=0
       nk = 0        
      endif          ! nk

      if (kth.eq.0) hithin(i,j) = hice
      if (kth.eq.1) hithick(i,j) = hice
      nkij(i,j,kth) = nk
      do k=0,nk
       tij(i,j,kth,k) = t(k)
      enddo

 200  continue

      call bound(hithin)
      call bound(hithick)
      call bound(hsthin)
      call bound(hsthick)
      call bound(compact)

      return
      end
c-----------------------------------------------------------------------
      subroutine euler(phi,trhs)
c.. incorporates the inhomogeneous terms in the transport equations,
c.. using a time-splitting method
      implicit none
#include "ice.h"
      integer i, j
      real phi(0:imt+1,0:jmt+1), trhs(imt,jmt), dtime

      dtime = 0.5*dt
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
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
#include "ice.h"
      integer i, j, k, ix, iy, iord
      real dive, divn, ph, xnu, ynu, betax, betay
      real donor, y1, y2, a, h, eps
      real phi(0:imt+1,0:jmt+1), phiavg(0:imt+1,0:jmt+1),
     1 uee(0:imt+1,0:jmt+1), vnn(0:imt+1,0:jmt+1),
     1 ua(0:imt+1,0:jmt+1), va(0:imt+1,0:jmt+1),
     1 fx(0:imt+1,0:jmt+1), fy(0:imt+1,0:jmt+1)

      data eps/1.e-15/
      donor(y1,y2,a,h)=0.5*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

c$DOACROSS LOCAL(i,j)
      do j=1,jmt
      do i=1,imt
       uee(i,j)=0.5*(u(i,j)+u(i,j-1))
       vnn(i,j)=0.5*(v(i,j)+v(i-1,j))
      enddo
      enddo
      call bound(uee)
      call bound(vnn)
 
c.. upwind 
c$DOACROSS LOCAL(i,j)
      do j=0,jmt
      do i=0,imt
       fx(i,j)=donor(phi(i,j),phi(i+1,j),uee(i,j),HTE(i,j))
       fy(i,j)=donor(phi(i,j),phi(i,j+1),vnn(i,j),HTN(i,j))
      enddo
      enddo

c$DOACROSS LOCAL(i,j)
      do j=1,jmt
      do i=1,imt
       phi(i,j)=phi(i,j)-(fx(i,j)-fx(i-1,j)+fy(i,j)-fy(i,j-1))
     1                          /tarea(i,j)
      enddo
      enddo
      call bound(phi)

      iord=3                             ! hard-wired MPDATA parameter 
c      if (iord.eq.0) goto 555           ! Donor Cell

      do 300 k=1,iord

c$DOACROSS LOCAL(i,j)
      do j=0,jmt
      do i=0,imt
       phiavg(i,j) = 0.25*(phi(i,j)+phi(i+1,j)+phi(i+1,j+1)+phi(i,j+1))
      enddo
      enddo

c.. antidiffusive velocities
c$DOACROSS LOCAL(i,j,dive,divn)
      do 100 j=1,jmt
      do 100 i=1,imt
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
      call bound(ua)
      call bound(va)

c.. upwind with antidiffusive velocities
c$DOACROSS LOCAL(i,j)
      do j=0,jmt
      do i=0,imt
       fx(i,j)=donor(phi(i,j),phi(i+1,j),ua(i,j),HTE(i,j))
       fy(i,j)=donor(phi(i,j),phi(i,j+1),va(i,j),HTN(i,j))
      enddo
      enddo

      ix=-1
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
      do i=1,imt
       phi(i,j)=phi(i,j)-(fx(i,j)-fx(i-1,j)+fy(i,j)-fy(i,j-1))
     1                          /tarea(i,j)
        if (phi(i,j).lt.-1.e-12) then 
c          print*, i,j,phi(i,j)
           ix=i
           iy=j
           ph=phi(i,j)
        elseif (phi(i,j).lt.0.) then
          phi(i,j) = 0. 
        endif
      xnu=dt*ua(i,j)/htn(i,j)
      ynu=dt*va(i,j)/hte(i,j)
c       if (xnu.ge.0.5) print*, i,j,'  xnu=',xnu
c       if (ynu.ge.0.5) print*, i,j,'  ynu=',ynu
      betax=(phiavg(i,j)-phiavg(i,j-1))/(phiavg(i,j)+phiavg(i,j-1)+eps)
      betay=(phiavg(i,j)-phiavg(i-1,j))/(phiavg(i,j)+phiavg(i-1,j)+eps)
c       if (betax.gt.1.00001) print*, i,j,'  betax=',betax
c       if (betay.gt.1.00001) print*, i,j,'  betay=',betay
      enddo
      enddo
      if (ix.ge.0) then
        print*, ix,iy,' transport unstable ',istep,k
      endif

      call bound(phi)

c.. update advective velocities
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
      do i=1,imt
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
      subroutine div_u
c.. calculates divergence of U
c.. this routine is provided for your convenience--it is not necessary
      implicit none
#include "ice.h"
      integer i, j
      real factor, UTE, UTW, VTN, VTS

c$DOACROSS LOCAL(i,j,factor,UTE,UTW,VTN,VTS)
      do j=1,jmt
       do i=1,imt

      factor = 0.5*dxtr(i,j)*dytr(i,j)
      UTE =  factor*HTE(i,j)
      UTW = -factor*HTE(i-1,j)
      VTN =  factor*HTN(i,j)
      VTS = -factor*HTN(i,j-1)
      
      divu(i,j) = UTE*(u(i,j)+u(i,j-1)) + UTW*(u(i-1,j) + u(i-1,j-1))
     1         + VTN*(v(i,j)+v(i-1,j)) + VTS*(v(i,j-1) + v(i-1,j-1))

       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
