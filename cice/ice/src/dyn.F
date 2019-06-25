c $Id: dyn.F,v 1.3 1998/04/14 14:59:48 eclare Exp $

c     program dyn.F
c.. Elastic-viscous-plastic sea ice dynamics model 
c.. Computes ice velocity
c
c.. author Elizabeth C. Hunke
c..        Fluid Dynamics Group, Los Alamos National Laboratory
c
c.. See E. C. Hunke and J. K. Dukowicz. An elastic-viscous-plastic model 
c..     for sea ice dynamics. J. Phys. Oceanogr., 1997.
c
c.. Copyright, 1997.  The Regents of the University of California.
c.. This software was produced under a U.S. Government contract 
c.. (W-7405-ENG-36) by Los Alamos National Laboratory, which is operated
c.. by the University of California for the U.S. Department of Energy.  
c.. The U.S. Government is licensed to use, reproduce, and distribute this 
c.. software.  Permission is granted to the public to copy and use this 
c.. software without charge, provided that this Notice and any statement
c.. of authorship are reproduced on all copies.  Neither the Government
c.. nor the University makes any warranty, express or implied, or assumes
c.. any liability or responsibility for the use of this software.

c-----------------------------------------------------------------------
      subroutine viscosity
c.. Computes the rates of strain, and the bulk and shear viscosities 
c.. zeta and eta.  The rates of strain (xi*) and viscosities are
c.. calculated for each of the four triangles in each grid cell 
c.. (north, south, east, west).
      implicit none
#include "ice.h"
      integer i, j
      real prs,pf,cc,dd,xi11n,xi12n,xi22n,xi11e,xi12e,xi22e,
     1 xi11s,xi12s,xi22s,xi11w,xi12w,xi22w,deltan,deltae,deltas,deltaw,
     1 zetamax

c$DOACROSS LOCAL(i,j,prs,pf,cc,dd,xi11n,xi12n,xi22n,xi11e,xi12e,xi22e,
c$&   xi11s,xi12s,xi22s,xi11w,xi12w,xi22w,deltan,deltae,deltas,deltaw,
c$&   zetamax)
      do j=1,jmt
       do i=1,imt
      prs = 0.5*prss(i,j)           ! g/s^2  (P/2, varies with c*H)
      pf = prs/tiny                ! initializes zeta to large values

c.. initialize zeta                                ! g/s
      zetan(i,j) = pf
      zetae(i,j) = pf
      zetas(i,j) = pf
      zetaw(i,j) = pf

c.. rates of strain                                ! 1/s
      cc = 0.5*(u(i,j) + u(i-1,j) - u(i,j-1) - u(i-1,j-1))*dytr(i,j) 
      dd = 0.5*(v(i,j) + v(i,j-1) - v(i-1,j) - v(i-1,j-1))*dxtr(i,j)

      xi11n = (u(i,j) - u(i-1,j))/HTN(i,j)
      xi12n = 0.5*((v(i,j) - v(i-1,j))/HTN(i,j) + cc)
      xi22n = 0.5*(v(i,j) + v(i-1,j) - v(i,j-1) - v(i-1,j-1))
     1 *dytr(i,j)

      xi11e = 0.5*(u(i,j) + u(i,j-1) - u(i-1,j) - u(i-1,j-1))
     1 *dxtr(i,j)
      xi12e = 0.5*((u(i,j) - u(i,j-1))/HTE(i,j) + dd)
      xi22e = (v(i,j) - v(i,j-1))/HTE(i,j)

      xi11s = (u(i,j-1) - u(i-1,j-1))/HTN(i,j-1)
      xi12s = 0.5*((v(i,j-1) - v(i-1,j-1))/HTN(i,j-1) + cc)
      xi22s = xi22n

      xi11w = xi11e
      xi12w = 0.5*((u(i-1,j) - u(i-1,j-1))/HTE(i-1,j) + dd)
      xi22w = (v(i-1,j) - v(i-1,j-1))/HTE(i-1,j)

c.. Delta (in the denominator of zeta, eta)        ! 1/s
      deltan = sqrt( (xi11n**2+xi22n**2)*ecc2p
     1  + 4.*xi12n**2*ecc2 + xi11n*xi22n*ecc2m)
      deltae = sqrt( (xi11e**2+xi22e**2)*ecc2p
     1  + 4.*xi12e**2*ecc2 + xi11e*xi22e*ecc2m)
      deltas = sqrt( (xi11s**2+xi22s**2)*ecc2p
     1  + 4.*xi12s**2*ecc2 + xi11s*xi22s*ecc2m)
      deltaw = sqrt( (xi11w**2+xi22w**2)*ecc2p
     1  + 4.*xi12w**2*ecc2 + xi11w*xi22w*ecc2m)

      deltan = max(1.e-20,deltan)
      deltae = max(1.e-20,deltae)
      deltas = max(1.e-20,deltas)
      deltaw = max(1.e-20,deltaw)

c.. bulk viscosity zeta, bounded by maximum, minimum values
      zetamax = 2.5e8*prss(i,j)  !Hibler pg 819    ! g/s

      zetan(i,j) = prs/deltan      
      zetan(i,j) = min(zetamax,zetan(i,j))                         
      zetan(i,j) = max(zetamin,zetan(i,j))

      zetae(i,j) = prs/deltae      
      zetae(i,j) = min(zetamax,zetae(i,j))                         
      zetae(i,j) = max(zetamin,zetae(i,j))

      zetas(i,j) = prs/deltas      
      zetas(i,j) = min(zetamax,zetas(i,j))                         
      zetas(i,j) = max(zetamin,zetas(i,j))        

      zetaw(i,j) = prs/deltaw      
      zetaw(i,j) = min(zetamax,zetaw(i,j))                         
      zetaw(i,j) = max(zetamin,zetaw(i,j)) 

c.. mask zeta (and therefore eta)
      zetan(i,j) = zetan(i,j)*hm(i,j)
      zetae(i,j) = zetae(i,j)*hm(i,j)
      zetas(i,j) = zetas(i,j)*hm(i,j)
      zetaw(i,j) = zetaw(i,j)*hm(i,j)

c.. shear viscosity eta                            ! g/s
      etan(i,j) = zetan(i,j)*ecc2
      etae(i,j) = zetae(i,j)*ecc2
      etas(i,j) = zetas(i,j)*ecc2
      etaw(i,j) = zetaw(i,j)*ecc2

       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine stressprep
c.. Computes quantities used in the subroutines stress and stepu.
c.. This subroutine is cryptic - I apologize - but it made the code
c.. faster by almost a factor of 2 on a Cray-YMP.  Here we compute
c.. quantities needed in the stress tensor (sigma) and momentum (u)
c.. equations, but which do not change during the subcycling. Many of
c.. these variables are grouped in parties of four, one for each 
c.. triangle (north, east, south, west).
      implicit none
#include "ice.h"
      integer i, j
      real Econst,ey,e2,en,ee,es,ew,zn,ze,zs,zw,c2n,c2e,c2s,c2w,
     1 c3n,c3e,c3s,c3w,d2n,d2e,d2s,d2w,a2n,a2e,a2s,a2w,prs
c.. used to compute E = Econst*c*H
      Econst = 2.*eyc*rhoi*xymin*dtei**2           ! g/cm s^2

c$DOACROSS LOCAL(i,j,ey,e2,en,ee,es,ew,zn,ze,zs,zw,c2n,c2e,c2s,c2w,
c$&   c3n,c3e,c3s,c3w,d2n,d2e,d2s,d2w,a2n,a2e,a2s,a2w,prs)
      do j=1,jmt
       do i=1,imt

      if (hm(i,j).gt.tiny) then

      ey = Econst*compact(i,j)*hithick(i,j)        ! E, g/s^2
      ey = max(ey,tiny)

      e2   = 0.5*ey                    ! E/2
      edy(i,j)  = e2*dytr(i,j)         ! E/(2*dy)
      edx(i,j)  = e2*dxtr(i,j)         ! E/(2*dx)
      eHN(i,j)  = e2/HTN(i,j)          ! E/(2*HTN)
      eHE(i,j)  = e2/HTE(i,j)          !   etc
      eHNm(i,j) = e2/HTN(i,j-1)
      eHEm(i,j) = e2/HTE(i-1,j)

      en = e2/etan(i,j)                !     E
      ee = e2/etae(i,j)                !  -------
      es = e2/etas(i,j)                !   2*eta
      ew = e2/etaw(i,j)

      zn = e2/zetan(i,j)               !     E
      ze = e2/zetae(i,j)               !  --------
      zs = e2/zetas(i,j)               !   2*zeta
      zw = e2/zetaw(i,j)

      c2n = dtei + en                  !   1      E
      c2e = dtei + ee                  !  --- + -----
      c2s = dtei + es                  !  dte   2*eta
      c2w = dtei + ew

      c3n = 0.5*(en - zn)              !  E      1     1
      c3e = 0.5*(ee - ze)              !  - * ( --- - ---- )
      c3s = 0.5*(es - zs)              !  4     eta   zeta
      c3w = 0.5*(ew - zw)              

      d2n = c2n - c3n                  !  1      E     E      1     1
      d2e = c2e - c3e                  ! --- + ----- - - * ( --- - ---- )
      d2s = c2s - c3s                  ! dte   2*eta   4     eta   zeta
      d2w = c2w - c3w

      h2n(i,j) = 1./c2n                ! this rapidly gets out of hand
      h2e(i,j) = 1./c2e
      h2s(i,j) = 1./c2s
      h2w(i,j) = 1./c2w

      a2n = h2n(i,j)/(d2n - c3n)
      a2e = h2e(i,j)/(d2e - c3e)
      a2s = h2s(i,j)/(d2s - c3s)
      a2w = h2w(i,j)/(d2w - c3w)

      b2n(i,j) = a2n*d2n
      b2e(i,j) = a2e*d2e
      b2s(i,j) = a2s*d2s
      b2w(i,j) = a2w*d2w

      a2na(i,j) = a2n*c3n
      a2ea(i,j) = a2e*c3e
      a2sa(i,j) = a2s*c3s
      a2wa(i,j) = a2w*c3w

      prs  = 0.5*prss(i,j)
      prssn(i,j) = prs*zn              !    P*E
      prsse(i,j) = prs*ze              !  --------
      prsss(i,j) = prs*zs              !   4*zeta
      prssw(i,j) = prs*zw

      endif

       enddo
      enddo

c$DOACROSS LOCAL(i,j)
      do j=1,jmt
      do i=1,imt
        HTN4(i,j) = 0.25/HTN(i,j)
        HTE4(i,j) = 0.25/HTE(i,j)
        dxt8(i,j) = 0.125*dxtr(i,j)
        dyt8(i,j) = 0.125*dytr(i,j)
      enddo
      enddo

c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
        fm(i,j) = fcor(i,j)*umass(i,j)        ! Coriolis * mass
c.. for water stress
        waterx(i,j) = uvm(i,j)*(gwatx(i,j)*cw - gwaty(i,j)*sw)
        watery(i,j) = uvm(i,j)*(gwaty(i,j)*cw + gwatx(i,j)*sw)
c.. air stress and -mg*gradH_o term (tilt)
c..   calculate tilt from geostrophic currents
        strairx(i,j) = uvm(i,j)*(strairx(i,j) - fm(i,j)*gwaty(i,j))
        strairy(i,j) = uvm(i,j)*(strairy(i,j) + fm(i,j)*gwatx(i,j))
c..   tilt from POP via flux coupler
c**        strairx(i,j) = uvm(i,j)*(strairx(i,j) - g*umass(i,j)*tiltx(i,j))
c**        strairy(i,j) = uvm(i,j)*(strairy(i,j) - g*umass(i,j)*tilty(i,j))
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine stress
c.. Calculates the internal stress components, sigma_ij, in the four
c.. triangles of each cell.
      implicit none
#include "ice.h"
      integer i, j
      real dun,dus,due,duw,dvn,dvs,dve,dvw,cc,dd,
     1 xi11n,xi12n,xi22n,xi11e,xi12e,xi22e,xi11s,xi12s,xi22s,xi11w,
     1 xi12w,xi22w,c4n,c4e,c4s,c4w,c5n,c5e,c5s,c5w

c$DOACROSS LOCAL(i,j,dun,dus,due,duw,dvn,dvs,dve,dvw,cc,dd,
c$&   xi11n,xi12n,xi22n,xi11e,xi12e,xi22e,xi11s,xi12s,xi22s,xi11w,
c$&   xi12w,xi22w,c4n,c4e,c4s,c4w,c5n,c5e,c5s,c5w)
      do j=1,jmt
       do i=1,imt

c**      if (hm(i,j).gt.tiny) then
      if (icehm(i,j)) then

      dun = u(i,j) - u(i-1,j)
      dus = u(i,j-1) - u(i-1,j-1)
      due = u(i,j) - u(i,j-1)
      duw = u(i-1,j) - u(i-1,j-1)

      dvn = v(i,j) - v(i-1,j)
      dvs = v(i,j-1) - v(i-1,j-1)
      dve = v(i,j) - v(i,j-1)
      dvw = v(i-1,j) - v(i-1,j-1)

      cc = 0.5*edy(i,j)*(due + duw)
      dd = 0.5*edx(i,j)*(dvn + dvs)

c.. NOTE these are rates of strain * E
      xi11n = 2.0*dun*eHN(i,j)
      xi12n = dvn*eHN(i,j) + cc
      xi22n = edy(i,j)*(dve + dvw)

      xi11e = edx(i,j)*(dun + dus)
      xi12e = due*eHE(i,j) + dd
      xi22e = 2.0*dve*eHE(i,j)

      xi11s = 2.0*dus*eHNm(i,j)
      xi12s = dvs*eHNm(i,j) + cc
      xi22s = xi22n

      xi11w = xi11e
      xi12w = duw*eHEm(i,j) + dd
      xi22w = 2.0*dvw*eHEm(i,j)

c.. solve for the three components of sigma in each triangle
      c4n = dtei*sig11n(i,j) + xi11n - prssn(i,j)
      c4e = dtei*sig11e(i,j) + xi11e - prsse(i,j)
      c4s = dtei*sig11s(i,j) + xi11s - prsss(i,j)
      c4w = dtei*sig11w(i,j) + xi11w - prssw(i,j)

      c5n = dtei*sig22n(i,j) + xi22n - prssn(i,j)
      c5e = dtei*sig22e(i,j) + xi22e - prsse(i,j)
      c5s = dtei*sig22s(i,j) + xi22s - prsss(i,j)
      c5w = dtei*sig22w(i,j) + xi22w - prssw(i,j)

      sig11n(i,j) = a2na(i,j)*c5n + c4n*b2n(i,j)        ! g/s^2
      sig11e(i,j) = a2ea(i,j)*c5e + c4e*b2e(i,j)
      sig11s(i,j) = a2sa(i,j)*c5s + c4s*b2s(i,j)
      sig11w(i,j) = a2wa(i,j)*c5w + c4w*b2w(i,j)

      sig22n(i,j) = a2na(i,j)*c4n + c5n*b2n(i,j)
      sig22e(i,j) = a2ea(i,j)*c4e + c5e*b2e(i,j)
      sig22s(i,j) = a2sa(i,j)*c4s + c5s*b2s(i,j)
      sig22w(i,j) = a2wa(i,j)*c4w + c5w*b2w(i,j)

      sig12n(i,j) = h2n(i,j)*(xi12n + dtei*sig12n(i,j))
      sig12e(i,j) = h2e(i,j)*(xi12e + dtei*sig12e(i,j))
      sig12s(i,j) = h2s(i,j)*(xi12s + dtei*sig12s(i,j))
      sig12w(i,j) = h2w(i,j)*(xi12w + dtei*sig12w(i,j))

      endif

       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mass_prss
c.. Computes ice mass, pressure (strength), and ice extent masks
      implicit none
#include "ice.h"
      integer i, j
      real tmass(0:imt+1,0:jmt+1)

c.. total mass of ice and snow, centered in T-cell
c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
c** if you, like Hibler (1979), have only thick ice
c        tmass(i,j) = rhoi*compact(i,j)*hithick(i,j)    ! g/cm^2
c** for two types of ice (thick and thin), plus snow on each type
        tmass(i,j) = rhoi*(compact(i,j)*hithick(i,j)
     1        + (1.-compact(i,j))*hithin(i,j))
     1           + rhos*(compact(i,j)*hsthick(i,j)
     1        + (1.-compact(i,j))*hsthin(i,j))

        tmass(i,j) = tmass(i,j)*hm(i,j)     ! mask
       enddo
      enddo

c.. mass centered at velocity nodes (U-cells)
c$DOACROSS LOCAL(i,j)
      do j=0,jmt
       do i=0,imt
        umass(i,j) = 0.25*(tmass(i,j) + tmass(i+1,j)
     1         + tmass(i,j+1) + tmass(i+1,j+1))        ! g/cm^2
       enddo
      enddo

c.. augmented masks (land + open ocean)
c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
       do i=0,imt+1
        ! ice extent mask (T-cells)
        icehm(i,j) = (hm(i,j).gt. tiny) .and. (tmass(i,j).gt.0.01)
       enddo
      enddo
      do j=0,jmt
       do i=0,imt
        ! ice extent mask (U-cells)
        iceuvm(i,j) = (uvm(i,j).gt.tiny) .and. (umass(i,j).gt.0.01)
       enddo
      enddo
      do j=1,jmt
       do i=1,imt
        ! ice extent mask + points north and east of ice pack (U-cells)
        iceuvmp(i,j) = iceuvm(i-1,j) .or. iceuvm(i,j-1) 
     1          .or. iceuvm(i-1,j-1) .or. iceuvm(i,j)
       enddo
      enddo

c.. pressure P 
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
        prss(i,j) = pst*hithick(i,j)*compact(i,j)
     1                *exp(-cst*(1.-compact(i,j)))     ! g/s^2
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine stepu
c.. Calculation of the surface stresses
c.. Integration of the momentum equation to find velocity (u,v)
      implicit none
#include "ice.h"
      integer i, j
      real uorel,vorel,vrel,taux,tauy,umassdtei,cca,ccb,ab2,
     1 s11,s12,s21,s22,c1,c2
      real s11ns(imt,jmt-1), s11ew(imt,jmt)
      real s22ew(imt-1,jmt), s22ns(imt,jmt)
      real s12ewi(imt-1,jmt), s12ns(imt,jmt)
      real s12nsj(imt,jmt-1), s12ew(imt,jmt)

c.. some cryptic but useful arrays
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
      do i=1,imt
       if (iceuvmp(i,j)) then
        s11ew(i,j) = dxt8(i,j)*(sig11e(i,j) + sig11w(i,j))
        s22ns(i,j) = dyt8(i,j)*(sig22n(i,j) + sig22s(i,j))
        s12ns(i,j) = dyt8(i,j)*(sig12n(i,j) + sig12s(i,j))
        s12ew(i,j) = dxt8(i,j)*(sig12e(i,j) + sig12w(i,j))
       endif
      enddo
      enddo

c$DOACROSS LOCAL(i,j)
      do j=1,jmt
      do i=1,imt-1
       if (iceuvmp(i,j)) then
        s22ew(i,j) = HTE4(i,j)*(sig22e(i,j) + sig22w(i+1,j))
        s12ewi(i,j) = HTE4(i,j)*(sig12e(i,j) + sig12w(i+1,j))
       endif
      enddo
      enddo

c$DOACROSS LOCAL(i,j)
      do j=1,jmt-1
      do i=1,imt
       if (iceuvmp(i,j)) then
        s11ns(i,j) = HTN4(i,j)*(sig11s(i,j+1) + sig11n(i,j))
        s12nsj(i,j) = HTN4(i,j)*(sig12s(i,j+1) + sig12n(i,j))
       endif
      enddo
      enddo

c.. integrate the momentum equation
c$DOACROSS LOCAL(i,j,uorel,vorel,vrel,taux,tauy,umassdtei,cca,ccb,ab2,
c$&   s11,s12,s21,s22,c1,c2)
      do j=1,jmt-1
      do i=1,imt-1

       if (iceuvm(i,j)) then
c.. geostrophic ocean currents relative to ice velocity
        uorel = gwatx(i,j) - u(i,j)
        vorel = gwaty(i,j) - v(i,j)
c.. (magnitude of relative geostrophic ocean current)*rhow*drag
        vrel = dragw*sqrt(uorel**2 + vorel**2)  ! cm/s
c.. ice/ocean stress
        strocnx(i,j) = vrel*waterx(i,j)
        strocny(i,j) = vrel*watery(i,j)

c**       if (uvm(i,j).gt.tiny.and.umass(i,j).gt.0.01) then     
c.. TOTAL surface stress--includes wind stress, ice/ocean stress, tilt
        taux = strairx(i,j) + strocnx(i,j)
        tauy = strairy(i,j) + strocny(i,j)

c.. alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
        umassdtei = umass(i,j)*dtei             ! m/dte, g/cm^2 s
        cca = umassdtei + vrel*cw               ! alpha, g/cm^2 s
        ccb = fm(i,j) + vrel*sw                 ! beta,  g/cm^2 s
        ab2 = cca**2 + ccb**2

c.. more cryptic stuff
        s11 = - s11ns(i,j) + s11ns(i+1,j)
     &   + s11ew(i+1,j+1) + s11ew(i+1,j) - s11ew(i,j+1) - s11ew(i,j)

        s12 = - s12ewi(i,j) + s12ewi(i,j+1)
     &   + s12ns(i+1,j+1) + s12ns(i,j+1) - s12ns(i+1,j) - s12ns(i,j)

        s21 = - s12nsj(i,j) + s12nsj(i+1,j)
     &   + s12ew(i+1,j+1) + s12ew(i+1,j) - s12ew(i,j+1) - s12ew(i,j)

        s22 = - s22ew(i,j) + s22ew(i,j+1)
     &   + s22ns(i+1,j+1) + s22ns(i,j+1) - s22ns(i+1,j) - s22ns(i,j)

c.. finally, the velocity components
        c1 = s11 + s12 + taux + umassdtei*u(i,j)
        c2 = s21 + s22 + tauy + umassdtei*v(i,j)
        u(i,j) = (cca*c1 + ccb*c2)/ab2          ! cm/s
        v(i,j) = (cca*c2 - ccb*c1)/ab2
      else
c.. set velocity and stress to zero on land and (nearly) open water
        u(i,j) = 0.
        v(i,j) = 0.
        strocnx(i,j) = 0.
        strocny(i,j) = 0.
      endif

      enddo
      enddo

      return
      end


