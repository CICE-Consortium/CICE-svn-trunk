c $Id: dyn.F,v 2.1 1999/12/10 22:26:52 eclare Exp $

c-----------------------------------------------------------------------
c.. Elastic-viscous-plastic sea ice dynamics model 
c.. Computes ice velocity
c
c.. author Elizabeth C. Hunke
c..        Fluid Dynamics Group, Los Alamos National Laboratory
c
c.. See E. C. Hunke and J. K. Dukowicz. An elastic-viscous-plastic model 
c..     for sea ice dynamics. J. Phys. Oceanogr., 1997.
c
c.. Copyright, 1999.  The Regents of the University of California.
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
      subroutine evp
c.. elastic-viscous-plastic dynamics driver
      implicit none
      include "ice.h"          ! ice code common blocks
      integer k

      call timer_start(2)

      call dynprep             ! preparation for dynamics

         do k=1,ndte           ! subcycling
      call stress              ! the stress tensor equation
      call stepu               ! total surface stress, momentum equation
      call bound_sw(u)         ! Periodic/Neumann boundary conditions
      call bound_sw(v)         ! Periodic/Neumann boundary conditions
c      call div_u              ! velocity divergence

         enddo

      call timer_stop(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine init_evp
c.. Initialize parameters and variables needed for the evp dynamics
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real dte             ! subcycling timestep for EVP dynamics, s
      real ecc, c2n, ecoef, d2n, dcoef, tdamp2

      ! elastic time step
      dte = dt/real(ndte)        ! s
      dtei = 1./dte              ! 1/s
      tdamp2 = 2.*eyc*dt         ! 2(wave damping time scale T)
      if (my_task.eq.master_task) 
     &  write(6,*) 'dt = ',dt,'  dte = ',dte,'  tdamp = ',0.5*tdamp2
      ! major/minor axis length ratio, squared
      ecc  = 4.
      ecci = 1.                  ! 4./ecc
      eccm = 2.-0.5*ecci
      eccp = 1.+ecci*0.25
      ! constants for stress equation
      Tdtei = dtei*tdamp2
      c2n = Tdtei + ecc
      ecoef = 0.5*(ecc-1.)
      d2n = c2n - ecoef
      ccoef = 1./c2n
      dcoef = ccoef/(Tdtei + 1.)
      bcoef = dcoef*d2n
      acoef = dcoef*ecoef

      ! Coriolis parameter
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
c        fcor(i,j) = 1.46e-4     ! Hibler 1979, Northern Hemisphere; 1/s
        fcor(i,j) = 2.*(7.292e-5)*sin(ULAT(i,j))  ! 1/s
       enddo
      enddo

c.....   velocity  ....................

c$OMP PARALLEL DO PRIVATE(i,j) 
      do j=jlo,jhi
       do i=ilo,ihi
        u(i,j) = u0*uvm(i,j)        ! m/s
        v(i,j) = v0*uvm(i,j)        ! m/s
       enddo
      enddo

      call bound(u)
      call bound(v)
c      call div_u

c.....   stress tensor  .................  ! kg/s^2
c..  3 components (sigma_11, sigma_12=sigma_21, sigma_22)
c..  4 overlapping triangles per grid cell (northeast, southeast, etc)
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
      sig11ne(i,j) = 0.
      sig11se(i,j) = 0.
      sig11sw(i,j) = 0.
      sig11nw(i,j) = 0.

      sig22ne(i,j) = 0.
      sig22se(i,j) = 0.
      sig22sw(i,j) = 0.
      sig22nw(i,j) = 0.

      sig12ne(i,j) = 0.
      sig12se(i,j) = 0.
      sig12sw(i,j) = 0.
      sig12nw(i,j) = 0.
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dynprep
c.. Computes quantities needed in the stress tensor (sigma) 
c.. and momentum (u) equations, but which do not change during 
c.. the thermodynamics/transport time step:
c..   ice mass and ice extent masks,
c..   pressure (strength), and part of the forcing stresses
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j, n
      real tmass(imt_local,jmt_local), chi_tmp, chs_tmp
      real umass(ilo:ihi,jlo:jhi)
      logical tmphm(imt_local,jmt_local)

c$OMP PARALLEL DO PRIVATE(i,j,n,chi_tmp,chs_tmp)
      do j=1,jmt_local
       do i=1,imt_local
c.. total mass of ice and snow, centered in T-cell
        chi_tmp = 0.           ! sum c * hice  over categories
        chs_tmp = 0.           ! sum c * hsnow over categories
        do n=1,ncat
          chi_tmp = chi_tmp + compact(i,j,n)*hice (i,j,n)
          chs_tmp = chs_tmp + compact(i,j,n)*hsnow(i,j,n)
        enddo
        tmass(i,j) = (rhoi*chi_tmp + rhos*chs_tmp)*hm(i,j) ! kg/m^2
       enddo
      enddo

c.. mass centered at velocity nodes (U-cells)
       call to_ugrid(tmass,umass)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        umassdtei(i,j) = umass(i,j)*dtei               ! m/dte, kg/m^2 s
       enddo
      enddo

c.. augmented masks (land + open ocean)
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        ! ice extent mask (T-cells)
        tmphm(i,j) = (hm(i,j).gt. tiny) .and. (tmass(i,j).gt.0.1)
     1         .and. (compact(i,j,ncat).gt.c_min)
       enddo
      enddo
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        ! extend ice extent mask (T-cells) to points around pack
        icehm(i,j) =  
     1   tmphm(i-1,j+1) .or. tmphm(i,j+1) .or. tmphm(i+1,j+1) .or.
     1   tmphm(i-1,j)   .or. tmphm(i,j)   .or. tmphm(i+1,j)   .or.
     1   tmphm(i-1,j-1) .or. tmphm(i,j-1) .or. tmphm(i+1,j-1)
        icehm(i,j) = icehm(i,j) .and. (hm(i,j).eq.1.) ! remask land points
        ! ice extent mask (U-cells)
        iceuvm(i,j) = (uvm(i,j).gt.tiny) .and. (umass(i,j).gt.0.1)
     1         .and. (compact(i,j,ncat).gt.c_min)
       enddo
      enddo

c.. pressure and forcing terms; set sigma=0 for no ice
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
       if (icehm(i,j)) then
        prss(i,j) = compact(i,j,ncat)*hice(i,j,ncat)
     &       * pst*exp(-cst*(1.-compact(i,j,ncat)))  ! pressure P, kg/s^2
        fm(i,j) = fcor(i,j)*umass(i,j)               ! Coriolis * mass
c.. for water stress
        waterx(i,j) = uvm(i,j)*(gwatx(i,j)*cosw - gwaty(i,j)*sinw)
        watery(i,j) = uvm(i,j)*(gwaty(i,j)*cosw + gwatx(i,j)*sinw)
c.. air stress and -mg*gradH_o term (tilt)
c..   calculate tilt from geostrophic currents
        strairx(i,j) = uvm(i,j)*(strairx(i,j) - fm(i,j)*gwaty(i,j))
        strairy(i,j) = uvm(i,j)*(strairy(i,j) + fm(i,j)*gwatx(i,j))
c..   tilt from POP via flux coupler
c**        strairx(i,j) = uvm(i,j)*(strairx(i,j) 
c**     &                           - grav*umass(i,j)*tiltx(i,j))
c**        strairy(i,j) = uvm(i,j)*(strairy(i,j) 
c**     &                           - grav*umass(i,j)*tilty(i,j))
       endif

       if (.not.tmphm(i,j)) then
      sig11ne(i,j) = 0.
      sig11sw(i,j) = 0.
      sig11nw(i,j) = 0.
      sig11se(i,j) = 0.

      sig22ne(i,j) = 0.
      sig22sw(i,j) = 0.
      sig22nw(i,j) = 0.
      sig22se(i,j) = 0.

      sig12ne(i,j) = 0.
      sig12sw(i,j) = 0.
      sig12nw(i,j) = 0.
      sig12se(i,j) = 0.
       endif

       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine stress
c.. Computes the rates of strain and internal stress components for each
c.. of the four triangles in each grid cell (north, south, east, west).
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real dun,due,dus,duw,dvn,dve,dvs,dvw,
     1 xi11ne,xi12ne,xi22ne,xi11sw,xi12sw,xi22sw,
     1 xi11nw,xi12nw,xi22nw,xi11se,xi12se,xi22se,
     1 deltane,deltase,deltasw,deltanw,
     1 c4ne,c4sw,c4nw,c4se,c5ne,c5sw,c5nw,c5se,
     1 prsdne, prsdsw, prsdnw, prsdse,
     1 prs2ne, prs2sw, prs2nw, prs2se
      real rcon
      data rcon / 4.e-9 /

      call timer_start(6)

c$OMP PARALLEL DO PRIVATE(i,j,
c$OMP&   dun,due,dus,duw,dvn,dve,dvs,dvw,
c$OMP&   xi11ne,xi12ne,xi22ne,xi11sw,xi12sw,xi22sw,xi11nw,xi12nw,xi22nw,
c$OMP&   xi11se,xi12se,xi22se,c4ne,c4sw,c4nw,c4se,c5ne,c5sw,c5nw,c5se,
c$OMP&   deltane,deltasw,deltanw,deltase,
c$OMP&   prsdne,prsdsw,prsdnw,prsdse,prs2ne,prs2sw,prs2nw,prs2se)
      do j=jlo,jhi
       do i=ilo,ihi

      if (icehm(i,j)) then

      dun = (u(i,j) - u(i-1,j))*dxtr(i,j)
      dus = (u(i,j-1) - u(i-1,j-1))*dxtr(i,j)
      due = (u(i,j) - u(i,j-1))*dytr(i,j)
      duw = (u(i-1,j) - u(i-1,j-1))*dytr(i,j)

      dvn = (v(i,j) - v(i-1,j))*dxtr(i,j)
      dvs = (v(i,j-1) - v(i-1,j-1))*dxtr(i,j)
      dve = (v(i,j) - v(i,j-1))*dytr(i,j)
      dvw = (v(i-1,j) - v(i-1,j-1))*dytr(i,j)

c.. rates of strain                                ! 1/s
      xi11ne = dun
      xi11sw = dus
      xi11nw = dun
      xi11se = dus

      xi12ne = 0.5*(due + dvn)
      xi12sw = 0.5*(duw + dvs)
      xi12nw = 0.5*(duw + dvn)
      xi12se = 0.5*(due + dvs)

      xi22ne = dve
      xi22sw = dvw
      xi22nw = dvw
      xi22se = dve

c.. Delta (in the denominator of zeta, eta)        ! 1/s
      deltane = sqrt( (xi11ne**2+xi22ne**2)*eccp
     1  + xi12ne**2*ecci + xi11ne*xi22ne*eccm)
      deltasw = sqrt( (xi11sw**2+xi22sw**2)*eccp
     1  + xi12sw**2*ecci + xi11sw*xi22sw*eccm)
      deltanw = sqrt( (xi11nw**2+xi22nw**2)*eccp
     1  + xi12nw**2*ecci + xi11nw*xi22nw*eccm)
      deltase = sqrt( (xi11se**2+xi22se**2)*eccp
     1  + xi12se**2*ecci + xi11se*xi22se*eccm)

c.. replacement pressure P/Delta
      prsdne = prss(i,j)/max(rcon,deltane)
      prsdsw = prss(i,j)/max(rcon,deltasw)
      prsdnw = prss(i,j)/max(rcon,deltanw)
      prsdse = prss(i,j)/max(rcon,deltase)

c.. P/2 
      prs2ne = 0.5*prsdne*deltane
      prs2sw = 0.5*prsdsw*deltasw
      prs2nw = 0.5*prsdnw*deltanw
      prs2se = 0.5*prsdse*deltase

c.. save replacement pressure for principal stress calculation
      prs_sig(i,j) = 2.*prs2ne  ! match triangle with that used in stress.F !
      if (deltane.le.rcon) prs_sig(i,j) = 0.   ! masks zetamax stresses

c.. cryptic stuff
      c4ne =  Tdtei*sig11ne(i,j) + prsdne*xi11ne - prs2ne
      c4sw =  Tdtei*sig11sw(i,j) + prsdsw*xi11sw - prs2sw
      c4nw =  Tdtei*sig11nw(i,j) + prsdnw*xi11nw - prs2nw
      c4se =  Tdtei*sig11se(i,j) + prsdse*xi11se - prs2se

      c5ne =  Tdtei*sig22ne(i,j) + prsdne*xi22ne - prs2ne
      c5sw =  Tdtei*sig22sw(i,j) + prsdsw*xi22sw - prs2sw
      c5nw =  Tdtei*sig22nw(i,j) + prsdnw*xi22nw - prs2nw
      c5se =  Tdtei*sig22se(i,j) + prsdse*xi22se - prs2se

c.. finally
      sig11ne(i,j) = acoef*c5ne + c4ne*bcoef             ! kg/s^2
      sig11sw(i,j) = acoef*c5sw + c4sw*bcoef
      sig11nw(i,j) = acoef*c5nw + c4nw*bcoef
      sig11se(i,j) = acoef*c5se + c4se*bcoef

      sig22ne(i,j) = acoef*c4ne + c5ne*bcoef
      sig22sw(i,j) = acoef*c4sw + c5sw*bcoef
      sig22nw(i,j) = acoef*c4nw + c5nw*bcoef
      sig22se(i,j) = acoef*c4se + c5se*bcoef

      sig12ne(i,j) = (Tdtei*sig12ne(i,j) + prsdne*xi12ne)*ccoef
      sig12sw(i,j) = (Tdtei*sig12sw(i,j) + prsdsw*xi12sw)*ccoef
      sig12nw(i,j) = (Tdtei*sig12nw(i,j) + prsdnw*xi12nw)*ccoef
      sig12se(i,j) = (Tdtei*sig12se(i,j) + prsdse*xi12se)*ccoef

      endif

       enddo
      enddo

      call timer_stop(6)

      return
      end
c-----------------------------------------------------------------------
      subroutine stepu
c.. Calculation of the surface stresses
c.. Integration of the momentum equation to find velocity (u,v)
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real vrel,cca,ccb,ab2,s11,s12,s21,s22,c1,c2
      real s1122(8,imt_local,jmt_local),work(imt_local,jmt_local)
c     s1122(1)=s11n(imt_local,jmt_local), s1122(2)=s11s(imt_local,jmt_local)
c     s1122(3)=s22e(imt_local,jmt_local), s1122(4)=s22w(imt_local,jmt_local)
c     s1122(5)=s12n(imt_local,jmt_local), s1122(6)=s12s(imt_local,jmt_local)
c     s1122(7)=s12e(imt_local,jmt_local), s1122(8)=s12w(imt_local,jmt_local)

      call timer_start(7)

c.. some cryptic but useful arrays
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
      if (icehm(i,j)) then
       s1122(1,i,j) = (sig11nw(i,j) + sig11ne(i,j))*dxtr4(i,j)
       s1122(2,i,j) = (sig11sw(i,j) + sig11se(i,j))*dxtr4(i,j)
       s1122(3,i,j) = (sig22se(i,j) + sig22ne(i,j))*dytr4(i,j)
       s1122(4,i,j) = (sig22sw(i,j) + sig22nw(i,j))*dytr4(i,j)
       s1122(5,i,j) = (sig12nw(i,j) + sig12ne(i,j))*dxtr4(i,j)
       s1122(6,i,j) = (sig12sw(i,j) + sig12se(i,j))*dxtr4(i,j)
       s1122(7,i,j) = (sig12se(i,j) + sig12ne(i,j))*dytr4(i,j)
       s1122(8,i,j) = (sig12sw(i,j) + sig12nw(i,j))*dytr4(i,j)
      else
       s1122(1,i,j) = 0.
       s1122(2,i,j) = 0.
       s1122(3,i,j) = 0.
       s1122(4,i,j) = 0.
       s1122(5,i,j) = 0.
       s1122(6,i,j) = 0.
       s1122(7,i,j) = 0.
       s1122(8,i,j) = 0.
      endif
      enddo
      enddo
      call bound8_ne(s1122)

c.. integrate the momentum equation
c$OMP PARALLEL DO PRIVATE(i,j,vrel,cca,ccb,ab2,
c$OMP&   s11,s12,s21,s22,c1,c2)
      do j=jlo,jhi
      do i=ilo,ihi

       if (iceuvm(i,j)) then
c.. (magnitude of relative geostrophic ocean current)*rhow*drag
        vrel = dragw*sqrt((gwatx(i,j) - u(i,j))**2 + 
     &                    (gwaty(i,j) - v(i,j))**2)  ! m/s
c.. ice/ocean stress
        strocnx(i,j) = vrel*waterx(i,j)
        strocny(i,j) = vrel*watery(i,j)

c.. alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
        cca = umassdtei(i,j) + vrel*cosw            ! alpha, kg/m^2 s
        ccb = fm(i,j) + vrel*sinw                   ! beta,  kg/m^2 s
        ab2 = cca**2 + ccb**2

c.. more cryptic stuff
        s11=s1122(2,i+1,j+1)-s1122(2,i,j+1)+s1122(1,i+1,j)-s1122(1,i,j)
        s22=s1122(4,i+1,j+1)-s1122(4,i+1,j)+s1122(3,i,j+1)-s1122(3,i,j)
        s21=s1122(6,i+1,j+1)-s1122(6,i,j+1)+s1122(5,i+1,j)-s1122(5,i,j)
        s12=s1122(8,i+1,j+1)-s1122(8,i+1,j)+s1122(7,i,j+1)-s1122(7,i,j)

c.. finally, the velocity components
        c1 = s11 + s12 
     &     + strairx(i,j) + strocnx(i,j) + umassdtei(i,j)*u(i,j)
        c2 = s21 + s22 
     &     + strairy(i,j) + strocny(i,j) + umassdtei(i,j)*v(i,j)
        u(i,j) = (cca*c1 + ccb*c2)/ab2              ! m/s
        v(i,j) = (cca*c2 - ccb*c1)/ab2

c.. ocean-ice stress for coupling  ...the sign will be reversed later...
        strocnx(i,j) = strocnx(i,j) - vrel*(u(i,j)*cosw - v(i,j)*sinw)
        strocny(i,j) = strocny(i,j) - vrel*(v(i,j)*cosw + u(i,j)*sinw)

      else
c.. set velocity and stress to zero on land and (nearly) open water
        u(i,j) = 0.
        v(i,j) = 0.
        strocnx(i,j) = 0.
        strocny(i,j) = 0.
      endif

      enddo
      enddo

      call timer_stop(7)

      return
      end

c-------------------------------------------------------------------------
      subroutine principal_stress(rec)
c.. computes principal stresses for comparison with the theoretical 
c.. yield curve

      implicit none
      include "ice.h"          ! ice code common blocks
      include "netcdf.inc"     ! include netcdf
      real sig11(ilo:ihi,jlo:jhi),sig12(ilo:ihi,jlo:jhi),
     &     sig22(ilo:ihi,jlo:jhi),diff
      integer i, j
      integer rec                              ! record number

c.. choose one set of stresses
c.. (averaging throws stresses off the yield curve, 
c..  especially at the corners)
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
        do i=ilo,ihi
          sig11(i,j) = sig11ne(i,j)
          sig12(i,j) = sig12ne(i,j)
          sig22(i,j) = sig22ne(i,j)
        enddo
      enddo

c.. principal stresses
c$OMP PARALLEL DO PRIVATE(i,j,diff)
      do j=jlo,jhi
        do i=ilo,ihi
          diff=sig11(i,j)-sig22(i,j)
          if(prs_sig(i,j).gt.tiny) then
            sig1(i,j)=(0.5*(sig11(i,j)+sig22(i,j)
     &              +sqrt(diff**2+4.*sig12(i,j)**2)))/prs_sig(i,j)
            sig2(i,j)=(0.5*(sig11(i,j)+sig22(i,j)
     &              -sqrt(diff**2+4.*sig12(i,j)**2)))/prs_sig(i,j)
          else 
            sig1(i,j)=1000.
            sig2(i,j)=1000.
          endif
        enddo
      enddo

      return
      end
c-------------------------------------------------------------------------
