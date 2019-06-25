c $Id: therm.F,v 1.4 1998/04/14 14:59:57 eclare Exp $

c     program therm.F      
c   Semtner 3-layer thermodynamics with implicit surface fluxes
c   .. note .. search for 'stand-alone' if running uncoupled
c-----------------------------------------------------------------------
      subroutine thermo_semtner
      implicit none
#include "ice.h"
c.. Local variables
      integer
     &  i, j, k,        ! generic loop indices
     &  nk,             ! nonintuitive code for the number of ice layers
     &  nkm1            ! nk - 1
      real
     &  aa,             ! constant in the saturation vapor pressure formula
     &  bb,             ! constant in the saturation vapor pressure formula
     &  cheat,          ! heat capacity - stand-alone only
     &  const,          ! temporary constant coefficient
     &  const1,         ! layer thickness
     &  const2,         ! layer thickness
     &  dd,             ! latent heat of sublimation
     &  dhb,            ! amount of bottom ablation or accretion
     &  dhi,            ! amount of melting on upper surface
     &  dhs,            ! amount of snowmelt
     &  Dl,             ! coefficient in the latent heat formula
     &  dmix,           ! depth of the mixed layer - stand-alone only
     &  Ds,             ! coefficient in the sensible heat formula
     &  enk,            ! real(nk)
     &  enkhi,          ! enk/hi
     &  enk2hi,         ! 2*enk/hi
     &  f(0:nkmax),     ! interfacial fluxes
     &  fa,             ! flux leaving ice/snow to atmosphere
     &  fs,             ! flux below top interface
     &  fslls,          ! sum sensible, latent, longwave, shortwave fluxes
     &  hi,             ! ice thickness
     &  hinew,          ! new ice thickness
     &  hs,             ! snow depth
     &  hsoldt,         ! old snow depth
     &  pntrate,        ! factor for penetrating radiation, including albedo
     &  Qcoef,          ! coefficient in surface specific humidity formula
     &  qmax,           ! maximum heat allowed in brine pockets
     &  qstor,          ! heat stored in brine pockets
     &  solpen,         ! fraction of penetrating solar radiation 
     &  t(0:nkmax),     ! internal temperatures of ice layers
     &  tb,             ! ice bottom temperature
     &  tk,             ! surface temperature in Kelvin
     &  tocean,         ! ocean temperature = sst except in stand-alone mode
     &  tp(nkmax),      ! internal ice layer temperature (temporary)
     &  tpkel,          ! tp in Kelvin
     &  tsfc,           ! surface temperature
     &  tsi,            ! snow-ice interface temperature
     &  ttmp,           ! temporary temperature
     &  z(0:nkmax),     ! below-surface depth for updating layer thicknesses
     &  zp(0:nkmax)     ! below-surface depth for updating layer thicknesses

c$DOACROSS LOCAL(aa,bb,cheat,const,const1,const2,dd,dhb,
c$&   dhi,dhs,Dl,dmix,Ds,enk,enkhi,enk2hi,f,fa,fs,
c$&   fslls,hi,hinew,hs,hsoldt,i,j,k,nk,nkm1,pntrate,Qcoef,
c$&   qmax,qstor,solpen,t,tb,tk,tocean,tp,tpkel,
c$&   tsfc,tsi,ttmp,z,zp),
c$&   SHARE(albice,albij,albsnow,albw,dt,fresh,fsensible,
c$&   flatent,flwout,fsw,hiold,hitherm,hmin,hsold,hstherm,kth,nkij,
c$&   qb,qi,qs,qstorij,rhocii,rhocsi,rhoi,rhos,salinity,fhnet,
c$&   sigma,sk,smin,snow,solarIo,tij,tsfcij,toceanij,tf,yk)
      do 130 j = 1,jmt
      do 130 i = 1,imt

      hi = hiold(i,j)
      hs = hsold(i,j)
      nk = nkij(i,j,kth)
      nkm1 = nk-1
      enk = real(nk)
      qstor = qstorij(i,j,kth)
      do 5 k=0,nk
        t(k) = tij(i,j,kth,k)
 5    continue
      tsfc = tsfcij(i,j,kth)
      tocean = sst(i,j)
cstand-alone      tocean = toceanij(i,j,kth)
      tb = tf(i,j)
      dhs = 0.                             ! initialize snowmelt 

c  constants needed for sensible and latent heat fluxes
c      rhoa=1.3e-3                         ! air density, g/cm^3
c      cc = 1.75e-3                        ! bulk xfer coef, P&W (1979)
c      cp = 0.2398                         ! air sp heat, cal/g K 
      Ds = 5.459e-7                        ! rhoa*cc*cp, cal/cm^3 K
      Qcoef = 0.622*6.11/1013.             ! PIPS code/Maykut (1986)
      aa = 21.8746                         ! values from P&W (1979) and
      bb = 265.5                           !   PIPS code
      dd = 677.326                         ! latent heat of sublim., cal/g
      Dl = dd*2.275e-6                     ! rhoa*cc*dd, cal/cm^3 

      qmax = 0.5*qi*(hi-hmin)              ! max heat in brine pockets
c-----------------------------------------------------------------------
c.. original thermodynamics code at a single gridpoint
c***********************************************************************
      if (nk.gt.0) then      ! one or more layers > hmin thick
c***********************************************************************
        enkhi = enk/hi
        enk2hi = 2.*enkhi
c-----------------------------------------------------------------------
        if (hs.lt.tiny) then     ! no snow
c-----------------------------------------------------------------------
c    allow penetrating radiation to be stored up to some maximum value
      solpen = solarIo
      if (qstor.gt.qmax) solpen = 0.
      qstor = qstor + dt*solpen*(1.-albice)*fsw(i,j)

c    determine ice surface temperature
      const = yk*enk2hi
      pntrate = (1.-solpen)*(1.-albice)
      ttmp = t(1)

      call tsurface(i,j,tsfc,pntrate,ttmp,const,tk,fslls,
     1              aa,bb,Qcoef,Dl,Ds)

      tsi = tk - 273.15                                                
      if (tsi.gt.-0.1) tsi = -0.1
      tsfc = tsi

c    determine fluxes
      fa = -fslls + sigma*tk*tk*tk*tk
      f(0) = enk2hi*yk*(t(1)-tsi)
      f(1) = yk*(t(2)-t(1))*enkhi

      if (tsi.lt.-0.1) then
c    allow snowfall on frozen ice surface
        if (snow(i,j).gt.0) hs=1.
        dhi = 0.
        if (nk.gt.1) then
c    update upper ice layer temperature
          tp(1) = t(1) + dt*(f(1)-fa)*enkhi*rhocii
          call brine(hi,qstor,enk,enkhi,nk,f,t,tp)
        endif
      else        
c    if ice is melting, update thickness and temperature
        dhi = dt*(fa-f(0))/qi
        if (qstor.ge.qmax) then
          qstor = qstor + dhi*qi
          if (qstor.lt.0) qstor = 0.
          dhi = 2.*dhi
        endif
        if (nk.gt.1) then
          tp(1) = t(1) + dt*rhocii*(f(1)-f(0))/(hi/enk+dhi)
          call brine(hi,qstor,enk,enkhi,nk,f,t,tp)
        else
          hi = hi+dhi
          enkhi = 1./hi
          enk2hi = 2.*enkhi
          dhi = 0.
        endif
       endif

c.. send fluxes to flux coupler
      albij(i,j,kth) = albice
      tk = tsfc + 273.15
      flwout(i,j,kth) = -sigma*tk*tk*tk*tk ! outgoing longwave 
c-----------------------------------------------------------------------
        else        ! snow (hs.ge.tiny)
c-----------------------------------------------------------------------
      if (hs.gt.smin) then
        const = 2.*sk/hs
      else
        const = sk/(hs+(hi*sk)/(yk*enk*2.))
        t(0) = t(1)
      endif

c.. determine snow surface temperature   ...........
      pntrate = 1.-albsnow
      ttmp = t(0)

      call tsurface(i,j,tsfc,pntrate,ttmp,const,tk,fslls,
     1              aa,bb,Qcoef,Dl,Ds)

      if (tsfc.gt.0) tsfc = 0.

c   find temperature at snow-ice interface
      tsi = (hi*sk*t(0) + enk*hs*yk*t(1))/(hi*sk + enk*hs*yk)

c   compute fluxes
      fa = -fslls + sigma*tk*tk*tk*tk
      fs = const*(t(0)-tsfc)
      f(0) = enk2hi*yk*(t(1)-tsi)
      f(1) = yk*(t(2)-t(1))*enkhi
      
      if (hs.lt.smin) then
        tsi = (hi*sk*tsfc + 2.*enk*hs*yk*t(1))
     1                 /(hi*sk + 2.*enk*hs*yk)
        fs = sk*(tsi-tsfc)/hs
        f(0) = fs
        t(0) = (tsi+tsfc)*0.5
      endif

c    when snow is not melting, update temperature and add snowfall
      if (tsfc.lt.0) then
        if (hs.ge.smin) t(0) = t(0) + dt*(f(0)-fa)*rhocsi/hs
        hs = hs + dt*snow(i,j)
      else
c    update temperature and thickness of melting snow
        hsoldt = hs
        hs = hs + dt*(fa-fs)/qs
        if (hsoldt.ge.smin) t(0) = t(0) + dt*(f(0)-fs)*rhocsi/hs
        if (hs.lt.0.1) hs = 0.
        if (hs.lt.tiny) t(0) = 0.
c    melting snow returned to ocean
        dhs = hs-hsoldt
      endif

      if (t(0).gt.0) t(0) = 0.

      dhi =0.

      if (nk.gt.1) then
c    update temperature in top layer of ice
        tp(1) = t(1) + dt*(f(1)-f(0))*enkhi*rhocii
        if (tp(1).gt.-0.1) tp(1) = -0.1
        call brine(hi,qstor,enk,enkhi,nk,f,t,tp)
      endif

c.. send fluxes to flux coupler
      albij(i,j,kth) = albsnow
      tk = tsfc + 273.15
      flwout(i,j,kth) = -sigma*tk*tk*tk*tk ! outgoing longwave 
c-----------------------------------------------------------------------
        endif       ! one or more ice layers, w/ or w/o snow
c-----------------------------------------------------------------------
c    determine amount of bottom ablation or accretion
      f(nk) = yk*(tb-t(nk))*enk2hi
      dhb = dt*(f(nk)+fw(i,j))/qb
      if (dhb.lt.0. .and. qstor.ge.qmax) then
        qstor = qstor + dhb*qb
        if (qstor.lt.0) qstor = 0.
        dhb = 2.*dhb
      endif

      hinew = hi + dhi + dhb

c        dhs < 0 : melting,    snow -> ocean  (snow surface)
c        dhi < 0 : melting,     ice -> ocean  (upper ice surface)
c        dhb > 0 : freezing,  ocean -> ice    (lower ice surface)
c        dhb < 0 : melting,     ice -> ocean
c    return heat to ocean--heat of fusion or with meltwater (cal/cm^2 s)
c      fhnet(i,j,kth) = abs(dhi*qi + dhb*qb + dhs*qs)/dt 
      fhnet(i,j,kth) = dhb*qb/dt 
c    fresh water flux into ocean  (g/cm^2 s)
      fresh(i,j,kth) = -(rhoi*(dhi+dhb)*(1.-salinity)
     1                 + rhos*dhs) / dt
c    salt flux into ocean (g/cm^2 s)
c      salt(i,j,kth) = -rhoi*(dhi+dhb)*salinity/dt

c    determine depths below surface of interim and final interfaces
      z(0) = 0.
      zp(0) = 0.
      const1 = hinew/enk
      do 200 k=1,nk
        z(k)=z(k-1)-const1
 200  continue

      const2 = hi/enk
      zp(1) = -const2 - dhi
      if (nk.gt.2) then
        do 210 k=2,nkm1
          zp(k)=zp(k-1)-const2
 210    continue
      endif

      if (dhb.gt.0) then
c    bottom accretion
        zp(nk)=zp(nkm1)-const2
        tp(nk) = t(nk) + dt*(f(nk)-f(nkm1))*enkhi*rhocii
        tp(nk+1)=tb

        do 220 k=1,nk
          t(k) = ((z(k-1)-zp(k))*tp(k) + (zp(k)-z(k))*tp(k+1))
     1                *enk/hinew
 220    continue

      else   ! (dhb.le.0)
c    bottom ablation
      zp(nk) = z(nk)  
      tp(nk) = t(nk) + dt*rhocii*(f(nk)-f(nkm1))/(hi/enk+dhb)

      do 230 k=1,nk
        if (z(k).lt.zp(k).and.zp(k).lt.z(k-1)) then
          t(k) = ((z(k-1)-zp(k))*tp(k) + (zp(k)-z(k))*tp(k+1))
     1                *enk/hinew
        elseif (z(k).lt.zp(k-1).and.zp(k-1).lt.z(k-1)) then
          t(k) = ((z(k-1)-zp(k-1))*tp(k-1) + (zp(k-1)-z(k))*tp(k))
     1                *enk/hinew
        else
          t(k) = tp(k)
        endif
 230  continue
      endif   ! (dhb)

      hi = hinew

c***********************************************************************
      elseif (nk.eq.0) then    ! 1 layer < hmin thick
c***********************************************************************
c-----------------------------------------------------------------------
        if (hs.lt.tiny) then    ! no snow
c-----------------------------------------------------------------------
      const = yk/hi
      solpen = solarIo*0.4             ! if nkmax = 0
      if (nkmax.ne.0) solpen = 0.

c.. determine ice surface temperature   ...........
      pntrate = (1.-solpen)*(1.-albice)
      ttmp = tb

      call tsurface(i,j,tsfc,pntrate,ttmp,const,tk,fslls,
     1              aa,bb,Qcoef,Dl,Ds)

      tsi = tk - 273.15                                                
      if (tsi.gt.-0.1) tsi = -0.1
      tsfc = tsi

c    determine fluxes
      fa = -fslls + sigma*tk*tk*tk*tk
      f(0) = yk*(tb-tsi)/hi

      if (tsi.lt.-0.1) then
c    allow snowfall on frozen ice surface
        if (snow(i,j).gt.0) hs=1.
        dhi = 0.
      else
        dhi = dt*(fa-f(0))/qi
      endif

c.. send fluxes to flux coupler
      albij(i,j,kth) = albice
      tk = tsfc + 273.15
      flwout(i,j,kth) = -sigma*tk*tk*tk*tk ! outgoing longwave 
c-------------------------------------------------------------
        else         ! snow  (hs.ge.tiny)
c-------------------------------------------------------------
      const = sk/(hs+(hi*sk/yk))

c.. determine snow surface temperature   ...........
      pntrate = 1.-albsnow
      ttmp = tb

      call tsurface(i,j,tsfc,pntrate,ttmp,const,tk,fslls,
     1              aa,bb,Qcoef,Dl,Ds)

      if (tsfc.gt.0) tsfc = 0.

c   find temperature at snow-ice interface
      tsi = (hi*sk*tsfc + hs*yk*tb) / (hi*sk + hs*yk)

c   compute fluxes
      fa = -fslls + sigma*tk*tk*tk*tk
      fs = sk*(tsi-tsfc)/hs
      f(0) = fs

c    when snow is not melting, update temperature and add snowfall
      if (tsfc.lt.0) then
        hs = hs + dt*snow(i,j)
      else
c    update temperature and thickness of melting snow
        hsoldt = hs
        hs = hs + dt*(fa-fs)/qs
        if (hs.lt.0.1) hs = 0.
c    melting snow returned to ocean
        dhs = hs-hsoldt
      endif

      dhi =0.

c.. send fluxes to flux coupler
      albij(i,j,kth) = albsnow
      tk = tsfc + 273.15
      flwout(i,j,kth) = -sigma*tk*tk*tk*tk ! outgoing longwave 
c--------------------------------------------------------------------
        endif       ! 1 ice layer < hmin thick, w/ or w/o snow
c--------------------------------------------------------------------
c    determine amount of bottom ablation or accretion
      dhb = dt*(f(0)+fw(i,j))/qb
      hinew = hi + dhi + dhb
      hi = hinew

c    return heat to ocean--heat of fusion or with meltwater (cal/cm^2 s)
      fhnet(i,j,kth) = dhb*qb/dt 
c    fresh water flux into ocean (g/cm^2 s)
      fresh(i,j,kth) = -(rhoi*(dhi+dhb)*(1.-salinity)
     1                 + rhos*dhs) / dt
c    salt flux into ocean (g/cm^2 s)
c      salt(i,j,kth) = -rhoi*(dhi+dhb)*salinity/dt


      if (hi.ge.hmin.and.nkmax.gt.0) then
        t(0) = 0.5*(tsfc + tsi)
        t(1) = 0.5*(tsi + tb)
      else                            ! heat/mass exchange with ocean ???
       if (hi.lt.1.) then
        if (hs.gt.0.) then
          hs = hs - (1.-hi)*0.9/0.33
          hi = 1.
          if (hs.lt.0.1) hs = 0.
        endif
        if (hs.lt.0.1) then
c    when ice disappears, predict ocean temperature until freezing again  
          hi = 0.
          tocean = tb
        endif
       endif
      endif
c***********************************************************************
      else       ! no ice  (nk=-1)
c***********************************************************************
      qstor = 0.
      t(0) = 0.
      albij(i,j,kth) = 0.                  ! albedo not included for flux cpl
cstand-alone      albij(i,j,kth) = albw            

cstand-alone:  when no ice is present, apply fluxes to upper ocean layer
c  These constants are used only in stand-alone mode
      aa = 17.2694                         ! values from P&W (1979) and
      bb = 237.3                           ! PIPS code
      dd = 597.5                           ! ltnt heat of vaporiz., cal/g
      Dl = dd*2.275e-6                     ! rhoa*cc*dd, cal/cm^3 
      dmix = 30.e2                         ! depth of the mixed layer, cm
      cheat = 1.                           ! heat capacity, cal/cm^3 C
      const = cheat*dmix/dt
      ttmp = tb
      pntrate = 1.-albw
cstand-alone      call tsurface(i,j,tocean,pntrate,ttmp,const,tk,fslls,
cstand-alone     1              aa,bb,Qcoef,Dl,Ds)

      if (tocean.le.tb) then
         hi=1.                                                             
c         albij(i,j,kth) = albice
      endif

      ! heat of fusion transferred to ocean upon freezing (cal/cm^2 s)
      fhnet(i,j,kth) = hi*qb/dt
      ! fresh water flux into ocean (g/cm^2 s)
      fresh(i,j,kth) = -rhoi*hi/dt
      ! salt flux into ocean (g/cm^2 s)
c      salt(i,j,kth) = -rhoi*hi*salinity/dt

c.. correct fluxes so that no ocean info goes to flux coupler
      fsensible(i,j,kth) = 0.              ! sensible heat flux
      flatent(i,j,kth) = 0.                ! latent heat flux
      flwout(i,j,kth) = 0.                 ! outgoing longwave 
c***********************************************************************
      endif
c***********************************************************************
c..  update variables on grid 
      nkij(i,j,kth) = nk
      qstorij(i,j,kth) = qstor
      do k=0,nk
        tij(i,j,kth,k) = t(k)
      enddo
      tsfcij(i,j,kth) = tsfc
      toceanij(i,j,kth) = tocean

      hitherm(i,j) = hi
      hstherm(i,j) = hs

 130  continue

c..  calculate rates of change of ice and snow thicknesses
c..  f = dh/dt due to thermodynamics 
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
      do i=1,imt
        fice(i,j) = (hitherm(i,j) - hiold(i,j))/dt 
        fsnow(i,j) = (hstherm(i,j) - hsold(i,j))/dt 
      enddo
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine tsurface(i,j,tsfc,pntrate,ttmp,const,tk,fslls,
     1              aa,bb,Qcoef,Dl,Ds)
c.. iteration for ice or snow surface temperature   ...........
c.. returns tsfc, tk, fslls
c.. This routine produces surface fluxes sent to the coupler 
      implicit none
#include "ice.h"
      integer i, j, k
      real tsfc, pntrate, ttmp, const, tk, fslls,
     1 aa, bb, Qcoef, Dl, Ds, tpkel, stp3, denom, arg,
     1 Qsfc, fltnt, fsns, fsl2, dtemp

      do k=1,20
        tpkel = tsfc + 273.15  
        stp3 = sigma*tpkel*tpkel*tpkel
        denom = tpkel - 273.16 + bb
        arg = aa*(tpkel-273.16)/denom
      ! specific humidity at the surface 
        Qsfc = Qcoef*exp(arg)
      ! latent heat flux - to first order
        fltnt = Dl*wind(i,j)*(Qa(i,j)-Qsfc)     ! cal/cm^2 s
      ! sensible heat flux - to first order
        fsns = Ds*wind(i,j)*(tair(i,j)-tpkel)   ! cal/cm^2 s
      ! sum sensible + latent + longwave + shortwave contribution
        fslls = fsns + fltnt + flw(i,j)  
     1   + pntrate*fsw(i,j)     ! cal/cm^2 s
      ! from second order term in Taylor series for sensible and latent heats
        fsl2 = wind(i,j)*(Ds + Dl*arg*bb*Qsfc/denom**2)
      ! dT, temperature change
        dtemp = (fslls - stp3*tpkel + (ttmp-tsfc)*const) 
     1   /(4.*stp3 + const + fsl2)
      ! new temperature
        tk = tpkel + dtemp  
      ! limit change allowed in iteration
        if (abs(dtemp).gt.10.) tk = tpkel + 10.*dtemp/abs(dtemp)
      ! reset
        tpkel = tk
        stp3 = sigma*tpkel*tpkel*tpkel
        tsfc = tpkel - 273.15              ! surface temperature, deg C
      ! jump out
        if (abs(dtemp).le.0.5) goto 666
      enddo       ! end iteration
 666  continue

c.. send fluxes to flux coupler
      fsensible(i,j,kth) = fsns            ! sensible heat flux
      flatent(i,j,kth) = fltnt             ! latent heat flux

      return
      end
c---------------------------------------------------------------------
      subroutine brine(hi,qstor,enk,enkhi,nk,f,t,tp)
c.. a bit of the thermodynamics code that is used repeatedly
      implicit none
#include "ice.h"
      real hi, qstor, enk, enkhi, qneed, 
     1 tp(nkmax), f(0:nkmax), t(0:nkmax), fact
      integer nk, k

        if (qstor.gt.0. .and. tp(1).lt.-0.1) then
      ! use stored heat in brine pockets to keep temperature at -0.1 
      ! until heat is used up
      qneed = -(tp(1)+0.1)*rhoci*hi/enk
      if (qstor.gt.qneed) then
        qstor = qstor - qneed
        tp(1) = -0.1
      else
        fact = qstor/qneed
        tp(1) = -0.1*fact + tp(1)*(1.-fact)
        qstor = 0.
      endif

      else

      if (tp(1).gt.-0.1) tp(1) = -0.1
      if (nk.gt.2) then
      ! update temperatures of internal ice layers
      do 300 k=2,nk-1
        f(k) = yk*(t(k+1)-t(k))*enkhi
        tp(k) = t(k) + dt*(f(k)-f(k-1))*enkhi*rhocii
 300  continue
      endif

      endif
      return
      end
c----------------------------------------------------------------------
