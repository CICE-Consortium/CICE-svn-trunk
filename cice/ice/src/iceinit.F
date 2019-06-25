c $Id: iceinit.F,v 1.6 1998/04/14 14:59:54 eclare Exp $

c     program iceinit.F
c   parameter and variable initializations
c-----------------------------------------------------------------------
      subroutine input_data
c.. Namelist variables, set to default values; may be altered
c.. at run time
      implicit none
#include "ice.h"

      namelist /inp/ hice, hsnow, hzero, 
     & u0, v0, cpct, pst, cst, eyc, dt, 
     & ndte, npt, iplot, ktherm, kstress, kread,
     & year
      namelist /fileinp/ mask_file, compact_file, thick_file,
     & vel_file, temps_file, fluxes_file, grid_file,
     & dump_file, restrt_file

c.. default values
      dt=7200.          ! time step dt_zeta, s 
      ndte=100          ! number of subcycles:  ndte=dt/dte
      npt=12            ! total number of time steps (dt) 
      iplot=1           ! how often output is written (10 = once/10 dt)
      cpct=0.9          ! initial ice compactness
      hice=200.         ! initial ice thickness, cm
      hsnow=10.         ! initial snow thickness, cm
      hzero=5.          ! thick/thin discriminating thickness, cm
      u0=0.             ! initial u velocity component, cm/s
      v0=0.             ! initial v velocity component, cm/s
      pst = 2.75e5      ! P* pressure coefficient, dyne/cm^2
      cst = 20.         ! c* pressure parameter
      eyc = 0.25        ! E_o coefficient in formula for E
      kread=0           ! type of read for restart
      ktherm = 0        ! type of thermodynamics data - stand-alone only 
      kstress = 0       ! type of dynamics data - stand-alone only
      year = 0000       ! initial year

      grid_file    = 'grid.d'
      mask_file    = 'ice_mask.year.nc'
      compact_file = 'ice_compact.year.nc'
      thick_file   = 'ice_thick.year.nc'
      vel_file     = 'ice_vel.year.nc'
      temps_file   = 'ice_temps.year.nc'
      fluxes_file  = 'ice_fluxes.year.nc'
      dump_file    = 'ice_restart.year.d'
      restrt_file  = 'ice_restart.year.d'

c.. read from input files
      open  (20,file='ice.inp',status='unknown')       ! parameters
      read  (20,inp)
      close (20)

 10   write (6,900) dt, ndte, npt, iplot
      write (6,920) ktherm, kstress, kread, year
      write (6,910) cpct, hice, hsnow, hzero
      write (6,905) u0, v0, cst, pst, eyc
c**      read (5,inp,err=10)

      open  (21,file='icefiles.inp',status='unknown')  ! file names
      read  (21,fileinp)
      close (21)

 900  format ('  dt=',f7.1,'  ndte=',i3,'  npt=',i6,
     1    '  iplot=',i3)
 920  format ('  ktherm=',i2,'  kstress=',i2,'  kread=',i2,
     1    '  year=',i8)
 910  format ('  cpct=',f5.1,'  hice=',f5.1,'  hsnow=',f5.1,
     1    '  hzero=',f5.1)
 905  format ('  u0=',f5.1,'  v0=',f5.1,'  cst=',f4.1,
     1    '  pst=',1pe8.2,'  Eyc=',1pe8.2)

      return
      end
c     ------------------------------------------------------      
      subroutine init_par
c.. Defines physical constants and numerical parameters
      implicit none
#include "ice.h"
      integer i, j
      real rhow, phiwater, rhocs

      istep0 = 0
c      secyr = 3.1104e7    ! 360-day year
      secyr = 3.1536e7     ! 365-day year
      monthl = 0           ! value of month at the last time step
      tiny = 1.e-11        ! a general-use small number
      max_comp = 0.999     ! maximum ice compactness

c.. elastic time step
      dte = dt/real(ndte)  ! s
      dtei = 1./dte        ! 1/s

c.. Transport parameters
      gamma = 2.e-7        ! snow aging timescale, 1/s

c.. Dynamics and Thermodynamics
      ! densities of ice, snow, seawater 
      rhoi = 0.91          ! Hibler 1979, g/cm^3
      rhos = 0.33          ! Semtner
      rhow = 1.03          ! Pond and Pickard, S=35., T=10.C, p=0

      rhoiw = rhow - rhoi
      rhosiw = rhos/(rhoi*rhow*dt)

      salinity = 0.004     ! salinity of ice, 4 ppt
      albice=0.64          ! albedo of the ice
c**      albice=0.                    ! if already included in shortwave
      albw = .1            ! albedo of the ocean
c**      albw = 0.                    ! if already included in shortwave
      solarIo=0.17         ! fraction of penetrating solar radiation

      ! minimum viscosity
      zetamin = 4.e11                 ! Hibler 1979, pg 823, g/s
      ! major/minor axis length ratio, squared
      ecc2 = 0.25
      ecc2m = 2.*(1.-ecc2)
      ecc2p = 1.+ecc2
      ! water drag coefficient * density  ! Hibler, 1979
      dragw = 0.0055*rhow             ! g/cm^3
      ! water turning angle
c      phiwater = 25.*acos(-1.)/180.  ! 25 deg converted to rad:  pi=acos(-1.)
c      sw = sin(phiwater)
c      cw = cos(phiwater)
      sw = 0.                         ! phiwater = 0 deg
      cw = 1.
      ! Coriolis parameter
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
c        fcor(i,j) = 1.46e-4          ! Hibler 1979, Northern Hemisphere; 1/s
        fcor(i,j) = 2.*(7.292e-5)*sin(ULAT(i,j))  ! 1/s
       enddo
      enddo

c.. Thermodynamics constants
      sk = 0.00074         ! conductivity of snow, cal/cm s C
      yk = .00486          ! conductivity of ice, cal/cm s C
      qs = 26.2            ! heat of fusion of snow, cal/cm^3
      qi = 72.0            ! heat of fusion of surface ice, cal/cm^3
      qb = 72.0            ! heat of fusion of bottom ice, cal/cm^3
      rhocs = .165         ! (density)*(specific heat) of snow, cal/cm^3 C
      rhocsi = 1./rhocs
      rhoci = .45          ! (density)*(specific heat) of ice, cal/cm^3 C
      rhocii = 1./rhoci
      hmin = sqrt(2.*yk*dt*rhocii) + 1.  ! min ice layer thickness, cm
      smin = sqrt(2.*sk*dt*rhocsi) + 1.  ! min snow layer thickness, cm
      sigma = 1.355e-12                  ! cal/s cm^2 K^4

      return
      end
c-----------------------------------------------------------------------
      subroutine init_var
c.. Initialize variables to be sent to coupler
      implicit none
#include "ice.h"
      integer i, j, k, kt
      real s0

c.....  thermodynamics  ...............

c$DOACROSS LOCAL(i,j,kt,k)
      do j = 1,jmt
      do i = 1,imt
        strocnx(i,j) = 0.           ! ice-ocean stress, x-direction (U-cell)
        strocny(i,j) = 0.           ! ice-ocean stress, x-direction (U-cell)
        do kt = 0,1                 ! 0 thin, 1 thick
          tsfcij(i,j,kt) = 0.       ! surface temperature of ice/snow, C
          albij(i,j,kt) = 0.        ! surface albedo
          flatent(i,j,kt) = 0.      ! latent heat flux, cal/s cm^2
          fsensible(i,j,kt) = 0.    ! sensible heat flux, cal/s cm^2
          flwout(i,j,kt) = 0.       ! outgoing longwave radiation, cal/s cm^2
          fresh(i,j,kt) = 0.        ! fresh water flux, g/cm^2 s
          fhnet(i,j,kt) = 0.        ! net heat flux to ocean, cal/s cm^2
        enddo
      enddo
      enddo

c.....  thicknesses and compactness  ...............

c$DOACROSS LOCAL(i,j)
      do j=0,jmt+1
      do i=0,imt+1
       hithick(i,j) = 0. 
       hithin(i,j) = 0.
       hsthick(i,j) = 0.
       hsthin(i,j) = 0.
       compact(i,j) = 0.
      enddo
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine init_more_var
c.. Completes variable initialization on the grid, given fluxes from coupler
      implicit none
#include "ice.h"
      integer i, j, k, kt
      real s0

c$DOACROSS LOCAL(i,j,kt,k)
      do j = 1,jmt
      do i = 1,imt
        do kt = 0,1                        !!! default: no ice 
          nkij(i,j,kt) = -1                ! no ice or snow layers
          tsfcij(i,j,kt) = tf(i,j)         ! ice/snow surface temperature
          qstorij(i,j,kt) = 0.             ! no heat stored in brine pockets
          toceanij(i,j,kt) = sst(i,j)      ! deg C
          do k=0,nkij(i,j,kt)
            tij(i,j,kt,k) = tf(i,j)        ! ice layer temperatures, deg C
          enddo
        enddo
        if (sst(i,j).lt.tf(i,j)) then            !!! cold sfc => ice exists
          tsfcij(i,j,0) = -20.             ! deg C       
          tsfcij(i,j,1) = -20.             ! deg C       
          nkij(i,j,1) = nkmax
          nkij(i,j,0) = 0
          compact(i,j) = cpct*hm(i,j)             
          hithick(i,j) = hice*hm(i,j)      ! cm
          hithin(i,j) = 10.*hm(i,j)        ! cm
          hsthick(i,j) = hsnow*hm(i,j)     ! cm 
          hsthin(i,j) = 1.5*hm(i,j)        ! cm
        endif
      enddo
      enddo

      call bound(hithin)
      call bound(hithick)
      call bound(hsthin)
      call bound(hsthick)
      call bound(compact)

c..  update parameters for thermodynamics ............
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
      do i=1,imt
       if (hithick(i,j).lt.tiny) nkij(i,j,1) = -1
       if (hithin(i,j).lt.tiny) nkij(i,j,0) = -1
       if (compact(i,j).lt.tiny) then
        hithick(i,j) = 0.
        hsthick(i,j) = 0.
        nkij(i,j,1) = -1
       elseif (compact(i,j).gt.max_comp) then
        hithin(i,j) = 0.
        hsthin(i,j) = 0.
        nkij(i,j,0) = -1
       elseif (hithin(i,j).le.hzero) then
        nkij(i,j,0) = 0
       endif
      enddo
      enddo

c.....   velocity  ....................

c$DOACROSS LOCAL(i,j) 
      do j=1,jmt
       do i=1,imt
        u(i,j) = u0*uvm(i,j)        ! cm/s
        v(i,j) = v0*uvm(i,j)        ! cm/s
       enddo
      enddo

      call boundu(u)
      call boundu(v)
c      call div_u

c.....   stress tensor  .................
c.. stress tensor, sigma                 ! g/s^2
c..  3 components (sigma_11, sigma_12=sigma_21, sigma_22)
c..  4 triangles per grid cell (north, east, south, west)
c..  initialized at steady state for grad(u)=0:
c..  sig11 = -P/2 = sig22, sig12 = 0.
c$DOACROSS LOCAL(i,j,s0)
      do j=1,jmt
       do i=1,imt
      s0 = -0.5*pst*compact(i,j)*hithick(i,j)
     1             *exp(-20.*(1.-compact(i,j))) 
      sig11n(i,j) = s0
      sig11e(i,j) = s0
      sig11s(i,j) = s0
      sig11w(i,j) = s0

      sig22n(i,j) = s0
      sig22e(i,j) = s0
      sig22s(i,j) = s0
      sig22w(i,j) = s0

      sig12n(i,j) = 0.
      sig12e(i,j) = 0.
      sig12s(i,j) = 0.
      sig12w(i,j) = 0.
       enddo
      enddo

      return
      end
