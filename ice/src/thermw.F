c $Id: thermw.F,v 2.1 2000/01/14 23:07:08 eclare Exp $

c-----------------------------------------------------------------------
c.. thermodynamics (vertical physics) based on M. Winton 3-layer model
c.. See Bitz, C. M. and W. H. Lipscomb, 1999:  "An energy-conserving 
c..       thermodynamic sea ice model for climate study."  J. Geophys. 
c..       Res., 104, 15669 - 15677.
c..     Winton, M., 1999:  "A reformulated three-layer sea ice model."  
c..       Submitted to J. Atmos. Ocean. Technol.  

c.. authors Elizabeth C. Hunke and William Lipscomb
c..         Fluid Dynamics Group, Los Alamos National Laboratory
c-----------------------------------------------------------------------
      subroutine thermo_winton(n,fice,fsnow,dqice,dTsfc)
c.. Compute temperature change using Winton model with 2 ice layers, of
c.. which only the top layer has a variable heat capacity.

      implicit none
      include "ice.h"          ! ice code common blocks
      include "thermw.h"       ! thermodynamics common blocks

c.. Arguments
      integer n                     ! category number
      real fice(ilo:ihi,jlo:jhi)    ! d(hi)/dt from thermodynamics (cm/s)
      real fsnow(ilo:ihi,jlo:jhi)   ! d(hs)/dt from thermodynamics (cm/s)
      real dqice(ilo:ihi,jlo:jhi,nlyr)   ! dq/dt from thermodynamics 
      real dTsfc(ilo:ihi,jlo:jhi)   ! d(Tsfc)/dt from thermodynamics (deg/s)
c
c Local variables
c
      integer i,j,k                 ! generic indices
      integer niter,nitermax        ! number of iterations
      data nitermax / 10 /

      real area
      real hi
      real hs

      real frsnow        ! fractional snow cover
      real fswabs        ! net SW down at surface (W m-2)
      real fswpen        ! SW penetrating beneath surface (W m-2)
      real fswlwdn       ! flwdn + SW absorbed at surface (W m-2)
      real fswint        ! SW absorbed in ice (W m-2)
      real fswocn        ! SW passed through ice to ocean (W m-2)

      real flwup         ! upward LW at surface (W m-2)
      real fsh           ! surface downward sensible heat (W m-2)
      real flh           ! surface downward latent heat (W m-2)
      real df0dT         ! deriv of f0 wrt Tsf (W m-2 deg-1)
      real f0            ! net surf flux, w/out conduction (W m-2)
      real fct           ! heat conducted to top surface
      real fcb           ! heat conducted to bottom surface

      real Tice(nlyr)    ! internal ice temperatures

      real k12, k32      ! thermal conductivity terms
      real a10, b10      ! coefficients in quadratic eqn for T1
      real a1, b1, c1    ! coefficients in quadratic eqn for T1
      real Tsf           ! surface (ice or snow) temperature
      real dTsf          ! change in Tsf
      real Tsf_start     ! old value of Tsf
 
      real etop          ! energy for top melting (J m-2)
      real ebot          ! energy for bottom melting (J m-2)
      real hnew(nlyr)    ! new ice layer thickness (m)
      real hlyr          ! individual ice layer thickness (m)
      real dhi           ! change in ice thickness
      real dhs           ! change in snow thickness
      real rq            ! rho * q for a layer
      real rqh           ! rho * q * h for a layer
      real qbot          ! q for new ice at bottom surf (J m-3)
      real b2, b3        ! terms in quadratic formula
      real dt2           ! 2*dt
      real dt4           ! 4*dt
      real dt6           ! 6*dt
      real esurp         ! surplus energy from melting (J m-2)
      real qicen (nlyr)  ! ice enthalpy (J m-3)
      real mwater0       ! fresh water mass gained/lost (kg/m^2)
      real msalt0        ! salt gained/lost  (kg/m^2)

      dt2 = 2.*dt
      dt4 = 4.*dt
      dt6 = 6.*dt


c$OMP PARALLEL DO PRIVATE(i,j,k,niter,area,hi,hs,
c$OMP&   fswabs,fswpen,fswlwdn,fswint,flwup,fsh,flh,df0dT,f0,
c$OMP&   fct,fcb,Tice,k12,k32,qicen,fswocn,mwater0,msalt0,
c$OMP&   a10,b10,a1,b1,c1,Tsf,dTsf,Tsf_start,etop,ebot,hnew,
c$OMP&   hlyr,dhi,dhs,rq,rqh,qbot,b2,b3,esurp,frsnow)
c$OMP&   SHARED(dt2,dt4,dt6)
      do 130 j = jlo,jhi
      do 130 i = ilo,ihi

      hi = 0.      ! initialize ice state
      hs = 0.
      area = 0.
      Tsf = 0.
      do k=1,nlyr
       qicen(k) = 0.
      enddo
      flh = 0.     ! initialize fluxes
      fsh = 0.
      flwup = 0.
      fswocn = 0.
      esurp = 0.   ! initialize energies
      etop = 0.
      ebot = 0.
      mwater0 = 0. ! initialize water, salt mass
      msalt0 = 0.
      frsnow = 0.

         if (hm(i,j).gt.tiny.and.hice(i,j,n).ge.himin) then

      hi = hice(i,j,n)
      hs = hsnow(i,j,n)
      area = compact(i,j,n)
      Tsf = Tsfc(i,j,n)
      do k=1,nlyr 
        qicen(k) = qice(i,j,n,k)
      enddo
      qicen(2) = qicen(2) + Lfresh   !!!! note !!!!
      hlyr = hi / rnlyr

      !  mass of fresh water and salt initially present in ice
      mwater0 = rhos*hs + rhoi*hi
      msalt0  = rhoi*hi*salinity

c fractional snow cover
      frsnow = 0.
      if (hs .gt. 0.) frsnow = 1.
c
c Compute SW flux absorbed at surface and penetrating to layer 1.
c
      fswabs  = fsw(i,j) * (1. - albedo(i,j,n))
      fswpen  = fswabs * (1. - frsnow) * i0
      fswocn = fswpen * exp(-ksolar*hi)
      fswint = fswpen - fswocn
      fswlwdn = fswabs - fswpen + flw(i,j)
c
c Compute conductivity terms at layer interfaces.
c
      k12 = kfresh4*ksnow / (ksnow*hi+kfresh4*hs)
      k32 = kfresh2 / hi
c
      a1 = cfresh
      b1 = qicen(1) + (cwater-cfresh)*Tmlt(1) - Lfresh
      c1 = Lfresh * Tmlt(1)
      Tice(1) = 0.5 * (-b1 - sqrt(b1*b1-4.*a1*c1)) / a1
      Tice(2) = (Lfresh-qicen(2)) / cfresh

      if (Tice(1).gt.0.or.Tice(2).gt.0) then
          write (6,*) 'Tice(1) = ',Tice(1),i,j,n
          write (6,*) 'Tice(2) = ',Tice(2),i,j,n
      endif
c
c Compute coefficients used in quadratic formula.
c     rhoi_cfresh = rhoi*cfresh, rhoi_Lfresh = rhoi*Lfresh
c     dt2 = 2*dt, etc.
c
      a10 = rhoi_cfresh*hi/dt2 +
     $     k32 * (dt4*k32 + rhoi_cfresh*hi) 
     $     / (dt6*k32 + rhoi_cfresh*hi) 
      b10 = -hi*
     $     (rhoi_cfresh*Tice(1)+rhoi_Lfresh*Tmlt(1)/Tice(1))/dt2 
     $     - k32 * (dt4*k32*Tf(i,j)+rhoi_cfresh*hi*Tice(2))
     $     / (dt6*k32 + rhoi_cfresh*hi) - fswint
      c1 = rhoi_Lfresh*hi*Tmlt(1) / dt2
c
c Compute new surface and internal temperatures; iterate until 
c Tsfc converges.    
c
      do 100 niter = 1,nitermax
c
c Save temperatures at start of iteration.
c
            Tsf_start = Tsf
c
c Compute top surface flux.
c
            call lwturb_flux (i,j,n,Tsf,flwup,flh,fsh,df0dT)
            f0 = fswlwdn + flwup + fsh + flh
c
c Compute new top layer and surface temperatures.
c If Tsfc is computed to be > 0 C, fix Tsfc = 0 and recompute T1
c with different coefficients. 
c           
         a1 = a10 - k12*df0dT / (k12-df0dT)
         b1 = b10 - k12*(f0-df0dT*Tsf) / (k12-df0dT)
         Tice(1) = -(b1 + sqrt(b1*b1-4.*a1*c1)) / (2.*a1)
         dTsf = (f0 + k12*(Tice(1)-Tsf)) / (k12-df0dT)
         Tsf = Tsf + dTsf
         if (Tsf .gt. 0.) then
            a1 = a10 + k12
            b1 = b10          ! note b1 = b10 - k12*Tf0
            Tice(1) = (-b1 - sqrt(b1*b1-4.*a1*c1)) / (2.*a1)
            Tsf = 0.
            call lwturb_flux (i,j,n,Tsf,flwup,flh,fsh,df0dT)
            f0 = fswlwdn + flwup + fsh + flh
            dTsf = 0.
         endif
c
c Check for convergence.  If no convergence, then repeat.
c               
c
c Convergence test: Make sure Tsfc has converged, within prescribed error.  
c (Energy conservation is guaranteed within machine roundoff, even
c if Tsfc has not converged.)
c If no convergence, then repeat.
c
         if (abs(dTsf).lt.Terrmax) then
            goto 150
         elseif (niter.eq.nitermax) then
            write (6,*) 'thermw conv err ',istep,i,j,dTsf
            write (6,*) 'thermw conv err ',hice(i,j,1),hice(i,j,2),n
            stop
         endif

100   continue  ! surface temperature iteration
150   continue
c
c Compute new bottom layer temperature.
c
         Tice(2) = (dt2*k32*(Tice(1)+2.*Tf(i,j))
     $           + rhoi_cfresh*hi*Tice(2))
     $           / (dt6*k32 + rhoi_cfresh*hi)
c
c Compute final flux values at surfaces.
c
      fct = k12*(Tsf-Tice(1))
      fcb = kfresh4*(Tice(2)-Tf(i,j))/hi
      if (Tsf.lt.0.) f0 = f0 + df0dT*dTsf
c
c Compute new enthalpy for each layer.
c
      qicen(1) = -cwater*Tmlt(1) + cfresh*(Tmlt(1)-Tice(1)) + 
     $              Lfresh*(1.-Tmlt(1)/Tice(1))
      qicen(2) = -cfresh*Tice(2) + Lfresh
c
c Make sure internal ice temperatures do not exceed Tmlt.
c (This should not happen for reasonable values of i0.)
c
c      do k = 1,2
c         if (Tice(k) .ge. Tmlt(k)) then 
c            write (6,*)  'Bug: Ice melt due to T > Tmlt'
c         end if
c      enddo

c......................................................................
c
c.. Compute growth and/or melting at the top and bottom surfaces.......
c......................................................................
c
c Compute enthalpy of new ice growing at bottom surface.
c
      qbot =  -cfresh*Tf(i,j) + Lfresh
c     
c Compute energy available for melting/growth.
c
      if (Tsf .eq. 0.) etop = (f0 - fct) * dt
      ebot = (fcb-fw(i,j)) * dt
c
c Initialize layer thicknesses.
c Make sure internal ice temperatures do not exceed Tmlt.
c If they do, then eliminate the layer.  (Dont think this will happen
c for reasonable values of i0.)
c
      do k = 1, nlyr
         hnew(k) = hlyr
      enddo
c
c Top melt: snow, then ice.
c
      if (etop .gt. 0.) then
         if (hs. gt. 0.) then
            rq =  rhos * qsnow
            rqh = rq * hs
            if (etop .lt. rqh) then
               hs = hs - etop/rq
               etop = 0.
            else
               etop = etop - rqh 
               hs = 0.
            endif
         endif
                  
         do k = 1, nlyr
            if (etop .gt. 0.) then
               rq =  rhoi * qicen(k)
               rqh = rq * hnew(k)
               if (etop .lt. rqh) then
                  hnew(k) = hnew(k) - etop / rq
                  etop = 0.
               else
                  etop = etop - rqh
                  hnew(k) = 0.
               endif
            endif
         enddo
      endif
c
c If ice is gone and melting energy remains
c      
      if (etop .gt. 0.) then
c         write (6,*)  'Ice melts from top, cat', n
         go to 200
      endif
c
c Bottom melt/growth. 
c
      if (ebot .lt. 0.) then
         dhi = -ebot / (qbot * rhoi)
         ebot = 0.
         k = nlyr
         qicen(k) = (hnew(k)*qicen(k)+dhi*qbot) / (hnew(k)+dhi)
         hnew(k) = hnew(k) + dhi
      else               
         do k = nlyr, 1, -1
            if (ebot .gt. 0. .and. hnew(k) .gt. 0.) then
               rq =  rhoi * qicen(k)
               rqh = rq * hnew(k)
               if (ebot .lt. rqh) then
                  hnew(k) = hnew(k) - ebot / rq
                  ebot = 0.
               else
                  ebot = ebot - rqh
                  hnew(k) = 0.
               endif
            endif
         enddo
c
c If ice melts completely and snow is left, remove the snow with 
c energy from the mixed layer
c         
         if (ebot .gt. 0. .and. hs .gt. 0.) then
            rq =  rhos * qsnow
            rqh = rq * hs
            if (ebot .lt. rqh) then
               hs = hs - ebot / rq
               esurp = esurp - rhos*qsnow*hs
               ebot = 0.
c              write (6,*)  'Ice melts from bottom; snow is left:'
c              write (6,*)  'Year', yn, ', Day', dn, ', Category', n
            else
               ebot = ebot - rqh
               hs = 0.
            endif
         endif
                  
         if (ebot .gt. 0.) then
c           write (6,*)  'Ice melts from bottom, day', dn,',cat', n
            go to 200
         endif

      endif
c
c Compute new total ice thickness.
c
      hi = 0.
      do k = 1, nlyr
         hi = hi + hnew(k)
      enddo

c If hi < himin, melt the ice. 
      if (hi .lt. himin) then
         esurp = esurp - rhos*qsnow*hs
         do k = 1, nlyr
            esurp = esurp - rhoi*qicen(k)*hnew(k)
         enddo
         hi = 0.
         hs = 0.
         Tsf=0.
         do k = 1, nlyr
           qicen(k) = 0.
         enddo
         go to 200
      endif

c
c Let it snow
c
        frsnow = 0.  ! temporarily indicates whether fresh snow accumulates
        if (Tsf.lt.0.) then
          frsnow = 1.
          hs = hs + snow(i,j)*dt
        endif

c If there is enough snow to lower the ice/snow interface below 
c freeboard, convert enough snow to ice to bring the interface back 
c to sea-level.  Adjust enthalpy of top ice layer accordingly.
c
      if ( hs .gt. hi*rhoiw/rhos ) then
c               write (6,*)  'Freeboard adjusts, cat', n, ', day', dn
         dhi = (hs * rhos - hi * rhoiw) / rhoml
         dhs = dhi * rhoi / rhos
         rqh = rhoi*qicen(1)*hnew(1) + rhos*qsnow*dhs
         hnew(1) = hnew(1) + dhi
         qicen(1) = rqh / (rhoi*hnew(1))
         hi = hi + dhi
         hs = hs - dhs
      end if

      hlyr = hi/rnlyr
      call new_layers_winton(hs,hi,hlyr,hnew,qicen)

200   continue
c     
c Compute surplus energy left over from melting.
c            
      esurp = esurp + etop + ebot 

       else   ! 0 < hi < himin
c If hi < himin, melt the ice.
         esurp = esurp - rhos*qsnow*hsnow(i,j,n)
         do k = 1, nlyr
            esurp = esurp - rhoi*qicen(k)*hice(i,j,n)*0.5
            qicen(k) = 0.
         enddo
         hi = 0.
         hs = 0.
         Tsf = 0.

       endif  ! hm and hi 

c.. store fluxes 
      flatent  (i,j,n) = area*flh                            ! W/m^2
      fsensible(i,j,n) = area*fsh                            ! W/m^2
      flwout   (i,j,n) = area*flwup                          ! W/m^2
      fswthru  (i,j,n) = area*fswocn                         ! W/m^2
      fresh(i,j,n) = area*((mwater0 - (rhos*hs + rhoi*hi))/dt
     &             + frsnow*snow(i,j))                         ! kg/m^2 s
      fsalt(i,j,n) = area*(msalt0 - rhoi*hi*salinity)*0.001/dt ! kg/m^2 s
      fhnet(i,j,n) = 0.
      if (frzmlt(i,j).le.0.) fhnet(i,j,n) = area*(fw(i,j)+esurp/dt) ! W/m^2

c..  calculate rates of change of ice and snow thicknesses
c..  f = dh/dt due to thermodynamics 
        fice(i,j) = (hi - hice(i,j,n))/dt 
        fsnow(i,j) = (hs - hsnow(i,j,n))/dt 
        dTsfc(i,j) = (Tsf - Tsfc(i,j,n))/dt 
        dqice(i,j,1) = (qicen(1) - qice(i,j,n,1))/dt
        dqice(i,j,2) = (qicen(2) - qice(i,j,n,2) - Lfresh)/dt

130   continue  ! loop over i,j

      return
      end       ! main thermo routine
c --------------------------------------------------------------------
      subroutine lwturb_flux (i,j,n,T,flwup,flh,fsh,df0dT)
c
c Compute net downward LW, SH, and LH flux at surface and derivative
c with respect to temperature for a given category.
c
      implicit none
      include "ice.h"          ! ice code common blocks
      include "thermw.h"       ! thermodynamics common blocks
c
c Input variables
c
      integer i,j,n
      real T               ! Surface temperature
c
c Output variables
c
      real flwup           ! upward LW at surface (W m-2)
      real fsh             ! surface downward sensible heat (W m-2)
      real flh             ! surface downward latent heat (W m-2)
      real df0dT           ! deriv of f0 wrt Tsf (W m-2 deg-1)
c
c Local variables
c
      real Tkel            ! temperature in Kelvin
      real Tkel2           ! Tkel^2
      real Tkel3           ! Tkel^3
      real Tkel4           ! Tkel^4
      real dflwdT          ! deriv of flw wrt Tsf (W m-2 deg-1)      
      real dfshdT          ! deriv of fsh wrt Tsf (W m-2 deg-1)
      real dflhdT          ! deriv of flh wrt Tsf (W m-2 deg-1)
      real Qsfc

      Tkel = T + Tf0kel
      Tkel2 = Tkel * Tkel
      Tkel3 = Tkel2 * Tkel
      Tkel4 = Tkel3 * Tkel
c
c Compute upward longwave.
c
c***  NOTE this is commented out for *net* flw data ***
      flwup = - sigma * Tkel4
c
c Compute downward latent and sensible heat fluxes over water and ice.
c Saturation vapor pressure over ice is given by
c     esat(T) = (6.11 mb) * exp [(Ls/Rv)*(T-Tf0)/(Tf0*T)]
c where Ls is latent heat of sublimation at 0 C
c     Ls/(Rv*Tf0) = 22.47
c
      Qsfc = Qcoef * exp(22.47*(1.-Tf0kel/Tkel))
      flh  = lhcoef(i,j,n)*(Qa(i,j) - Qsfc)
      fsh  = shcoef(i,j,n)*(potT(i,j) - Tkel)
c
c Over ice, compute derivative of fluxes wrt surface temp.
c
      if (n .ge. 1) then
         dflhdT = - lhcoef(i,j,n)*lvrrv*Qsfc/Tkel2
         dfshdT = - shcoef(i,j,n)
         dflwdT = -4. * sigma * Tkel3
         df0dT  = dflwdT + dfshdT + dflhdT
      endif

      return
      end
c -----------------------------------------------------------------------
      subroutine new_layers_winton(hs,hi,hlyr,hnew,qicen)
c
c Repartition into equal-thickness layers, conserving energy.
c This is the 2-layer version.
c
      implicit none
      include "ice.h"    ! ice code common blocks
      include "thermw.h" ! thermodynamics common blocks

      real hnew(nlyr)    ! new ice layer thickness (m)
      real hs
      real hi
      real hlyr          ! individual ice layer thickness (m)
      real qicen (nlyr)  ! ice enthalpy (J m-3)
c
c Local variables
c
      real f1            ! Fraction of upper layer ice in new layer
      real qh1, qh2      ! qice*h for layers 1 and 2
      real qhtot         ! qh1 + qh2
      real q2tmp         ! Temporary value of qice for layer 2

      if (hnew(1).gt.hnew(2)) then  ! Layer 1 gives ice to layer 2 
         f1 = (hnew(1)-hlyr)/hlyr
         q2tmp = f1*qicen(1) + (1.-f1)*qicen(2)
         if (q2tmp.gt.Lfresh) then
            qicen(2) = q2tmp
         else            ! Keep q2 fixed to avoid q2<Lfresh and T2>0
            qh2 = hlyr*qicen(2)
            qhtot = hnew(1)*qicen(1) + hnew(2)*qicen(2)
            qh1 = qhtot - qh2
            qicen(1) = qh1/hlyr
         endif
      else               ! Layer 2 gives ice to layer 1
         f1 = hnew(1)/hlyr
         qicen(1) = f1*qicen(1) + (1.-f1)*qicen(2)
      endif

      return
      end
c --------------------------------------------------------------------
      subroutine sfc_albedo
c.. Compute surface albedo 
      implicit none
      include "ice.h"          ! ice code common blocks
      include "thermw.h"       ! thermodynamics common blocks
      integer n                ! category number
      integer i,j
      real frsnow              ! fractional snow cover

      do n=1,ncat

c$OMP PARALLEL DO PRIVATE(i,j,frsnow)
      do j = jlo,jhi
      do i = ilo,ihi

      frsnow = 0.
      if (hsnow(i,j,n) .gt. 0.) frsnow = 1.

      if (Tsfc(i,j,n) .lt. 0.) then
         albedo(i,j,n) = frsnow*albsnodry + 
     $     (1.-frsnow)*(albicemin + (albicemax - albicemin) 
     $                *(1.-exp(-hice(i,j,n)/halb)))
      else             
         albedo(i,j,n) = frsnow*albsnowet + 
     $     (1.-frsnow)*(albicemin + (albicemax - albicemin) 
     $                *(1.-exp(-hice(i,j,n)/halb)))
      endif

      enddo  
      enddo
      enddo
 
      return
      end
c --------------------------------------------------------------------
      subroutine frazil(vice)
c.. frazil ice growth in leads, polynyas or open ocean
      implicit none
      include "ice.h"          ! ice code common blocks
      include "thermw.h"       ! thermodynamics common blocks
      integer i,j
      real vice(ilo:ihi,jlo:jhi)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        if (frzmlt(i,j).gt.0.) then
         qice(i,j,0,2) = -cfresh*Tf(i,j)
         qice(i,j,0,1) = cfresh*(Tmlt(1) - Tf(i,j)) - cwater*Tmlt(1)
     &               + Lfresh*(1.-Tmlt(1)/Tf(i,j))
         vice(i,j) = tarea(i,j) * 2.*frzmlt(i,j)*dt
     &       /(rhoi*(qice(i,j,0,1) + qice(i,j,0,2) + Lfresh))
         Tsfc(i,j,0) = Tair(i,j) - Tf0kel
        else
         vice(i,j) = 0.
        endif
       enddo
      enddo

      return
      end
c --------------------------------------------------------------------
