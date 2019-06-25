c $Id: ice_dyn_evp.F,v 1.2 2001/02/23 21:06:53 schramm Exp $
c=======================================================================
!---! Elastic-viscous-plastic sea ice dynamics model 
!---! Computes ice velocity
!---!
!---! author Elizabeth C. Hunke
!---!        Fluid Dynamics Group, Los Alamos National Laboratory
!---!
!---! See E. C. Hunke and J. K. Dukowicz. An elastic-viscous-plastic 
!---!     model for sea ice dynamics. J. Phys. Oceanogr., 1997.
!---!
!---! Copyright, 1999.  The Regents of the University of California.
!---! This software was produced under a U.S. Government contract 
!---! (W-7405-ENG-36) by Los Alamos National Laboratory, which is 
!---! operated by the University of California for the U.S. Department 
!---! of Energy.  The U.S. Government is licensed to use, reproduce, and 
!---! distribute this software.  Permission is granted to the public to 
!---! copy and use this software without charge, provided that this 
!---! Notice and any statement of authorship are reproduced on all 
!---! copies.  Neither the Government nor the University makes any 
!---! warranty, express or implied, or assumes any liability or 
!---! responsibility for the use of this software.
c=======================================================================

      module ice_dyn_evp

      use ice_kinds_mod
      use ice_domain
      use ice_grid
      use ice_constants

      implicit none

      integer (kind=int_kind) ::
     &   kdyn         ! type of dynamics ( 1 = evp )
     &,  ndte         ! number of subcycles:  ndte=dt/dte

      logical (kind=log_kind) ::
     &   evp_damping  ! if true, use evp damping procedure

      real (kind=dbl_kind), parameter ::
     &   dragw = 0.00536_dbl_kind * rhow 
                      ! drag coefficient for water on ice *rhow (kg/m^3)
     &,  eyc = 0.36_dbl_kind
                      ! coefficient for calculating the parameter E
     &,  cosw = c1    ! cos(ocean turning angle)  ! turning angle = 0
     &,  sinw = c0    ! sin(ocean turning angle)  ! turning angle = 0

      real (kind=dbl_kind) ::
     &   u0           ! constant coefficient for initial u field (m/s)
     &,  v0           ! constant coefficient for initial v field (m/s)
     &,  ecci         ! 1/e^2 
     &,  dte2T        ! 
     &,  denom1       ! 
     &,  denom2       ! 

      real (kind=dbl_kind), dimension (imt_local,jmt_local) ::
     &   u            ! x-component of velocity (m/s)
     &,  v            ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   shear    ! strain rate II component  (1/s)
     &,  divu     ! strain rate I component, velocity divergence  (1/s)
     &,  Delta

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   waterx       ! for ocean stress calculation, x (m/s)
     &,  watery       ! for ocean stress calculation, y (m/s)
     &,  forcex       ! work array: combined atm stress and ocn tilt, x
     &,  forcey       ! work array: combined atm stress and ocn tilt, y
     &,  umassdtei    ! mass of U-cell/dte (kg/m^2 s)
     &,  fcor         ! Coriolis parameter (1/s)
     &,  fm           ! Coriolis param. * mass in U-cell (kg/s)
     &,  rcon         ! for damping criterion (kg/s)
     &,  prss         ! pressure P (centered in T-cell) (kg/s)
     &,  prs_sig      ! replacement pressure, for stress calc
     &,  strintx        ! divergence of internal ice stress, x (N/m^2)
     &,  strinty        ! divergence of internal ice stress, y (N/m^2)
     &,  sig1         ! principal stress component (diagnostic)
     &,  sig2         ! principal stress component (diagnostic)

      ! ice stress tensor in each corner of T cell (kg/s^2)
      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   stressp (4,ilo:ihi,jlo:jhi) ! sigma11+sigma22
     &,  stressm (4,ilo:ihi,jlo:jhi) ! sigma11-sigma22
     &,  stress12(4,ilo:ihi,jlo:jhi) ! sigma12

c=======================================================================

      contains

c=======================================================================

      subroutine evp(kstrngth)

!---!-------------------------------------------------------------------
!---! elastic-viscous-plastic dynamics driver
!---!-------------------------------------------------------------------

      use ice_timers

      integer (kind=int_kind), intent(in) ::
     &   kstrngth    
      integer (kind=int_kind) :: k

      call ice_timer_start(2)  ! dynamics

      call evp_prep(kstrngth)  ! preparation for dynamics

      do k=1,ndte              ! subcycling
        call stress(k)         ! the stress tensor equation
        call stepu             ! total surface stress, momentum equation
        call bound_sw(u)       ! Periodic/Neumann boundary conditions
        call bound_sw(v)       ! Periodic/Neumann boundary conditions
      enddo

      call evp_finish          ! ice-ocean stress

      call ice_timer_stop(2)   ! dynamics

      end subroutine evp

c=======================================================================

      subroutine init_evp

!---!-------------------------------------------------------------------
!---! Initialize parameters and variables needed for the evp dynamics
!---!-------------------------------------------------------------------

      use ice_calendar

      integer (kind=int_kind) :: i, j, k
      real (kind=dbl_kind) ::  
     &   dte             ! subcycling timestep for EVP dynamics, s
     &,  ecc             ! (ratio of major to minor ellipse axes)^2
     &,  tdamp2          ! 2(wave damping time scale T)

      ! elastic time step
      dte = dt/real(ndte)        ! s
      dtei = c1/dte              ! 1/s
      tdamp2 = c2*eyc*dt         ! s
      if (my_task.eq.master_task) 
     &  write(6,*) 'dt = ',dt,'  dte = ',dte,'  tdamp = ',p5*tdamp2
      ! major/minor axis length ratio, squared
      ecc  = c4
      ecci = p25                  ! 1/ecc
      ! constants for stress equation
      dte2T = dte/tdamp2                    ! unitless
      denom1 = c1/(c1+dte2T)
      denom2 = c1/(c1+dte2T*ecc)
      rcon = 1230._dbl_kind*eyc*dt*dtei**2  ! kg/s  

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        ! Coriolis parameter
c        fcor(i,j) = 1.46e-4_dbl_kind ! Hibler 1979, Northern Hemisphere; 1/s
        fcor(i,j) = c2*(7.292e-5_dbl_kind)*sin(ULAT(i,j))  ! 1/s

        ! velocity
        u(i,j) = c0       ! m/s
        v(i,j) = c0       ! m/s

        ! stress tensor,  kg/s^2
        do k=1,4
         stressp (k,i,j) = c0
         stressm (k,i,j) = c0
         stress12(k,i,j) = c0
        enddo
       enddo
      enddo

      call bound(u)
      call bound(v)

      end subroutine init_evp

c=======================================================================

      subroutine evp_prep(kstrngth)

!---!-------------------------------------------------------------------
!---! Computes quantities needed in the stress tensor (sigma) 
!---! and momentum (u) equations, but which do not change during 
!---! the thermodynamics/transport time step:
!---!   wind stress shift to U grid,
!---!   ice mass and ice extent masks,
!---!   pressure (strength), and part of the forcing stresses
!---! initializes ice velocity for new points to ocean sfc current
!---!-------------------------------------------------------------------

      use ice_flux
      use ice_calendar
      use ice_state
      use ice_mechred_cice

      integer (kind=int_kind), intent(in) ::
     &   kstrngth    
      integer (kind=int_kind) :: i, j, k, n
      real (kind=dbl_kind) :: 
     &   umass(ilo:ihi,jlo:jhi)
     &,  aiu(ilo:ihi,jlo:jhi)   ! ice fraction on u-grid
      logical (kind=log_kind) :: tmphm(imt_local,jmt_local)
      logical (kind=log_kind) :: iceumask_old(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! total mass of ice and snow, centered in T-cell
      !-----------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(i,j,n)
      do j=1,jmt_local
       do i=1,imt_local
        if (tmask(i,j)) then
         tmass(i,j) = (rhoi*vice(i,j) + rhos*vsno(i,j)) ! kg/m^2
        endif
       enddo
      enddo

      !-----------------------------------------------------------------
      ! convert dynamics variables to U grid
      !-----------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        strairx(i,j) = strairxT(i,j) ! prep to convert to U grid
        strairy(i,j) = strairyT(i,j)
       enddo
      enddo
      call t2ugrid(strairx)
      call t2ugrid(strairy)
      call to_ugrid(tmass,umass)
      call to_ugrid(aice,aiu)

      !-----------------------------------------------------------------
      ! convenient variable for evp
      !-----------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        umassdtei(i,j) = umass(i,j)*dtei               ! m/dte, kg/m^2 s
       enddo
      enddo

      !-----------------------------------------------------------------
      ! augmented masks (land + open ocean)
      !-----------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        ! ice extent mask (T-cells)
        tmphm(i,j) = tmask(i,j) .and. (aice(i,j).gt.puny)
       enddo
      enddo
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        ! extend ice extent mask (T-cells) to points around pack
        icetmask(i,j) =  
     &   tmphm(i-1,j+1) .or. tmphm(i,j+1) .or. tmphm(i+1,j+1) .or.
     &   tmphm(i-1,j)   .or. tmphm(i,j)   .or. tmphm(i+1,j)   .or.
     &   tmphm(i-1,j-1) .or. tmphm(i,j-1) .or. tmphm(i+1,j-1)
        icetmask(i,j) = icetmask(i,j) .and. tmask(i,j) ! remask land points
        ! ice extent mask (U-cells)
         iceumask_old(i,j) = iceumask(i,j)  ! save
        iceumask(i,j) = (umask(i,j)) .and. (aiu(i,j).gt.puny)   
       enddo
      enddo

      !-----------------------------------------------------------------
      ! pressure and forcing terms; set sigma=0 for no ice;
      ! initialize ice velocity in cells previously empty to ocn current
      !-----------------------------------------------------------------
      call ice_strength(kstrngth)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
       if (icetmask(i,j)) then
        prss(i,j) = strength(i,j)
        fm(i,j) = fcor(i,j)*umass(i,j)               ! Coriolis * mass

        if (umask(i,j)) then
        ! for ocean stress
        waterx(i,j) = uocn(i,j)*cosw - vocn(i,j)*sinw
        watery(i,j) = vocn(i,j)*cosw + uocn(i,j)*sinw
#ifndef coupled
        ! calculate tilt from geostrophic currents if needed
c        strtltx(i,j) = -fm(i,j)*vocn(i,j)
c        strtlty(i,j) =  fm(i,j)*uocn(i,j)
        ! combine tilt with wind stress
        forcex(i,j) = strairx(i,j) - fm(i,j)*vocn(i,j)
        forcey(i,j) = strairy(i,j) + fm(i,j)*uocn(i,j)
#else
        forcex(i,j) = strairx(i,j) - gravit*umass(i,j)*strtltx(i,j)
        forcey(i,j) = strairy(i,j) - gravit*umass(i,j)*strtlty(i,j)
#endif
        endif  ! umask

       else    ! .not. icetmask
        divu(i,j) = c0
        Delta(i,j) = c0
        shear(i,j) = c0

        do k=1,4
         stressp (k,i,j) = c0
         stressm (k,i,j) = c0
         stress12(k,i,j) = c0
        enddo

       endif  ! icetmask

        ! initialize velocity for new ice points to ocean sfc current
        if( iceumask(i,j) .and. (.not. iceumask_old(i,j))) then
          u(i,j) = uocn(i,j)
          v(i,j) = vocn(i,j)
        endif 

       enddo
      enddo

      end subroutine evp_prep

c=======================================================================

      subroutine stress(ksub)

!---!-------------------------------------------------------------------
!---! Computes the rates of strain and internal stress components for 
!---! each of the four corners on each T-grid cell
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: ksub  ! subcycling step
      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: 
     &  divune, divunw, divuse, divusw             ! divergence
     &, tensionne, tensionnw, tensionse, tensionsw ! tension
     &, shearne, shearnw, shearse, shearsw         ! shearing
     &, Deltane, Deltanw, Deltase, Deltasw         ! Delta 
     &, c0ne, c0nw, c0se, c0sw                     ! useful combinations
     &, c1ne, c1nw, c1se, c1sw

      do j=jlo,jhi
       do i=ilo,ihi

      if (icetmask(i,j)) then

      !-----------------------------------------------------------------
      ! strain rates            ! 1/s
      !-----------------------------------------------------------------
      ! divergence  =  e_11 + e_22
      divune    = cyp(i,j)*u(i  ,j  ) - dyt(i,j)*u(i-1,j  )
     &          + cxp(i,j)*v(i  ,j  ) - dxt(i,j)*v(i  ,j-1)
      divunw    = cym(i,j)*u(i-1,j  ) + dyt(i,j)*u(i  ,j  )
     &          + cxp(i,j)*v(i-1,j  ) - dxt(i,j)*v(i-1,j-1)
      divusw    = cym(i,j)*u(i-1,j-1) + dyt(i,j)*u(i  ,j-1)
     &          + cxm(i,j)*v(i-1,j-1) + dxt(i,j)*v(i-1,j  )
      divuse    = cyp(i,j)*u(i  ,j-1) - dyt(i,j)*u(i-1,j-1)
     &          + cxm(i,j)*v(i  ,j-1) + dxt(i,j)*v(i  ,j  )

      ! tension strain rate  =  e_11 - e_22
      tensionne = -cym(i,j)*u(i  ,j  ) - dyt(i,j)*u(i-1,j  )
     &          +  cxm(i,j)*v(i  ,j  ) + dxt(i,j)*v(i  ,j-1)
      tensionnw = -cyp(i,j)*u(i-1,j  ) + dyt(i,j)*u(i  ,j  )
     &          +  cxm(i,j)*v(i-1,j  ) + dxt(i,j)*v(i-1,j-1)
      tensionsw = -cyp(i,j)*u(i-1,j-1) + dyt(i,j)*u(i  ,j-1)
     &          +  cxp(i,j)*v(i-1,j-1) - dxt(i,j)*v(i-1,j  )
      tensionse = -cym(i,j)*u(i  ,j-1) - dyt(i,j)*u(i-1,j-1)
     &          +  cxp(i,j)*v(i  ,j-1) - dxt(i,j)*v(i  ,j  )

      ! shearing strain rate  =  e_12
      shearne = -cym(i,j)*v(i  ,j  ) - dyt(i,j)*v(i-1,j  )
     &        -  cxm(i,j)*u(i  ,j  ) - dxt(i,j)*u(i  ,j-1)
      shearnw = -cyp(i,j)*v(i-1,j  ) + dyt(i,j)*v(i  ,j  )
     &        -  cxm(i,j)*u(i-1,j  ) - dxt(i,j)*u(i-1,j-1)
      shearsw = -cyp(i,j)*v(i-1,j-1) + dyt(i,j)*v(i  ,j-1)
     &        -  cxp(i,j)*u(i-1,j-1) + dxt(i,j)*u(i-1,j  )
      shearse = -cym(i,j)*v(i  ,j-1) - dyt(i,j)*v(i-1,j-1)
     &        -  cxp(i,j)*u(i  ,j-1) + dxt(i,j)*u(i  ,j  )

      ! Delta (in the denominator of zeta, eta)
      Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
      Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
      Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
      Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))

      !-----------------------------------------------------------------
      ! save quantities for mechanical redistribution
      !-----------------------------------------------------------------
      if (ksub.eq.ndte) then
      divu(i,j) = p25*(divune + divunw + divuse + divusw)*tarear(i,j)
      Delta(i,j) = p25*(Deltane+Deltanw+Deltase+Deltasw)*tarear(i,j)
      ! diagnostic only
      ! shear = sqrt(tension**2 + shearing**2) 
      shear(i,j) = p25*tarear(i,j)*sqrt(
     &  (tensionne + tensionnw + tensionse + tensionsw)**2
     & +(  shearne +   shearnw +   shearse +   shearsw)**2)
      endif

      !-----------------------------------------------------------------
      ! replacement pressure/Delta                   ! kg/s
      ! save replacement pressure for principal stress calculation
      !-----------------------------------------------------------------
      if (evp_damping) then
        ! enforce damping criterion
        c0ne = min(prss(i,j)/max(Deltane,c4*tinyarea(i,j)),rcon(i,j))
        c0nw = min(prss(i,j)/max(Deltanw,c4*tinyarea(i,j)),rcon(i,j))
        c0sw = min(prss(i,j)/max(Deltasw,c4*tinyarea(i,j)),rcon(i,j))
        c0se = min(prss(i,j)/max(Deltase,c4*tinyarea(i,j)),rcon(i,j))
c!!        prs_sig(i,j) = prss(i,j)*Deltane/max(Deltane,puny) ! northeast
        prs_sig(i,j) = prss(i,j)*Deltane/max(Deltane,c4*tinyarea(i,j)) ! ne
      else
        ! original version
        c0ne = prss(i,j)/max(Deltane,tinyarea(i,j))
        c0nw = prss(i,j)/max(Deltanw,tinyarea(i,j))
        c0sw = prss(i,j)/max(Deltasw,tinyarea(i,j))
        c0se = prss(i,j)/max(Deltase,tinyarea(i,j))
        prs_sig(i,j) = c0ne*Deltane ! northeast
      endif

      c1ne = c0ne*dte2T
      c1nw = c0nw*dte2T
      c1sw = c0sw*dte2T
      c1se = c0se*dte2T

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

      stressp(1,i,j) = (stressp(1,i,j) + c1ne*(divune - Deltane))*denom1
      stressp(2,i,j) = (stressp(2,i,j) + c1nw*(divunw - Deltanw))*denom1
      stressp(3,i,j) = (stressp(3,i,j) + c1sw*(divusw - Deltasw))*denom1
      stressp(4,i,j) = (stressp(4,i,j) + c1se*(divuse - Deltase))*denom1

      stressm(1,i,j) = (stressm(1,i,j) + c1ne*tensionne)*denom2
      stressm(2,i,j) = (stressm(2,i,j) + c1nw*tensionnw)*denom2
      stressm(3,i,j) = (stressm(3,i,j) + c1sw*tensionsw)*denom2
      stressm(4,i,j) = (stressm(4,i,j) + c1se*tensionse)*denom2

      stress12(1,i,j) = (stress12(1,i,j) + c1ne*shearne*p5)*denom2
      stress12(2,i,j) = (stress12(2,i,j) + c1nw*shearnw*p5)*denom2
      stress12(3,i,j) = (stress12(3,i,j) + c1sw*shearsw*p5)*denom2
      stress12(4,i,j) = (stress12(4,i,j) + c1se*shearse*p5)*denom2

      endif

       enddo
      enddo

      end subroutine stress

c=======================================================================

      subroutine stepu

!---!-------------------------------------------------------------------
!---! Calculation of the surface stresses
!---! Integration of the momentum equation to find velocity (u,v)
!---!-------------------------------------------------------------------

      use ice_flux

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: 
     &   vrel,cca,ccb,ab2,cc1,cc2, taux,tauy
     &,  str(8,imt_local,jmt_local)
     &,  ssigpn, ssigps, ssigpe, ssigpw
     &,  ssigmn, ssigms, ssigme, ssigmw
     &,  ssig12n, ssig12s, ssig12e, ssig12w
     &,  ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122
     &,  csigpne, csigpnw, csigpse, csigpsw
     &,  csigmne, csigmnw, csigmse, csigmsw
     &,  csig12ne, csig12nw, csig12se, csig12sw
     &,  str12ew, str12we, str12ns, str12sn
     &,  strp_tmp, strm_tmp, str12_tmp

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
       if (icetmask(i,j)) then
        ssigpn  = stressp(1,i,j)  + stressp(2,i,j)
        ssigps  = stressp(3,i,j)  + stressp(4,i,j)
        ssigpe  = stressp(1,i,j)  + stressp(4,i,j)
        ssigpw  = stressp(2,i,j)  + stressp(3,i,j)
        ssigp1  =(stressp(1,i,j)  + stressp(3,i,j))*p055
        ssigp2  =(stressp(2,i,j)  + stressp(4,i,j))*p055

        ssigmn  = stressm(1,i,j)  + stressm(2,i,j)
        ssigms  = stressm(3,i,j)  + stressm(4,i,j)
        ssigme  = stressm(1,i,j)  + stressm(4,i,j)
        ssigmw  = stressm(2,i,j)  + stressm(3,i,j)
        ssigm1  =(stressm(1,i,j)  + stressm(3,i,j))*p055
        ssigm2  =(stressm(2,i,j)  + stressm(4,i,j))*p055

        ssig12n = stress12(1,i,j) + stress12(2,i,j)
        ssig12s = stress12(3,i,j) + stress12(4,i,j)
        ssig12e = stress12(1,i,j) + stress12(4,i,j)
        ssig12w = stress12(2,i,j) + stress12(3,i,j)
        ssig121 =(stress12(1,i,j) + stress12(3,i,j))*p111
        ssig122 =(stress12(2,i,j) + stress12(4,i,j))*p111

        csigpne = p111*stressp(1,i,j) + ssigp2 + p027*stressp(3,i,j)
        csigpnw = p111*stressp(2,i,j) + ssigp1 + p027*stressp(4,i,j)
        csigpsw = p111*stressp(3,i,j) + ssigp2 + p027*stressp(1,i,j)
        csigpse = p111*stressp(4,i,j) + ssigp1 + p027*stressp(2,i,j)

        csigmne = p111*stressm(1,i,j) + ssigm2 + p027*stressm(3,i,j)
        csigmnw = p111*stressm(2,i,j) + ssigm1 + p027*stressm(4,i,j)
        csigmsw = p111*stressm(3,i,j) + ssigm2 + p027*stressm(1,i,j)
        csigmse = p111*stressm(4,i,j) + ssigm1 + p027*stressm(2,i,j)

        csig12ne = p222*stress12(1,i,j) + ssig122 + p055*stress12(3,i,j)
        csig12nw = p222*stress12(2,i,j) + ssig121 + p055*stress12(4,i,j)
        csig12sw = p222*stress12(3,i,j) + ssig122 + p055*stress12(1,i,j)
        csig12se = p222*stress12(4,i,j) + ssig121 + p055*stress12(2,i,j) 

        str12ew = dxt2(i,j)*(p333*ssig12e + p166*ssig12w)
        str12we = dxt2(i,j)*(p333*ssig12w + p166*ssig12e)
        str12ns = dyt2(i,j)*(p333*ssig12n + p166*ssig12s)
        str12sn = dyt2(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
        strp_tmp  = dyt4(i,j)*(p333*ssigpn  + p166*ssigps)
        strm_tmp  = dyt4(i,j)*(p333*ssigmn  + p166*ssigms)
        ! northeast (i,j)
        str(1,i,j) = -strp_tmp - strm_tmp - str12ew 
     &  + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

        ! northwest (i+1,j)
        str(2,i,j) = strp_tmp + strm_tmp - str12we 
     &  + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

        strp_tmp  = dyt4(i,j)*(p333*ssigps  + p166*ssigpn)
        strm_tmp  = dyt4(i,j)*(p333*ssigms  + p166*ssigmn)
        ! southeast (i,j+1)
        str(3,i,j) = -strp_tmp - strm_tmp + str12ew 
     &  + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

        ! southwest (i+1,j+1)
        str(4,i,j) = strp_tmp + strm_tmp + str12we 
     &  + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw


      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
        strp_tmp  = dxt4(i,j)*(p333*ssigpe  + p166*ssigpw)
        strm_tmp  = dxt4(i,j)*(p333*ssigme  + p166*ssigmw)
        ! northeast (i,j)
        str(5,i,j) = -strp_tmp + strm_tmp - str12ns 
     &  - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

        ! southeast (i,j+1)
        str(6,i,j) = strp_tmp - strm_tmp - str12sn 
     &  - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

        strp_tmp  = dxt4(i,j)*(p333*ssigpw  + p166*ssigpe)
        strm_tmp  = dxt4(i,j)*(p333*ssigmw  + p166*ssigme)
        ! northwest (i+1,j)
        str(7,i,j) = -strp_tmp + strm_tmp + str12ns
     &  - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

        ! southwest (i+1,j+1)
        str(8,i,j) = strp_tmp - strm_tmp + str12sn 
     &  - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

       else
        str(1,i,j) = c0
        str(2,i,j) = c0
        str(3,i,j) = c0
        str(4,i,j) = c0
        str(5,i,j) = c0
        str(6,i,j) = c0
        str(7,i,j) = c0
        str(8,i,j) = c0
       endif
      enddo
      enddo
      call bound_narr_ne(8,str)

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(i,j,vrel,cca,ccb,ab2,
c$OMP&   s11,s12,s21,s22,cc1,cc2,taux,tauy)
      do j=jlo,jhi
      do i=ilo,ihi

       if (iceumask(i,j)) then
        ! (magnitude of relative ocean current)*rhow*drag
        vrel = dragw*sqrt((uocn(i,j) - u(i,j))**2 + 
     &                    (vocn(i,j) - v(i,j))**2)  ! m/s
        ! ice/ocean stress
        taux = vrel*waterx(i,j) ! NOTE this is not the entire
        tauy = vrel*watery(i,j) ! ocn stress term

        ! alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
        cca = umassdtei(i,j) + vrel         ! alpha, kg/m^2 s
        ccb = fm(i,j)                       ! beta,  kg/m^2 s
        ab2 = cca**2 + ccb**2

        ! divergence of the internal stress tensor
        strintx(i,j) = uarear(i,j)*
     &      (str(1,i,j) + str(2,i+1,j) + str(3,i,j+1) + str(4,i+1,j+1))
        strinty(i,j) = uarear(i,j)*
     &      (str(5,i,j) + str(6,i,j+1) + str(7,i+1,j) + str(8,i+1,j+1))

        ! finally, the velocity components
        cc1 = strintx(i,j) + forcex(i,j) + taux + umassdtei(i,j)*u(i,j)
        cc2 = strinty(i,j) + forcey(i,j) + tauy + umassdtei(i,j)*v(i,j)

        u(i,j) = (cca*cc1 + ccb*cc2)/ab2              ! m/s
        v(i,j) = (cca*cc2 - ccb*cc1)/ab2

      !-----------------------------------------------------------------
      ! ocean-ice stress for coupling 
      !-----------------------------------------------------------------
        strocnx(i,j) = taux  
        strocny(i,j) = tauy

      else
      !-----------------------------------------------------------------
      ! set velocity and stress to zero on land and (nearly) open water
      !-----------------------------------------------------------------
        u(i,j) = c0
        v(i,j) = c0
        strocnx(i,j) = c0
        strocny(i,j) = c0
        strintx(i,j) = c0
        strinty(i,j) = c0
      endif

      enddo
      enddo

      end subroutine stepu

c=======================================================================

      subroutine evp_finish

!---!-------------------------------------------------------------------
!---! calculation of the ice-ocean stress
!---! ...the sign will be reversed later...
!---!-------------------------------------------------------------------

      use ice_flux

      real (kind=dbl_kind) :: vrel
      integer (kind=int_kind) :: i, j

c$OMP PARALLEL DO PRIVATE(i,j,vrel,factor,UTE,UTW,VTN,VTS)
      do j=jlo,jhi
      do i=ilo,ihi
      if (iceumask(i,j)) then
        ! ocean-ice stress for coupling
        vrel = dragw*sqrt((uocn(i,j) - u(i,j))**2 + 
     &                    (vocn(i,j) - v(i,j))**2)  ! m/s
        strocnx(i,j) = strocnx(i,j) - vrel*(u(i,j)*cosw - v(i,j)*sinw)
        strocny(i,j) = strocny(i,j) - vrel*(v(i,j)*cosw + u(i,j)*sinw)
      else
        ! set stress to zero on land and (nearly) open water
        strocnx(i,j) = c0
        strocny(i,j) = c0
      endif

      !-----------------------------------------------------------------
      ! convert strocn to T grid
      !-----------------------------------------------------------------
      strocnxT(i,j) = strocnx(i,j)  ! prepare to shift
      strocnyT(i,j) = strocny(i,j)

      enddo
      enddo

      call u2tgrid(strocnxT)        ! shift
      call u2tgrid(strocnyT)

      end subroutine evp_finish

c=======================================================================

      subroutine principal_stress

!---!-------------------------------------------------------------------
!---! computes principal stresses for comparison with the theoretical 
!---! yield curve
!---! northeast values
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
        do i=ilo,ihi
          if(prs_sig(i,j).gt.puny) then
           sig1(i,j)=(p5*(stressp(1,i,j)
     &              +sqrt(stressm(1,i,j)**2+c4*stress12(1,i,j)**2)))
     &              /prs_sig(i,j)
           sig2(i,j)=(p5*(stressp(1,i,j)
     &              -sqrt(stressm(1,i,j)**2+c4*stress12(1,i,j)**2)))
     &              /prs_sig(i,j)
          else 
            sig1(i,j)=1000._dbl_kind
            sig2(i,j)=1000._dbl_kind
          endif
        enddo
      enddo

      end subroutine principal_stress

c=======================================================================

      end module ice_dyn_evp

c=======================================================================
