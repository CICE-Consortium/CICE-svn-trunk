! $Id: ice_dyn_evp.F,v 1.15 2004/02/25 17:35:13 eclare Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_dyn_evp - elastic-viscous-plastic sea ice dynamics model 
!
! !DESCRIPTION:
!
! Elastic-viscous-plastic sea ice dynamics model \\
! Computes ice velocity and deformation \\
!
! See: \\
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. {\em J. Phys. Oceanogr.}, {\bf 27}, 1849--1867. \\
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. {\em Journal of Computational Physics}, {\bf 170},
! 18--38. \\
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere---Incorporation of Metric Terms. {\em Monthly Weather Review},
! {\bf 130}, 1848--1865.\\
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.\\
!
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! Summer 2003:  Vectorized by Clifford Chen (Fujitsu) and
!               William H. Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_dyn_evp
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain
      use ice_grid
      use ice_constants
      use ice_state
      use ice_work, only:  worka, workb
!
!EOP
!
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
     &   ecci         ! 1/e^2 
     &,  dte2T        ! dte/2T
     &,  denom1       ! constants for stress equation
     &,  denom2       !

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   shear        ! strain rate II component (1/s)
     &,  divu         ! strain rate I component, velocity divergence (1/s)
     &,  Delta        ! function of strain rates (1/s)

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
     &,  strintx      ! divergence of internal ice stress, x (N/m^2)
     &,  strinty      ! divergence of internal ice stress, y (N/m^2)
     &,  sig1         ! principal stress component (diagnostic)
     &,  sig2         ! principal stress component (diagnostic)
     &,  aiu          ! ice fraction on u-grid

      ! ice stress tensor in each corner of T cell (kg/s^2)
      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   stressp_1, stressp_2, stressp_3, stressp_4   ! sigma11+sigma22
     &,  stressm_1, stressm_2, stressm_3, stressm_4   ! sigma11-sigma22
     &,  stress12_1,stress12_2,stress12_3,stress12_4  ! sigma12

      integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) ::
     &   indxti       ! compressed index in i-direction
     &,  indxtj       ! compressed index in j-direction
     &,  indxui       ! compressed index in i-direction
     &,  indxuj       ! compressed index in j-direction

      integer (kind=int_kind) ::
     &   icellt       ! total count when icetmask is true
     &,  icellu       ! total count when iceumask is true

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: evp - elastic-viscous-plastic dynamics driver
!
! !INTERFACE:
!
      subroutine evp(kstrngth)
!
! !DESCRIPTION:
!
! Elastic-viscous-plastic dynamics driver
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
      use ice_timers
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   kstrngth    
!
!EOP
!
      integer (kind=int_kind) :: k

      call ice_timer_start(2)  ! dynamics

      call evp_prep(kstrngth)  ! preparation for dynamics

      do k=1,ndte              ! subcycling
        call stress(k)         ! the stress tensor equation
        call stepu             ! total surface stress, momentum equation
        call bound_sw(uvel)    ! Periodic/Neumann boundary conditions
        call bound_sw(vvel)    ! Periodic/Neumann boundary conditions
      enddo

      call evp_finish          ! ice-ocean stress

      call ice_timer_stop(2)   ! dynamics

      end subroutine evp

!=======================================================================
!BOP
!
! !IROUTINE: init_evp - initialize parameters needed for evp dynamics
!
! !INTERFACE:
!
      subroutine init_evp
!
! !DESCRIPTION:
!
! Initialize parameters and variables needed for the evp dynamics
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
      use ice_calendar
      use ice_fileunits
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
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
     & write(nu_diag,*) 'dt = ',dt,'  dte = ',dte,  
     &                  '  tdamp = ', p5*tdamp2
      ! major/minor axis length ratio, squared
      ecc  = c4
      ecci = p25                  ! 1/ecc
      ! constants for stress equation
      dte2T = dte/tdamp2                    ! unitless
      denom1 = c1/(c1+dte2T)
      denom2 = c1/(c1+dte2T*ecc)
      rcon = 1230._dbl_kind*eyc*dt*dtei**2  ! kg/s  

      do j=jlo,jhi
       do i=ilo,ihi
        ! Coriolis parameter
c        fcor(i,j) = 1.46e-4_dbl_kind ! Hibler 1979, Northern Hemisphere; 1/s
        fcor(i,j) = c2*omega*sin(ULAT(i,j))  ! 1/s

        ! velocity
        uvel(i,j) = c0       ! m/s
        vvel(i,j) = c0       ! m/s

       enddo
      enddo

      ! stress tensor,  kg/s^2
      do j=jlo,jhi
       do i=ilo,ihi
         stressp_1 (i,j) = c0
         stressp_2 (i,j) = c0
         stressp_3 (i,j) = c0
         stressp_4 (i,j) = c0
         stressm_1 (i,j) = c0
         stressm_2 (i,j) = c0
         stressm_3 (i,j) = c0
         stressm_4 (i,j) = c0
         stress12_1(i,j) = c0
         stress12_2(i,j) = c0
         stress12_3(i,j) = c0
         stress12_4(i,j) = c0
       enddo
      enddo

      call bound(uvel)
      call bound(vvel)

      end subroutine init_evp

!=======================================================================
!BOP
!
! !IROUTINE: evp_prep - compute quantities needed for stress tensor and mom eqns
!
! !INTERFACE:
!
      subroutine evp_prep(kstrngth)
!
! !DESCRIPTION:
!
! Computes quantities needed in the stress tensor (sigma)  \\
! and momentum (u) equations, but which do not change during  \\
! the thermodynamics/transport time step: \\
! --wind stress shift to U grid, \\
! --ice mass and ice extent masks, \\
! --pressure (strength), and part of the forcing stresses \\
! initializes ice velocity for new points to ocean sfc current  \\
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
      use ice_flux
      use ice_calendar
      use ice_mechred_cice
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   kstrngth    
!
!EOP
!
      integer (kind=int_kind) :: i, j, k, n

      real (kind=dbl_kind) :: 
     &   umass(ilo:ihi,jlo:jhi) ! ice mass on u-grid

      real (kind=dbl_kind), parameter ::
     &   a_min = p001  ! minimum ice area
     &,  m_min = p01   ! minimum ice mass

      logical (kind=log_kind) ::
     &   tmphm(imt_local,jmt_local)
     &,  iceumask_old(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! total mass of ice and snow, centered in T-cell
      ! NOTE: vice and vsno must be up to date in all grid cells,
      !       including ghost cells
      !-----------------------------------------------------------------
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
      do j=jlo,jhi
       do i=ilo,ihi
        ! Factor of aice needed for correct treatment of free drift
        strairx(i,j) = strairxT(i,j)*aice(i,j) ! prep to convert to U grid
        strairy(i,j) = strairyT(i,j)*aice(i,j)
       enddo
      enddo
      call t2ugrid(strairx)
      call t2ugrid(strairy)
      call to_ugrid(tmass,umass)
      call to_ugrid(aice,aiu)

      !-----------------------------------------------------------------
      ! convenient variable for evp
      !-----------------------------------------------------------------
      do j=jlo,jhi
       do i=ilo,ihi
        umassdtei(i,j) = umass(i,j)*dtei               ! m/dte, kg/m^2 s
       enddo
      enddo

      !-----------------------------------------------------------------
      ! augmented masks (land + open ocean)
      !-----------------------------------------------------------------
      do j=1,jmt_local
       do i=1,imt_local
        ! ice extent mask (T-cells)
        tmphm(i,j) = tmask(i,j) .and. (aice (i,j).gt.a_min)
     &                          .and. (tmass(i,j).gt.m_min)
       enddo
      enddo
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
         iceumask(i,j) = (umask(i,j)) .and. (aiu  (i,j).gt.a_min)   
     &                                .and. (umass(i,j).gt.m_min)
       enddo
      enddo

      !-----------------------------------------------------------------
      ! pressure and forcing terms; set sigma=0 for no ice;
      ! initialize ice velocity in cells previously empty to ocn current
      !-----------------------------------------------------------------
      call ice_strength(kstrngth)

      do j=jlo,jhi
       do i=ilo,ihi
       if (icetmask(i,j)) then
        prss(i,j) = strength(i,j)
        fm(i,j) = fcor(i,j)*umass(i,j)               ! Coriolis * mass

        if (umask(i,j)) then
        ! for ocean stress
        waterx(i,j) = uocn(i,j)*cosw - vocn(i,j)*sinw
        watery(i,j) = vocn(i,j)*cosw + uocn(i,j)*sinw
        ! combine tilt with wind stress
#ifndef coupled
        ! calculate tilt from geostrophic currents if needed
        strtltx(i,j) = -fm(i,j)*vocn(i,j)
        strtlty(i,j) =  fm(i,j)*uocn(i,j)
#else
        strtltx(i,j) = -gravit*umass(i,j)*ss_tltx(i,j)
        strtlty(i,j) = -gravit*umass(i,j)*ss_tlty(i,j)
#endif
        forcex(i,j) = strairx(i,j) + strtltx(i,j)
        forcey(i,j) = strairy(i,j) + strtlty(i,j)
        endif  ! umask

       else    ! .not. icetmask
        stressp_1 (i,j) = c0
        stressp_2 (i,j) = c0
        stressp_3 (i,j) = c0
        stressp_4 (i,j) = c0
        stressm_1 (i,j) = c0
        stressm_2 (i,j) = c0
        stressm_3 (i,j) = c0
        stressm_4 (i,j) = c0
        stress12_1(i,j) = c0
        stress12_2(i,j) = c0
        stress12_3(i,j) = c0
        stress12_4(i,j) = c0

        divu(i,j) = c0
        Delta(i,j) = c0
        shear(i,j) = c0

       endif  ! icetmask

        ! initialize velocity for new ice points to ocean sfc current
        if( iceumask(i,j) .and. (.not. iceumask_old(i,j))) then
          uvel(i,j) = uocn(i,j)
          vvel(i,j) = vocn(i,j)
        endif 

       enddo
      enddo

      icellt = 0
      do j=jlo,jhi
        do i=ilo,ihi
          if (icetmask(i,j)) then
            icellt = icellt + 1
            indxti(icellt) = i
            indxtj(icellt) = j
          endif
        enddo
      enddo

      icellu = 0
      do j=jlo,jhi
        do i=ilo,ihi
          if (iceumask(i,j)) then
            icellu = icellu + 1
            indxui(icellu) = i
            indxuj(icellu) = j
          endif
        enddo
      enddo

      end subroutine evp_prep

!=======================================================================
!BOP
!
! !IROUTINE: stress - computes strain rates and internal stress components
!
! !INTERFACE:
!
      subroutine stress(ksub)
!
! !DESCRIPTION:
!
! Computes the rates of strain and internal stress components for \\
! each of the four corners on each T-grid cell
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: ksub  ! subcycling step
!
!EOP
!
      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) :: 
     &  divune, divunw, divuse, divusw             ! divergence
     &, tensionne, tensionnw, tensionse, tensionsw ! tension
     &, shearne, shearnw, shearse, shearsw         ! shearing
     &, Deltane, Deltanw, Deltase, Deltasw         ! Delta 
     &, c0ne, c0nw, c0se, c0sw                     ! useful combinations
     &, c1ne, c1nw, c1se, c1sw

      integer (kind=int_kind) ::
     &  ij      ! loop index, combination of i and j loops

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij=1,icellt
       i = indxti(ij)
       j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
      ! divergence  =  e_11 + e_22
      divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  )
     &          + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
      divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  )
     &          + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
      divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1)
     &          + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
      divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1)
     &          + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )

      ! tension strain rate  =  e_11 - e_22
      tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  )
     &          +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
      tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  )
     &          +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
      tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1)
     &          +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
      tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1)
     &          +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )

      ! shearing strain rate  =  e_12
      shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  )
     &        -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
      shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  )
     &        -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
      shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1)
     &        -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
      shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1)
     &        -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )

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

      stressp_1(i,j) = (stressp_1(i,j) + c1ne*(divune - Deltane))*denom1
      stressp_2(i,j) = (stressp_2(i,j) + c1nw*(divunw - Deltanw))*denom1
      stressp_3(i,j) = (stressp_3(i,j) + c1sw*(divusw - Deltasw))*denom1
      stressp_4(i,j) = (stressp_4(i,j) + c1se*(divuse - Deltase))*denom1

      stressm_1(i,j) = (stressm_1(i,j) + c1ne*tensionne)*denom2
      stressm_2(i,j) = (stressm_2(i,j) + c1nw*tensionnw)*denom2
      stressm_3(i,j) = (stressm_3(i,j) + c1sw*tensionsw)*denom2
      stressm_4(i,j) = (stressm_4(i,j) + c1se*tensionse)*denom2

      stress12_1(i,j) = (stress12_1(i,j) + c1ne*shearne*p5)*denom2
      stress12_2(i,j) = (stress12_2(i,j) + c1nw*shearnw*p5)*denom2
      stress12_3(i,j) = (stress12_3(i,j) + c1sw*shearsw*p5)*denom2
      stress12_4(i,j) = (stress12_4(i,j) + c1se*shearse*p5)*denom2

      enddo                     ! ij

      end subroutine stress

!=======================================================================
!BOP
!
! !IROUTINE: stepu - surface stresses and integrates mom eqn for u,v
!
! !INTERFACE:
!
      subroutine stepu
!
! !DESCRIPTION:
!
! Calculation of the surface stresses \\
! Integration of the momentum equation to find velocity (u,v)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) :: 
     &   vrel,cca,ccb,ab2,cc1,cc2, taux,tauy
     &,  str(imt_local,jmt_local,8)
     &,  ssigpn, ssigps, ssigpe, ssigpw
     &,  ssigmn, ssigms, ssigme, ssigmw
     &,  ssig12n, ssig12s, ssig12e, ssig12w
     &,  ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122
     &,  csigpne, csigpnw, csigpse, csigpsw
     &,  csigmne, csigmnw, csigmse, csigmsw
     &,  csig12ne, csig12nw, csig12se, csig12sw
     &,  str12ew, str12we, str12ns, str12sn
     &,  strp_tmp, strm_tmp, str12_tmp

      integer (kind=int_kind) ::
     &   ij      ! loop index, combination of i and j loops

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

      str(:,:,:) = c0

      do ij=1,icellt
        i = indxti(ij)
        j = indxtj(ij)

        ssigpn  = stressp_1(i,j) + stressp_2(i,j)
        ssigps  = stressp_3(i,j) + stressp_4(i,j)
        ssigpe  = stressp_1(i,j) + stressp_4(i,j)
        ssigpw  = stressp_2(i,j) + stressp_3(i,j)
        ssigp1  =(stressp_1(i,j) + stressp_3(i,j))*p055
        ssigp2  =(stressp_2(i,j) + stressp_4(i,j))*p055

        ssigmn  = stressm_1(i,j) + stressm_2(i,j)
        ssigms  = stressm_3(i,j) + stressm_4(i,j)
        ssigme  = stressm_1(i,j) + stressm_4(i,j)
        ssigmw  = stressm_2(i,j) + stressm_3(i,j)
        ssigm1  =(stressm_1(i,j) + stressm_3(i,j))*p055
        ssigm2  =(stressm_2(i,j) + stressm_4(i,j))*p055

        ssig12n = stress12_1(i,j) + stress12_2(i,j)
        ssig12s = stress12_3(i,j) + stress12_4(i,j)
        ssig12e = stress12_1(i,j) + stress12_4(i,j)
        ssig12w = stress12_2(i,j) + stress12_3(i,j)
        ssig121 =(stress12_1(i,j) + stress12_3(i,j))*p111
        ssig122 =(stress12_2(i,j) + stress12_4(i,j))*p111

        csigpne = p111*stressp_1(i,j) + ssigp2 + p027*stressp_3(i,j)
        csigpnw = p111*stressp_2(i,j) + ssigp1 + p027*stressp_4(i,j)
        csigpsw = p111*stressp_3(i,j) + ssigp2 + p027*stressp_1(i,j)
        csigpse = p111*stressp_4(i,j) + ssigp1 + p027*stressp_2(i,j)

        csigmne = p111*stressm_1(i,j) + ssigm2 + p027*stressm_3(i,j)
        csigmnw = p111*stressm_2(i,j) + ssigm1 + p027*stressm_4(i,j)
        csigmsw = p111*stressm_3(i,j) + ssigm2 + p027*stressm_1(i,j)
        csigmse = p111*stressm_4(i,j) + ssigm1 + p027*stressm_2(i,j)

        csig12ne = p222*stress12_1(i,j) + ssig122 + p055*stress12_3(i,j)
        csig12nw = p222*stress12_2(i,j) + ssig121 + p055*stress12_4(i,j)
        csig12sw = p222*stress12_3(i,j) + ssig122 + p055*stress12_1(i,j)
        csig12se = p222*stress12_4(i,j) + ssig121 + p055*stress12_2(i,j)

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
        str(i,j,1) = -strp_tmp - strm_tmp - str12ew 
     &  + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

        ! northwest (i+1,j)
        str(i,j,2) = strp_tmp + strm_tmp - str12we 
     &  + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

        strp_tmp  = dyt4(i,j)*(p333*ssigps  + p166*ssigpn)
        strm_tmp  = dyt4(i,j)*(p333*ssigms  + p166*ssigmn)

        ! southeast (i,j+1)
        str(i,j,3) = -strp_tmp - strm_tmp + str12ew 
     &  + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

        ! southwest (i+1,j+1)
        str(i,j,4) = strp_tmp + strm_tmp + str12we 
     &  + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw


      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
        strp_tmp  = dxt4(i,j)*(p333*ssigpe  + p166*ssigpw)
        strm_tmp  = dxt4(i,j)*(p333*ssigme  + p166*ssigmw)

        ! northeast (i,j)
        str(i,j,5) = -strp_tmp + strm_tmp - str12ns 
     &  - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

        ! southeast (i,j+1)
        str(i,j,6) = strp_tmp - strm_tmp - str12sn 
     &  - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

        strp_tmp  = dxt4(i,j)*(p333*ssigpw  + p166*ssigpe)
        strm_tmp  = dxt4(i,j)*(p333*ssigmw  + p166*ssigme)

        ! northwest (i+1,j)
        str(i,j,7) = -strp_tmp + strm_tmp + str12ns
     &  - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

        ! southwest (i+1,j+1)
        str(i,j,8) = strp_tmp - strm_tmp + str12sn 
     &  - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

      enddo                     ! ij

      call bound_narr_ne(8,str)

      !-----------------------------------------------------------------
      ! set velocity and stress to zero on land and (nearly) open water
      ! use working arrays to avoid u and v vector dependency
      !_________________________________________________________________

      do j=jlo,jhi
      do i=ilo,ihi
         if (.not.iceumask(i,j)) then
            worka(i,j) = c0
            workb(i,j) = c0
            strocnx(i,j) = c0
            strocny(i,j) = c0
            strintx(i,j) = c0
            strinty(i,j) = c0
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      do ij=1,icellu
        i=indxui(ij)
        j=indxuj(ij)

        ! (magnitude of relative ocean current)*rhow*drag*aice
        vrel = aiu(i,j)*dragw*sqrt((uocn(i,j) - uvel(i,j))**2 + 
     &                             (vocn(i,j) - vvel(i,j))**2)  ! m/s
        ! ice/ocean stress
        taux = vrel*waterx(i,j) ! NOTE this is not the entire
        tauy = vrel*watery(i,j) ! ocn stress term

        ! alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
        cca = umassdtei(i,j) + vrel * cosw         ! alpha, kg/m^2 s
        ccb = fm(i,j)        + vrel * sinw         ! beta,  kg/m^2 s
        ab2 = cca**2 + ccb**2

        ! divergence of the internal stress tensor
        strintx(i,j) = uarear(i,j)*
     &      (str(i,j,1) + str(i+1,j,2) + str(i,j+1,3) + str(i+1,j+1,4))
        strinty(i,j) = uarear(i,j)*
     &      (str(i,j,5) + str(i,j+1,6) + str(i+1,j,7) + str(i+1,j+1,8))

        ! finally, the velocity components
        cc1 = strintx(i,j) + forcex(i,j) + taux 
     &       + umassdtei(i,j)*uvel(i,j)
        cc2 = strinty(i,j) + forcey(i,j) + tauy 
     &       + umassdtei(i,j)*vvel(i,j)

        worka(i,j) = (cca*cc1 + ccb*cc2)/ab2              ! m/s
        workb(i,j) = (cca*cc2 - ccb*cc1)/ab2

      !-----------------------------------------------------------------
      ! ocean-ice stress for coupling
      ! scale to full grid cell
      !-----------------------------------------------------------------
        strocnx(i,j) = taux/aiu(i,j)  
        strocny(i,j) = tauy/aiu(i,j)  

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! recover u and v arrays
      !-----------------------------------------------------------------
      do j=jlo,jhi
      do i=ilo,ihi
         uvel(i,j) = worka(i,j)
         vvel(i,j) = workb(i,j)
      enddo
      enddo

      end subroutine stepu

!=======================================================================
!BOP
!
! !IROUTINE: evp_finish - calculates ice-ocean stress 
!
! !INTERFACE:
!
      subroutine evp_finish
!
! !DESCRIPTION:
!
! Calculation of the ice-ocean stress \\
! ...the sign will be reversed later...
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) :: vrel

      do j=jlo,jhi
      do i=ilo,ihi
      if (iceumask(i,j)) then
        ! ocean-ice stress for coupling
        vrel = dragw*sqrt((uocn(i,j) - uvel(i,j))**2 + 
     &                    (vocn(i,j) - vvel(i,j))**2)  ! m/s
        strocnx(i,j) = strocnx(i,j) 
     &               - vrel*(uvel(i,j)*cosw - vvel(i,j)*sinw)
        strocny(i,j) = strocny(i,j) 
     &               - vrel*(vvel(i,j)*cosw + uvel(i,j)*sinw)
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

!=======================================================================
!BOP
!
! !IROUTINE: principal_stress - computes principal stress for yield curve
!
! !INTERFACE:
!
      subroutine principal_stress
!
! !DESCRIPTION:
!
! Computes principal stresses for comparison with the theoretical 
! yield curve; northeast values
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
        do i=ilo,ihi
          if(prs_sig(i,j).gt.puny) then
           sig1(i,j)=(p5*(stressp_1(i,j)
     &              +sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2)))
     &              /prs_sig(i,j)
           sig2(i,j)=(p5*(stressp_1(i,j)
     &              -sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2)))
     &              /prs_sig(i,j)
          else 
            sig1(i,j)=1000._dbl_kind
            sig2(i,j)=1000._dbl_kind
          endif
        enddo
      enddo

      end subroutine principal_stress

!=======================================================================

      end module ice_dyn_evp

!=======================================================================
