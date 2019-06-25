c $Id: $
c=======================================================================
!---! Driver for mechanical redistribution and strength computations
!---!
!---! authors William H. Lipscomb
!---!         Elizabeth C. Hunke 
c=======================================================================

      module ice_mechred_cice

      use ice_model_size
      use ice_constants
      use ice_state
      use ice_itd
      use ice_grid
      use ice_fileunits

      implicit none

      integer (kind=int_kind) ::
     &   kstrength        ! 1 for Rothrock 1975 pressure formulation

      real (kind=dbl_kind), parameter ::     
     &   Gstar = 0.15_dbl_kind ! used to compute participation function
     &,  Hstar = 25._dbl_kind  ! determines mean thickness of ridged ice (m)
                               ! Flato & Hibler have Hstar = 100
     &,  Cf = 17._dbl_kind     ! ratio of ridging work to PE change in ridging
     &,  Cs = p5               ! fraction of shear energy contrbtng to ridging
     &,  fsnowrdg = p5         ! snow fraction that survives in ridging
     &,  Gstari = c1/Gstar     ! 1/Gstar

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::     
     &   divu_adv         ! divergence rate due to advection (1/s)
     &,  aksum            ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,0:ncat) ::     
     &   athorn           ! participation function; fraction of ridging/
                          !  closing associated w/ category n
      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ncat) ::     
     &   hrmin            ! minimum ridge thickness
     &,  hrmax            ! maximum ridge thickness
     &,  krdg             ! mean ridge thickness/thickness of ridging ice 
                                      ! 
      real (kind=dbl_kind), dimension (0:ncat) ::     
     &   athornn          ! participation function; fraction of ridging/
                          !  closing associated w/ category n
      real (kind=dbl_kind), dimension (ncat) ::     
     &   hrminn           ! minimum ridge thickness
     &,  hrmaxn           ! maximum ridge thickness
     &,  krdgn            ! mean ridge thickness / thickness of ridging ice
     &,  ardg1            ! area of ridging ice per category
     &,  vrdg1            ! volume of ridging ice per category
     &,  ardg2            ! area of new ridged ice per category
     &,  vrdg2            ! volume of new ridged ice per category

      real (kind=dbl_kind) ::
     &   aksumn           ! ratio of area removed to area ridged
     &,  ardgtot          ! total ice area ridged
     &,  anewtot          ! total area of new ridges
     &,  vrdgtot          ! total ice volume ridged

c=======================================================================

      contains

c=======================================================================

      subroutine init_mechred

!---!-------------------------------------------------------------------
!---! Initialize constants for ridging
!---!
!---! note... This definition does not belong here but we've kept it
!---! for compatibility with the CCSM version.  It ought to be in
!---! init_itd.    -- ech
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i,j,layer

      !-----------------------------------------------------------------
      ! vectors identifying first and last layer in each bin
      !-----------------------------------------------------------------
      ilyr1(1)=1                      ! if nilyr  = 4
      ilyrn(1)=nilyr                  !   ilyr1 = { 1,5,9 } 
      do layer=2,ncat                 !   ilyrn = { 4,8,12} etc
         ilyr1(layer)=ilyrn(layer-1)+1
         ilyrn(layer)=ilyrn(layer-1)+nilyr
      enddo

      end subroutine init_mechred

c=======================================================================

      subroutine ridge_prep
!---!-------------------------------------------------------------------
!---! preparation for ridging and strength calculations
!---!-------------------------------------------------------------------

      integer (kind=int_kind) ::
     &   i,j      
     &,  n                ! thickness category index

      real (kind=dbl_kind) ::     
     &   Gsum(-1:ncat)    ! Gsum(n) = sum of areas in categories 0 to n
     &,  dPE              ! change in potential energy due to ridging
     &,  hin(ncat)        ! ice thickness for each cat        (m)

      real (kind=dbl_kind), parameter ::     
     &   Cp = p5 * gravit * (rhow-rhoi)*rhoi/rhow  ! proprtnl to pot. energy

      call aggregate      ! total ice area

      do j=jlo,jhi
      do i=ilo,ihi
      call to_column(i,j) ! load column variables from state variables

      if (ai0 .lt. Gstar) then  

      !-----------------------------------------------------------------
      ! Compute the participation function athornn(n); this is analogous to
      ! a(h) = b(h)g(h) as defined in Thorndike et al. (1975).
      !
      !               area lost from category n due to ridging/closing
      !  athornn(n) = ------------------------------------------------
      !                  total area lost due to ridging/closing
      !
      ! Assume b(h) = (2/Gstar) * (1 - G(h)/Gstar). 
      ! The expressions for athornn are found by integrating b(h)g(h) between
      ! the category boundaries.
      !-----------------------------------------------------------------
        Gsum(-1) = c0
        Gsum(0) = Gsum(-1) + ai0
        do n = 1, ncat
          Gsum(n) = Gsum(n-1) + ain(n)
        enddo
        if (aice(i,j) .gt. c1) then  ! renormalize Gsum
          do n = 0, ncat
            Gsum(n) = Gsum(n) / aice(i,j)
          enddo
        endif

        do n = 0, ncat
          if (Gsum(n) .lt. Gstar) then
            athornn(n) = Gstari * (Gsum(n)-Gsum(n-1)) * 
     &              (c2 - (Gsum(n-1)+Gsum(n))*Gstari)
          elseif (Gsum(n-1) .lt. Gstar) then
            athornn(n) = Gstari * (Gstar-Gsum(n-1)) * 
     &              (c2 - (Gsum(n-1)+Gstar)*Gstari)
          else
            athornn(n) = c0
          endif
        enddo
         
      !-----------------------------------------------------------------
      ! Compute max and min ridged ice thickness for each ridging category.
      ! Assume ice is uniformly distributed between hrmin and hrmax.
      ! aksumn = net/total = weighted sum of net area loss terms:
      ! total area removed = area of ice that ridges
      !   net area removed = total - area of new ridges
      !-----------------------------------------------------------------
        aksumn = athornn(0)
        do n = 1, ncat
          if (athornn(n) .gt. puny) then
            hin(n) = vin(n)/ain(n)
            hrminn(n) = c2*hin(n) ! Hibler (1980)
            hrmaxn(n) = c2*sqrt(Hstar*hin(n)) ! Hibler (1980)
            if (hrmaxn(n) .lt. hrminn(n)) hrmaxn(n) = hrminn(n) + 0.01
                      ! This can happen for large hin and small Hstar
            krdgn(n) = p5 * (hrmaxn(n) + hrminn(n)) / hin(n)
            aksumn = aksumn + athornn(n) * (c1 - c1/krdgn(n))
          else                ! neglect ridging for this category
            hrminn(n) = c0
            hrmaxn(n) = c0
            krdgn(n) = c1
          endif
        enddo

        athorn(i,j,0) = athornn(0)    ! save for routine ridge_ice
        do n = 1, ncat
         athorn(i,j,n) = athornn(n) ! save for routine ridge_ice
         hrmin(i,j,n) = hrminn(n)   ! save for routine ridge_ice
         hrmax(i,j,n) = hrmaxn(n)   ! save for routine ridge_ice
         krdg(i,j,n) = krdgn(n)     ! save for routine ridge_ice
        enddo
        aksum(i,j) = aksumn           ! save for routine ridge_ice
         
      else  ! ai0 > Gstar
      !-----------------------------------------------------------------
      ! ice is dispersed and therefore weak
      !-----------------------------------------------------------------

        athorn(i,j,0) = c0
        do n = 1, ncat
          athorn(i,j,n) = c0
          hrmin(i,j,n) = c0
          hrmax(i,j,n) = c0
          krdg(i,j,n) = c0
        enddo
        aksum(i,j) = c1

      endif ! ai0

      enddo
      enddo

      return
      end subroutine ridge_prep

c=======================================================================

      subroutine ridge_ice(Delta, divu)

!---!-------------------------------------------------------------------
!---! Computes changes in the thickness distribution due to divergence
!---!  and shear.  
!---!
!---! From Flato and Hibler (1995), the rate of energy dissipation is
!---!     M = P * [Cs*0.5*(Delta-eps1) - min(eps1,0)]
!---! where P is the ice strength
!---!       Cs accounts for energy lost in shear without ridging 
!---!       eps1 = divergence rate
!---!       eps2 = shear rate
!---!       Delta = sqrt (eps1^2 + eps2^2/e^2)   
!---!       e = aspect ratio of elliptical yield curve
!---! See Flato and Hibler (1995) and Stern et al. (1995) for details.
!---!
!---! In the first term of this equation, which represents ridging
!---! due to shear, we set eps1 = divu, the divergence rate computed 
!---! from the velocity field in the dynamics module.  But in the 
!---! second term, we set eps1 = divu_adv, the divergence rate computed 
!---! from the advection/transport module.  Thus we ensure that the area 
!---! of ice/open water removed by convergence is equal to the net area
!---! that enters during transport.
!---! 
!---! Note: The expression in brackets is the NET area of ice and water
!---! removed; i.e. the open water area removed, plus the area of ice
!---! ridged, minus the area of new ridges.  The NET area removed, armvnet, 
!---! is less than the TOTAL area removed, armvtot: 
!---!              armvnet = armvtot * aksumn, where aksumn <=1.
!---!-------------------------------------------------------------------

      use ice_timers
      use ice_calendar
      use ice_flux
      use ice_domain

      integer (kind=int_kind) ::
     &   i,j
     &,  n, n1, n2         ! thickness category indices
     &,  k                 ! ice layer index

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(in) ::
     &   Delta    ! (in the denominator of zeta, eta)             (1/s)
     &,  divu     ! strain rate I component, velocity divergence  (1/s)

      real (kind=dbl_kind) ::     
     &   vice1, vice2      ! ice volume, initial and final
     &,  eice1, eice2      ! ice vol * qice, initial and final
     &,  vlyr              ! ice layer volume
     &,  armvtot           ! area removed due to ridging and closing
     &,  armvnet           ! net area removed (armvtot - area of new ridges)
     &,  ardg(ncat)        ! area of ice ridged
     &,  anew(ncat)        ! area of new ridge
     &,  afrac             ! fraction of category area ridged
     &,  virdg(ncat)       ! ice volume of ridging ice
     &,  vsrdg(ncat)       ! snow volume of ridging ice
     &,  eirdg(nilyr,ncat) ! ice energy of ridging ice
     &,  esrdg(ncat)       ! snow energy of ridging ice
     &,  msnow_ocn         ! mass of snow added to ocean (kg m-2)
     &,  esnow_ocn         ! energy needed to melt msnow_ocn (J m-2)
     &,  farea             ! fraction of new ridge area going to n2
     &,  fvol              ! fraction of new ridge volume going to n2
     &,  hL, hR            ! left and right limits of integration
     &,  dhr               ! hrmax - hrmin
     &,  dhr2              ! hrmax^2 - hrmin^2
     &,  Hbtmp             ! temporary value of hin_max(ncat)
     &,  asurp(ncat)       ! excess ridging area when ardg > ain

      real (kind=dbl_kind), parameter ::     
     &   big = 1.0e+8_dbl_kind

      call ice_timer_start(6)  ! ridging 

      call ridge_prep      ! variables needed for both ridging and strength

      do j=jlo,jhi
      do i=ilo,ihi

      if (aice0(i,j) .lt. Gstar) then   ! otherwise no ridging occurs

        call to_column(i,j) ! load column variables from state variables
        athornn(0) = athorn(i,j,0)
        do n = 1, ncat
          athornn(n) = athorn(i,j,n)
          hrminn(n) = hrmin(i,j,n) 
          hrmaxn(n) = hrmax(i,j,n)
          krdgn(n) = krdg(i,j,n)
        enddo
        aksumn = aksum(i,j) 

        ! quantities that ridging should conserve
        vice1 = c0
        eice1 = c0
        do n = 1, ncat
          vice1 = vice1 + vin(n)
          do k = 1, nilyr
            eice1 = eice1 + ein(k,n)
          enddo
        enddo
         
        ! net ice/water area removed, from Flato and Hibler (1995)
        armvnet = (Cs*p5*(Delta(i,j)-divu(i,j)) - 
     &             min(divu_adv(i,j),c0)) * dt

        ! Make sure enough ice will be removed that the new fractional 
        ! area of ice satisfies aice <= 1

        if (aice(i,j).gt.c1+puny .and. armvnet.lt.aice(i,j)-c1) then
c          write(nu_diag,*) ''
c 	  write(nu_diag,*) 'Correcting armvnet'
c          write(nu_diag,*) 'istep1,my_task,i,j:',istep1,my_task,i,j
c          write(nu_diag,*) 'armvnet =',armvnet
c          write(nu_diag,*) 'aice - 1 =', aice(i,j) - c1
          armvnet = aice(i,j) - c1
        endif

        ! total ice/water area removed
        armvtot = armvnet / aksumn

        ! Reset hin_max(ncat) temporarily to ensure that all ridged 
        ! ice is thinner than hin_max(ncat)
        Hbtmp = hin_max(ncat)
        hin_max(ncat) = big

        ! area to be ridged in each category
        do n = 1, ncat
          asurp(n) = c0       !  Initialize
        enddo

        ! Initialize diagnostics
        ardgtot = c0           ! area of ridging ice
        anewtot = c0           ! area of new ridge
        vrdgtot = c0           ! ridged ice volume

      !-------------------------------------------------------------------
      ! Compute the area, volume, and energy of ice ridging in each
      ! category, along with the area of the resulting ridge.
      !-------------------------------------------------------------------

        do n1 = 1, ncat
          if (athornn(n1) .gt. puny) then
            ardg(n1) = armvtot * athornn(n1) + asurp(n1) 
                                          ! area of ridging ice
            if (ardg(n1) .gt. ain(n1)) then
              write(nu_diag,*) my_task,i,j,'ardg > ain, category', n1
              write(nu_diag,*) 'ardg =',ardg(n1),'ain =',ain(n1)
              if (n1 .lt. ncat) then
                asurp(n1+1) = (ardg(n1)-ain(n1))*(c1-c1/krdgn(n1))
     &               / (c1-c1/krdgn(n1+1)) ! increase ridging in next cat
              else
                write(nu_diag,*) my_task,':',i,j,
     &              'Excess ridging in category ncat'
                     stop
              endif
              ardg(n1) = ain(n1)
            endif
            anew(n1) = ardg(n1) / krdgn(n1) ! area of new ridge
            afrac = ardg(n1) / ain(n1) ! fraction of cat 1 ridging 
            virdg(n1) = vin(n1) * afrac
            vsrdg(n1) = vsn(n1) * afrac
            do k = 1, nilyr
               eirdg(k,n1) = ein(k,n1) * afrac
            enddo
            esrdg(n1) = esn(n1) * afrac
            ardgtot = ardgtot + ardg(n1)
            anewtot = anewtot + anew(n1)
            vrdgtot = vrdgtot + virdg(n1)
          else                ! athornn(n1) = 0
            ardg(n1) = c0
          endif
        enddo

      !-----------------------------------------------------------------
      ! For each ridging category n1:
      ! (1) Remove area, volume and energy from n1.
      ! (2) Compute the fraction of the area and volume going to each
      !     category n2.  Assume Hibler (1980) redistribution function.
      ! (3) Add area, volume, and energy to n2.
      !-----------------------------------------------------------------

        do n1 = 1, ncat
          if (ardg(n1) .gt. puny) then   

            ! Remove area, volume, and energy from ridging category.
            ain(n1) = ain(n1) - ardg(n1)
            vin(n1) = vin(n1) - virdg(n1)
            vsn(n1) = vsn(n1) - vsrdg(n1)
            do k = 1, nilyr
              ein(k,n1) = ein(k,n1) - eirdg(k,n1)
            enddo
            esn(n1) = esn(n1) - esrdg(n1)

            ! Place part of the snow lost by ridging into the ocean.
            msnow_ocn = rhos*vsrdg(n1)*(c1-fsnowrdg)
            esnow_ocn =      esrdg(n1)*(c1-fsnowrdg)
            Fresh(i,j) = Fresh(i,j) + msnow_ocn/dt 
                      ! fresh water source for ocean
            Fhnet(i,j) = Fhnet(i,j) - esnow_ocn/dt 
                      ! esnow_ocn > 0, so ocean loses energy to melt snow
            
            ! Compute the fraction of ridged ice area and volume going to 
            ! each thickness category, and transfer accordingly.
            dhr = hrmaxn(n1) - hrminn(n1)
            dhr2 = hrmaxn(n1)*hrmaxn(n1) - hrminn(n1)*hrminn(n1)
               
            do n2 = 1, ncat
                  
              if (hrminn(n1).ge.hin_max(n2) .or. 
     $           hrmaxn(n1).le.hin_max(n2-1)) then
                hL = c0
                hR = c0
              else
                hL = max(hrminn(n1),hin_max(n2-1))
                hR = min(hrmaxn(n1),hin_max(n2))
              endif
               
              if (hR-hL .gt. puny) then
                farea = (hR-hL) / dhr 
                ! fraction of ridged ice in cat n1 that moves to cat n2
                ain(n2) = ain(n2) + farea*anew(n1)  

                fvol = (hR*hR - hL*hL) / dhr2 ! associated ice volume
                vin(n2) = vin(n2) + fvol*virdg(n1)
                vsn(n2) = vsn(n2) + fvol*vsrdg(n1)*fsnowrdg ! move snow 
                do k = 1,nilyr ! and energy too
                  ein(k,n2) = ein(k,n2) + fvol*eirdg(k,n1)
                enddo
                esn(n2) = esn(n2) + fvol*esrdg(n1)*fsnowrdg
              endif         ! hR-hL > puny

            enddo            ! loop over n2 (new ridges)            
          endif               ! ardg(n1) > puny
        enddo                  ! loop over n1 (ridging categories)

        hin_max(ncat) = Hbtmp   ! restore original value

      !-----------------------------------------------------------------
      ! check volume and energy conservation
      !-----------------------------------------------------------------
        vice2 = c0
        eice2 = c0
        do n = 1, ncat
          vice2 = vice2 + vin(n)
          do k = 1, nilyr
            eice2 = eice2 + ein(k,n)
          enddo
        enddo

        if (abs(vice2-vice1) .gt. eps04) then 
          write(nu_diag,*) my_task,':',i,j,'Ice volume error, ridging'
          write(nu_diag,*) my_task,':',i,j,'vice1 =',vice1,
     &                                     'vice2 =',vice2
          write(nu_diag,*) my_task,':',i,j,'Difference =', vice1 - vice2
          do n=1,ncat
          write(nu_diag,*) my_task,':',i,j,n,'ain =',ain(n),
     &                                  'vin =',vin(n)
          enddo
c          stop
         endif
        if (abs(eice2-eice1) .gt. eps04) then 
          write(nu_diag,*) my_task,':',i,j,'Ice energy error, ridging'
          write(nu_diag,*) my_task,':',i,j,'eice1 =', eice1
          write(nu_diag,*) my_task,':',i,j,'eice2 =', eice2
          write(nu_diag,*) my_task,':',i,j,'Difference =', eice1 - eice2
c          stop
        endif

        call from_column(i,j)  ! reload state variables

      endif                     ! ai0 < Gstar 

      enddo
      enddo

      !-----------------------------------------------------------------
      ! check that total ice area <= 1
      !-----------------------------------------------------------------

      call aggregate
      do j=jlo,jhi
        do i=ilo,ihi
          if (aice(i,j) .gt. c1) then
            if (aice(i,j) .gt. c1+puny) then
              write(nu_diag,*) i,j, 'After ridging, aice =', aice(i,j)
              stop
            else
              aice(i,j) = c1
            endif
          endif
          aice0(i,j) = c1 - aice(i,j)
        enddo
      enddo

      call ice_timer_stop(6)  ! ridging 

      end subroutine ridge_ice

c=======================================================================

      subroutine ice_strength(kstrngth)

!---!-------------------------------------------------------------------
!---! Compute strength of ice pack, prss, defined as the energy (J m-2) 
!---! dissipated per unit area removed from the ice pack under compression,
!---! and assumed proportional to the change in potential energy caused
!---! by ridging.  The procedure is similar to that of Hibler (1980).
!---! We compute the area of ice that ridges in each category when a unit
!---! area of ice plus open water is removed, then compute the change
!---! in potential energy resulting from this ridging, and finally
!---! multiply by a constant to account for frictional losses.
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) ::
     &   kstrngth    

      integer (kind=int_kind) ::
     &   i,j      
     &,  n                   ! thickness category index
      real (kind=dbl_kind) ::     
     &   dPE                 ! change in potential energy due to ridging
      real (kind=dbl_kind), parameter ::     
     &   Cp = p5 * gravit * (rhow-rhoi)*rhoi/rhow  ! g = 9.8...
     &,  Pstar = 2.75e4_dbl_kind ! for empirical strength formula

      call ridge_prep        ! compute ridging variables

      do j=jlo,jhi
      do i=ilo,ihi

      strength(i,j)=c0
      if (aice0(i,j).lt.puny) aice0(i,j) = c0

      if (kstrngth.eq.1) then
      if (aice(i,j).ge.c1-Gstar) then  

        call to_column(i,j)  ! load column variables from state variables
        athornn(0) = athorn(i,j,0)
        do n = 1, ncat
          athornn(n) = athorn(i,j,n)
          hrminn(n) = hrmin(i,j,n) 
          hrmaxn(n) = hrmax(i,j,n)
          krdgn(n) = krdg(i,j,n)
          aksumn = aksum(i,j) 
        enddo

        dPE = c0             ! change in potential energy due to ridging 
        do n = 1, ncat
          if (athornn(n) .gt. puny) then
            dPE = dPE - athornn(n)*(vin(n)/ain(n))**2 
                                          ! PE loss from ridging ice
            dPE = dPE + athornn(n)/(3.*krdgn(n)) * 
     &           (hrmaxn(n)**3 - hrminn(n)**3) / (hrmaxn(n)-hrminn(n))
                                          ! PE gain from new ridge
          endif
        enddo
        dPE = Cp * dPE / aksumn  ! Cp = (g/2) * (rhow-rhoi) * (rhoi/rhow)

        strength(i,j) = Cf * dPE ! ice strength P, kg/s^2
                                 ! Cf accounts for frictional dissipation
                                 ! See Flato and Hibler (1995)
      else  ! aice < 1-Gstar     ! ice is dispersed and therefore weak
        strength(i,j) = c0
      endif
      else  ! kstrngth ne 1:  Hibler (1979) form
        strength(i,j) = Pstar*vice(i,j)*exp(-c20*(c1-aice(i,j)))
      endif

      enddo
      enddo

      end subroutine ice_strength

c=======================================================================

      end module ice_mechred_cice

c=======================================================================
