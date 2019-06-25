c $Id: $
c=======================================================================
!---! Compute changes in ice/snow thickness and internal temperatures.
!---! Multiple ice layers, with heat capacity in each layer depending 
!---! on T and S.
!---!     
!---! author W. H. Lipscomb
!---!     
!---! See Bitz, C.M., and W.H. Lipscomb, 1999: 
!---! An energy-conserving thermodynamic model of sea ice,
!---! J. Geophys. Res., 104, 15,669-15,677. 
c=======================================================================

      module ice_therm_cice

      use ice_model_size
      use ice_domain 
      use ice_constants
      use ice_state
      use ice_history
      use ice_fileunits

      implicit none

      real (kind=dbl_kind), save ::
     &   rnilyr             !  real(nilyr)

      ! Error limits for temperature calculation. 
      ! Larger values make code faster. 
      ! Do not increase either value above 1.e-2, or results may change.
      real (kind=dbl_kind), parameter ::
     &   ferrmax = 1.e-3    ! max allowed energy flux error (W m-2)
     &,  Tsf_errmax = 5.e-4 ! max allowed error in Tsfc
     &,  hsnowmin = 1.0e-6  ! min thickness for which Tsnow computed (m)

      real (kind=dbl_kind) ::
     &   salin(nilyr+1)     ! salinity (ppt)
     &,  Tmlt(nilyr+1)      ! melting temp, -mu * salinity
     &,  ai                 ! ice area, aice
     &,  hi                 ! ice thickness, hice
     &,  hs                 ! snow thickness, hsnow
     &,  Tsf                ! surface temperature, Tsfc
     &,  qi(nilyr)          ! ice enthalpy, qice
     &,  qs                 ! snow enthalpy, qsnow
     &,  Tice(nilyr)        ! internal ice layer temperatures
     &,  Tsnow              ! internal snow temperature
     &,  fwnet              ! fw, with correction for surplus energy
     &,  fsh                ! surface downward sensible heat (W m-2)
     &,  flh                ! surface downward latent heat (W m-2)
     &,  flwup              ! upward LW at surface (W m-2)
     &,  fswsfc             ! SW absorbed at ice/snow surface (W m-2)
     &,  fswint             ! SW absorbed in ice (W m-2)
     &,  f0                 ! net flux to top of ice/snow, not inc. fct
     &,  fswocn             ! SW through ice to ocean (W m-2)
     &,  fswopen            ! SW absorbed in open water (W m-2)
     &,  fct                ! downward cond flux at top surface (W m-2)
     &,  fcb                ! downward cond flux at bottom surface (W m-2)
     &,  msub               ! ice mass sublimated/condensed (kg m-2) 
     &,  dfshdT             ! deriv of fsh wrt Tsfc (W m-2 deg-1)
     &,  dflhdT             ! deriv of flh wrt Tsfc (W m-2 deg-1)
     &,  dflwdT             ! deriv of flwup wrt Tsfc (W m-2 deg-1)

c=======================================================================

      contains

c=======================================================================

      subroutine init_therm_cice

!---!-------------------------------------------------------------------
!---! setup layers for each category
!---! setup the salinity profile and the melting temperature
!---! for each layer
!---!-------------------------------------------------------------------

      use ice_itd

      integer (kind=int_kind) :: layer,nc, nsum
      real (kind=dbl_kind) ::  zn       ! normalized sea ice thickness
      real (kind=dbl_kind), parameter ::
     &   nsal = 0.407_dbl_kind
     &,  msal = 0.573_dbl_kind

      rnilyr = real(nilyr)

      ! salinity and melting temperature profile
      do layer=1,nilyr
        zn=(real(layer)-p5)/rnilyr
        salin(layer)=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
c        salin(layer)=saltmax ! for isosaline ice
      enddo
      salin(nilyr+1)=saltmax
      do layer=1,nilyr+1
        Tmlt(layer)=-salin(layer)*depressT
      enddo

      end subroutine init_therm_cice

c=======================================================================

      subroutine thermo_cice(i,j,Fbot,
     &         dhi1,dhin,dhsn,hfocn,qin,qsn)

!---!-------------------------------------------------------------------
!---! heat budget over open water and ice 
!---! NOTE the wind stress is also computed here for later use
!---!-------------------------------------------------------------------

      use ice_itd
      use ice_flux
      use ice_calendar
      use ice_atmo
      use ice_diagnostics !! debugging

      integer (kind=int_kind), intent(in) :: i, j
      real (kind=dbl_kind), intent(in) ::
     &   Fbot

      real (kind=dbl_kind), intent(out) ::
     &   dhi1            ! melt at bottom and top plus sublimation  (m)
     &,  dhin(ncat)      ! ice thickness change (m)
     &,  dhsn(ncat)      ! snow thickness change (m)
     &,  hfocn           ! heat flx  to ice-snow from ocn  (J/m^2)
     &,  qin(nilyr,ncat) ! enthalpy per unit volume     (J/m^3)
     &,  qsn(ncat)       ! snow enthalpy per unit volume     (J/m^3)

      integer (kind=int_kind) :: 
     &   n               ! thickness category index
     &,  k               ! ice layer index

      real (kind=dbl_kind) ::
     &   hlyr            ! ice layer thickness
     &,  fwav            ! area-weighted average fw
     &,  dhs             ! new snow depth (m)
     &,  dvs             ! new snow volume (m)
     &,  qsnew           ! enthalpy of new snow (J kg-1) 
     &,  hsnew           ! new depth of snow (J kg-1) 
     &,  fsnowtot        ! snowfall rate (kg m-2 s-1)
     &,  einit           ! initial energy of melting (J m-2)
     &,  efinal          ! final energy of melting (J m-2)
     &,  einp            ! energy input during timestep (J m-2)
     &,  ferr            ! energy conservation error (W m-2)
     &,  ch(ncat)        ! heat transfer coefficient
     &,  Tfb             ! ocean freezing temperature (deg C)      
     &,  fsnon           ! flux of snow, atmos to ice (kg m-2 s-1)
     &,  frshn           ! flux of water, ice to ocean (kg m-2 s-1)
     &,  fsaltn          ! flux of salt, ice to ocean (kg m-2 s-1)
     &,  fsubn           ! flux of vapor, atmos to ice (kg m-2 s-1)
     &,   Trefn          ! air tmp rfrnc level                  (K)
     &,  strxn           ! air/ice zonal  strss,           (N/m^2)
     &,  stryn           ! air/ice merdnl strss,           (N/m^2)
     &,  rdn_in          ! initial value for rdn in stability routine
     &,  hi_init         ! initial value for hi, for each category
     &,  hs_init         ! initial value for hs, for each category
     &,  dummy           ! placeholder for stability call

      Tref(i,j) = c0 
      dhi1 = c0

      do n = 1,ncat

        dhin(n) = c0
        dhsn(n) = c0
        lhcoef(i,j,n) = c0
        shcoef(i,j,n) = c0
        fsnon = c0
        dhs = c0

        if (ain(n) .gt. puny) then

      !-----------------------------------------------------------------
      ! prep for air to ice heat, momentum, radiative and water fluxes 
      !-----------------------------------------------------------------

           rdn_in = vonkar/log(zref/iceruf)
           call stability( i,j,n, Tsfn(n), rdn_in,
     $         strxn   ,stryn   ,Trefn, dummy,dummy,dummy)

           fwnet = Fbot    ! ocean energy used by the ice < 0

      !-----------------------------------------------------------------
      ! Load variables used for thermo calculations.
      ! Note: ain(n) does not change and is not needed.
      !-----------------------------------------------------------------

           hi = vin(n) / ain(n)
           hs = vsn(n) / ain(n)
           hi_init = hi  ! save values
           hs_init = hs
           Tsf = Tsfn(n)
           Tfb = Tf(i,j)
           do k = 1, nilyr
             qi(k) = -ein(k,n)*rnilyr/vin(n) ! qi > 0
           enddo
           if (hs .gt. hsnowmin) then
             qs = -esn(n) / vsn(n)
           else
             qs = rhos * Lfresh
           endif

           einit = -hs*qs                  ! Initial energy per unit 
           hlyr = hi/nilyr                 ! ice area, relative to liquid 
           do k = 1, nilyr                 ! water at 0 C
             einit = einit - hlyr*qi(k) 
           enddo

      !-----------------------------------------------------------------
      ! New internal temperatures
      !-----------------------------------------------------------------

           call new_temps(i,j,n,Tfb)

      !-----------------------------------------------------------------
      ! Thickness change due to growth and/or melting
      ! Repartition into layers of equal thickness
      !-----------------------------------------------------------------

           call new_thickness(i,j,n,Tfb,Fbot,ain(n))

           ! final energy, including energy of sublimated/condensed ice
           efinal = -hs*qs
           hlyr = hi/nilyr
           do k = 1, nilyr
             efinal = efinal - hlyr*qi(k)
           enddo
           efinal = efinal - msub*Lvap

      !-----------------------------------------------------------------
      ! Check for energy conservation by comparing the change in energy 
      ! to the net energy input
      !-----------------------------------------------------------------
      ! Note that f0 - flh = fsw + flw + fsh; i.e., the latent heat
      ! flux is not included in the energy input;  de is the energy 
      ! change in the system ice + vapor, and the latent heat lost by 
      ! the ice is equal to that gained by the vapor.
      !-----------------------------------------------------------------

           einp = (f0 - flh + fswint - fwnet) * dt
           ferr = abs(efinal-einit-einp) / dt
           if (ferr .gt. ferrmax) then
             write(nu_diag,*) 'Energy error, cat',n,
     &                        'my_task,i,j',my_task,i,j
             write(nu_diag,*) 'Ferr =', ferr
             write(nu_diag,*) 'Eerr =', ferr*dt
             write(nu_diag,*) 'einit =', einit
             write(nu_diag,*) 'efinal =', efinal
             write(nu_diag,*) 'de =', efinal-einit
             write(nu_diag,*) 'einp =', einp
             write(nu_diag,*) 'f0 =', f0
             write(nu_diag,*) 'Fsw =', Fsw(i,j)
             write(nu_diag,*) 'Flw =', Flw(i,j)
             write(nu_diag,*) 'fsh =', fsh
             write(nu_diag,*) 'flh =', flh
             write(nu_diag,*) 'msub*Lvap =', msub*Lvap
             write(nu_diag,*) 'fswint =', fswint
             write(nu_diag,*) 'fwnet =', fwnet          
             write(nu_diag,*) 'Fbot =', Fbot
             call print_state('ice state at stop',i,j)
             stop
           endif

      !-----------------------------------------------------------------
      ! Let it snow
      !-----------------------------------------------------------------
           if (Tsfn(n).lt.-puny.and.Fsnow(i,j).gt.c0) then
c           if (Fsnow(i,j).gt.c0) then
             dhs = Fsnow(i,j)/rhos * dt
             qsnew = rhos * (Lfresh - cp_ice*min(Tair(i,j),c0))
             hsnew = hs + dhs
             if (hsnew.gt.puny) then
                qs = (hs*qs + dhs*qsnew)/hsnew
                qs = max(qs,rhos*Lfresh) ! avert roundoff errors if hs=0
                hs = hsnew
             else
                qs = c0
                hs = c0
             endif
           else ! convert snow to fresh water runoff
             Fresh(i,j) = Fresh(i,j) + Fsnow(i,j)*ain(n)
           endif

      !-----------------------------------------------------------------
      ! change in thicknesses, new enthalpy and surface temperature
      !-----------------------------------------------------------------
           dhin(n) = hi - hi_init
           dhsn(n) = hs - hs_init
           dhi1 = dhin(1)  ! used only for ncat=1 case

           Tsfn(n) = Tsf
           do k = 1, nilyr
             qin(k,n) = -qi(k)  ! note sign for consistency with ein
           enddo
           qsn(n) = -qs

      !-----------------------------------------------------------------
      ! flux of water and salt from ice to ocean for each category
      ! Note: These are fluxes per unit area of grid cell, not per
      ! unit area of ice
      !-----------------------------------------------------------------

        fsubn = -msub/dt ! msub>0 => sublimation, msub<0 => condensation
        frshn = fsubn - (rhoi*dhin(n) + rhos*(dhsn(n)-dhs))/dt 
        fsaltn =  - rhoi*dhin(n)*ice_ref_salinity*p001/dt 

        call merge_fluxes(i,j,n,frshn,fsaltn,fsubn,
     &      fsh,flh,flwup,fswocn,fwnet,Trefn,strxn,stryn)

        hfocn = fwnet*ain(n)

      !-----------------------------------------------------------------
      ! remove categories with very small areas
      !-----------------------------------------------------------------

        elseif (ain(n).le.puny) then
          ai0 = ai0 + ain(n)
          call zerocat(Tfb,ain(n),vin(n),vsn(n),ein(:,n),esn(n),
     &            Tsfn(n),hs,hi)

        endif  ! ain(n)

      enddo     ! Thermodynamic loop over n categories

      end subroutine thermo_cice

c=======================================================================

      subroutine new_temps(i,j,n,Tfb)

!---!-------------------------------------------------------------------
!---! Compute changes in ice, snow, and surface temperature for each category.
!---! Include effects of salinity on sea ice heat capacity and 
!---! conductivity, in a way that conserves energy.
!---! Finite differencing is backward implicit, 2nd order accurate in 
!---! space and 1st order in time.
!---!-------------------------------------------------------------------

      use ice_albedo
      use ice_calendar
      use ice_flux
      use ice_diagnostics
      use ice_atmo
      use ice_itd

      integer (kind=int_kind), intent(in) :: 
     &   i, j            ! grid cell indices
     &,  n               ! ice category number
      real (kind=dbl_kind), intent(in) ::
     &   Tfb             ! ocean freezing temperature (deg C)      

      integer (kind=int_kind), parameter :: 
     &   nitermax = 50   ! max number of iterations allowed

      integer (kind=int_kind) :: 
     &   k               ! ice level index
     &,  nmin            ! index of first matrix row
     &,  nmat            ! matrix dimension
     &,  niter           ! number of iterations to compute temps

      real (kind=dbl_kind), parameter ::
     &   alph = 3.       ! constant used to get 2nd order accurate fluxes
     &,  bet = -1./3.    ! constant used to get 2nd order accurate fluxes
     &,  alph_bet = alph+bet
     &,  betak  = 0.13   ! constant in formula for k (W m-1 ppt-1)
     &,  kimin  = 0.10   ! min conductivity of saline ice (W m-1 deg-1)
     &,  rhos_cp_ice = rhos*cp_ice ! convenient constant
     &,  i0vis = 0.7     ! fraction of penetrating solar rad (visible)

      real (kind=dbl_kind) ::
     &   fswabsv         ! fswabs in vis (wvlngth < 700nm)  (W/m^2)
     &,  fswabsi         ! fswabs in nir (wvlngth > 700nm)  (W/m^2)
     &,  flwdabs         ! down lw absorbed heat flx       (W/m^2)
     &,  hsn(ncat)       ! snow thickness for each cat       (m)
     &,  hlyr            ! ice layer thickness (m)
     &,  aa1, bb1, cc1   ! terms in quadratic formula
     &,  fswabs          ! net SW down at surface (W m-2)
     &,  fswpen          ! SW penetrating beneath surface (W m-2)
     &,  tran(0:nilyr)   ! transmitted frac of penetrating SW
     &,  Ipen(nilyr)     ! SW absorbed in particular layer (W m-2) 
     &,  frsnow          ! fractional snow coverage
     &,  df0dT           ! deriv of f0 wrt Tsf
     &,  ci              ! specific heat of ice (J kg-1 deg-1)
     &,  ki(nilyr+1)     ! thermal cond of ice (W m-1 deg-1)
     &,  kh(nilyr+1)     ! ki / hlyr
     &,  khs             ! ksno / hs
     &,  etas            ! dt / (rhos * cp_ice * hs)
     &,  eta(nilyr)      ! dt / (rhoi * ci * hlyr)
     &,  dt_rhoi_hlyr    ! dt/(rhoi*hlyr)
     &,  LTT(nilyr)      ! Lfresh*Tmlt/Tice_old (used to compute ci)
     &,  Tsnow_old       ! Tsnow at beginning of time step
     &,  Tice_old(nilyr) ! Tice at beginning of time step
     &,  Tsf_start       ! Tsf at start of iteration
     &,  Tsnow_start     ! Tsnow at start of iteration
     &,  Tice_start(nilyr)! Tice at start of iteration
     &,  diag(-1:nilyr)  ! diagonal matrix elements
     &,  sbdiag(-1:nilyr)! sub-diagonal matrix elements
     &,  spdiag(-1:nilyr)! super-diagonal matrix elements
     &,  rhs(-1:nilyr)   ! rhs of tri-diagonal matrix eqn.
     &,  spdiag2         ! term to right of superdiag on top row
     &,  Tmat(-1:nilyr)  ! matrix output temperatures
     &,  dTsf            ! Tsf - Tsf_start
     &,  dTsf_prev       ! dTsf_uncor from previous iteration
     &,  e1              ! total energy at beginning of time step
     &,  e2              ! total energy at end of iteration
     &,  ferr            ! energy error; if > ferrmax, keep iterating
     &,  dummy           ! place holder for subroutine call

      logical (kind=log_kind) ::
     &   converge        ! = true when acceptable solution found

      !-----------------------------------------------------------------
      ! total energy per area of ice (relative to 0 C) at start of 
      ! time step
      !-----------------------------------------------------------------

      e1 = -hs * qs          !!!! already computed above
      hlyr = hi/rnilyr
      do k = 1, nilyr
         e1 = e1 - hlyr * qi(k)
      enddo

      !-----------------------------------------------------------------
      ! Compute temperatures from enthalpies
      !-----------------------------------------------------------------

      Tsnow = (Lfresh - qs/rhos) / cp_ice
      do k = 1, nilyr
        if (Tmlt(k) .le. -puny) then
          aa1 = cp_ice
          bb1 = (cp_ocn-cp_ice)*Tmlt(k) + qi(k)/rhoi - Lfresh 
          cc1 = Lfresh * Tmlt(k)
          Tice(k) =  (-bb1 - sqrt(bb1*bb1 - 4.*aa1*cc1)) / (2.*aa1)
        else
          Tice(k) = (Lfresh - qi(k)/rhoi) / cp_ice
        endif
      end do

      !-----------------------------------------------------------------
      ! Bug check: Make sure initial temperatures do not exceed Tmlt and
      ! are not unreasonably low
      !-----------------------------------------------------------------

      if (Tsnow .gt. c0) then 
        if (Tsnow .gt. puny) then 
          write(nu_diag,*) ''
          write(nu_diag,*) 'Starting thermo, Ts > 0, cat',n
          write(nu_diag,*) 'Tsnow=',Tsnow
          write(nu_diag,*) 'istep1, my_task, i,j',istep1,my_task, i,j
          write(nu_diag,*) 'qs',qs
          call print_state('ice state at stop',i,j)
          stop        
        else  ! roundoff error
          Tsnow = c0
          qs = Lfresh*rhos
        endif
      endif
      do k = 1,nilyr
        if (Tice(k) .gt. Tmlt(k)) then 
          write(nu_diag,*) ''
          write(nu_diag,*) 'Starting thermo, T > Tmlt, cat',n,',layer',k
          write(nu_diag,*) 'Tice=',Tice(k),', Tmlt=',Tmlt(k)
          write(nu_diag,*) 'istep1, my_task, i,j',istep1,my_task, i,j
          write(nu_diag,*) 'qi',qi(k)
          call print_state('ice state at stop',i,j)
          stop
        endif
        if (Tice(k) .lt. -273._dbl_kind) then
          write(nu_diag,*) ''
          write(nu_diag,*) 'Starting thermo T < -273, cat',n,', layer',k
          write(nu_diag,*) 'Tice =', Tice(k)
          write(nu_diag,*) 'istep1, my_task, i,j',istep1,my_task, i,j
          call print_state('ice state at stop',i,j)
          stop
        endif
      enddo

      !-----------------------------------------------------------------
      ! Compute thermal conductivity at interfaces (held fixed during 
      ! the subsequent iteration);   conductivity must be >= kimin
      !-----------------------------------------------------------------

      ki(1) = kice + betak*salin(1) / min(-puny,Tice(1))
      do k = 2, nilyr
        ki(k) = kice + betak*p5*(salin(k-1)+salin(k)) /
     &        min (-puny, p5*(Tice(k-1)+Tice(k)))
      enddo
      ki(nilyr+1) = kice + betak*salin(nilyr+1)/Tfb 
      do k = 1, nilyr+1
        ki(k) = max (ki(k), kimin)
        kh(k) = ki(k) / hlyr
      enddo

      if (hs .gt. hsnowmin) then
        khs = ksno / hs
        kh(1) = 2.*kh(1)*khs / (kh(1)+khs)  ! snow-ice interface
      endif

      !-----------------------------------------------------------------
      ! transmittance for each ice layer 
      ! two-band, but transmit only visible (no infrared)
      !-----------------------------------------------------------------

      tran(0) = c1
      do k = 1, nilyr
        tran(k) = exp(-kappav * hlyr * real(k))
      enddo

      !-----------------------------------------------------------------
      ! snow thickness and fractional snow cover
      ! The snow fraction declines linearly from 1 to 0 as hsnow 
      ! declines from snowpatch to 0.
      !-----------------------------------------------------------------

      hsn(n) = vsn(n) / ain(n) 
      frsnow = hsn(n) / (hsn(n) + snowpatch)

      !-----------------------------------------------------------------
      ! shortwave flux absorbed at surface, absorbed internally,
      ! and penetrating to mixed layer.  
      !-----------------------------------------------------------------

      fswabsv  = swvdr(i,j)*(c1-alvdrn(i,j,n)) 
     &         + swvdf(i,j)*(c1-alvdfn(i,j,n))
      fswabsi  = swidr(i,j)*(c1-alidrn(i,j,n))
     &         + swidf(i,j)*(c1-alidfn(i,j,n))
      fswabs   = fswabsv + fswabsi

      fswpen = fswabsv * (c1-frsnow) * i0vis
c    &       + fswabsi * (c1-frsnow) * i0nir  ! i0nir = 0
      fswsfc = fswabs - fswpen
      do k = 1,nilyr
        Ipen(k) = fswpen * (tran(k-1)-tran(k))   
      enddo
      fswocn = fswpen * tran(nilyr)
      fswint  = fswpen - fswocn
      sabs(i,j) = sabs(i,j) + (fswsfc+fswint)*ain(n) ! for history diagnostics

      !-----------------------------------------------------------------
      ! Save initial temperatures
      !-----------------------------------------------------------------

      Tsnow_old = Tsnow
      do k = 1, nilyr
        Tice_old(k) = Tice(k)
      enddo

      !-----------------------------------------------------------------
      ! quantities that are used repeatedly in the matrix loop
      !-----------------------------------------------------------------

      dt_rhoi_hlyr = dt / (rhoi*hlyr)
      do k = 1, nilyr
        LTT(k) = Lfresh*Tmlt(k)/Tice_old(k)
      enddo
      if (hs .gt. hsnowmin) etas = dt / (rhos_cp_ice*hs)

      !-----------------------------------------------------------------
      ! Iterate until new temperatures converge
      !-----------------------------------------------------------------

      do niter = 1, nitermax      

        ! Update radiative/turbulent flux, f0, and derivative wrt Tsf.  
        call ice_sfc_flux(i,j,n,Tsf,flwdabs,flwup,
     &    fsh,flh,dflhdT,dfshdT,dflwdT)

        f0 = fswsfc + flwdabs + flwup + fsh + flh
        df0dT = dflwdT + dfshdT + dflhdT

        ! If Tsf = 0, make sure f0 > fct.  If not, then reset Tsf to 
        ! slightly less than zero (but not greater than -puny).
        if (Tsf .gt. -puny) then
          if (hs .gt. hsnowmin) then
            fct = c2 * khs * (Tsf-Tsnow)
          else
            fct = kh(1) * (alph*(Tsf-Tice(1)) + bet*(Tsf-Tice(2)))
          endif
          if (f0 .lt. fct) Tsf = -puny
        endif

        ! Save starting temperatures
        Tsf_start = Tsf
        Tsnow_start = Tsnow
        do k = 1, nilyr
          Tice_start(k) = Tice(k)
        enddo

        ! For each internal ice layer compute the specific heat ci, 
        ! a function of the starting temperature and of the current 
        ! guess for the final temperature.
        do k = 1, nilyr
          ci = cp_ice - LTT(k)/Tice(k)
          eta(k) = dt_rhoi_hlyr / ci
        enddo

      !-----------------------------------------------------------------
      !    Set up system of linear equations.
      !-----------------------------------------------------------------
      !    Four cases:
      !      (1) Cold surface (Tsf<0), snow present, (nilyr+2) matrix rows
      !      (2) Melting surface (Tsf=0), snow present, (nilyr+1) matrix rows
      !      (3) Cold surface (Tsf<0), no snow, (nilyr+1) matrix rows
      !      (4) Melting surface (Tsf=0), no snow, nilyr matrix rows
      !-----------------------------------------------------------------

        if (hs .gt. hsnowmin) then
          if (Tsf. le. -puny) then      ! case 1
            nmin = -1
               
            sbdiag(-1)  = c0
            spdiag(-1)  = c2 * khs                       
            diag(-1)    = df0dt - spdiag(-1)
            rhs(-1)     = df0dT*Tsf - f0

            sbdiag(0)   = -etas * c2 * khs
            spdiag(0)   = -etas * kh(1)
            diag(0)     = c1 + etas * (c2*khs + kh(1))
            rhs(0)      = Tsnow_old

          else                ! Tsf = 0, case 2
            nmin = 0
               
            sbdiag(0)   = c0
            spdiag(0)   = -etas * kh(1)
            diag(0)     = c1 + etas * (c2*khs + kh(1))
            rhs(0)      = Tsnow_old + etas*c2*khs*Tsf

          end if

          sbdiag(1) = -eta(1) * kh(1)
          spdiag(1) = -eta(1) * kh(2)
          diag(1)   = c1 - sbdiag(1) - spdiag(1)
          rhs(1)    = Tice_old(1) + eta(1)*Ipen(1)
            
        else                   ! no snow
            
          if (Tsf .le. -puny) then     ! case 3
            nmin = 0

            sbdiag(0)  = c0
            spdiag(0)  = alph * kh(1)                          
            spdiag2    = bet * kh(1)                              
            diag(0)    = df0dT -spdiag(0) - spdiag2
            rhs(0)     = df0dT*Tsf - f0
               
            sbdiag(1)  = -(alph+bet) * eta(1) * kh(1)
            spdiag(1)  = -eta(1) * (kh(2) - bet*kh(1))
            diag(1)    = c1 + eta(1) * (kh(2) + alph*kh(1))
            rhs(1)     = Tice_old(1) + eta(1)*Ipen(1)

            diag(0)   = spdiag(1)*diag(0) - spdiag2*sbdiag(1)
            spdiag(0) = spdiag(1)*spdiag(0) - spdiag2*diag(1)
            rhs(0)    = spdiag(1)*rhs(0) - spdiag2*rhs(1)  

          else                ! Tsf = 0, case 4
            nmin = 1
               
            sbdiag(1)  = c0
            spdiag(1)  = -eta(1) * (kh(2) - bet*kh(1))
            diag(1)    = c1 + eta(1) * (kh(2) + alph*kh(1))
            rhs(1)     = Tice_old(1) + eta(1)*Ipen(1) 
c     $    + (alph+bet)*eta(1)*kh(1)*Tsf ! necessary if Tsf were in Kelvin

          endif
        endif

        ! Ice interior
        if (nilyr .gt. 2) then
          do k = 2, nilyr-1
            sbdiag(k) =  -eta(k)*kh(k)
            spdiag(k) =  -eta(k)*kh(k+1)
            diag(k)   =  c1 - sbdiag(k) - spdiag(k)
            rhs(k)    =  Tice_old(k) + eta(k)*Ipen(k)
           enddo
        endif

        ! Bottom surface
        k = nilyr
        sbdiag(k)   = -eta(k) * (kh(k) - bet*kh(k+1)) 
        spdiag(k)   = c0
        diag(k)     = c1 + eta(k) * (kh(k) + alph*kh(k+1))
        rhs(k)      = Tice_old(k) + eta(k)*Ipen(k) + 
     $        (alph+bet)*eta(k)*kh(k+1)*Tfb

        nmat = nilyr - nmin + 1

      !-----------------------------------------------------------------
      ! Solve tridiagonal matrix to obtain the new temperatures.
      !-----------------------------------------------------------------
        call tridiag (sbdiag(nmin:nilyr), diag(nmin:nilyr), 
     &     spdiag(nmin:nilyr), rhs(nmin:nilyr), Tmat(nmin:nilyr), nmat)

        ! Reload temperatures
        if (hs .gt. hsnowmin) then
          if (Tsf .le. -puny) Tsf = Tmat(-1)
          Tsnow = Tmat(0)
        else   ! hs <= hsnowmin
          if (Tsf .le. -puny) Tsf = Tmat(0)
        endif
        do k = 1, nilyr
          Tice(k) = Tmat(k)
        end do

      !-----------------------------------------------------------------
      !    Convergence test.  Five conditions must hold.
      !    (1) Tsf <= 0 C.
      !    (2) Tsf is not oscillating; i.e., if both dTsf(niter) and
      !        dTsf(niter-1) have magnitudes greater than puny, then
      !        dTsf(niter)/dTsf(niter-1) cannot be a negative number
      !        with magnitude greater than 0.5.  
      !    (3) abs(dTsf) < Tsf_errmax
      !    (4) If Tsf = 0 C, then the downward turbulent/radiative 
      !        flux, f0, must be greater than or equal to the downward
      !        conductive flux, fct.
      !    (5) The net energy added to the ice per unit time must equal 
      !        the net change in internal ice energy per unit time,
      !        within the prescribed error ferrmax.
      !    When all five conditions are satisfied, leave the loop.
      !-----------------------------------------------------------------

        converge = .true.
        dTsf = Tsf - Tsf_start

        ! If Tsf > 0, set Tsf = 0 and reset the ice and snow 
        ! temperatures to midway between their starting and
        ! new (or melting) values.

        if (Tsf .gt. puny) then                   ! condition (1)
          Tsf = c0
          dTsf = -Tsf_start
          Tsnow = p5 * (Tsnow_start + min (Tsnow, c0))
          do k = 1, nilyr
            Tice(k) = p5 * (Tice_start(k) + min(Tice(k),Tmlt(k)))
          enddo
          converge = .false.
          go to 500
        endif

        ! If Tsf is oscillating, reset the temperatures to midway
        ! between their starting and new values.

        if (niter.gt.1 .and. Tsf_start.le.-puny .and. 
     &       abs(dTsf).gt.puny .and. abs(dTsf_prev).gt.puny) then 
           if (-dTsf/dTsf_prev .gt. p5) then      ! condition (2)
              Tsf = p5 * (Tsf_start + Tsf)
              dTsf = p5 * dTsf
              Tsnow = p5 * (Tsnow_start + Tsnow)
              do k = 1, nilyr
                 Tice(k) = p5 * (Tice_start(k) + Tice(k))
              enddo
              converge = .false.
              go to 500
           endif
        endif

        ! Check remaining conditions

        if (abs(dTsf) .gt. Tsf_errmax) then       ! condition (3)
          converge = .false.
          go to 500
        endif

        f0 = f0 + dTsf*df0dT
        if (hs .gt. hsnowmin) then
          fct = c2 * khs * (Tsf-Tsnow)
        else
          fct = kh(1) * (alph*(Tsf-Tice(1)) + bet*(Tsf-Tice(2)))
        endif
        if (Tsf.gt.-puny .and. f0.lt.fct) then    ! condition (4)
           converge = .false.
           go to 500
        endif

        fcb = kh(nilyr+1) * 
     $        (alph*(Tice(nilyr)-Tfb) + bet*(Tice(nilyr-1)-Tfb))
        qs = rhos * (Lfresh - cp_ice*Tsnow)
        e2 = -hs * qs
        do k = 1, nilyr
          qi(k) = rhoi * (cp_ice*(Tmlt(k)-Tice(k)) 
     $           + Lfresh*(c1-Tmlt(k)/Tice(k)) - cp_ocn*Tmlt(k)) 
          e2 = e2 - hlyr*qi(k)
        enddo
        ferr = abs( (e2-e1)/dt - (fct - fcb + fswint) )
        if (ferr .gt. 0.9*ferrmax) then           ! condition (5)
          ! factor of 0.9 allows for further errors by other subroutines
          converge = .false.
        endif
         
 500    if (converge) then
           exit                 ! solution is acceptable; leave loop
        else
           dTsf_prev = dTsf
        endif

      enddo  ! temperature iteration

      if (.not.converge) then
      write(nu_diag,*) 'Finished loop without converging,',' cat',n
      write(nu_diag,*) 'dTsf, Tsf_errmax',dTsf, Tsf_errmax
      write(nu_diag,*) 'Tsf, f0, fct',Tsf, f0, fct
      write(nu_diag,*) 'Flux conservation error =', ferr
      write(nu_diag,*) 'snow and ice layer temperatures'
      write(nu_diag,*) Tsnow,(Tice(k),k=1,nilyr)
      write(nu_diag,*) 'starting temperatures'
      write(nu_diag,*) Tsnow_old,(Tice_old(k),k=1,nilyr)
      call print_state('ice state at stop',i,j)
      stop
      endif

      !-----------------------------------------------------------------
      ! Update fluxes that depend on Tsf.
      !-----------------------------------------------------------------

      flwup = flwup + dTsf*dflwdT
      fsh = fsh + dTsf*dfshdT
      flh = flh + dTsf*dflhdT

      end subroutine new_temps

c=======================================================================

      subroutine tridiag (a, b, c, r, x, N)

!---!-------------------------------------------------------------------
!---! Tridiagonal matrix solver.
!---! The expression to be solved is of the form Ax = r, where A is an 
!---! N x N matrix and x and r are vectors of length N. We have 'a' denoting
!---! the subdiagonal of the matrix, 'b' the main diagonal, and 'c' the 
!---! superdiagonal.  The unknown vector is x.
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: 
     &   N               ! dimension of matrix
      real (kind=dbl_kind), intent(in) ::
     &   a(nilyr+2)      ! subdiagonal
     &,  b(nilyr+2)      ! diagonal
     &,  c(nilyr+2)      ! superdiagonal
     &,  r(nilyr+2)      ! right-hand side
      real (kind=dbl_kind), intent(inout) ::
     &   x(nilyr+2)      ! unknown vector

      integer (kind=int_kind) :: 
     &   k               ! counting index
      real (kind=dbl_kind) ::
     &   bet             ! temporary variable
     &,  gam(nilyr+2)    ! temporary variable

      bet = b(1)
      x(1) = r(1) / bet
      do k = 2, N
        gam(k) = c(k-1)/bet
        bet = b(k) - a(k)*gam(k)
        x(k) = (r(k) - a(k)*x(k-1)) / bet
      enddo
      do k = N-1, 1, -1
        x(k) = x(k) - gam(k+1)*x(k+1)
      enddo

      end subroutine tridiag

c=======================================================================

      subroutine new_thickness(i,j,n,Tfb,Fbot,ai)

!---!-------------------------------------------------------------------
!---! Compute growth and/or melting at the top and bottom surfaces.
!---! Then repartition into equal-thickness layers, conserving energy.
!---!-------------------------------------------------------------------

      use ice_calendar

      integer (kind=int_kind), intent(in) :: 
     &   n             ! ice category number
     &,  i,j           ! spatial indices
      real (kind=dbl_kind), intent(in) ::
     &   Fbot
     &,  Tfb           ! ocean freezing temperature (deg C)      
     &,  ai            ! category area

      integer (kind=int_kind) :: 
     &   n2            ! thickness category index
     &,  k             ! ice level index
      real (kind=dbl_kind) ::
     &   hlyr          ! ice layer thickness (m)
     &,  hnew(nilyr)   ! new ice layer thicknesses (m)
     &,  esub          ! energy for sublimation/condensation (J m-2)
     &,  etop          ! energy for top melting (J m-2)
     &,  ebot          ! energy for bottom melting (J m-2)
     &,  esurp         ! energy left over when all ice is melted (J m-2)
     &,  hssub         ! depth of condensed/sublimated snow (m)
     &,  hisub         ! depth of condensed/sublimated ice (m)
     &,  emax          ! energy needed to melt ice or snow layer (J m-2)
     &,  qbot          ! enthalpy of ice growing at bottom surface (J m-3)
     &,  qsub          ! energy needed to sublimate ice/snow (J m-3) 
     &,  dhi           ! change in ice thickness
     &,  dhs           ! change in snow thickness
     &,  hsmlt         ! reduction of hs from top melting
     &,  hqtot         ! sum of hq for two layers

      !-----------------------------------------------------------------
      ! Initialize 
      !-----------------------------------------------------------------

      hsmlt = c0
      msub = c0
      hlyr = hi/rnilyr
      do k = 1, nilyr
        hnew(k) = hlyr
      enddo

      !-----------------------------------------------------------------
      ! energy available for melting/growth
      !-----------------------------------------------------------------

      esub = -flh * dt   

      if (abs(Tsf) .lt. puny) then
        etop = (f0 - fct) * dt
      else
        etop = c0
      endif

      if (etop .lt. -puny) then
      write(nu_diag,*) ''
      write(nu_diag,*) 'etop < 0, cat',n
      write(nu_diag,*) 'etop =', etop
      write(nu_diag,*) 'f0 =', f0
      write(nu_diag,*) 'fct =', fct
      stop
      endif
      
      ebot = (fcb - Fbot) * dt

      !-----------------------------------------------------------------
      ! Sublimation/condensation of ice/snow
      !-----------------------------------------------------------------

      if (esub .gt. puny) then   ! Ice or snow sublimates
                                 ! Note msub < 0 for sublimation
        if (hs .gt. puny) then
          qsub = qs + rhos*Lvap
          emax = qsub * hs
          if (esub .lt. emax) then
            hssub = esub/qsub
            hs = hs - hssub
            esub = c0
          else
            esub = esub - emax
            hssub = hs
            hs = c0
          endif
          msub = msub - hssub * rhos
        endif

        do k = 1, nilyr
          if (esub .gt. puny) then
            qsub = qi(k) + rhoi*Lvap
            emax = qsub * hnew(k)
            if (esub .lt. emax) then
              hisub = esub/qsub
              hnew(k) = hnew(k) - hisub
              esub = c0
            else
              esub = esub - emax
              hisub = hnew(k)
              hnew(k) = c0
            endif
            msub = msub - hisub * rhoi
          else
            exit
          endif
        enddo

        if (esub .gt. puny) then
          ! ice and snow are gone and melting energy remains
          hi = c0
          go to 200
        endif

      elseif (esub .lt. -puny) then   ! Water deposited onto ice

        if (hs .gt. puny) then   ! Add snow with enthalpy qs
          hssub = -esub / (qs + rhos*Lvap)     
          hs = hs + hssub
          msub = msub + hssub * rhos
        else                     ! Add ice with enthalpy qi(1)
          hisub = -esub / (qi(1) + rhoi*Lvap)  
          hnew(1) = hnew(1) + hisub
          msub = msub + hisub * rhoi
        endif
        esub = c0

      endif  ! initial esub > puny

      !-----------------------------------------------------------------
      ! Top melt
      !-----------------------------------------------------------------

      if (etop .gt. puny) then

        if (hs .gt. puny) then
          emax = qs * hs
          if (etop .lt. emax) then
            hsmlt = etop/qs                  
            hs = hs - hsmlt
            etop = c0
          else
            etop = etop - emax
            hsmlt = hs
            hs = c0
            qs = c0
          endif
        endif
      
        do k = 1, nilyr
          if (etop .gt. puny) then
            emax = qi(k) * hnew(k)
            if (etop .lt. emax) then
              dhi = etop/qi(k)
              meltt(i,j) = meltt(i,j) + dhi*ai ! for history diagnostics
              hnew(k) = hnew(k) - dhi
              etop = c0
            else
              etop = etop - emax
              meltt(i,j) = meltt(i,j) + hnew(k)*ai ! for history diagnostics
              hnew(k) = c0
            endif
          else
            exit
          endif
        enddo

        if (etop .gt. puny) then
          ! ice and snow are gone and melting energy remains
          hi = c0
          go to 200
        endif

      endif   ! initial etop > puny


      !-----------------------------------------------------------------
      ! Bottom melt/growth
      !-----------------------------------------------------------------

      if (ebot .lt. -puny) then

        ! enthalpy of new ice growing at bottom surface
        qbot = rhoi * (cp_ice*(Tmlt(nilyr+1)-Tfb)
     $        + Lfresh*(c1-Tmlt(nilyr+1)/Tfb) - cp_ocn*Tmlt(nilyr+1)) 
        dhi = -ebot / qbot
        ebot = c0
        k = nilyr
        hqtot = hnew(k)*qi(k) + dhi*qbot
        hnew(k) = hnew(k) + dhi
        qi(k) = hqtot / hnew(k)
        congel(i,j) = congel(i,j) + dhi*ai ! for history diagnostics

      elseif (ebot .gt. puny) then   ! Ice melts at bottom surface

        do k = nilyr, 1, -1
          if (ebot.gt.puny .and. hnew(k).gt.puny) then
            emax = qi(k) * hnew(k)
            if (ebot .lt. emax) then
              dhi = ebot/qi(k)
              meltb(i,j) = meltb(i,j) + dhi*ai ! for history diagnostics
              hnew(k) = hnew(k) - dhi
              ebot = c0
            else
              ebot = ebot - emax
              meltb(i,j) = meltb(i,j) + hnew(k)*ai ! for history diagnostics
              hnew(k) = c0
            endif
          else
            exit
          endif
        enddo

        if (ebot .gt. puny) then  ! ice is gone but snow may be left

           if (hs .gt. puny) then
              emax = qs * hs
              if (ebot .lt. emax) then
                 hs = hs - ebot/qs
                 ebot = c0
                 ! convert remaining snow to ice
                 hnew(1) = hs * rhos/rhoi ! reduce thickness
                 qi(1)   = qs * rhoi/rhos ! increase q to conserve energy
                 hs = c0
                 qs = c0
              else
                 ebot = ebot - emax
                 hs = c0 
                 ! ice and snow are gone and melting energy remains
                 hi = c0
                 go to 200  
              endif
           else  ! no snow
              ! ice and snow are gone and melting energy remains
              hi = c0
              go to 200  
           endif

        endif

      endif   ! initial ebot < -puny


      !-----------------------------------------------------------------
      ! new total ice thickness
      !-----------------------------------------------------------------

      hi = c0
      do k = 1, nilyr
        hi = hi + hnew(k)
      enddo

      !-----------------------------------------------------------------
      ! Repartition into layers of equal thickness
      !-----------------------------------------------------------------

      if (hi .gt. puny) then
         call new_layers (hnew)
      else
         hi = c0
      endif

      !-----------------------------------------------------------------
      ! surplus energy left over from melting
      ! amount used is sent to ocean model
      !-----------------------------------------------------------------
            
 200  esurp = esub + etop + ebot    ! esurp is non-negative
      fwnet = fwnet + esurp/dt      ! amount of heat used by ice
          ! positive fwnet means energy goes from ice to ocean
          ! negative fwnet means ice used energy from ocean

      end subroutine new_thickness

c=======================================================================

      subroutine new_layers (hnew)

!---!-------------------------------------------------------------------
!---! Switch to ice layers of uniform thickness, conserving enthalpy.
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::
     &   hnew(nilyr)    ! new ice layer thicknesses (m)

      integer (kind=int_kind) :: 
     &   k,m,mstart     ! layer indices
      real (kind=dbl_kind) ::
     &   hlyr           ! ice layer thickness (m)
     &,  hq(nilyr)      ! h * q for a layer
     &,  ztmp(0:nilyr)  ! depth of layer boundaries (m)
     &,  zadj(0:nilyr)  ! adjusted depths, with equal hlyr (m)
     &,  hovlp          ! overlap between old and new layers (m)

      hlyr = hi/rnilyr

      ztmp(0) = c0
      zadj(0) = c0
      do k = 1, nilyr-1
        ztmp(k) = ztmp(k-1) + hnew(k)
        zadj(k) = zadj(k-1) + hlyr
      end do
      ztmp(nilyr) = hi
      zadj(nilyr) = hi

      ! Compute overlap between old layers (m) and new layers (k).
      ! Start with the first m layer of non-zero thickness.
      mstart = 1
 100  if (ztmp(mstart) .lt. puny) then
        mstart = mstart + 1
        go to 100
      endif
      
      do k = 1, nilyr
        hq(k) = c0
        do m = mstart, nilyr
          hovlp = min(ztmp(m),zadj(k)) - max(ztmp(m-1),zadj(k-1)) 
          if (hovlp .gt. c0) then
            hq(k) = hq(k) + hovlp * qi(m)
          else
            mstart = m-1
            go to 150
          endif
          if (m .eq. nilyr) mstart = nilyr
        enddo
 150  enddo

      ! enthalpy of new layers
      do k = 1, nilyr
        qi(k) = hq(k) / hlyr
      enddo

      end subroutine new_layers

c=======================================================================

      end module ice_therm_cice

c=======================================================================

