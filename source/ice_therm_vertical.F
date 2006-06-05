!=========================================================================
!BOP
!
! !MODULE: ice_therm_vertical - thermo calculations before call to coupler
!
! !DESCRIPTION:
!
! Update ice and snow internal temperatures and compute
! thermodynamic growth rates and atmospheric fluxes.
!
! NOTE: The thermodynamic calculation is split in two for load balancing.
!       First ice_therm_vertical computes vertical growth rates and coupler
!       fluxes.  Then ice_therm_itd does thermodynamic calculations not
!       needed for coupling.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW
!          Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
!
! !INTERFACE:
!
      module ice_therm_vertical
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size, only: ncat, nilyr, nslyr, ntilyr, ntslyr
      use ice_constants
      use ice_fileunits, only: nu_diag
!
!EOP
!
      implicit none
      save

      real (kind=dbl_kind) ::
     &   salin(nilyr+1)  ! salinity (ppt)
     &,  Tmlt(nilyr+1)   ! melting temp, -mu * salinity

      real (kind=dbl_kind), parameter, private ::
     &   ferrmax = 1.0e-3_dbl_kind  ! max allowed energy flux error (W m-2)
                                    ! recommend ferrmax < 0.01 W m-2
     &,  hsnomin = 1.0e-6_dbl_kind  ! min thickness for which Tsno computed (m)

      character (char_len) :: stoplabel

      logical (kind=log_kind) ::
     &   l_brine         ! if true, treat brine pocket effects

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: thermo_vertical - driver for pre-coupler thermodynamics
!
! !DESCRIPTION:
!
! Driver for updating ice and snow internal temperatures and
! computing thermodynamic growth rates and atmospheric fluxes.
!
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW
!
! !INTERFACE:


      subroutine thermo_vertical (nx_block,    ny_block,
     &                            dt,          icells,
     &                            indxi,       indxj,
     &                            aicen,       Tsfcn,
     &                            vicen,       vsnon,
     &                            eicen,       esnon,
     &                            flw,         potT,
     &                            Qa,          rhoa,
     &                            fsnow,
     &                            fbot,        Tbot,
     &                            lhcoef,      shcoef,
     &                            fswsfc,      fswint,
     &                            fswthrun,    Iswabs,
     &                            fsensn,      flatn,
     &                            fswabsn,     flwoutn,
     &                            evapn,       freshn,
     &                            fsaltn,      fhocnn,
     &                            meltt,       meltb,
     &                            congel,      snoice,
     &                            mlt_onset,   frz_onset,
     &                            yday,        abort_flag)

! 
! !USES:
!
      use ice_communicate, only: my_task, master_task
      use ice_calendar, only: istep1
      use ice_diagnostics, only: print_state
      use ice_exit
      use ice_ocean
      use ice_itd, only: ilyr1, slyr1, ilyrn, slyrn
      use ice_work, only: worka, workb
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  icells              ! number of cells with ice present

      integer (kind=int_kind), dimension (nx_block*ny_block),
     &   intent(in) ::
     &   indxi, indxj     ! compressed indices for cells with ice

      real (kind=dbl_kind), intent(in) ::
     &   dt      ! time step

      ! ice state variables
      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   aicen     ! concentration of ice
     &,  vicen     ! volume per unit area of ice          (m)
     &,  vsnon     ! volume per unit area of snow         (m)
     &,  Tsfcn     ! temperature of ice/snow top surface  (C)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr),
     &   intent(inout) ::
     &   eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nslyr),
     &   intent(inout) ::
     &   esnon     ! energy of melting for each snow layer (J/m^2)

      ! input from atmosphere
      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   flw         ! incoming longwave radiation (W/m^2)
     &,  potT        ! air potential temperature  (K) 
     &,  Qa          ! specific humidity (kg/kg) 
     &,  rhoa        ! air density (kg/m^3) 
     &,  fsnow       ! snowfall rate (kg m-2 s-1)
     &,  shcoef      ! transfer coefficient for sensible heat
     &,  lhcoef      ! transfer coefficient for latent heat
     &,  fswsfc      ! SW absorbed at ice/snow surface (W m-2)
     &,  fswint      ! SW absorbed in ice interior, below surface (W m-2)
     &,  fswthrun    ! SW through ice to ocean         (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(in) ::
     &   Iswabs      ! SW radiation absorbed in ice layers (W m-2)

      ! input from ocean
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   fbot        ! ice-ocean heat flux at bottom surface (W/m^2)
     &,  Tbot        ! ice bottom surface temperature (deg C)

      ! coupler fluxes to atmosphere
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out)::
     &   fsensn      ! sensible heat flux (W/m^2) 
     &,  flatn       ! latent heat flux   (W/m^2) 
     &,  fswabsn     ! shortwave flux absorbed in ice and ocean (W/m^2) 
     &,  flwoutn     ! outgoing longwave radiation (W/m^2) 
     &,  evapn       ! evaporative water flux (kg/m^2/s) 

      ! coupler fluxes to ocean
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out)::
     &   freshn      ! fresh water flux to ocean (kg/m^2/s)
     &,  fsaltn      ! salt flux to ocean (kg/m^2/s)
     &,  fhocnn      ! net heat flux to ocean (W/m^2) 

      ! diagnostic fields
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out)::
     &   meltt     ! top ice melt             (m/step-->cm/day) 
     &,  meltb     ! basal ice melt           (m/step-->cm/day) 
     &,  congel    ! basal ice growth         (m/step-->cm/day) 
     &,  snoice    ! snow-ice formation       (m/step-->cm/day) 
     &,  mlt_onset ! day of year that sfc melting begins 
     &,  frz_onset ! day of year that freezing begins (congel or frazil) 

      real (kind=dbl_kind), intent(in) ::
     &   yday      ! day of year

      logical (kind=log_kind), intent(out) ::
     &   abort_flag  ! if true, abort on return


!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  ilo,ihi,jlo,jhi ! beginning and end of physical domain
     &,  n               ! thickness category index
     &,  k               ! ice layer index
     &,  il1, il2        ! ice layer indices for eice
     &,  sl1, sl2        ! snow layer indices for esno

      real (kind=dbl_kind) ::
     &   dhi             ! change in ice thickness
     &,  dhs             ! change in snow thickness

! 2D state variables (thickness, temperature, enthalpy)

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   hin             ! ice thickness (m)
     &,  hsn             ! snow thickness (m)
     &,  hilyr           ! ice layer thickness
     &,  hslyr           ! snow layer thickness
     &,  hsn_new         ! thickness of new snow (m)
     &,  Tsf             ! ice/snow top surface temp, same as Tsfcn (deg C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr) ::
     &   qin             ! ice layer enthalpy, qin < 0 (J m-3)
     &,  Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) ::
     &   qsn             ! snow layer enthalpy, qsn < 0 (J m-3)
     &,  Tsn             ! internal snow layer temperatures

! other 2D flux and energy variables

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   fsurf           ! net flux to top surface, not including fcondtop
     &,  fcondtop        ! downward cond flux at top surface (W m-2)
     &,  fcondbot        ! downward cond flux at bottom surface (W m-2)
     &,  einit           ! initial energy of melting (J m-2)
     &,  efinal          ! final energy of melting (J m-2)

      logical (kind=log_kind) ::
     &   l_stop          ! if true, print diagnostics and abort model

      integer (kind=int_kind) ::
     &   istop, jstop    ! i and j indices of cell where model fails


      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.
            
      fsensn (:,:) = c0
      flatn  (:,:) = c0
      fswabsn(:,:) = c0
      flwoutn(:,:) = c0
      evapn  (:,:) = c0

      freshn (:,:) = c0
      fsaltn (:,:) = c0
      fhocnn (:,:) = c0

      !-----------------------------------------------------------------
      ! Compute variables needed for vertical thermo calculation
      !-----------------------------------------------------------------

      call init_vertical_profile (nx_block,  ny_block,
     &                            my_task,   istep1,
     &                            icells, 
     &                            indxi,     indxj,
     &                            aicen,
     &                            vicen,     vsnon,
     &                            Tsfcn,
     &                            eicen,     esnon,
     &                            hin,       hilyr,
     &                            hsn,       hslyr,
     &                            qin,       Tin,
     &                            qsn,       Tsn,
     &                            Tsf,       einit,
     &                            l_stop,
     &                            istop,     jstop)

      if (l_stop) then
         i = istop
         j = jstop
!!            stoplabel = 'ice state at stop'
!!            call print_state (stoplabel,i,j)
         call abort_ice('ice thermo: bad temperature')
      endif

      ! Save initial ice and snow thickness (for fresh and fsalt)
      worka(:,:) = hin(:,:)
      workb(:,:) = hsn(:,:)
            
      !-----------------------------------------------------------------
      ! Compute new surface temperature and internal ice and snow
      !  temperatures.
      !-----------------------------------------------------------------

      call temperature_changes (nx_block,      ny_block,
     &                          my_task,       istep1,
     &                          dt,            icells,   
     &                          indxi,         indxj,
     &                          rhoa,          flw,
     &                          potT,          Qa,
     &                          shcoef,        lhcoef,
     &                          fswsfc,      fswint,
     &                          fswthrun,    Iswabs,
     &                          hilyr,         hslyr,
     &                          qin,           Tin,
     &                          qsn,           Tsn,
     &                          Tsf,           Tbot,
     &                          fsensn,        flatn,
     &                          fswabsn,       flwoutn,
     &                          fsurf,
     &                          fcondtop,      fcondbot,
     &                          einit,         l_stop,
     &                          istop,         jstop)

      if (l_stop) then
         i = istop
         j = jstop
!!            stoplabel = 'ice state at stop'
!!            call print_state (stoplabel,i,j)
         call abort_ice('ice thermo: convergence error')
      endif

      !-----------------------------------------------------------------
      ! Compute growth and/or melting at the top and bottom surfaces.
      ! Add new snowfall.
      ! Repartition ice into equal-thickness layers, conserving energy.
      !-----------------------------------------------------------------

      call thickness_changes(nx_block,     ny_block,
     &                       dt,
     &                       yday,         icells,  
     &                       indxi,        indxj,
     &                       aicen,        efinal,
     &                       hin,          hilyr,
     &                       hsn,          hslyr,
     &                       qin,          qsn,
     &                       fbot,         Tbot,
     &                       flatn,        fsurf,
     &                       fcondtop,     fcondbot,
     &                       fsnow,        hsn_new,
     &                       fhocnn,       evapn,
     &                       meltt,        meltb,
     &                       congel,       snoice,
     &                       mlt_onset,    frz_onset)

      !-----------------------------------------------------------------
      ! Check for energy conservation by comparing the change in energy
      ! to the net energy input
      !-----------------------------------------------------------------

      call conservation_check_vthermo(nx_block, ny_block,
     &                                my_task,  istep1,
     &                                dt,       icells, 
     &                                indxi,    indxj,
     &                                fsurf,    flatn,
     &                                fhocnn,   fswint,
     &                                fsnow,
     &                                einit,    efinal,
     &                                l_stop,
     &                                istop,    jstop)

      if (l_stop) then
         i = istop
         j = jstop
!            stoplabel = 'ice state at stop'
!            call print_state (stoplabel,i,j)
         call abort_ice('ice thermo: energy conservation error')
      endif


      !-----------------------------------------------------------------
      ! Compute fluxes of water and salt from ice to ocean.
      ! evapn < 0 => sublimation, evapn > 0 => condensation
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
            
         dhi = hin(i,j) - worka(i,j)
         dhs = hsn(i,j) - workb(i,j)
               
         freshn(i,j) = evapn(i,j) -
     &                 (rhoi*dhi + rhos*(dhs-hsn_new(i,j))) / dt
         fsaltn(i,j) = -rhoi*dhi*ice_ref_salinity*p001/dt

      enddo                     ! ij


      !-----------------------------------------------------------------
      !  Given the vertical thermo state variables (hin, hsn, Tsf,
      !   qin, qsn,), compute the new ice state variables (vicen, vsnon,
      !   Tsfcn, eicen, esnon).
      !-----------------------------------------------------------------

      call update_state_vthermo(nx_block, ny_block,
     &                             icells, 
     &                             indxi,     indxj,
     &                             Tbot,      Tsf,         
     &                             hin,       hsn,
     &                             qin,       qsn,
     &                             aicen,  
     &                             vicen,     vsnon, 
     &                             Tsfcn,
     &                             eicen,     esnon)


      end subroutine thermo_vertical

!=======================================================================
!BOP
!
! !ROUTINE: init_thermo_vertical - initialize salinity and melting temp
!
! !DESCRIPTION:
!
! Initialize the vertical profile of ice salinity and melting temperature.
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine init_thermo_vertical
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind), parameter ::
     &   nsal    = 0.407_dbl_kind
     &,  msal    = 0.573_dbl_kind
     &,  saltmax = 3.2_dbl_kind    ! max salinity at ice base (ppt)
     &,  min_salin = 0.1_dbl_kind  ! threshold for brine pocket treatment 

      integer (kind=int_kind) :: k        ! ice layer index
      real (kind=dbl_kind)    :: zn       ! normalized ice thickness

      !-----------------------------------------------------------------
      ! Determine l_brine based on saltmax.
      ! Thermodynamic solver will not converge if l_brine is true and
      !  saltmax is close to zero.
      !-----------------------------------------------------------------

      if (saltmax > min_salin) then
         l_brine = .true.
      else
         l_brine = .false.
      endif

      !-----------------------------------------------------------------
      ! Prescibe vertical profile of salinity and melting temperature.
      !-----------------------------------------------------------------

      if (l_brine) then
         do k = 1, nilyr
            zn = (real(k)-p5) / real(nilyr)
            salin(k)=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
c            salin(k)=saltmax ! for isosaline ice
         enddo
         salin(nilyr+1) = saltmax
         do k = 1, nilyr+1
            Tmlt(k) = -salin(k)*depressT
         enddo
      else
         do k = 1, nilyr+1
            salin(k) = c0
            Tmlt(k) = c0
         enddo
      endif

      end subroutine init_thermo_vertical

!=======================================================================
!BOP
!
! !ROUTINE: frzmlt_bottom_lateral - bottom and lateral heat fluxes
!
! !DESCRIPTION:
!
! Adjust frzmlt to account for changes in fhocn since from_coupler.
! Compute heat flux to bottom surface.
! Compute fraction of ice that melts laterally.
!
! !REVISION HISTORY:
!
! authors C. M. Bitz, UW
!         William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine frzmlt_bottom_lateral (nx_block, ny_block,
     &                                  nghost,   dt,
     &                                  aice,     frzmlt,
     &                                  eicen,    esnon, 
     &                                  sst,      Tf,  
     &                                  strocnxT, strocnyT, 
     &                                  Tbot,     fbot,   
     &                                  rside)
!
! !USES:
!
      use ice_itd, only: ilyr1, slyr1

! !INPUT/OUTPUT PARAMETERS:

      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), intent(in) ::
     &   dt                  ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::
     &   aice        ! ice concentration
     &,  frzmlt      ! freezing/melting potential (W/m^2)
     &,  sst         ! sea surface temperature (C)
     &,  Tf          ! freezing temperature (C)
     &,  strocnxT    ! ice-ocean stress, x-direction
     &,  strocnyT    ! ice-ocean stress, y-direction

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntilyr),
     &   intent(in) ::
     &   eicen       ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntslyr),
     &   intent(in) ::
     &   esnon       ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block),
     &   intent(out) ::
     &   Tbot        ! ice bottom surface temperature (deg C)
     &,  fbot        ! heat flux to ice bottom  (W/m^2)
     &,  rside       ! fraction of ice that melts laterally
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j             ! horizontal indices
     &,  n                ! thickness category index
     &,  k                ! layer index
     &,  ilo,ihi,jlo,jhi  ! beginning and end of physical domain
     &,  ij               ! horizontal index, combines i and j loops
     &,  icells           ! number of cells with ice melting

      integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &   indxi, indxj     ! compressed indices for cells with ice melting

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   etot        ! total energy in column
     &,  fside       ! lateral heat flux (W/m^2)

      real (kind=dbl_kind) ::
     &   deltaT      ! SST - Tbot >= 0
     &,  ustar       ! skin friction velocity for fbot (m/s)
     &,  fhocn_init  ! initial value of fhocn (W/m^2)
     &,  wlat        ! lateral melt rate (m/s)
     &,  xtmp        ! temporarty variable

      ! Parameters for bottom melting

      ! 0.006 = unitless param for basal heat flx ala McPhee and Maykut

      real (kind=dbl_kind), parameter ::
     &   cpchr = -cp_ocn*rhow*0.006_dbl_kind

#ifdef coupled
      real (kind=dbl_kind), parameter ::
     &   ustar_min = 5.e-3_dbl_kind
#endif


#ifndef coupled
      real (kind=dbl_kind), parameter ::
     &   ustar_min = 5.e-2_dbl_kind   ! for zero currents
#endif

      ! Parameters for lateral melting

      real (kind=dbl_kind), parameter ::
     &   floediam = 300.0_dbl_kind    ! effective floe diameter (m)
     &,  alpha    = 0.66_dbl_kind     ! constant from Steele (unitless)
     &,  m1 = 1.6e-6_dbl_kind         ! constant from Maykut & Perovich
                                      ! (m/s/deg^(-m2))
     &,  m2 = 1.36_dbl_kind           ! constant from Maykut & Perovich
                                      ! (unitless)

      rside(:,:) = c0
      Tbot (:,:) = Tf(:,:)
      fbot (:,:) = c0

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt.
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (aice(i,j) > puny .and. frzmlt(i,j) < c0) then ! ice can melt
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo                     ! i
      enddo                     ! j

      do ij = 1, icells  ! cells where ice can melt
         i = indxi(ij)
         j = indxj(ij)

      !-----------------------------------------------------------------
      ! Use boundary layer theory for fbot.
      ! See Maykut and McPhee (1995): JGR, 100, 24,691-24,703.
      !-----------------------------------------------------------------

         deltaT = max((sst(i,j)-Tbot(i,j)),c0)

         ! strocnx has units N/m^2 so strocnx/rho has units m^2/s^2
         ustar = sqrt (sqrt(strocnxT(i,j)**2+strocnyT(i,j)**2)/rhow)
         ustar = max (ustar,ustar_min)

         fbot(i,j) = cpchr * deltaT * ustar ! < 0
         fbot(i,j) = max (fbot(i,j), frzmlt(i,j)) ! frzmlt < fbot < 0

!!! uncomment to use all frzmlt for standalone runs
!!!         fbot(i,j) = min (c0, frzmlt(i,j))

      !-----------------------------------------------------------------
      ! Compute rside.  See these references:
      !    Maykut and Perovich (1987): JGR, 92, 7032-7044
      !    Steele (1992): JGR, 97, 17,729-17,738
      !-----------------------------------------------------------------

         wlat = m1 * deltaT**m2 ! Maykut & Perovich
         rside(i,j) = wlat*dt*pi/(alpha*floediam) ! Steele
         rside(i,j) = max(c0,min(rside(i,j),c1))

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Compute heat flux associated with this value of rside.
      !-----------------------------------------------------------------

      fside(:,:) = c0

      do n = 1, ncat

         etot(:,:) = c0

         ! melting energy/unit area in each column, etot < 0

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               etot(i,j) = etot(i,j) + esnon(i,j,slyr1(n)+k-1)
            enddo               ! ij
         enddo

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               etot(i,j) = etot(i,j) + eicen(i,j,ilyr1(n)+k-1)
            enddo               ! ij
         enddo                  ! nilyr

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            ! lateral heat flux
            fside(i,j) = fside(i,j) + rside(i,j)*etot(i,j)/dt ! fside < 0
         enddo                  ! ij

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Limit bottom and lateral heat fluxes if necessary.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu


      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         xtmp = (fbot(i,j) + fside(i,j)) / (frzmlt(i,j)-puny)
         xtmp = c1/ max(xtmp, c1)
         fbot(i,j)  = fbot(i,j)  * xtmp
         rside(i,j) = rside(i,j) * xtmp
      enddo                     ! ij

      end subroutine frzmlt_bottom_lateral

!=======================================================================
!BOP
!
! !ROUTINE: init_vertical_profile - initial thickness, enthalpy, temperature
!
! !DESCRIPTION:
!
! Given the state variables (vicen, vsnon, eicen, esnon, Tsfcn),
! compute variables needed for the vertical thermodynamics
! (hin, hsn, qin, qsn, Tin, Tsn, Tsf).
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine init_vertical_profile(nx_block, ny_block,
     &                                 my_task,  istep1,
     &                                 icells,
     &                                 indxi,    indxj,
     &                                 aicen,    vicen,
     &                                 vsnon,    Tsfcn,
     &                                 eicen,    esnon,
     &                                 hin,      hilyr,
     &                                 hsn,      hslyr,
     &                                 qin,      Tin,
     &                                 qsn,      Tsn,
     &                                 Tsf,
     &                                 einit,    l_stop,
     &                                 istop,    jstop)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  my_task             ! task number (diagnostic only)
     &,  istep1              ! time step index (diagnostic only)
     &,  icells              ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block),
     &   intent(in) ::
     &   indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aicen     ! concentration of ice
     &,  vicen     ! volume per unit area of ice          (m)
     &,  vsnon     ! volume per unit area of snow         (m)
     &,  Tsfcn     ! temperature of ice/snow top surface  (C)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr),
     &   intent(in) ::
     &   eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nslyr),
     &   intent(in) ::
     &   esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(out)::
     &   hin             ! ice thickness (m)
     &,  hsn             ! snow thickness (m)
     &,  hilyr           ! ice layer thickness
     &,  hslyr           ! snow layer thickness
     &,  Tsf             ! ice/snow surface temperature, Tsfcn
     &,  einit           ! initial energy of melting (J m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(out) ::
     &   qin             ! ice layer enthalpy (J m-3)
     &,  Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr),
     &   intent(out) ::
     &   qsn             ! snow enthalpy
     &,  Tsn             ! snow temperature

      logical (kind=log_kind), intent(inout) ::
     &   l_stop          ! if true, print diagnostics and abort model

      integer (kind=int_kind), intent(out) ::
     &   istop, jstop    ! i and j indices of cell where model fails
!
!EOP
!
      real (kind=dbl_kind), parameter ::
     &   Tmin = -100._dbl_kind ! min allowed internal temperature (deg C)

      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! ice layer index

      real (kind=dbl_kind) ::
     &   aa1, bb1, cc1   ! terms in quadratic formula
     &,  Tmax            ! maximum allowed snow/ice temperature (deg C)

      logical (kind=log_kind) ::   ! for vector-friendly error checks
     &   tsno_high       ! flag for Tsn > Tmax
     &,  tice_high       ! flag for Tin > Tmlt
     &,  tsno_low        ! flag for Tsn < Tmin
     &,  tice_low        ! flag for Tin < Tmin

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      hin     (:,:) = c0
      hsn     (:,:) = c0
      hilyr   (:,:) = c0
      hslyr   (:,:) = c0
      Tsf     (:,:) = c0
      einit   (:,:) = c0

      do k = 1, nslyr
         qsn(:,:,k) = c0
         Tsn(:,:,k) = c0
      enddo
      
      do k = 1, nilyr
         qin(:,:,k) = c0
         Tin(:,:,k) = c0
      enddo

      tsno_high = .false.
      tice_high = .false.
      tsno_low  = .false.
      tice_low  = .false.

      !-----------------------------------------------------------------
      ! Load arrays for vertical thermo calculation.
      !-----------------------------------------------------------------
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !-----------------------------------------------------------------
      ! Surface temperature, ice and snow thickness
      ! Initialize internal energy
      !-----------------------------------------------------------------

         Tsf(i,j)    = Tsfcn(i,j)
         hin(i,j)    = vicen(i,j) / aicen(i,j)
         hsn(i,j)    = vsnon(i,j) / aicen(i,j)
         hilyr(i,j)    = hin(i,j) / real(nilyr)
         hslyr(i,j)    = hsn(i,j) / real(nslyr)

      enddo                     ! ij

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Snow enthalpy and maximum allowed snow temperature
      !
      ! Tmax based on the idea that dT ~ dq / (rhos*cp_ice)
      !                             dq ~ q dv / v
      !                             dv ~ puny = eps11
      ! where 'd' denotes an error due to roundoff.
      !-----------------------------------------------------------------

            if (hslyr(i,j) > hsnomin) then
               ! qsn, esnon < 0              
               qsn  (i,j,k) = esnon(i,j,k)*real(nslyr)/vsnon(i,j) 
               Tmax = -qsn(i,j,k)*puny /
     &                 (rhos*cp_ice*vsnon(i,j))
            else
               qsn  (i,j,k) = -rhos * Lfresh
               Tmax = puny
            endif

      !-----------------------------------------------------------------
      ! Compute snow temperatures from enthalpies.
      ! Note: qsn <= -rhos*Lfresh, so Tsn <= 0.
      !-----------------------------------------------------------------
            Tsn(i,j,k) = (Lfresh + qsn(i,j,k)/rhos)/cp_ice
  
      !-----------------------------------------------------------------
      ! Check for Tsn > Tmax (allowing for roundoff error) and Tsn < Tmin.
      !-----------------------------------------------------------------
            if (Tsn(i,j,k) > Tmax) then
               tsno_high = .true.
            elseif (Tsn(i,j,k) < Tmin) then
               tsno_low  = .true.
            endif

         enddo                  ! ij
      enddo                     ! nslyr

      !-----------------------------------------------------------------
      ! If Tsn is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------
      if (tsno_high) then
         do k = 1, nslyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               if (hslyr(i,j) > hsnomin) then
                  Tmax = -qsn(i,j,k)*puny /
     &                     (rhos*cp_ice*vsnon(i,j))
               else
                  Tmax = puny
               endif

               if (Tsn(i,j,k) > Tmax) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'Starting thermo, Tsn > Tmax'
                  write(nu_diag,*) 'Tsn=',Tsn(i,j,k)
                  write(nu_diag,*) 'Tmax=',Tmax
                  write(nu_diag,*) 'istep1, my_task, i, j:',
     &                              istep1, my_task, i, j
                  write(nu_diag,*) 'qsn',qsn(i,j,k)
                  l_stop = .true.
                  istop = i
                  jstop = j
                  return
               endif

            enddo               ! ij
         enddo                  ! nslyr
      endif                     ! tsno_high

      if (tsno_low) then
         do k = 1, nslyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               if (Tsn(i,j,k) < Tmin) then ! allowing for roundoff error
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'Starting thermo, Tsn < Tmin'
                  write(nu_diag,*) 'Tsn=', Tsn(i,j,k)
                  write(nu_diag,*) 'Tmin=', Tmin
                  write(nu_diag,*) 'istep1, my_task, i, j:',
     &                              istep1, my_task, i, j
                  write(nu_diag,*) 'qsn', qsn(i,j,k)
                  l_stop = .true.
                  istop = i
                  jstop = j
                  return
               endif

            enddo               ! ij
         enddo                  ! nslyr
      endif                     ! tsno_low

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (Tsn(i,j,k) > c0) then   ! correct roundoff error
               Tsn(i,j,k) = c0
               qsn(i,j,k) = -rhos*Lfresh
            endif

      !-----------------------------------------------------------------
      ! initial energy per unit area of ice/snow, relative to 0 C
      !-----------------------------------------------------------------
            einit(i,j) = einit(i,j) + hslyr(i,j)*qsn(i,j,k)

         enddo                  ! ij
      enddo                     ! nslyr

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Compute ice enthalpy
      !-----------------------------------------------------------------
            ! qin, eicen < 0
            qin(i,j,k) = eicen(i,j,k)*real(nilyr)/vicen(i,j)  

      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------

            if (l_brine) then
               aa1 = cp_ice
               bb1 = (cp_ocn-cp_ice)*Tmlt(k) - qin(i,j,k)/rhoi - Lfresh 
               cc1 = Lfresh * Tmlt(k)
               Tin(i,j,k) =  (-bb1 - sqrt(bb1*bb1 - 4.*aa1*cc1)) / 
     &                       (2.*aa1)
               Tmax = Tmlt(k)

            else                ! fresh ice
               Tin(i,j,k) = (Lfresh + qin(i,j,k)/rhoi) / cp_ice
               Tmax = -qin(i,j,k)*puny/(rhos*cp_ice*vicen(i,j))
                         ! as above for snow
            endif

      !-----------------------------------------------------------------
      ! Check for Tin > Tmax and Tin < Tmin
      !-----------------------------------------------------------------
            if (Tin(i,j,k) > Tmax) then
               tice_high = .true.
            elseif (Tin(i,j,k) < Tmin) then
               tice_low  = .true.
            endif

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! If Tin is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------

         if (tice_high) then
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               if (l_brine) then
                  Tmax = Tmlt(k)
               else             ! fresh ice
                  Tmax = -qin(i,j,k)*puny/(rhos*cp_ice*vicen(i,j))
               endif

               if (Tin(i,j,k) > Tmax) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'Starting thermo, T > Tmax, layer', k
                  write(nu_diag,*) 'Tin=',Tin(i,j,k),', Tmax=',Tmax
                  write(nu_diag,*) 'istep1, my_task, i, j:',
     &                              istep1, my_task, i, j
                  write(nu_diag,*) 'qin',qin(i,j,k)
                  l_stop = .true.
                  istop = i
                  jstop = j
                  return
               endif
            enddo               ! ij
         endif                  ! tice_high

         if (tice_low) then
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               if (Tin(i,j,k) < Tmin) then
                  write(nu_diag,*) ' '
                  write(nu_diag,*) 'Starting thermo T < Tmin, layer', k
                  write(nu_diag,*) 'Tin =', Tin(i,j,k)
                  write(nu_diag,*) 'Tmin =', Tmin
                  write(nu_diag,*) 'istep1, my_task, i, j:',
     &                              istep1, my_task, i, j
                  l_stop = .true.
                  istop = i
                  jstop = j
                  return
               endif
            enddo               ! ij
         endif                  ! tice_low

      !-----------------------------------------------------------------
      ! initial energy per unit area of ice/snow, relative to 0 C
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (Tin(i,j,k) > c0) then ! correct roundoff error
               Tin(i,j,k) = c0
               qin(i,j,k) = -rhoi*Lfresh
            endif
            
            einit(i,j) = einit(i,j) + hilyr(i,j)*qin(i,j,k) 
            
         enddo                  ! ij


      enddo                     ! nilyr

      end subroutine init_vertical_profile

!=======================================================================
!BOP
!
! !ROUTINE: temperature_changes  - new vertical temperature profile
!
! !DESCRIPTION:
!
! Compute new surface temperature and internal ice and snow
! temperatures.  Include effects of salinity on sea ice heat
! capacity in a way that conserves energy (Bitz and Lipscomb, 1999).
!
! New temperatures are computed iteratively by solving a tridiagonal
! system of equations; heat capacity is updated with each iteration.
! Finite differencing is backward implicit.
!
! See Bitz, C.M., and W.H. Lipscomb, 1999:
! An energy-conserving thermodynamic model of sea ice,
! J. Geophys. Res., 104, 15,669-15,677.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine temperature_changes (nx_block, ny_block,
     &                                my_task,  istep1,
     &                                dt,       icells,   
     &                                indxi,    indxj,
     &                                rhoa,     flw,    
     &                                potT,     Qa,
     &                                shcoef,   lhcoef,
     &                                fswsfc,   fswint,
     &                                fswthrun, Iswabs,
     &                                hilyr,    hslyr,
     &                                qin,      Tin,
     &                                qsn,      Tsn,
     &                                Tsf,      Tbot,
     &                                fsensn,   flatn,
     &                                fswabsn,  flwoutn,
     &                                fsurf,  
     &                                fcondtop, fcondbot, 
     &                                einit,    l_stop,
     &                                istop,    jstop)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  my_task         ! task number (diagnostic only)
     &,  istep1          ! time step index (diagnostic only)
     &,  icells          ! number of cells with aicen > puny

      real (kind=dbl_kind), intent(in) ::
     &   dt              ! time step

      integer (kind=int_kind), dimension(nx_block*ny_block),
     &   intent(in) ::
     &   indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   rhoa            ! air density (kg/m^3)
     &,  flw             ! incoming longwave radiation (W/m^2)
     &,  potT            ! air potential temperature  (K)
     &,  Qa              ! specific humidity (kg/kg)
     &,  shcoef          ! transfer coefficient for sensible heat
     &,  lhcoef          ! transfer coefficient for latent heat
     &,  hilyr           ! ice layer thickness (m)
     &,  hslyr           ! snow layer thickness (m)
     &,  Tbot            ! ice bottom surface temperature (deg C)
     &,  einit           ! initial energy of melting (J m-2)
     &,  fswsfc          ! SW absorbed at ice/snow surface (W m-2)
     &,  fswint          ! SW absorbed in ice interior below surface (W m-2)
     &,  fswthrun        ! SW through ice to ocean         (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(in) ::
     &   Iswabs          ! SW radiation absorbed in ice layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out)::
     &   fsensn          ! surface downward sensible heat (W m-2)
     &,  flatn           ! surface downward latent heat (W m-2)
     &,  fswabsn         ! shortwave absorbed by ice (W m-2)
     &,  flwoutn         ! upward LW at surface (W m-2)
     &,  fsurf           ! net flux to top surface, not including fcondtop
     &,  fcondtop        ! downward cond flux at top surface (W m-2)
     &,  fcondbot        ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout)::
     &   Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(inout) ::
     &   qin             ! ice layer enthalpy (J m-3)
     &,  Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr),
     &   intent(inout) ::
     &   qsn             ! snow layer enthalpy (J m-3)
     &,  Tsn             ! internal snow layer temperatures

      logical (kind=log_kind), intent(inout) ::
     &   l_stop          ! if true, print diagnostics and abort model

      integer (kind=int_kind), intent(out) ::
     &   istop, jstop    ! i and j indices of cell where model fails
!
!EOP
!
      integer (kind=int_kind), parameter ::
     &   nitermax = 50    ! max number of iterations in temperature solver
     &,  nmat = nslyr + nilyr + 1  ! matrix dimension

      real (kind=dbl_kind), parameter ::
     &   Tsf_errmax = 5.e-4_dbl_kind ! max allowed error in Tsf
                                     ! recommend Tsf_errmax < 0.01 K

      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! ice layer index
     &,  niter           ! iteration counter in temperature solver

      integer (kind=int_kind) ::
     &   isolve          ! number of cells with temps not converged

      integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &   indxii, indxjj  ! compressed indices for cells not converged

      logical (kind=log_kind), dimension (nx_block,ny_block) ::
     &   l_snow          ! true if snow temperatures are computed
     &,  l_cold          ! true if surface temperature is computed

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   Tsf_start       ! Tsf at start of iteration
     &,  dTsf            ! Tsf - Tsf_start
     &,  dTsf_prev       ! dTsf from previous iteration
     &,  dfsurf_dT       ! derivative of fsurf wrt Tsf
     &,  dfsens_dT       ! deriv of fsens wrt Tsf (W m-2 deg-1)
     &,  dflat_dT        ! deriv of flat wrt Tsf (W m-2 deg-1)
     &,  dflwout_dT      ! deriv of flwout wrt Tsf (W m-2 deg-1)
     &,  dt_rhoi_hlyr    ! dt/(rhoi*hilyr)
     &,  avg_Tsi         ! = 1. if new snow/ice temps avg'd w/starting temps
     &,  enew            ! new energy of melting after temp change (J m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr) ::
     &   Tin_init        ! Tin at beginning of time step
     &,  Tin_start       ! Tin at start of iteration
     &,  etai            ! dt / (rho * cp * h) for ice layers

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) ::
     &   Tsn_init        ! Tsn at beginning of time step
     &,  Tsn_start       ! Tsn at start of iteration
     &,  etas            ! dt / (rho * cp * h) for snow layers

      real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr+nslyr+1)::
     &   kh              ! effective conductivity at interfaces (W m-2 deg-1)
     &,  diag            ! diagonal matrix elements
     &,  sbdiag          ! sub-diagonal matrix elements
     &,  spdiag          ! super-diagonal matrix elements
     &,  rhs             ! rhs of tri-diagonal matrix equation
     &,  Tmat            ! matrix output temperatures

      real (kind=dbl_kind) ::
     &   ci              ! specific heat of sea ice (J kg-1 deg-1)
     &,  avg_Tsf         ! = 1. if Tsf averaged w/Tsf_start, else = 0.
     &,  ferr            ! energy conservation error (W m-2)

      logical (kind=log_kind), dimension (nx_block,ny_block) ::
     &   converged      ! = true when local solution has converged

      logical (kind=log_kind) ::
     &   all_converged  ! = true when all cells have converged

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      fsensn  (:,:) = c0
      flatn   (:,:) = c0
      fswabsn (:,:) = c0
      flwoutn (:,:) = c0
      fsurf   (:,:) = c0
      fcondtop(:,:) = c0
      fcondbot(:,:) = c0

      all_converged   = .false.

      converged (:,:) = .false.
      l_snow    (:,:) = .false.
      l_cold    (:,:) = .true.

      dTsf      (:,:) = c0
      dTsf_prev (:,:) = c0
      Tsf_start (:,:) = c0
      enew      (:,:) = c0
      avg_Tsi   (:,:) = c0 
      dt_rhoi_hlyr(:,:) = c0
      dfsurf_dT (:,:) = c0
      dfsens_dt (:,:) = c0
      dflat_dt  (:,:) = c0
      dflwout_dt(:,:) = c0

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         dt_rhoi_hlyr(i,j) = dt / (rhoi*hilyr(i,j))  ! hilyr > 0
         if (hslyr(i,j) > hsnomin) l_snow(i,j) = .true.
      enddo                     ! ij

      do k = 1, nslyr

         etas     (:,:,k) = c0
         Tsn_init (:,:,k) = Tsn(:,:,k) ! beginning of time step
         Tsn_start(:,:,k) = Tsn(:,:,k) ! beginning of iteration

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            if (l_snow(i,j)) etas(i,j,k) = dt/(rhos*cp_ice*hslyr(i,j))
         enddo                  ! ij
      enddo                     ! k

      do k = 1, nilyr
         etai     (:,:,k) = c0
         Tin_init (:,:,k) = Tin(:,:,k)   ! beginning of time step
         Tin_start(:,:,k) = Tin(:,:,k)   ! beginning of iteration
      enddo

      !-----------------------------------------------------------------
      ! Initialize matrix elements.
      ! Note: When we do not need to solve for the surface or snow
      !       temperature, we solve dummy equations with solution T = 0.
      !-----------------------------------------------------------------

      do k = 1, nilyr+nslyr+1
         sbdiag(:,:,k) = c0
         diag  (:,:,k) = c1
         spdiag(:,:,k) = c0
         rhs   (:,:,k) = c0
         Tmat  (:,:,k) = c0
      enddo

      !-----------------------------------------------------------------
      ! Compute thermal conductivity at interfaces (held fixed during
      !  subsequent iterations).
      ! Ice and snow interfaces are combined into one array (kh) to
      !  simplify the logic.
      !-----------------------------------------------------------------

      call conductivity (nx_block, ny_block,
     &                   l_snow,   icells, 
     &                   indxi,    indxj,
     &                   hilyr,    hslyr,
     &                   Tin,      kh )


      !-----------------------------------------------------------------
      ! Solve for new temperatures.
      ! Iterate until temperatures converge with minimal energy error.
      !-----------------------------------------------------------------

      do niter = 1, nitermax

         if (all_converged) then  ! thermo calculation is done
            exit
         else                     ! identify cells not yet converged
            isolve = 0
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (.not.converged(i,j)) then
                  isolve = isolve + 1
                  indxii(isolve) = i
                  indxjj(isolve) = j
               endif
            enddo               ! ij
         endif

      !-----------------------------------------------------------------
      ! Update radiative and turbulent fluxes and their derivatives
      ! with respect to Tsf.
      !-----------------------------------------------------------------

         call surface_fluxes (nx_block,    ny_block,
     &                        isolve,     
     &                        indxii,      indxjj,
     &                        Tsf,         fswsfc,
     &                        rhoa,        flw,
     &                        potT,        Qa,
     &                        shcoef,      lhcoef,
     &                        flwoutn,     fsensn,
     &                        flatn,       fsurf,
     &                        dflwout_dT,  dfsens_dT,
     &                        dflat_dT,    dfsurf_dT)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            i = indxii(ij)   ! NOTE: not indxi and indxj
            j = indxjj(ij)

      !-----------------------------------------------------------------
      ! Compute conductive flux at top surface, fcondtop.
      ! If fsurf < fcondtop and Tsf = 0, then reset Tsf to slightly less
      !  than zero (but not less than -puny).
      !-----------------------------------------------------------------

            if (l_snow(i,j)) then
               fcondtop(i,j) = kh(i,j,1) * (Tsf(i,j) - Tsn(i,j,1))
            else
               fcondtop(i,j) = kh(i,j,1+nslyr) * (Tsf(i,j) - Tin(i,j,1))
            endif

            if (fsurf(i,j) < fcondtop(i,j))
     &           Tsf(i,j) = min (Tsf(i,j), -puny)

      !-----------------------------------------------------------------
      ! Save surface temperature at start of iteration
      !-----------------------------------------------------------------

            Tsf_start(i,j) = Tsf(i,j)

            if (Tsf(i,j) <= -puny) then
               l_cold(i,j) = .true.
            else
               l_cold(i,j) = .false.
            endif
         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Update specific heat of ice layers.
      ! To ensure energy conservation, the specific heat is a function of
      ! both the starting temperature and the (latest guess for) the
      ! final temperature.
      !-----------------------------------------------------------------

         do k = 1, nilyr
            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)

               if (l_brine) then
                  ci = cp_ice - Lfresh*Tmlt(k) / 
     &                          (Tin(i,j,k)*Tin_init(i,j,k))
               else
                  ci = cp_ice
               endif
               etai(i,j,k) = dt_rhoi_hlyr(i,j) / ci

            enddo
         enddo

         call get_matrix_elements (nx_block, ny_block,
     &                             isolve,     
     &                             indxii,   indxjj,   
     &                             l_snow,   l_cold,
     &                             Tsf,      Tbot,
     &                             fsurf,    dfsurf_dT,
     &                             Tin_init, Tsn_init, 
     &                             kh,       Iswabs,
     &                             etai,     etas,
     &                             sbdiag,   diag,
     &                             spdiag,   rhs)

      !-----------------------------------------------------------------
      ! Solve tridiagonal matrix to obtain the new temperatures.
      !-----------------------------------------------------------------

         call tridiag_solver (nx_block, ny_block,
     &                        isolve,     
     &                        indxii,   indxjj,
     &                        nmat,     sbdiag,
     &                        diag,     spdiag, 
     &                        rhs,      Tmat)

      !-----------------------------------------------------------------
      ! Determine whether the computation has converged to an acceptable
      ! solution.  Five conditions must be satisfied:
      !
      !    (1) Tsf <= 0 C.
      !    (2) Tsf is not oscillating; i.e., if both dTsf(niter) and
      !        dTsf(niter-1) have magnitudes greater than puny, then
      !        dTsf(niter)/dTsf(niter-1) cannot be a negative number
      !        with magnitude greater than 0.5.  
      !    (3) abs(dTsf) < Tsf_errmax
      !    (4) If Tsf = 0 C, then the downward turbulent/radiative 
      !        flux, fsurf, must be greater than or equal to the downward
      !        conductive flux, fcondtop.
      !    (5) The net energy added to the ice per unit time must equal 
      !        the net change in internal ice energy per unit time,
      !        within the prescribed error ferrmax.
      !
      ! For briny ice (the standard case), Tsn and Tin are limited
      !  to prevent them from exceeding their melting temperatures.
      !  (Note that the specific heat formula for briny ice assumes
      !  that T < Tmlt.)  
      ! For fresh ice there is no limiting, since there are cases
      !  when the only convergent solution has Tsn > 0 and/or Tin > 0.
      !  Above-zero temperatures are then reset to zero (with melting 
      !  to conserve energy) in the thickness_changes subroutine.
      !-----------------------------------------------------------------

         ! initialize global convergence flag
         all_converged = .true.

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

      !-----------------------------------------------------------------
      ! Reload Tsf from matrix solution
      !-----------------------------------------------------------------

            if (l_cold(i,j)) then
               if (l_snow(i,j)) then
                  Tsf(i,j) = Tmat(i,j,1)
               else
                  Tsf(i,j) = Tmat(i,j,1+nslyr)
               endif
            else                ! melting surface
               Tsf(i,j) = c0
            endif

      !-----------------------------------------------------------------
      ! Initialize convergence flag (true until proven false), dTsf,
      !  and temperature-averaging coefficients.
      ! Average only if test 1 or 2 fails.
      ! Initialize energy.
      !-----------------------------------------------------------------

            converged(i,j) = .true.
            dTsf(i,j) = Tsf(i,j) - Tsf_start(i,j)
            avg_Tsf      = c0
            avg_Tsi(i,j) = c0
            enew(i,j) = c0

      !-----------------------------------------------------------------
      ! Condition 1: check for Tsf > 0
      ! If Tsf > 0, set Tsf = 0, then average Tsn and Tin to force
      ! internal temps below their melting temps.
      !-----------------------------------------------------------------

            if (Tsf(i,j) > puny) then
               Tsf(i,j) = c0
               dTsf(i,j) = -Tsf_start(i,j)
               if (l_brine) avg_Tsi(i,j) = c1   ! avg with starting temp
               converged(i,j) = .false.
               all_converged = .false.

      !-----------------------------------------------------------------
      ! Condition 2: check for oscillating Tsf
      ! If oscillating, average all temps to increase rate of convergence.
      !-----------------------------------------------------------------

            elseif (niter > 1                  ! condition (2)
     &        .and. Tsf_start(i,j) <= -puny
     &        .and. abs(dTsf(i,j)) > puny
     &        .and. abs(dTsf_prev(i,j)) > puny
     &        .and. -dTsf(i,j)/(dTsf_prev(i,j)+puny*puny) > p5) then

               if (l_brine) then ! average with starting temp
                  avg_Tsf  = c1    
                  avg_Tsi(i,j) = c1
               endif
               dTsf(i,j) = p5 * dTsf(i,j)
               converged(i,j) = .false.
               all_converged = .false.
            endif

      !-----------------------------------------------------------------
      ! If condition 2 failed, average new surface temperature with
      !  starting value.
      !-----------------------------------------------------------------
            Tsf(i,j)  = Tsf(i,j)
     &                + avg_Tsf * p5 * (Tsf_start(i,j) - Tsf(i,j))

         enddo  ! ij

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)

      !-----------------------------------------------------------------
      ! Reload Tsn from matrix solution
      !-----------------------------------------------------------------

               if (l_snow(i,j)) then
                  Tsn(i,j,k) = Tmat(i,j,k+1)
               else
                  Tsn(i,j,k) = c0
               endif
               if (l_brine) Tsn(i,j,k) = min(Tsn(i,j,k), c0)

      !-----------------------------------------------------------------
      ! If condition 1 or 2 failed, average new snow layer
      !  temperatures with their starting values.
      !-----------------------------------------------------------------
               Tsn(i,j,k) = Tsn(i,j,k)
     &                   + avg_Tsi(i,j)*p5*(Tsn_start(i,j,k)-Tsn(i,j,k))

      !-----------------------------------------------------------------
      ! Compute qsn and increment new energy.
      !-----------------------------------------------------------------
               qsn(i,j,k) = -rhos * (Lfresh - cp_ice*Tsn(i,j,k))
               enew(i,j) = enew(i,j) + hslyr(i,j) * qsn(i,j,k)

               Tsn_start(i,j,k) = Tsn(i,j,k) ! for next iteration

            enddo               ! ij
         enddo                  ! nslyr

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)

      !-----------------------------------------------------------------
      ! Reload Tin from matrix solution
      !-----------------------------------------------------------------
               Tin(i,j,k) = Tmat(i,j,k+1+nslyr)
               if (l_brine) Tin(i,j,k) = min(Tin(i,j,k), Tmlt(k))

      !-----------------------------------------------------------------
      ! If condition 1 or 2 failed, average new ice layer
      !  temperatures with their starting values.
      !-----------------------------------------------------------------
               Tin(i,j,k) = Tin(i,j,k)
     &                   + avg_Tsi(i,j)*p5*(Tin_start(i,j,k)-Tin(i,j,k))

      !-----------------------------------------------------------------
      ! Compute qin and increment new energy.
      !-----------------------------------------------------------------
               qin(i,j,k) = -rhoi * (cp_ice*(Tmlt(k)-Tin(i,j,k))
     &                             + Lfresh*(c1-Tmlt(k)/Tin(i,j,k))
     &                             - cp_ocn*Tmlt(k))
               enew(i,j) = enew(i,j) + hilyr(i,j) * qin(i,j,k)

               Tin_start(i,j,k) = Tin(i,j,k) ! for next iteration

            enddo               ! ij
         enddo                  ! nilyr

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

      !-----------------------------------------------------------------
      ! Condition 3: check for large change in Tsf
      !-----------------------------------------------------------------

            if (abs(dTsf(i,j)) > Tsf_errmax) then
               converged(i,j) = .false.
               all_converged = .false.
            endif

      !-----------------------------------------------------------------
      ! Condition 4: check for fsurf < fcondtop with Tsf > 0
      !-----------------------------------------------------------------

            fsurf(i,j) = fsurf(i,j) + dTsf(i,j)*dfsurf_dT(i,j)
            if (l_snow(i,j)) then
               fcondtop(i,j) = kh(i,j,1) * (Tsf(i,j)-Tsn(i,j,1))
            else
               fcondtop(i,j) = kh(i,j,1+nslyr) * (Tsf(i,j)-Tin(i,j,1))
            endif

            if (Tsf(i,j) > -puny .and. fsurf(i,j) < fcondtop(i,j)) then
               converged(i,j) = .false.
               all_converged = .false.
            endif

      !-----------------------------------------------------------------
      ! Condition 5: check for energy conservation error
      ! Change in internal ice energy should equal net energy input.
      !-----------------------------------------------------------------

            fcondbot(i,j) = kh(i,j,1+nslyr+nilyr) *
     &                     (Tin(i,j,nilyr)   - Tbot(i,j))

            ferr = abs( (enew(i,j)-einit(i,j))/dt
     &           - (fcondtop(i,j) - fcondbot(i,j) + fswint(i,j)) )

            ! factor of 0.9 allows for roundoff errors later
            if (ferr > 0.9*ferrmax) then               ! condition (5)
               converged(i,j) = .false.
               all_converged = .false.
            endif

            dTsf_prev(i,j) = dTsf(i,j)

         enddo                  ! ij

      enddo                     ! temperature iteration

      if (.not.all_converged) then

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Check for convergence failures.
      !-----------------------------------------------------------------
            if (.not.converged(i,j)) then
               write(nu_diag,*) 'Thermo iteration does not converge,',
     &                          'istep1, my_task, i, j:',
     &                           istep1, my_task, i, j
               write(nu_diag,*) 'Ice thickness:',  hilyr(i,j)*nilyr
               write(nu_diag,*) 'Snow thickness:', hslyr(i,j)*nslyr
               write(nu_diag,*) 'dTsf, Tsf_errmax:',dTsf(i,j),Tsf_errmax
               write(nu_diag,*) 'Tsf:', Tsf(i,j)
               write(nu_diag,*) 'fsurf:', fsurf(i,j)
               write(nu_diag,*) 'fcondtop, fcondbot, fswint',
     &                         fcondtop(i,j), fcondbot(i,j), fswint(i,j)
               write(nu_diag,*) 'Flux conservation error =', ferr
               write(nu_diag,*) 'Initial snow temperatures:'
               write(nu_diag,*) (Tsn_init(i,j,k),k=1,nslyr)
               write(nu_diag,*) 'Initial ice temperatures:'
               write(nu_diag,*) (Tin_init(i,j,k),k=1,nilyr)
               write(nu_diag,*) 'Final snow temperatures:'
               write(nu_diag,*) (Tsn(i,j,k),k=1,nslyr)
               write(nu_diag,*) 'Final ice temperatures:'
               write(nu_diag,*) (Tin(i,j,k),k=1,nilyr)
               l_stop = .true.
               istop = i
               jstop = j
               return
            endif
         enddo                  ! ij
      endif                     ! all_converged

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! update fluxes that depend on Tsf
         flwoutn(i,j) = flwoutn(i,j) + dTsf(i,j) * dflwout_dT(i,j)
         fsensn(i,j)  = fsensn(i,j)  + dTsf(i,j) * dfsens_dT(i,j)
         flatn(i,j)   = flatn(i,j)   + dTsf(i,j) * dflat_dT(i,j)

         ! absorbed shortwave flux for coupler
         fswabsn(i,j) = fswsfc(i,j) + fswint(i,j) + fswthrun(i,j)

      enddo                     ! ij

      end subroutine temperature_changes

!=======================================================================
!BOP
!
! !ROUTINE: conductivity - compute ice thermal conductivity
!
! !DESCRIPTION:
!
! Compute thermal conductivity at interfaces (held fixed during
!  the subsequent iteration).
!
! NOTE: Ice conductivity must be >= kimin
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine conductivity (nx_block, ny_block,
     &                         l_snow,   icells, 
     &                         indxi,    indxj,
     &                         hilyr,    hslyr, 
     &                         Tin,      kh)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  icells          ! number of cells with aicen > puny

      logical (kind=log_kind), dimension(nx_block,ny_block),
     &   intent(in) ::
     &   l_snow          ! true if snow temperatures are computed

      integer (kind=int_kind), dimension(nx_block*ny_block),
     &   intent(in) ::
     &   indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   hilyr           ! ice layer thickness (same for all ice layers)
     &,  hslyr           ! snow layer thickness (same for all snow layers)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(in) ::
     &   Tin             ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr+nslyr+1),
     &   intent(out) ::
     &   kh              ! effective conductivity at interfaces (W m-2 deg-1)
!
!EOP
!
      real (kind=dbl_kind), parameter ::
     &   betak   = 0.13_dbl_kind ! constant in formula for k (W m-1 ppt-1)
     &,  kimin   = 0.10_dbl_kind ! min conductivity of saline ice (W m-1 deg-1)

      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! vertical index

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr) ::
     &   kilyr           ! thermal cond at ice layer midpoints (W m-1 deg-1)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) ::
     &   kslyr           ! thermal cond at snow layer midpoints (W m-1 deg-1)

      ! initialize
      do k = 1, nilyr+nslyr+1
         kh(:,:,k) = c0
      enddo

      ! interior snow layers (simple for now, but may be fancier later)
      do k = 1, nslyr
         kslyr(:,:,k) = c0

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            kslyr(i,j,k) = ksno
         enddo
      enddo                     ! nslyr

      ! interior ice layers
      do k = 1, nilyr
         kilyr(:,:,k) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            kilyr(i,j,k) = kice + betak*salin(k)/min(-puny,Tin(i,j,k))
            kilyr(i,j,k) = max (kilyr(i,j,k), kimin)
         enddo
      enddo                     ! nilyr

      ! top snow interface, top and bottom ice interfaces
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! top of snow layer; top surface of top ice layer
         if (l_snow(i,j)) then
            kh(i,j,1)       = c2 * kslyr(i,j,1) / hslyr(i,j)
            kh(i,j,1+nslyr) = c2 * kslyr(i,j,nslyr) * kilyr(i,j,1) /
     &                       ( kslyr(i,j,nslyr)*hilyr(i,j) + 
     &                         kilyr(i,j,1    )*hslyr(i,j) )
         else
            kh(i,j,1)       = c0
            kh(i,j,1+nslyr) = c2 * kilyr(i,j,1) / hilyr(i,j)
         endif

         ! bottom surface of bottom ice layer
         kh(i,j,1+nslyr+nilyr) = c2 * kilyr(i,j,nilyr) / hilyr(i,j)

      enddo                     ! ij

      ! interior snow interfaces
      do k = 2, nslyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (l_snow(i,j)) then
               kh(i,j,k) = c2 * kslyr(i,j,k-1) * kslyr(i,j,k) /
     &                      ((kslyr(i,j,k-1) + kslyr(i,j,k))*hslyr(i,j))
            else
               kh(i,j,k) = c0
            endif
         enddo                  ! ij
      enddo                     ! nilyr

      ! interior ice interfaces
      do k = 2, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            kh(i,j,k+nslyr) = c2 * kilyr(i,j,k-1) * kilyr(i,j,k) /
     &                      ((kilyr(i,j,k-1) + kilyr(i,j,k))*hilyr(i,j))
         enddo                  ! ij
      enddo                     ! nilyr

      end subroutine conductivity

!=======================================================================
!BOP
!
! !ROUTINE: surface_fluxes - surface radiative and turbulent fluxes
!
! !DESCRIPTION:
!
! Compute radiative and turbulent fluxes and their derivatives
! with respect to Tsf.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine surface_fluxes (nx_block,   ny_block,
     &                           isolve,     
     &                           indxii,     indxjj,
     &                           Tsf,        fswsfc,
     &                           rhoa,       flw,
     &                           potT,       Qa,
     &                           shcoef,     lhcoef,
     &                           flwoutn,    fsensn,
     &                           flatn,      fsurf,
     &                           dflwout_dT, dfsens_dT,
     &                           dflat_dT,   dfsurf_dT)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  isolve              ! number of cells with temps not converged

      integer (kind=int_kind), dimension(nx_block*ny_block),
     &   intent(in) ::
     &   indxii, indxjj  ! compressed indices for cells not converged

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   Tsf             ! ice/snow surface temperature, Tsfcn
     &,  fswsfc          ! SW absorbed at ice/snow surface (W m-2)
     &,  rhoa            ! air density (kg/m^3)
     &,  flw             ! incoming longwave radiation (W/m^2)
     &,  potT            ! air potential temperature  (K)
     &,  Qa              ! specific humidity (kg/kg)
     &,  shcoef          ! transfer coefficient for sensible heat
     &,  lhcoef          ! transfer coefficient for latent heat

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   fsensn          ! surface downward sensible heat (W m-2)
     &,  flatn           ! surface downward latent heat (W m-2)
     &,  flwoutn         ! upward LW at surface (W m-2)
     &,  fsurf           ! net flux to top surface, not including fcondtop

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   dfsens_dT       ! deriv of fsens wrt Tsf (W m-2 deg-1)
     &,  dflat_dT        ! deriv of flat wrt Tsf (W m-2 deg-1)
     &,  dflwout_dT      ! deriv of flwout wrt Tsf (W m-2 deg-1)
     &,  dfsurf_dT       ! derivative of fsurf wrt Tsf
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! ice layer index

      real (kind=dbl_kind) ::
     &   TsfK            ! ice/snow surface temperature (K)
     &,  Qsfc            ! saturated surface specific humidity (kg/kg)
     &,  dqsatdt         ! derivative of qsat wrt surface temperature
     &,  qsat            ! the saturation humidity of air (kg/m^3)
     &,  flwdabs         ! downward longwave absorbed heat flx (W/m^2)
     &,  tmpvar          ! 1/TsfK

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, isolve
         i = indxii(ij)         ! NOTE: not indxi and indxj
         j = indxjj(ij)

         ! ice surface temperature in Kelvin
         TsfK = Tsf(i,j) + Tffresh
         tmpvar = c1/TsfK

         ! saturation humidity
         qsat    = qqqice * exp(-TTTice*tmpvar)
         Qsfc    = qsat / rhoa(i,j)
         dqsatdt = TTTice * tmpvar*tmpvar * qsat

         ! longwave radiative flux
         flwdabs =  emissivity * flw(i,j)
         flwoutn(i,j) = -emissivity * stefan_boltzmann * TsfK**4

         ! downward latent and sensible heat fluxes
         fsensn(i,j) = shcoef(i,j) * (potT(i,j) - TsfK)
         flatn(i,j)  = lhcoef(i,j) * (Qa(i,j) - Qsfc)

         ! derivatives wrt surface temp
         dflwout_dT(i,j) = - emissivity*stefan_boltzmann * c4*TsfK**3
         dfsens_dT(i,j)  = - shcoef(i,j)
         dflat_dT(i,j)   = - lhcoef(i,j) * dqsatdt

         fsurf(i,j) = fswsfc(i,j) + flwdabs + flwoutn(i,j)
     &              + fsensn(i,j) + flatn(i,j)
         dfsurf_dT(i,j) = dflwout_dT(i,j)
     &                   + dfsens_dT(i,j) + dflat_dT(i,j)

      enddo                     ! ij

      end subroutine surface_fluxes

!=======================================================================
!BOP
!
! !ROUTINE: get_matrix_elements - compute tridiagonal matrix elements
!
! !DESCRIPTION:
!
! Compute terms in tridiagonal matrix that will be solved to find
!  the new vertical temperature profile
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
!
! revised March 2004 by William H. Lipscomb for multiple snow layers
!
! !INTERFACE:
!
      subroutine get_matrix_elements (nx_block, ny_block,
     &                                isolve,     
     &                                indxii,   indxjj,   
     &                                l_snow,   l_cold,
     &                                Tsf,      Tbot,
     &                                fsurf,    dfsurf_dT,
     &                                Tin_init, Tsn_init, 
     &                                kh,       Iswabs,
     &                                etai,     etas,
     &                                sbdiag,   diag,
     &                                spdiag,   rhs)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  isolve          ! number of cells with temps not converged

      integer (kind=int_kind), dimension(nx_block*ny_block),
     &   intent(in) ::
     &   indxii, indxjj  ! compressed indices for cells not converged

      logical (kind=log_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   l_snow          ! true if snow temperatures are computed
     &,  l_cold          ! true if surface temperature is computed

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   Tsf             ! ice/snow top surface temp (deg C)
     &,  Tbot            ! ice bottom surface temperature (deg C)
     &,  fsurf           ! net flux to top surface, not including fcondtop

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   dfsurf_dT       ! derivative of fsurf wrt Tsf

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(in) ::
     &   etai            ! dt / (rho*cp*h) for ice layers
     &,  Tin_init        ! ice temp at beginning of time step
     &,  Iswabs          ! absorbed SW flux in ice layers

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr),
     &   intent(in) ::
     &   etas            ! dt / (rho*cp*h) for snow layers
     &,  Tsn_init        ! snow temp at beginning of time step
                         ! Note: no absorbed SW in snow layers

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr+nilyr+1),
     &   intent(in) ::
     &   kh              ! effective conductivity at layer interfaces

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr+nilyr+1),
     &   intent(inout) ::
     &   sbdiag          ! sub-diagonal matrix elements
     &,  diag            ! diagonal matrix elements
     &,  spdiag          ! super-diagonal matrix elements
     &,  rhs             ! rhs of tri-diagonal matrix eqn.
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k, ks, ki, kr   ! vertical indices and row counters

      !-----------------------------------------------------------------
      ! Compute matrix elements
      !
      ! Four possible cases to solve:
      !   (1) Cold surface (Tsf < 0), snow present
      !   (2) Melting surface (Tsf = 0), snow present
      !   (3) Cold surface (Tsf < 0), no snow
      !   (4) Melting surface (Tsf = 0), no snow
      !-----------------------------------------------------------------

      do ij = 1, isolve
         i = indxii(ij)
         j = indxjj(ij)

      !-----------------------------------------------------------------
      ! Tsf equation for case of cold surface (with or without snow)
      !-----------------------------------------------------------------
         if (l_cold(i,j)) then
            if (l_snow(i,j)) then
               k = 1
            else                ! no snow
               k = 1 + nslyr
            endif
            kr = k
            
            sbdiag(i,j,kr) = c0
            diag  (i,j,kr) = dfsurf_dT(i,j) - kh(i,j,k)
            spdiag(i,j,kr) = kh(i,j,k)
            rhs   (i,j,kr) = dfsurf_dT(i,j)*Tsf(i,j) - fsurf(i,j)
            
         endif                  ! l_cold
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! top snow layer
      !-----------------------------------------------------------------
      k = 1
      kr = 2

      do ij = 1, isolve
         i = indxii(ij)
         j = indxjj(ij)

         if (l_snow(i,j)) then
            if (l_cold(i,j)) then
               sbdiag(i,j,kr) = -etas(i,j,k) * kh(i,j,k)
               spdiag(i,j,kr) = -etas(i,j,k) * kh(i,j,k+1)
               diag  (i,j,kr) = c1
     &                        + etas(i,j,k) * (kh(i,j,k) + kh(i,j,k+1))
               rhs   (i,j,kr) = Tsn_init(i,j,k)
            else                ! melting surface
               sbdiag(i,j,kr) = c0
               spdiag(i,j,kr) = -etas(i,j,k) * kh(i,j,k+1)
               diag  (i,j,kr) = c1
     &                        + etas(i,j,k) * (kh(i,j,k) + kh(i,j,k+1))
               rhs   (i,j,kr) = Tsn_init(i,j,k)
     &                        + etas(i,j,k)*kh(i,j,k)*Tsf(i,j)
            endif               ! l_cold
         endif                  ! l_snow
      enddo                     ! ij
      
      !-----------------------------------------------------------------
      ! remaining snow layers
      !-----------------------------------------------------------------

      if (nslyr > 1) then
         do k = 2, nslyr
            kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)

               if (l_snow(i,j)) then
                  sbdiag(i,j,kr) = -etas(i,j,k) * kh(i,j,k)
                  spdiag(i,j,kr) = -etas(i,j,k) * kh(i,j,k+1)
                  diag  (i,j,kr) = c1
     &                         + etas(i,j,k) * (kh(i,j,k) + kh(i,j,k+1))
                  rhs   (i,j,kr) = Tsn_init(i,j,k)
               endif
            enddo               ! ij
         enddo                  ! nslyr
      endif                     ! nslyr > 1

      !-----------------------------------------------------------------
      ! top ice layer
      !-----------------------------------------------------------------

      ki = 1
      k  = ki + nslyr
      kr = k + 1

      do ij = 1, isolve
         i = indxii(ij)
         j = indxjj(ij)

         if (.not.l_snow(i,j) .and. .not.l_cold(i,j)) then
            sbdiag(i,j,kr) = c0
            spdiag(i,j,kr) = -etai(i,j,ki) * kh(i,j,k+1)
            diag  (i,j,kr) = c1
     &                     + etai(i,j,ki) * (kh(i,j,k) + kh(i,j,k+1))
            rhs   (i,j,kr) = Tin_init(i,j,ki)
     &                     + etai(i,j,ki)*Iswabs(i,j,ki)
     &                     + etai(i,j,ki)*kh(i,j,k)*Tsf(i,j)
         else
            sbdiag(i,j,kr) = -etai(i,j,ki) * kh(i,j,k)
            spdiag(i,j,kr) = -etai(i,j,ki) * kh(i,j,k+1)
            diag  (i,j,kr) = c1
     &                     + etai(i,j,ki) * (kh(i,j,k) + kh(i,j,k+1))
            rhs   (i,j,kr) = Tin_init(i,j,ki)
     &                     + etai(i,j,ki)*Iswabs(i,j,ki)
         endif

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! interior ice layers
      !-----------------------------------------------------------------
      do ki = 2, nilyr-1
           
         k  = ki + nslyr
         kr = k + 1
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

            sbdiag(i,j,kr) = -etai(i,j,ki) * kh(i,j,k)
            spdiag(i,j,kr) = -etai(i,j,ki) * kh(i,j,k+1)
            diag  (i,j,kr) = c1
     &                     + etai(i,j,ki) * (kh(i,j,k) + kh(i,j,k+1))
            rhs   (i,j,kr) = Tin_init(i,j,ki)
     &                     + etai(i,j,ki)*Iswabs(i,j,ki)

         enddo                  ! ij
      enddo                     ! nilyr

      !-----------------------------------------------------------------
      ! bottom ice layer
      !-----------------------------------------------------------------

      ki = nilyr
      k  = ki + nslyr
      kr = k + 1
      
      do ij = 1, isolve
         i = indxii(ij)
         j = indxjj(ij)
         
         sbdiag(i,j,kr) = -etai(i,j,ki) * kh(i,j,k)
         spdiag(i,j,kr) = c0
         diag  (i,j,kr) = c1 
     &                  + etai(i,j,ki) * (kh(i,j,k) + kh(i,j,k+1))
         rhs   (i,j,kr) = Tin_init(i,j,ki)
     &                  + etai(i,j,ki)*Iswabs(i,j,ki)
     &                  + etai(i,j,ki)*kh(i,j,k+1)*Tbot(i,j)
           
      enddo                     ! ij

      end subroutine get_matrix_elements

!=======================================================================
!BOP
!
! !ROUTINE: tridiag_solver - tridiagonal matrix solver
!
! !DESCRIPTION:
!
! Tridiagonal matrix solver--used to solve the implicit vertical heat
! equation in ice and snow
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine tridiag_solver (nx_block, ny_block,
     &                           isolve,     
     &                           indxii,   indxjj,
     &                           nmat,     sbdiag,
     &                           diag,     spdiag, 
     &                           rhs,      xout)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  isolve          ! number of cells with temps not converged

      integer (kind=int_kind), dimension(nx_block*ny_block),
     &   intent(in) ::
     &   indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), intent(in) ::
     &   nmat            ! matrix dimension

      real (kind=dbl_kind), dimension (nx_block,ny_block,nmat),
     &     intent(in) ::
     &   sbdiag          ! sub-diagonal matrix elements
     &,  diag            ! diagonal matrix elements
     &,  spdiag          ! super-diagonal matrix elements
     &,  rhs             ! rhs of tri-diagonal matrix eqn.

      real (kind=dbl_kind), dimension (nx_block,ny_block,nmat),
     &     intent(inout) ::
     &   xout            ! solution vector
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! row counter

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   wbeta           ! temporary matrix variable

      real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr+nslyr+1)::
     &   wgamma          ! temporary matrix variable

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, isolve
         i = indxii(ij)
         j = indxjj(ij)

         wbeta(i,j) = diag(i,j,1)
         xout(i,j,1) = rhs(i,j,1) / wbeta(i,j)

      enddo                     ! ij

      do k = 2, nmat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

            wgamma(i,j,k) = spdiag(i,j,k-1) / wbeta(i,j)
            wbeta(i,j) = diag(i,j,k) - sbdiag(i,j,k)*wgamma(i,j,k)
            xout(i,j,k) = (rhs(i,j,k) - sbdiag(i,j,k)*xout(i,j,k-1))
     &                   / wbeta(i,j)

         enddo                  ! ij
      enddo                     ! k

      do k = nmat-1, 1, -1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)

            xout(i,j,k) = xout(i,j,k) - wgamma(i,j,k+1)*xout(i,j,k+1)

         enddo                  ! ij
      enddo                     ! k

      end subroutine tridiag_solver

!=======================================================================
!BOP
!
! !ROUTINE: thickness changes - top and bottom growth/melting
!
! !DESCRIPTION:
!
! Compute growth and/or melting at the top and bottom surfaces.
! Convert snow to ice if necessary.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine thickness_changes (nx_block,  ny_block,
     &                              dt,
     &                              yday,      icells,  
     &                              indxi,     indxj,     
     &                              aicen,     efinal, 
     &                              hin,       hilyr,   
     &                              hsn,       hslyr,
     &                              qin,       qsn,
     &                              fbot,      Tbot,
     &                              flatn,     fsurf,
     &                              fcondtop,  fcondbot,
     &                              fsnow,     hsn_new, 
     &                              fhocnn,    evapn,
     &                              meltt,     meltb,
     &                              congel,    snoice,   
     &                              mlt_onset, frz_onset)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block),
     &   intent(in) ::
     &   indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), intent(in) ::
     &   dt              ! time step
     &,  yday            ! day of the year

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aicen           ! fractional concentration of ice
     &,  fbot            ! ice-ocean heat flux at bottom surface (W/m^2)
     &,  Tbot            ! ice bottom surface temperature (deg C)
     &,  fsnow           ! snowfall rate (kg m-2 s-1)
     &,  flatn           ! surface downward latent heat (W m-2)
     &,  fsurf           ! net flux to top surface, not including fcondtop
     &,  fcondtop        ! downward cond flux at top surface (W m-2)
     &,  fcondbot        ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(inout) ::
     &   qin             ! ice layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr),
     &   intent(inout) ::
     &   qsn             ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   hilyr           ! ice layer thickness (m)
     &,  hslyr           ! snow layer thickness (m)
     &,  hin             ! total ice thickness (m)
     &,  hsn             ! total snow thickness (m)
     &,  meltt           ! top ice melt             (m/step-->cm/day)
     &,  meltb           ! basal ice melt           (m/step-->cm/day)
     &,  congel          ! basal ice growth         (m/step-->cm/day)
     &,  snoice          ! snow-ice formation       (m/step-->cm/day)
     &,  mlt_onset       ! day of year that sfc melting begins
     &,  frz_onset       ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out)::
     &   efinal          ! final energy of melting (J m-2)
     &,  hsn_new         ! thickness of new snow (m)
     &,  fhocnn          ! fbot, corrected for any surplus energy (W m-2)
     &,  evapn           ! ice/snow mass sublimated/condensed (kg m-2 s-1)
!
!EOP
!
      real (kind=dbl_kind), parameter ::
     &   qbotmax = -p5*rhoi*Lfresh  ! max enthalpy of ice growing at bottom

      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! vertical index

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   esub            ! energy for sublimation, > 0    (J m-2)
     &,  econ            ! energy for condensation, < 0   (J m-2)
     &,  etop_mlt        ! energy for top melting, > 0    (J m-2)
     &,  ebot_mlt        ! energy for bottom melting, > 0 (J m-2)
     &,  ebot_gro        ! energy for bottom growth, < 0  (J m-2)

      real (kind=dbl_kind) ::
     &   dhi             ! change in ice thickness
     &,  dhs             ! change in snow thickness
     &,  Ti              ! ice temperature
     &,  Ts              ! snow temperature
     &,  qbot            ! enthalpy of ice growing at bottom surface (J m-3)
     &,  qsub            ! energy/unit volume to sublimate ice/snow (J m-3)
     &,  hqtot           ! sum of h*q for two layers
     &,  wk1             ! temporary variable
     &,  qsnew           ! enthalpy of new snow (J m-3)
     &,  hstot           ! snow thickness including new snow (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr+1) ::
     &   zi1             ! depth of ice layer boundaries (m)
     &,  zi2             ! adjusted depths, with equal hilyr (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr+1) ::
     &   zs1             ! depth of snow layer boundaries (m)
     &,  zs2             ! adjusted depths, with equal hslyr (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr) ::
     &   dzi             ! ice layer thickness after growth/melting

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) ::
     &   dzs             ! snow layer thickness after growth/melting

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      efinal  (:,:) = c0
      hsn_new (:,:) = c0
      fhocnn  (:,:) = c0
      evapn   (:,:) = c0

      esub    (:,:) = c0
      econ    (:,:) = c0
      etop_mlt(:,:) = c0
      ebot_mlt(:,:) = c0
      ebot_gro(:,:) = c0

      do k = 1, nilyr+1
         zi1(:,:,k) = c0
         zi2(:,:,k) = c0
      enddo

      do k = 1, nslyr+1
         zs1(:,:,k) = c0
         zs2(:,:,k) = c0
      enddo

      do k = 1, nilyr
         dzi(:,:,k) = c0
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            dzi(i,j,k) = hilyr(i,j)
         enddo
      enddo

      do k = 1, nslyr
         dzs(:,:,k) = c0
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            dzs(i,j,k) = hslyr(i,j)
         enddo
      enddo

      !-----------------------------------------------------------------
      ! For l_brine = false (fresh ice), check for temperatures > 0.
      !  Melt ice or snow as needed to bring temperatures back to 0.
      ! For l_brine = true, this should not be necessary.
      !-----------------------------------------------------------------

      if (.not. l_brine) then 

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               Ts = (Lfresh + qsn(i,j,k)/rhos) / cp_ice
               if (Ts > c0) then
                  dhs = cp_ice*Ts*dzs(i,j,k) / Lfresh
                  dzs(i,j,k) = dzs(i,j,k) - dhs
                  qsn(i,j,k) = -rhos*Lfresh
               endif
            enddo
         enddo

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               Ti = (Lfresh + qin(i,j,k)/rhoi) / cp_ice
               if (Ti > c0) then
                  dhi = cp_ice*Ti*dzi(i,j,k) / Lfresh
                  dzi(i,j,k) = dzi(i,j,k) - dhi
                  qin(i,j,k) = -rhoi*Lfresh
               endif
            enddo               ! ij
         enddo                  ! k

      endif                     ! .not. l_brine

      !-----------------------------------------------------------------
      ! Compute energy available for sublimation/condensation, top melt,
      ! and bottom growth/melt.
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         wk1 = -flatn(i,j) * dt
         esub(i,j) = max(wk1, c0)     ! energy for sublimation, > 0
         econ(i,j) = min(wk1, c0)     ! energy for condensation, < 0

         wk1 = (fsurf(i,j) - fcondtop(i,j)) * dt
         etop_mlt(i,j) = max(wk1, c0)           ! etop_mlt > 0

         wk1 = (fcondbot(i,j) - fbot(i,j)) * dt
         ebot_mlt(i,j) = max(wk1, c0)           ! ebot_mlt > 0
         ebot_gro(i,j) = min(wk1, c0)           ! ebot_gro < 0

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Sublimation/condensation of ice/snow
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         !--------------------------------------------------------------
         ! Condensation (evapn > 0)
         ! Note: evapn here has unit of kg/m^2.  Divide by dt later.
         !--------------------------------------------------------------

         if (hsn(i,j) > puny) then   ! add snow with enthalpy qsn(i,j,1)
            dhs = econ(i,j) / (qsn(i,j,1) - rhos*Lvap) ! econ < 0, dhs > 0
            dzs(i,j,1) = dzs(i,j,1) + dhs
            evapn(i,j) = evapn(i,j) + dhs*rhos
         else                         ! add ice with enthalpy qin(i,j,1)
            dhi = econ(i,j) / (qin(i,j,1) - rhoi*Lvap) ! econ < 0, dhi > 0
            dzi(i,j,1) = dzi(i,j,1) + dhi
            evapn(i,j) = evapn(i,j) + dhi*rhoi
         endif
      enddo                     ! ij

         !--------------------------------------------------------------
         ! Sublimation of snow (evapn < 0)
         !--------------------------------------------------------------

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            qsub = qsn(i,j,k) - rhos*Lvap ! qsub < 0
            dhs  = max (-dzs(i,j,k), esub(i,j)/qsub)  ! esub > 0, dhs < 0
            dzs(i,j,k) = dzs(i,j,k) + dhs
            esub(i,j) = esub(i,j) - dhs*qsub
            esub(i,j) = max(esub(i,j), c0)   ! in case of roundoff error
            evapn(i,j) = evapn(i,j) + dhs*rhos
         enddo                  ! ij
      enddo                     ! nslyr

         !--------------------------------------------------------------
         ! Sublimation of ice (evapn < 0)
         !--------------------------------------------------------------

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            qsub = qin(i,j,k) - rhoi*Lvap              ! qsub < 0
            dhi  = max (-dzi(i,j,k), esub(i,j)/qsub) ! esub < 0, dhi < 0
            dzi(i,j,k) = dzi(i,j,k) + dhi
            esub(i,j) = esub(i,j) - dhi*qsub
            esub(i,j) = max(esub(i,j), c0)
            evapn(i,j) = evapn(i,j) + dhi*rhoi
         enddo                  ! ij
      enddo                     ! nilyr

      !-----------------------------------------------------------------
      ! Top melt
      ! (There is no top growth because there is no liquid water to
      !  freeze at the top.)
      !-----------------------------------------------------------------

         !--------------------------------------------------------------
         ! Melt snow
         !--------------------------------------------------------------

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            dhs = max(-dzs(i,j,k), etop_mlt(i,j)/qsn(i,j,k))
            dzs(i,j,k) = dzs(i,j,k) + dhs         ! qsn < 0, dhs < 0
            etop_mlt(i,j) = etop_mlt(i,j) - dhs*qsn(i,j,k)
            etop_mlt(i,j) = max(etop_mlt(i,j), c0) ! in case of roundoff error

            ! history diagnostics
            if (dhs < -puny .and. mlt_onset(i,j) < puny)
     &         mlt_onset(i,j) = yday

         enddo                  ! ij
      enddo                     ! nslyr

         !--------------------------------------------------------------
         ! Melt ice
         !--------------------------------------------------------------

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            dhi = max(-dzi(i,j,k), etop_mlt(i,j)/qin(i,j,k))
            dzi(i,j,k) = dzi(i,j,k) + dhi         ! qin < 0, dhi < 0
            etop_mlt(i,j) = etop_mlt(i,j) - dhi*qin(i,j,k)
            etop_mlt(i,j) = max(etop_mlt(i,j), c0)

            ! history diagnostics
            if (dhi < -puny .and. mlt_onset(i,j) < puny)
     &           mlt_onset(i,j) = yday
            meltt(i,j) = meltt(i,j) - dhi*aicen(i,j)

         enddo                  ! ij
      enddo                     ! k

       !----------------------------------------------------------------
       ! Bottom growth and melt
       !----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         !--------------------------------------------------------------
         ! Grow ice
         !--------------------------------------------------------------

         ! enthalpy of new ice growing at bottom surface
         qbot = -rhoi * (cp_ice * (Tmlt(nilyr+1)-Tbot(i,j))
     &                 + Lfresh * (c1-Tmlt(nilyr+1)/Tbot(i,j))
     &                 - cp_ocn * Tmlt(nilyr+1))
         qbot = min (qbot, qbotmax)      ! in case Tbot is close to Tmlt
         dhi  = ebot_gro(i,j) / qbot     ! dhi > 0

         hqtot = dzi(i,j,nilyr)*qin(i,j,nilyr) + dhi*qbot
         dzi(i,j,nilyr) = dzi(i,j,nilyr) + dhi

         if (dzi(i,j,nilyr) > puny)
     &        qin(i,j,nilyr) = hqtot / dzi(i,j,nilyr)

         ! history diagnostics
         congel(i,j) = congel(i,j) + dhi*aicen(i,j)
         if (dhi > puny .and. frz_onset(i,j) < puny)
     &           frz_onset(i,j) = yday

      enddo                     ! ij

         !--------------------------------------------------------------
         ! Melt ice
         !--------------------------------------------------------------

      do k = nilyr, 1, -1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            dhi = max(-dzi(i,j,k), ebot_mlt(i,j)/qin(i,j,k))
            dzi(i,j,k) = dzi(i,j,k) + dhi         ! qin < 0, dhi < 0
            ebot_mlt(i,j) = ebot_mlt(i,j) - dhi*qin(i,j,k)
            ebot_mlt(i,j) = max(ebot_mlt(i,j), c0)

            ! history diagnostics
            meltb(i,j) = meltb(i,j) - dhi*aicen(i,j)

         enddo                  ! ij
      enddo                     ! nilyr

         !--------------------------------------------------------------
         ! Melt snow (only if all the ice has melted)
         !--------------------------------------------------------------

      do k = nslyr, 1, -1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            dhs = max(-dzs(i,j,k), ebot_mlt(i,j)/qsn(i,j,k))
            dzs(i,j,k) = dzs(i,j,k) + dhs         ! qsn < 0, dhs < 0
            ebot_mlt(i,j) = ebot_mlt(i,j) - dhs*qsn(i,j,k)
            ebot_mlt(i,j) = max(ebot_mlt(i,j), c0)

         enddo                  ! ij
      enddo                     ! nslyr


      !-----------------------------------------------------------------
      ! Compute heat flux used by the ice (<=0).
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         fhocnn(i,j) = fbot(i,j)
     &               + (esub(i,j) + etop_mlt(i,j) + ebot_mlt(i,j))/dt
      enddo

!---!-----------------------------------------------------------------
!---! Add new snowfall at top surface.
!---!-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !----------------------------------------------------------------
      ! NOTE: If heat flux diagnostics are to work, new snow should
      !       have T = 0 (i.e. q = -rhos*Lfresh) and should not be
      !       converted to rain.
      !----------------------------------------------------------------

         if (fsnow(i,j) > c0) then

            hsn_new(i,j) = fsnow(i,j)/rhos * dt
            qsnew = -rhos*Lfresh
            hstot = dzs(i,j,1) + hsn_new(i,j)

            if (hstot > c0) then
               qsn(i,j,1) =  (dzs(i,j,1) * qsn(i,j,1)
     &                    + hsn_new(i,j) * qsnew) / hstot
               ! avoid roundoff errors
               qsn(i,j,1) = min(qsn(i,j,1), -rhos*Lfresh)

               dzs(i,j,1) = hstot
            endif
         endif

      enddo                     ! ij

    !-----------------------------------------------------------------
    ! Find the new ice and snow thicknesses.
    !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         hin(i,j) = c0
         hsn(i,j) = c0
      enddo

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            hin(i,j) = hin(i,j) + dzi(i,j,k)
         enddo                  ! ij
      enddo                     ! k

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            hsn(i,j) = hsn(i,j) + dzs(i,j,k)
         enddo                  ! ij
      enddo                     ! k

    !-------------------------------------------------------------------
    ! Convert snow to ice if snow lies below freeboard.
    !-------------------------------------------------------------------

      call freeboard (nx_block, ny_block,
     &                icells, 
     &                indxi,    indxj,
     &                aicen,    snoice,
     &                hin,      hsn,
     &                qin,      qsn,
     &                dzi,      dzs)

!---!-------------------------------------------------------------------
!---! Repartition the ice and snow into equal-thickness layers,
!---! conserving energy.
!---!-------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Compute desired layer thicknesses.
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (hin(i,j) > c0) then
            hilyr(i,j) = hin(i,j) / real(nilyr)
         else
            hin(i,j) = c0
            hilyr(i,j) = c0
         endif

         if (hsn(i,j) > c0) then
            hslyr(i,j) = hsn(i,j) / real(nslyr)
         else
            hsn(i,j) = c0
            hslyr(i,j) = c0
         endif
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Compute depths zi1 of old layers (unequal thickness).
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         zi1(i,j,1) = c0
         zi1(i,j,1+nilyr) = hin(i,j)
      enddo

      do k = 1, nilyr-1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            zi1(i,j,k+1) = zi1(i,j,k) + dzi(i,j,k)
         end do
      enddo

      !-----------------------------------------------------------------
      ! Compute depths zi2 of new layers (equal thickness).
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         zi2(i,j,1) = c0
         zi2(i,j,1+nilyr) = hin(i,j)
      enddo

      do k = 1, nilyr-1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            zi2(i,j,k+1) = zi2(i,j,k) + hilyr(i,j)
         end do
      enddo

      !-----------------------------------------------------------------
      ! Conserving energy, compute the enthalpy of the new equal layers.
      !-----------------------------------------------------------------

      call adjust_enthalpy (nx_block, ny_block,
     &                      nilyr,    icells, 
     &                      indxi,    indxj,
     &                      zi1,      zi2,   
     &                      hilyr,    hin,
     &                      qin)

      if (nslyr > 1) then

      !-----------------------------------------------------------------
      ! Compute depths zs1 of old layers (unequal thickness).
      !-----------------------------------------------------------------

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            zs1(i,j,1) = c0
            zs1(i,j,1+nslyr) = hsn(i,j)
         enddo

         do k = 1, nslyr-1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               zs1(i,j,k+1) = zs1(i,j,k) + dzs(i,j,k)
            end do
         enddo

      !-----------------------------------------------------------------
      ! Compute depths zs2 of new layers (equal thickness).
      !-----------------------------------------------------------------

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            zs2(i,j,1) = c0
            zs2(i,j,1+nslyr) = hsn(i,j)
         enddo

         do k = 1, nslyr-1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               zs2(i,j,k+1) = zs2(i,j,k) + hslyr(i,j)
            end do
         enddo

      !-----------------------------------------------------------------
      ! Conserving energy, compute the enthalpy of the new equal layers.
      !-----------------------------------------------------------------

         call adjust_enthalpy (nx_block, ny_block,
     &                         nslyr,    icells, 
     &                         indxi,    indxj,
     &                         zs1,      zs2,
     &                         hslyr,    hsn,
     &                         qsn)

      endif   ! nslyr > 1

      !-----------------------------------------------------------------
      ! Compute final ice-snow energy, including the energy of
      !  sublimated/condensed ice.
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         efinal(i,j) = -evapn(i,j)*Lvap
         evapn(i,j)  =  evapn(i,j)/dt
      enddo

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            efinal(i,j) = efinal(i,j) + hslyr(i,j)*qsn(i,j,k)
         enddo                  ! ij
      enddo

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            efinal(i,j) = efinal(i,j) + hilyr(i,j)*qin(i,j,k)
         enddo                  ! ij
      enddo                     ! k

      end subroutine thickness_changes

!=======================================================================
!BOP
!
! !ROUTINE: freeboard - snow-ice conversion
!
! !DESCRIPTION:
!
! If there is enough snow to lower the ice/snow interface below
! sea level, convert enough snow to ice to bring the interface back
! to sea level.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine freeboard (nx_block, ny_block,
     &                      icells, 
     &                      indxi,    indxj,
     &                      aicen,    snoice,
     &                      hin,      hsn,
     &                      qin,      qsn,
     &                      dzi,      dzs)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  icells              ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block),
     &   intent(in) ::
     &   indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aicen           ! fractional ice area

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   snoice      ! snow-ice formation       (m/step-->cm/day)
     &,  hin         ! ice thickness (m)
     &,  hsn         ! snow thickness (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(inout) ::
     &   qin         ! ice layer enthalpy (J m-3)
     &,  dzi         ! ice layer thicknesses (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr),
     &   intent(in) ::
     &   qsn         ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr),
     &   intent(inout) ::
     &   dzs         ! snow layer thicknesses (m)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! vertical index

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   dhin            ! change in ice thickness (m)
     &,  dhsn            ! change in snow thickness (m)
     &,  hqs             ! sum of h*q for snow (J m-2)

      real (kind=dbl_kind) ::
     &   wk1             ! temporary variable
     &,  dhs             ! snow to remove from layer (m)

      !-----------------------------------------------------------------
      ! Determine whether snow lies below freeboard.
      !-----------------------------------------------------------------

      dhin(:,:) = c0
      dhsn(:,:) = c0
      hqs (:,:) = c0

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         wk1 = hsn(i,j) - hin(i,j)*(rhow-rhoi)/rhos

         if (wk1 > puny .and. hsn(i,j) > puny) then  ! snow below freeboard
            dhsn(i,j) = min(wk1*rhoi/rhow, hsn(i,j)) ! snow to remove
            dhin(i,j) = dhsn(i,j) * rhos/rhoi        ! ice to add
         endif
      enddo

      !-----------------------------------------------------------------
      ! Adjust snow layer thickness.
      ! Compute energy to transfer from snow to ice.
      !-----------------------------------------------------------------

      do k = nslyr, 1, -1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (dhin(i,j) > puny) then
               dhs = min(dhsn(i,j), dzs(i,j,k)) ! snow to remove from layer
               hsn(i,j) = hsn(i,j) - dhs
               dzs(i,j,k) = dzs(i,j,k) - dhs
               dhsn(i,j) = dhsn(i,j) - dhs
               dhsn(i,j) = max(dhsn(i,j),c0)
               hqs(i,j) = hqs(i,j) + dhs * qsn(i,j,k)
            endif               ! dhin > puny
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Transfer volume and energy from snow to top ice layer.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (dhin(i,j) > puny) then
            wk1 = dzi(i,j,1) + dhin(i,j)
            hin(i,j) = hin(i,j) + dhin(i,j)
            qin(i,j,1) = (dzi(i,j,1)*qin(i,j,1) + hqs(i,j)) / wk1
            dzi(i,j,1) = wk1

            ! history diagnostic
            snoice(i,j) = snoice(i,j) + dhin(i,j)*aicen(i,j)
         endif               ! dhin > puny

      enddo                  ! ij

      end subroutine freeboard

!=======================================================================
!BOP
!
! !ROUTINE: adjust_enthalpy -- enthalpy of new layers
!
! !DESCRIPTION:
!
! Conserving energy, compute the new enthalpy of equal-thickness ice
! or snow layers.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine adjust_enthalpy (nx_block, ny_block,
     &                            nlyr,     icells, 
     &                            indxi,    indxj,
     &                            z1,       z2,   
     &                            hlyr,     hn,
     &                            qn)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nlyr                ! number of layers (nilyr or nslyr)
     &,  icells              ! number of cells with aicen > puny

      integer (kind=int_kind), dimension (nx_block*ny_block),
     &   intent(in) ::
     &   indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block,nlyr+1),
     &   intent(in) ::
     &   z1              ! interface depth for old, unequal layers (m)
     &,  z2              ! interface depth for new, equal layers (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   hlyr            ! new layer thickness (m)
     &,  hn              ! total thickness (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nlyr),
     &   intent(inout) ::
     &   qn              ! layer enthalpy (J m-3)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k, k1, k2       ! vertical indices

      real (kind=dbl_kind) ::
     &   hovlp           ! overlap between old and new layers (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   rhlyr           ! 1./hlyr

      real (kind=dbl_kind), dimension (nx_block,ny_block,nlyr) ::
     &   hq              ! h * q for a layer

      rhlyr(:,:) = c0
      do k = 1, nlyr
         hq(:,:,k) = c0
      enddo

      !-----------------------------------------------------------------
      ! Compute reciprocal layer thickness.
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (hn(i,j) > puny) rhlyr(i,j) = c1 / hlyr(i,j)
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Compute h*q for new layers (k2) given overlap with old layers (k1)
      !-----------------------------------------------------------------

      do k2 = 1, nlyr

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            hq(i,j,k2) = c0
         enddo

         do k1 = 1, nlyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               hovlp = min (z1(i,j,k1+1), z2(i,j,k2+1))
     &               - max (z1(i,j,k1),   z2(i,j,k2))
               hovlp = max (hovlp, c0)

               hq(i,j,k2) = hq(i,j,k2) + hovlp*qn(i,j,k1)
            enddo               ! ij
         enddo                  ! kold
      enddo                     ! k

      !-----------------------------------------------------------------
      ! Compute new enthalpies.
      !-----------------------------------------------------------------

      do k = 1, nlyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            qn(i,j,k) = hq(i,j,k) * rhlyr(i,j)
         enddo                  ! ij
      enddo                     ! k

      end subroutine adjust_enthalpy

!=======================================================================
!BOP
!
! !ROUTINE: conservation_check_vthermo - energy conservation check
!
! !DESCRIPTION:
!
! Check for energy conservation by comparing the change in energy
! to the net energy input.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine conservation_check_vthermo(nx_block, ny_block,
     &                                      my_task,  istep1,
     &                                      dt,       icells, 
     &                                      indxi,    indxj,
     &                                      fsurf,    flatn,
     &                                      fhocnn,   fswint,
     &                                      fsnow,
     &                                      einit,    efinal,
     &                                      l_stop,
     &                                      istop,    jstop)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  my_task             ! task number (diagnostic only)
     &,  istep1              ! time step index (diagnostic only)
     &,  icells              ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block),
     &   intent(in) ::
     &   indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), intent(in) ::
     &   dt              ! time step

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   fsurf           ! net flux to top surface, not including fcondtop
     &,  flatn           ! surface downward latent heat (W m-2)
     &,  fhocnn          ! fbot, corrected for any surplus energy
     &,  fswint          ! SW absorbed in ice interior, below surface (W m-2)
     &,  fsnow           ! snowfall rate (kg m-2 s-1)
     &,  einit           ! initial energy of melting (J m-2)
     &,  efinal          ! final energy of melting (J m-2)

      logical (kind=log_kind), intent(inout) ::
     &   l_stop          ! if true, print diagnostics and abort model

      integer (kind=int_kind), intent(out) ::
     &   istop, jstop    ! i and j indices of cell where model fails
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops

      real (kind=dbl_kind) ::
     &   einp            ! energy input during timestep (J m-2)
     &,  ferr            ! energy conservation error (W m-2)

      logical (kind=log_kind) ::   ! for vector-friendly error checks
     &   ferr_flag       ! flag for energy error, ferr > ferrmax

      ferr_flag = .false.     ! ferr <= ferrmax, initialization

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         einp = (fsurf(i,j) - flatn(i,j) + fswint(i,j) - fhocnn(i,j)
     &         - fsnow(i,j)*Lfresh) * dt
         ferr = abs(efinal(i,j)-einit(i,j)-einp) / dt
         if (ferr > ferrmax) ferr_flag = .true.
      enddo

      !----------------------------------------------------------------
      ! If energy is not conserved, print diagnostics and exit.
      !----------------------------------------------------------------
      if (ferr_flag) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Note that fsurf - flat = fsw + flw + fsens; i.e., the latent
      ! heat is not included in the energy input, since (efinal - einit)
      ! is the energy change in the system ice + vapor, and the latent
      ! heat lost by the ice is equal to that gained by the vapor.
      !-----------------------------------------------------------------

            einp = (fsurf(i,j) - flatn(i,j) + fswint(i,j) - fhocnn(i,j)
     &             -fsnow(i,j)*Lfresh) * dt
            ferr = abs(efinal(i,j)-einit(i,j)-einp) / dt

            if (ferr > ferrmax) then
               write(nu_diag,*) 'Thermo energy conservation error'
               write(nu_diag,*) 'istep1, my_task, i, j:',
     &                           istep1, my_task, i, j
               write(nu_diag,*) 'Flux error (W/m^2) =', ferr
               write(nu_diag,*) 'Energy error (J) =', ferr*dt
               write(nu_diag,*) 'Initial energy =', einit(i,j)
               write(nu_diag,*) 'Final energy =', efinal(i,j)
               write(nu_diag,*) 'efinal - einit =',
     &                           efinal(i,j)-einit(i,j)
               write(nu_diag,*) 'Input energy =', einp
               l_stop = .true.
               istop = i
               jstop = j
               return
            endif
         enddo
      endif

      end subroutine conservation_check_vthermo

!=======================================================================
!BOP
!
! !ROUTINE: update_state_vthermo - new state variables
!
! !DESCRIPTION:
!
! Given the vertical thermo state variables (hin, hsn, qin,
!  qsn, Tsf), compute the new ice state variables (vicen, vsnon,
!  eicen, esnon, Tsfcn).
! Zero out state variables if ice has melted entirely.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine update_state_vthermo (nx_block, ny_block,
     &                                 icells, 
     &                                 indxi,    indxj,
     &                                 Tf,       Tsf,
     &                                 hin,      hsn,
     &                                 qin,      qsn,
     &                                 aicen,    vicen,
     &                                 vsnon,    Tsfcn,
     &                                 eicen,    esnon)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  icells              ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block),
     &   intent(in) ::
     &   indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::
     &   Tf              ! freezing temperature (C)
     &,  Tsf             ! ice/snow surface temperature, Tsfcn
     &,  hin             ! ice thickness (m)
     &,  hsn             ! snow thickness (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(in) ::
     &   qin             ! ice layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr),
     &   intent(in) ::
     &   qsn             ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   aicen           ! concentration of ice
     &,  vicen           ! volume per unit area of ice          (m)
     &,  vsnon           ! volume per unit area of snow         (m)
     &,  Tsfcn           ! temperature of ice/snow top surface  (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(inout) ::
     &   eicen           ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr),
     &   intent(inout) ::
     &   esnon           ! energy of melting for each snow layer (J/m^2)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! ice layer index

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (hin(i,j) > c0) then
            ! aicen is already up to date
            vicen(i,j) = aicen(i,j) * hin(i,j)
            vsnon(i,j) = aicen(i,j) * hsn(i,j)
            Tsfcn(i,j) = Tsf(i,j)
         else  ! (hin(i,j) == c0)
            aicen(i,j) = c0
            vicen(i,j) = c0
            vsnon(i,j) = c0
            Tsfcn(i,j) = Tf(i,j)
         endif

      enddo                     ! ij

      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (hin(i,j) > c0) then
               eicen(i,j,k) = qin(i,j,k) * vicen(i,j)/real(nilyr)
            else
               eicen(i,j,k) = c0
            endif

         enddo                  ! ij
      enddo                     ! nilyr

      do k = 1, nslyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (hin(i,j) > c0) then
               esnon(i,j,k) = qsn(i,j,k) * vsnon(i,j)/real(nslyr)
            else
               esnon(i,j,k) = c0
            endif

         enddo                  ! ij
      enddo                     ! nslyr

      end subroutine update_state_vthermo



!=======================================================================

      end module ice_therm_vertical

!=======================================================================
