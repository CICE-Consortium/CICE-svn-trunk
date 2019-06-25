c $Id: ice_flux.F,v 1.20 2006/06/08 22:16:07 eclare Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_flux - flux variable declarations: coupler, diagnostic and internal
!
! !DESCRIPTION:
!
! Flux variable declarations; these include fields sent from the coupler
! ("in"), sent to the coupler ("out"), written to diagnostic history files
! ("diagnostic"), and used internally ("internal").
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      module ice_flux
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain
      use ice_constants
!
!EOP
!
      implicit none
      save

      !-----------------------------------------------------------------
      ! Dynamics component
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::

       ! in from ocean
     &   ss_tltx     ! sea surface slope, x-direction (m/m)
     &,  ss_tlty     ! sea surface slope, y-direction
     &,  uocn        ! ocean current, x-direction (m/s)
     &,  vocn        ! ocean current, y-direction (m/s)

       ! out to atmosphere

     &,  strairxT    ! stress on ice by air, x-direction
     &,  strairyT    ! stress on ice by air, y-direction

       ! out to ocean          T-cell (kg/m s^2)
     &,  strocnxT    ! ice-ocean stress, x-direction
     &,  strocnyT    ! ice-ocean stress, y-direction

       ! diagnostic
     &,  daidtd    ! ice area tendency due to transport   (s^-1)
     &,  dvidtd    ! ice volume tendency due to transport (m/s)

       ! internal              U-cell (kg/m s^2)
     &,  strocnx     ! ice-ocean stress, x-direction
     &,  strocny     ! ice-ocean stress, y-direction
     &,  strairx     ! stress on ice by air, x-direction
     &,  strairy     ! stress on ice by air, y-direction
     &,  strtltx     ! stress due to sea surface slope, x-direction
     &,  strtlty     ! stress due to sea surface slope, y-direction

      !-----------------------------------------------------------------
      ! Thermodynamics component
      !-----------------------------------------------------------------

       ! in from atmosphere

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   zlvl        ! atm level height (m)
     &,  uatm        ! wind speed (m/s)
     &,  vatm 
     &,  potT        ! air potential temperature  (K)
     &,  Tair        ! air temperature  (K)
     &,  Qa          ! specific humidity (kg/kg)
     &,  rhoa        ! air density (kg/m^3)
     &,  swvdr       ! sw down, visible, direct  (W/m^2)
     &,  swvdf       ! sw down, visible, diffuse (W/m^2)
     &,  swidr       ! sw down, near IR, direct  (W/m^2)
     &,  swidf       ! sw down, near IR, diffuse (W/m^2)
     &,  flw         ! incoming longwave radiation (W/m^2)
     &,  frain       ! rainfall rate (kg/m^2 s)
     &,  fsnow       ! snowfall rate (kg/m^2 s)

       ! in from ocean

     &,  frzmlt      ! freezing/melting potential (W/m^2)
     &,  sss         ! sea surface salinity (ppt)
     &,  sst         ! sea surface temperature (C)
     &,  Tf          ! freezing temperature (C)
     &,  qdp         ! deep ocean heat flux (W/m^2)
     &,  hmix        ! mixed layer depth (m)

       ! out to atmosphere 
       ! note albedos are in ice_albedo.F, Tsfc in ice_state.F

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   fsens       ! sensible heat flux (W/m^2)
     &,  flat        ! latent heat flux   (W/m^2)
     &,  fswabs      ! shortwave flux absorbed in ice and ocean (W/m^2)
     &,  flwout      ! outgoing longwave radiation (W/m^2)
     &,  evap        ! evaporative water flux (kg/m^2/s)
     &,  Tref        ! 2m atm reference temperature (K)
     &,  Qref        ! 2m atm reference sp humidity (kg/kg)

       ! out to ocean

     &,  fresh       ! fresh water flux to ocean (kg/m^2/s)
     &,  fsalt       ! salt flux to ocean (kg/m^2/s)
     &,  fhnet       ! net heat flux to ocean (W/m^2)
     &,  fswthru     ! shortwave penetrating to ocean (W/m^2)

       ! diagnostic 

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::       
     &   congel    ! basal ice growth         (m/step-->cm/day)
     &,  frazil    ! frazil ice growth        (m/step-->cm/day)
     &,  snoice    ! snow-ice formation       (m/step-->cm/day)
     &,  meltt     ! top ice melt             (m/step-->cm/day)
     &,  meltb     ! basal ice melt           (m/step-->cm/day)
     &,  meltl     ! lateral ice melt         (m/step-->cm/day)
     &,  daidtt    ! ice area tendency thermo.   (s^-1)
     &,  dvidtt    ! ice volume tendency thermo. (m/s)
     &,  mlt_onset ! day of year that sfc melting begins
     &,  frz_onset ! day of year that freezing begins (congel or frazil)

       ! NOTE: The following ocean diagnostic fluxes measure
       ! the same thing as their coupler counterparts but over
       ! different time intervals.  The coupler variables are 
       ! computed from one to_coupler call to the next, whereas
       ! the diagnostic variables are computed over a time step.
     &,  fresh_hist  ! fresh water flux to ocean (kg/m^2/s)
     &,  fsalt_hist  ! salt flux to ocean (kg/m^2/s)
     &,  fhnet_hist  ! net heat flux to ocean (W/m^2)
     &,  fswthru_hist! shortwave penetrating to ocean (W/m^2)

       ! internal

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   fsw         ! incoming shortwave radiation (W/m^2)
     &,  wind        ! wind speed (m/s)
     &,  shcoef      ! transfer coefficient for sensible heat
     &,  lhcoef      ! transfer coefficient for latent heat

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_flux_atm - initialize all atmospheric fluxes sent to coupler
!
! !INTERFACE:
!
      subroutine init_flux_atm
!
! !DESCRIPTION:
!
! Initialize all fluxes sent to coupler for use by the atm model
! and a few state quantities
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
      use ice_state, only : aice
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j    ! horizontal indices

      !-----------------------------------------------------------------
      ! fluxes sent
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi

         strairxT(i,j) = c0        ! wind stress, T grid
         strairyT(i,j) = c0
         fsens   (i,j) = c0
         flat    (i,j) = c0
         fswabs  (i,j) = c0
         flwout  (i,j) = c0
         evap    (i,j) = c0
         Tref    (i,j) = c0
         Qref    (i,j) = c0

      !-----------------------------------------------------------------
      ! other miscellaneous fields
      !-----------------------------------------------------------------
         strairx(i,j)  = c0        ! wind stress, U grid
         strairy(i,j)  = c0

      enddo
      enddo

      end subroutine init_flux_atm

!=======================================================================
!BOP
!
! !IROUTINE: init_flux_ocn - initialize ocean fluxes sent to coupler
!
! !INTERFACE:
!
      subroutine init_flux_ocn

!cdir$ inlinenever init_flux_ocn
!
! !DESCRIPTION:
!
! Initialize fluxes sent to coupler for use by the ocean model
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j    ! horizontal indices

      !-----------------------------------------------------------------
      ! fluxes sent
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         fresh(i,j)    = c0
         fsalt(i,j)    = c0
         fhnet(i,j)    = c0
         fswthru(i,j)  = c0

         qdp(i,j)      = c0
         hmix(i,j)     = c20
      enddo
      enddo

      end subroutine init_flux_ocn

!=======================================================================
!BOP
!
! !IROUTINE: init_diagnostics - initialize diagnostic fields
!
! !INTERFACE:
!
      subroutine init_diagnostics
!
! !DESCRIPTION:
!
! Initialize diagnostic fields written to history files.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      use ice_state, only: aice, vice

      integer (kind=int_kind) :: 
     &   i, j                   ! horizontal indices

      do j = jlo, jhi
      do i = ilo, ihi
         congel (i,j) = c0
         frazil (i,j) = c0
         snoice (i,j) = c0
         meltt  (i,j) = c0
         meltb  (i,j) = c0
         meltl  (i,j) = c0
         daidtt (i,j) = aice(i,j) ! temporarily used for initial area
         dvidtt (i,j) = vice(i,j) ! temporarily used for initial volume
         daidtd (i,j) = c0
         dvidtd (i,j) = c0
         fresh_hist(i,j) = c0
         fsalt_hist(i,j) = c0
         fhnet_hist(i,j) = c0
         fswthru_hist(i,j) = c0
      enddo                     ! i
      enddo                     ! j

      end subroutine init_diagnostics

!=======================================================================
!BOP
!
! !IROUTINE: merge_fluxes - aggregate flux information over ITD
!
! !INTERFACE:
!
      subroutine merge_fluxes (n,
     &     strxn,   stryn,  fsensn, flatn, fswabsn, 
     &     flwoutn, evapn,  Trefn,  Qrefn,
     &     freshn,  fsaltn, fhnetn, fswthrun)
!
! !DESCRIPTION:
!
! Aggregates flux information from all ice thickness categories
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &    n           ! thickness category index

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(in) ::
     &    strxn       ! air/ice zonal  strss,           (N/m**2)
     &,   stryn       ! air/ice merdnl strss,           (N/m**2)
     &,   fsensn      ! sensible heat flx               (W/m**2)
     &,   flatn       ! latent   heat flx               (W/m**2)
     &,   fswabsn     ! shortwave absorbed heat flx     (W/m**2)
     &,   flwoutn     ! upwd lw emitted heat flx        (W/m**2)
     &,   evapn       ! evaporation                     (kg/m2/s)
     &,   Trefn       ! air tmp reference level         (K)
     &,   Qrefn       ! air sp hum reference level      (kg/kg)
     &,   freshn      ! fresh water flux to ocean       (kg/m2/s)
     &,   fsaltn      ! salt flux to ocean              (kg/m2/s)
     &,   fhnetn      ! actual ocn/ice heat flx         (W/m**2)
     &,   fswthrun    ! sw radiation through ice bot    (W/m**2)
!
!EOP
!
      integer (kind=int_kind) ::
     &    i, j        ! horizontal indices

      do j = jlo,jhi
      do i = ilo,ihi

      ! atmo fluxes

         strairxT (i,j)  = strairxT (i,j)  + strxn  (i,j) * aicen(i,j,n)
         strairyT (i,j)  = strairyT (i,j)  + stryn  (i,j) * aicen(i,j,n)
         fsens    (i,j)  = fsens    (i,j)  + fsensn (i,j) * aicen(i,j,n)
         flat     (i,j)  = flat     (i,j)  + flatn  (i,j) * aicen(i,j,n)
         fswabs   (i,j)  = fswabs   (i,j)  + fswabsn(i,j) * aicen(i,j,n)
         flwout   (i,j)  = flwout   (i,j)  + 
     &           (flwoutn(i,j) - (c1-emissivity)*flw(i,j))* aicen(i,j,n)
         evap     (i,j)  = evap     (i,j)  + evapn  (i,j) * aicen(i,j,n)
         Tref     (i,j)  = Tref     (i,j)  + Trefn  (i,j) * aicen(i,j,n)
         Qref     (i,j)  = Qref     (i,j)  + Qrefn  (i,j) * aicen(i,j,n)

      ! ocean fluxes: update both coupler and history variables

         fresh     (i,j) = fresh     (i,j) + freshn(i,j)  * aicen(i,j,n)
         fresh_hist(i,j) = fresh_hist(i,j) + freshn(i,j)  * aicen(i,j,n)
         fsalt     (i,j) = fsalt     (i,j) + fsaltn(i,j)  * aicen(i,j,n)
         fsalt_hist(i,j) = fsalt_hist(i,j) + fsaltn(i,j)  * aicen(i,j,n)
         fhnet     (i,j) = fhnet     (i,j) + fhnetn(i,j)  * aicen(i,j,n)
         fhnet_hist(i,j) = fhnet_hist(i,j) + fhnetn(i,j)  * aicen(i,j,n)
         fswthru   (i,j) = fswthru(i,j)  + fswthrun(i,j)  * aicen(i,j,n)
         fswthru_hist(i,j) = fswthru_hist(i,j) 
     &                                   + fswthrun(i,j)  * aicen(i,j,n)
      enddo                     ! i
      enddo                     ! j

      end subroutine merge_fluxes

!=======================================================================

      end module ice_flux

!=======================================================================
