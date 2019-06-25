c $Id: $
c=======================================================================
!---! flux variable declarations
!---! these include fields sent from the coupler ("in"),
!---!                      sent to the coupler ("out"),
!---!                      and used internally ("internal")
!---!
!---! author Elizabeth C. Hunke, LANL
c=======================================================================

      module ice_flux

      use ice_kinds_mod
      use ice_domain
      use ice_constants

      implicit none

      !-----------------------------------------------------------------
      ! Dynamics component
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
       ! in
     &   strtltx     ! sea surface slope, x-direction (m/m)
     &,  strtlty     ! sea surface slope, y-direction
     &,  uocn        ! ocean current, x-direction (m/s)
     &,  vocn        ! ocean current, y-direction (m/s)
       ! out          T-cell (kg/m s^2)  
     &,  strocnxT    ! ice-ocean stress, x-direction
     &,  strocnyT    ! ice-ocean stress, y-direction
     &,  strairxT    ! stress on ice by air, x-direction
     &,  strairyT    ! stress on ice by air, y-direction
       ! internal     U-cell (kg/m s^2)  
     &,  strocnx     ! ice-ocean stress, x-direction
     &,  strocny     ! ice-ocean stress, y-direction
     &,  strairx     ! stress on ice by air, x-direction
     &,  strairy     ! stress on ice by air, y-direction

      !-----------------------------------------------------------------
      ! Thermodynamics component
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
       ! input from atmo
     &   zlvl        ! atm level height (m)
     &,  uatm        ! wind velocity    (m/s)
     &,  vatm 
     &,  potT        ! air potential temperature  (K)
     &,  Tair        ! air temperature  (K)
     &,  Qa          ! specific humidity (kg/kg)
     &,  rhoa        ! air density (kg/m^3)
     &,  swvdr       ! sw down, visible, direct  (W/m^2)
     &,  swvdf       ! sw down, visible, diffuse (W/m^2)
     &,  swidr       ! sw down, near IR, direct  (W/m^2)
     &,  swidf       ! sw down, near IR, diffuse (W/m^2)
     &,  Flw         ! incoming longwave radiation (W/m^2)
     &,  Frain       ! rainfall rate (kg/m^2 s)
     &,  Fsnow       ! snowfall rate (kg/m^2 s liquid equivalent)
       ! input from ocean
     &,  frzmlt      ! freezing/melting potential (W/m^2)
     &,  sss         ! sea surface salinity (ppt)
     &,  sst         ! sea surface temperature (C)
     &,  Tf          ! freezing temperature (C)
       ! out !! note albedos are in ice_albedo.H, Tsfc in ice_state.H
     &,  Fsensible   ! sensible heat flux (W/m^2)
     &,  Flatent     ! latent heat flux   (W/m^2)
     &,  Flwout      ! outgoing longwave radiation (W/m^2)
     &,  evap        ! evaporativewater flux (m/s)
     &,  Tref        ! 2m atm reference temperature (K)
     &,  Fresh       ! fresh water flux (kg/m^2 s)
     &,  Fsalt       ! net salt flux to ocean
     &,  Fhnet       ! net heat flux to ocean
     &,  Fswthru     ! shortwave penetrating to ocean
       ! internal
     &,  Fsw         ! incoming shortwave radiation (W/m^2)
     &,  wind        ! wind speed (m/s)
         ! diagnostic only
     &,  Frsh_hist   ! fresh water flux (kg/m^2 s) (Fresh for current dt)
     &,  Foht_hist   ! net heat flux to ocean (Fhnet for current dt)
     &,  Fswthru_hist! shortwave penetrating to ocean (Fswthru for current dt)
     &,  Fsalt_hist  ! salt flux to ocean (Fsalt for current dt)
         ! used by mixed layer model for restarts
     &,  Fhocn       ! net ocean heat flux towards ice

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ncat) ::
       ! internal
     &   shcoef      ! transfer coefficient for sensible heat
     &,  lhcoef      ! transfer coefficient for latent heat

c=======================================================================

      contains

c=======================================================================

      subroutine init_flux_atm(i,j)

!---!-------------------------------------------------------------------
!---! initialize all fluxes sent to coupler for use by the atm model
!---! and a few state quantities
!---!-------------------------------------------------------------------

      integer i,j

      !-----------------------------------------------------------------
      ! fluxes sent
      !-----------------------------------------------------------------
        Fsensible(i,j) = c0
        Flatent(i,j)  = c0
        Flwout(i,j)   = c0
        evap(i,j)     = c0
        strairxT(i,j) = c0
        strairyT(i,j) = c0
        Tref(i,j) = Tair(i,j)  ! over land and ice-free ocean

      !-----------------------------------------------------------------
      ! other miscellaneous fields
      !-----------------------------------------------------------------
        strairx(i,j) = c0
        strairy(i,j) = c0

      end subroutine init_flux_atm

c=======================================================================

      subroutine init_flux_ocn(i,j)

!---!-------------------------------------------------------------------
!---! initialize fluxes sent to coupler for use by the ocean model
!---!-------------------------------------------------------------------

      integer i,j

      !-----------------------------------------------------------------
      ! fluxes sent
      !-----------------------------------------------------------------
        Fresh(i,j)    = c0
        Fsalt(i,j)    = c0
        Fhnet(i,j)    = c0
        Fswthru(i,j)  = c0

      end subroutine init_flux_ocn

c=======================================================================

      subroutine merge_fluxes(i,j,nc,frshn,fsaltn,evapn,
     &      fshn,flhn,flwupn,swbot,Focn,Trefn,strxn,stryn)

!---!-------------------------------------------------------------------
!---! aggregates flux information from all ice thickness categories
!---!-------------------------------------------------------------------

      use ice_state

      integer (kind=int_kind), intent(in) :: i,j,nc
      real (kind=dbl_kind), intent(in) ::
     &    fshn    ! sensible heat flux            (W/m^2)
     &,   flhn    ! latent   heat flux            (W/m^2)
     &,   flwupn  ! upwd lw emitted heat flux     (W/m^2)
     &,   swbot   ! sw radiation through ice bot  (W/m^2)
     &,   Focn    ! actual ocn/ice heat flux      (W/m^2)
     &,   Trefn   ! air tmp rfrnc level               (K)
     &,   strxn   ! air/ice zonal  strss,         (N/m^2)
     &,   stryn   ! air/ice merdnl strss,         (N/m^2)
     &,   frshn   ! fresh water flux to ocean  (kg/m^2/s)
     &,   evapn   ! evaporation                (kg/m^2/s)
     &,  fsaltn   ! flux of salt, ice to ocean (kg/m^2/s)

       Fresh    (i,j) = Fresh    (i,j) + frshn *ain(nc)
       Fsalt    (i,j) = Fsalt    (i,j) + fsaltn*ain(nc)
       evap     (i,j) = evap     (i,j) + evapn *ain(nc)
       Fsensible(i,j) = Fsensible(i,j) + fshn  *ain(nc)
       Flatent  (i,j) = Flatent  (i,j) + flhn  *ain(nc)
       Flwout   (i,j) = Flwout   (i,j) + 
     &        (flwupn-(c1-emissivity)*Flw(i,j))*ain(nc)
       Fswthru  (i,j) = Fswthru  (i,j) + swbot *ain(nc)
       Fhnet    (i,j) = Fhnet    (i,j) + Focn  *ain(nc) 
       Tref     (i,j) = Tref     (i,j) + Trefn *ain(nc)

       strairxT (i,j) = strairxT (i,j) + strxn *ain(nc)
       strairyT (i,j) = strairyT (i,j) + stryn *ain(nc)

      end subroutine merge_fluxes

c=======================================================================

      end module ice_flux

c=======================================================================
