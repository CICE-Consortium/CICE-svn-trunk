c $Id: ice_albedo.F,v 1.8 2004/02/09 17:58:46 lipscomb Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_shortwave - reflected, transmitted, and absorbed solar rad
!
! !DESCRIPTION:
!
! The albedo and absorbed/transmitted flux parameterizations for
! snow over ice, bare ice and ponded ice.
!
! Presently, two methods are included:
!   (1) CCSM3 
!   (2) Delta-Eddington 
! as two distinct routines.
! Either can be called from the ice driver.
!
! The delta-Eddington method is described here:
!
! Briegleb, B. P., and B. Light (2006): A delta-Eddington parameterization
!    for sea ice shortwave radiation in the Community Climate System
!    Model, NCAR Technical Note, NCAR/TN-?+STR, National Center for 
!    Atmospheric Research, Boulder, Colo., ?? pp.
! [Still in draft form as of May 2006 - WHL]
!
! !REVISION HISTORY:
!
! name: originally ice_albedo
!
! authors:  Bruce P. Briegleb, NCAR
!           Elizabeth C. Hunke and William H. Lipscomb, LANL
! 2005, WHL: Moved absorbed_solar from ice_therm_vertical to this 
!            module and changed name from ice_albedo
! 2006, WHL: Added delta Eddington routines from Bruce Briegleb
!
! !INTERFACE:
!
      module ice_shortwave
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size
      use ice_constants
!
!EOP
!
      implicit none
      save

      character (len=char_len) ::
     &   shortwave     ! shortwave method, 'ccsm3' or 'dEdd'

      ! baseline albedos, set in namelist
      real (kind=dbl_kind) ::
     &   albicev   ! visible ice albedo for h > ahmax
     &,  albicei   ! near-ir ice albedo for h > ahmax
     &,  albsnowv  ! cold snow albedo, visible
     &,  albsnowi  ! cold snow albedo, near IR

      ! storage for approximate exponential for delta-Eddington
      integer (kind=int_kind), parameter :: 
     &   nmbexp = 1000000  ! number of exponential values in lookup table

      real (kind=dbl_kind), parameter :: 
     &   argmax = c10      ! maximum argument of exponential

      real (kind=dbl_kind) ::
     &   app_exp(0:nmbexp) ! table for exponentials
     &,  dx_exp            ! change in argument between table values
     &,  exp_min           ! minimum exponential value

      integer (kind=int_kind) ::
     &   indx_exp          ! index for exponential lookup table

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: shortwave_ccsm3 - driver for CCSM3 shortwave radiation
!
! !INTERFACE:
!
      subroutine shortwave_ccsm3 (nx_block, ny_block,
     &                            icells,
     &                            indxi,    indxj,
     &                            aicen,    vicen,
     &                            vsnon,    Tsfcn,
     &                            swvdr,    swvdf,
     &                            swidr,    swidf,
     &                            alvdrn,   alidrn,
     &                            alvdfn,   alidfn,
     &                            fswsfc,   fswint,
     &                            fswthru,  Iswabs)
!
! !DESCRIPTION:
!
! Driver for basic solar radiation from CCSM3.  Albedos and absorbed solar.
!
! !REVISION HISTORY:
!
! authors:  same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  icells              ! number of ice-covered grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block),
     &   intent(in) ::
     &   indxi        ! indices for ice-covered cells
     &,  indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   aicen        ! concentration of ice per category
     &,  vicen        ! volume of ice per category
     &,  vsnon        ! volume of ice per category
     &,  Tsfcn        ! surface temperature
     &,  swvdr        ! sw down, visible, direct  (W/m^2)
     &,  swvdf        ! sw down, visible, diffuse (W/m^2)
     &,  swidr        ! sw down, near IR, direct  (W/m^2)
     &,  swidf        ! sw down, near IR, diffuse (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(out) ::
     &   alvdrn       ! visible, direct, avg   (fraction)
     &,  alidrn       ! near-ir, direct, avg   (fraction)
     &,  alvdfn       ! visible, diffuse, avg  (fraction)
     &,  alidfn       ! near-ir, diffuse, avg  (fraction)
     &,  fswsfc       ! SW absorbed at ice/snow surface (W m-2)
     &,  fswint       ! SW absorbed in ice interior, below surface (W m-2)
     &,  fswthru      ! SW through ice to ocean (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &     intent(out) ::
     &   Iswabs       ! SW absorbed in particular layer (W m-2)
!
!EOP
!
      ! ice and snow albedo for each category

      real (kind=dbl_kind), dimension (nx_block,ny_block)::
     &   alvdrni    ! visible, direct, ice    (fraction)
     &,  alidrni    ! near-ir, direct, ice    (fraction)
     &,  alvdfni    ! visible, diffuse, ice   (fraction)
     &,  alidfni    ! near-ir, diffuse, ice   (fraction)
     &,  alvdrns    ! visible, direct, snow   (fraction)
     &,  alidrns    ! near-ir, direct, snow   (fraction)
     &,  alvdfns    ! visible, diffuse, snow  (fraction)
     &,  alidfns    ! near-ir, diffuse, snow  (fraction)

      !-----------------------------------------------------------------
      ! Compute albedos for ice and snow.
      !-----------------------------------------------------------------

      call compute_albedos (nx_block,   ny_block,
     &                      icells,
     &                      indxi,      indxj,
     &                      aicen,      vicen,
     &                      vsnon,      Tsfcn,
     &                      alvdrni,    alidrni,
     &                      alvdfni,    alidfni,
     &                      alvdrns,    alidrns,
     &                      alvdfns,    alidfns,
     &                      alvdrn,     alidrn,
     &                      alvdfn,     alidfn)

      !-----------------------------------------------------------------
      ! Compute solar radiation absorbed in ice and penetrating to ocean.
      !-----------------------------------------------------------------

      call absorbed_solar  (nx_block,   ny_block,
     &                      icells, 
     &                      indxi,      indxj,
     &                      aicen,
     &                      vicen,      vsnon,
     &                      swvdr,      swvdf,
     &                      swidr,      swidf,
     &                      alvdrni,    alvdfni,
     &                      alidrni,    alidfni,
     &                      alvdrns,    alvdfns,
     &                      alidrns,    alidfns,
     &                      fswsfc,     fswint,
     &                      fswthru,    Iswabs)

      end subroutine shortwave_ccsm3

!=======================================================================
!BOP
!
! !IROUTINE: compute_albedos - compute albedos for each thickness ategory
!
! !INTERFACE:
!
      subroutine compute_albedos (nx_block, ny_block,
     &                            icells,
     &                            indxi,    indxj,
     &                            aicen,    vicen,
     &                            vsnon,    Tsfcn,
     &                            alvdrni,  alidrni,
     &                            alvdfni,  alidfni,
     &                            alvdrns,  alidrns,
     &                            alvdfns,  alidfns,
     &                            alvdrn,   alidrn,
     &                            alvdfn,   alidfn)
!
! !DESCRIPTION:
!
! Compute albedos for each thickness category
!
! !REVISION HISTORY:
!
! authors:  same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  icells              ! number of ice-covered grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block),
     &   intent(in) ::
     &   indxi   ! compressed indices for ice-covered cells
     &,  indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   aicen       ! concentration of ice per category
     &,  vicen       ! volume of ice per category
     &,  vsnon       ! volume of ice per category
     &,  Tsfcn       ! surface temperature

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(out) ::
     &   alvdrni      ! visible, direct, ice   (fraction)
     &,  alidrni      ! near-ir, direct, ice   (fraction)
     &,  alvdfni      ! visible, diffuse, ice  (fraction)
     &,  alidfni      ! near-ir, diffuse, ice  (fraction)
     &,  alvdrns      ! visible, direct, snow  (fraction)
     &,  alidrns      ! near-ir, direct, snow  (fraction)
     &,  alvdfns      ! visible, diffuse, snow (fraction)
     &,  alidfns      ! near-ir, diffuse, snow (fraction)
     &,  alvdrn       ! visible, direct, avg   (fraction)
     &,  alidrn       ! near-ir, direct, avg   (fraction)
     &,  alvdfn       ! visible, diffuse, avg  (fraction)
     &,  alidfn       ! near-ir, diffuse, avg  (fraction)
!
!EOP
!
      real (kind=dbl_kind), parameter ::
     &   ahmax     = p5              ! thickness above which ice albedo 
                                     ! is constant (m)
     &,  dT_mlt    = c1              ! change in temp to give dalb_mlt 
                                     ! albedo change
     &,  dalb_mlt  = -0.075_dbl_kind ! albedo change per dT_mlt change
                                     ! in temp for ice
     &,  dalb_mltv = -p1             ! albedo vis change per dT_mlt change
                                     ! in temp for snow
     &,  dalb_mlti = -p15            ! albedo nir change per dT_mlt change
                                     ! in temp for snow

      integer (kind=int_kind) ::
     &   i, j, n

      real (kind=dbl_kind) ::
     &   hi      ! ice thickness  (m)
     &,  hs      ! snow thickness  (m)
     &,  albo    ! effective ocean albedo, function of ice thickness
     &,  fh      ! piecewise linear function of thickness
     &,  fT      ! piecewise linear function of surface temperature
     &,  dTs     ! difference of Tsfc and Timelt
     &,  fhtan   ! factor used in albedo dependence on ice thickness
     &,  asnow   ! fractional area of snow cover

      integer (kind=int_kind) ::
     &   ij      ! horizontal index, combines i and j loops


      fhtan = atan(ahmax*c4)

      alvdrni(:,:) = albocn
      alidrni(:,:) = albocn
      alvdfni(:,:) = albocn
      alidfni(:,:) = albocn

      alvdrns(:,:) = albocn
      alidrns(:,:) = albocn
      alvdfns(:,:) = albocn
      alidfns(:,:) = albocn

      alvdrn(:,:) = albocn
      alidrn(:,:) = albocn
      alvdfn(:,:) = albocn
      alidfn(:,:) = albocn
         
      !-----------------------------------------------------------------
      ! Compute albedo for each thickness category.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hi = vicen(i,j) / aicen(i,j)
         hs = vsnon(i,j) / aicen(i,j)            

         ! bare ice, thickness dependence
         fh = min(atan(hi*c4)/fhtan,c1)
         albo = albocn*(c1-fh)
         alvdfni(i,j) = albicev*fh + albo
         alidfni(i,j) = albicei*fh + albo

         ! bare ice, temperature dependence
         dTs = Timelt - Tsfcn(i,j)
         fT = min(dTs/dT_mlt-c1,c0)
         alvdfni(i,j) = alvdfni(i,j) - dalb_mlt*fT
         alidfni(i,j) = alidfni(i,j) - dalb_mlt*fT

         ! avoid negative albedos for thin, bare, melting ice
         alvdfni(i,j) = max (alvdfni(i,j), albocn)
         alidfni(i,j) = max (alidfni(i,j), albocn)

         if (hs > puny) then

            alvdfns(i,j) = albsnowv
            alidfns(i,j) = albsnowi

            ! snow on ice, temperature dependence
            alvdfns(i,j) = alvdfns(i,j) - dalb_mltv*fT
            alidfns(i,j) = alidfns(i,j) - dalb_mlti*fT

         endif                  ! hs > puny

         ! direct albedos (same as diffuse for now)
         alvdrni(i,j) = alvdfni(i,j)
         alidrni(i,j) = alidfni(i,j)
         alvdrns(i,j) = alvdfns(i,j)
         alidrns(i,j) = alidfns(i,j)


         ! fractional area of snow cover
         if (hs > puny) then
            asnow = hs / (hs + snowpatch)
         else
            asnow = c0
         endif

         ! combine ice and snow albedos (for coupler)
         alvdfn(i,j) = alvdfni(i,j)*(c1-asnow) +
     &                 alvdfns(i,j)*asnow
         alidfn(i,j) = alidfni(i,j)*(c1-asnow) +
     &                 alidfns(i,j)*asnow
         alvdrn(i,j) = alvdrni(i,j)*(c1-asnow) +
     &                 alvdrns(i,j)*asnow
         alidrn(i,j) = alidrni(i,j)*(c1-asnow) +
     &                 alidrns(i,j)*asnow

      enddo                     ! ij

      end subroutine compute_albedos

!=======================================================================
!BOP
!
! !ROUTINE: absorbed_solar - shortwave radiation absorbed by ice, ocean
!
! !DESCRIPTION:
!
! Compute solar radiation absorbed in ice and penetrating to ocean
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine absorbed_solar (nx_block, ny_block,
     &                           icells, 
     &                           indxi,    indxj,
     &                           aicen,
     &                           vicen,    vsnon,
     &                           swvdr,    swvdf,
     &                           swidr,    swidf,
     &                           alvdrni,  alvdfni,
     &                           alidrni,  alidfni,
     &                           alvdrns,  alvdfns,
     &                           alidrns,  alidfns,
     &                           fswsfc,   fswint,
     &                           fswthru,  Iswabs)
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

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aicen           ! fractional ice area
     &,  vicen           ! ice volume
     &,  vsnon           ! snow volume
     &,  swvdr           ! sw down, visible, direct  (W/m^2)
     &,  swvdf           ! sw down, visible, diffuse (W/m^2)
     &,  swidr           ! sw down, near IR, direct  (W/m^2)
     &,  swidf           ! sw down, near IR, diffuse (W/m^2)
     &,  alvdrni         ! visible, direct albedo,ice
     &,  alidrni         ! near-ir, direct albedo,ice
     &,  alvdfni         ! visible, diffuse albedo,ice
     &,  alidfni         ! near-ir, diffuse albedo,ice
     &,  alvdrns         ! visible, direct albedo, snow
     &,  alidrns         ! near-ir, direct albedo, snow
     &,  alvdfns         ! visible, diffuse albedo, snow
     &,  alidfns         ! near-ir, diffuse albedo, snow

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(out)::
     &   fswsfc          ! SW absorbed at ice/snow surface (W m-2)
     &,  fswint          ! SW absorbed in ice interior, below surface (W m-2)
     &,  fswthru         ! SW through ice to ocean (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(out) ::
     &   Iswabs          ! SW absorbed in particular layer (W m-2)
!
!EOP
!
      real (kind=dbl_kind), parameter ::
     &   i0vis = 0.70_dbl_kind  ! fraction of penetrating solar rad (visible)

      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! ice layer index

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   fswpen          ! SW penetrating beneath surface (W m-2)
     &,  trantop         ! transmitted frac of penetrating SW at layer top
     &,  tranbot         ! transmitted frac of penetrating SW at layer bot

      real (kind=dbl_kind) ::
     &   swabs           ! net SW down at surface (W m-2)
     &,  swabsv          ! swabs in vis (wvlngth < 700nm)  (W/m^2)
     &,  swabsi          ! swabs in nir (wvlngth > 700nm)  (W/m^2)
     &,  fswpenvdr       ! penetrating SW, vis direct
     &,  fswpenvdf       ! penetrating SW, vis diffuse
     &,  hi              ! ice thickness (m)
     &,  hs              ! snow thickness (m)
     &,  hilyr           ! ice layer thickness
     &,  asnow           ! fractional area of snow cover

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      fswsfc (:,:) = c0
      fswint (:,:) = c0
      fswthru(:,:) = c0
      Iswabs (:,:,:) = c0
      fswpen (:,:) = c0
      trantop(:,:) = c0
      tranbot(:,:) = c0

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hi  = vicen(i,j) / aicen(i,j)
         hs  = vsnon(i,j) / aicen(i,j)
         hilyr = hi / real(nilyr)

      !-----------------------------------------------------------------
      ! Fractional snow cover
      !-----------------------------------------------------------------
         if (hs > puny) then
            asnow = hs / (hs + snowpatch)
         else
            asnow = c0
         endif

      !-----------------------------------------------------------------
      ! Shortwave flux absorbed at surface, absorbed internally,
      !  and penetrating to mixed layer.
      ! This parameterization assumes that all IR is absorbed at the
      !  surface; only visible is absorbed in the ice interior or
      !  transmitted to the ocean.
      !-----------------------------------------------------------------

         swabsv  = swvdr(i,j) * ( (c1-alvdrni(i,j))*(c1-asnow)
     &                          + (c1-alvdrns(i,j))*asnow )
     &           + swvdf(i,j) * ( (c1-alvdfni(i,j))*(c1-asnow)
     &                          + (c1-alvdfns(i,j))*asnow )

         swabsi  = swidr(i,j) * ( (c1-alidrni(i,j))*(c1-asnow)
     &                          + (c1-alidrns(i,j))*asnow )
     &           + swidf(i,j) * ( (c1-alidfni(i,j))*(c1-asnow)
     &                          + (c1-alidfns(i,j))*asnow )

         swabs   = swabsv + swabsi

         fswpenvdr = swvdr(i,j) * (c1-alvdrni(i,j)) * (c1-asnow) * i0vis
         fswpenvdf = swvdf(i,j) * (c1-alvdfni(i,j)) * (c1-asnow) * i0vis

          ! no penetrating radiation in near IR
!         fswpenidr = swidr(i,j) * (c1-alidrni(i,j)) * (c1-asnow) * i0nir
!         fswpenidf = swidf(i,j) * (c1-alidfni(i,j)) * (c1-asnow) * i0nir  

         fswpen(i,j) = fswpenvdr + fswpenvdf
                      
         fswsfc(i,j) = swabs - fswpen(i,j)

         trantop(i,j) = c1  ! transmittance at top of ice

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! penetrating SW absorbed in each ice layer
      !-----------------------------------------------------------------

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            tranbot(i,j) = exp (-kappav * hilyr * real(k,kind=dbl_kind))
            Iswabs(i,j,k) = fswpen(i,j) * (trantop(i,j)-tranbot(i,j))

            ! bottom of layer k = top of layer k+1
            trantop(i,j) = tranbot(i,j)

         enddo                  ! ij
      enddo                     ! nilyr

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! SW penetrating thru ice into ocean
         fswthru(i,j) = fswpen(i,j) * tranbot(i,j)

         ! SW absorbed in ice interior
         fswint(i,j)  = fswpen(i,j) - fswthru(i,j)

      enddo                     ! ij

      end subroutine absorbed_solar

! End ccsm3 shortwave method
!=======================================================================
! Begin delta-Eddington shortwave method
!
!BOP
!
! !IROUTINE: init_dEdd - initialize delta-Eddington parameters
!
! !INTERFACE:
!
      subroutine init_dEdd
!
! !DESCRIPTION:
!
! Compute initial data for delta-Eddington method, specifically, 
! the approximate exponential look-up table.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS: none
!
!EOP
!     
      real (kind=dbl_kind) ::
     &   x  ! abscissa for exp evaluation

      ! evaluate approximate exponential table
      dx_exp = argmax / real(nmbexp)
      x          = c0
      app_exp(0) = c1
      do indx_exp=1,nmbexp
        x = x + dx_exp
        app_exp(indx_exp) = exp(-x)
      enddo 
      exp_min = exp(-argmax)
 
      end subroutine init_dEdd
 
!=======================================================================
!BOP
!
! !IROUTINE: shortwave_dEdd - compute shortwave using delta-Eddington method
!
! !INTERFACE:
!
      subroutine shortwave_dEdd  (nx_block, ny_block,
     &                            icells,
     &                            indxi,    indxj,
     &                            tlat,     coszen,
     &                            aicen,    vicen,
     &                            vsnon,    Tsfcn,
     &                            swvdr,    swvdf,
     &                            swidr,    swidf,
     &                            alvdrn,   alidrn,
     &                            alvdfn,   alidfn,
     &                            fswsfc,   fswint,
     &                            fswthru,  Iswabs)
!
!
! !DESCRIPTION:
!
!   Compute snow/ice shortwave, including albedos, absorbed and transmitted flux
!   For each category, compute shortwave albedos and fluxes for two surface types:
!   snow over ice (including bare ice) and ponded ice. For large snow depths, no 
!   ponds are assumed present, so a single snow-over-ice calculation is done. For
!   small snow (or no snow) depths, ponds are present, for which a weighted average 
!   over a prescribed pond fraction is done between bare ice and ponded ice.
!   Albedos and fluxes are stored for later use by thermodynamic routines. 
!   Invokes at most two calls to compute_dEdd, which sets inherent optical
!   properties appropriate for the surface type, and then calls solution_dEdd
!   to evaluate the delta-Eddington solution.
!
!   Two important prescribed quantities are the snow age and the pond fraction/depth.
!   These are empirically determined at the present from snow/ice state (snow depth
!   and surface temperature in the former, and calendar date in the latter).
!   The snow age is a linear weight between -1 for fresh, 0 for small, and +1 for 
!   large grain snow inherent optical properties. The pond fraction is simply
!   prescribed based on a calendar date, and depth is the pond fraction unless the
!   depth is greater than half the ice thickness, for which it takes half ice 
!   thickness value. Ideally, these quantities should be made prognostic in the
!   future.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR 
!
! !USES:
!
      use ice_calendar
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  icells              ! number of ice-covered grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block),
     &   intent(in) ::
     &   indxi   ! compressed indices for ice-covered cells
     &,  indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   tlat         ! latitude (radians)
     &,  coszen       ! cosine of solar zenith angle 
     &,  aicen        ! concentration of ice per category
     &,  vicen        ! volume of ice per category
     &,  vsnon        ! volume of ice per category
     &,  Tsfcn        ! surface temperature
     &,  swvdr        ! sw down, visible, direct  (W/m^2)
     &,  swvdf        ! sw down, visible, diffuse (W/m^2)
     &,  swidr        ! sw down, near IR, direct  (W/m^2)
     &,  swidf        ! sw down, near IR, diffuse (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(out) ::
     &   alvdrn       ! visible, direct, avg   (fraction)
     &,  alidrn       ! near-ir, direct, avg   (fraction)
     &,  alvdfn       ! visible, diffuse, avg  (fraction)
     &,  alidfn       ! near-ir, diffuse, avg  (fraction)
     &,  fswsfc       ! SW absorbed at ice/snow surface (W m-2)
     &,  fswint       ! SW absorbed in ice interior, below surface (W m-2)
     &,  fswthru      ! SW through ice to ocean (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &     intent(out) ::
     &   Iswabs         ! SW absorbed in particular layer (W m-2)
!      
!EOP
! 
! !LOCAL PARAMETERS:
!
      real (kind=dbl_kind),dimension (nx_block,ny_block) ::
     &   fnidr         ! fraction of direct to total down surface flux in nir

      ! albedos and fluxes for a given surface type
      real (kind=dbl_kind),dimension (nx_block,ny_block) ::
     &   avdr          ! vis direct albedo
     &,  avdf          ! vis diffuse albedo
     &,  aidr          ! nir direct albedo
     &,  aidf          ! nir direct albedo
     &,  fsfc          ! SW absorbed at ice/snow surface (W m-2)
     &,  fint          ! SW absorbed in ice interior, below surface (W m-2)
     &,  fthru         ! SW through ice to ocean (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr) ::
     &   Iabs          ! SW absorbed in particular layer (W m-2)
 
      real (kind=dbl_kind), dimension(nx_block,ny_block) ::
     &   hsn           ! snow thickness (m)
     &,  sage          ! snow age (-1 to +1)
     &,  hin           ! ice thickness (m)
     &,  hpn           ! pond depth (m)
     &,  fi            ! snow/bare ice fractional coverage (0 to 1)
     &,  fs            ! snow fractional coverage (0 to 1)
     &,  fp            ! pond fractional coverage (0 to 1)

      integer (kind=int_kind), dimension(nx_block,ny_block) ::
     &   srftyp        ! surface type over ice: (0=air, 1=snow, 2=pond)
 
      real (kind=dbl_kind), parameter ::
     &   hspond = 0.025_dbl_kind  ! snow depth (m) below which ponded ice exists
     &,  cp40   = 0.400_dbl_kind  ! snow age numerical parameter
     &,  cp50   = 0.500_dbl_kind  ! snow age numerical parameter
     &,  fpmin  = 0.100_dbl_kind  ! ponded ice minimum fractional coverage
     &,  fpmax  = 0.450_dbl_kind  ! ponded ice maximum fractional coverage
     &,  c165   = 165.0_dbl_kind  ! day offset for ponded ice fraction
     &,  c167   = 167.0_dbl_kind  ! normalization day for ponded ice fraction
 
      real (kind=dbl_kind) ::      
     &   dTs     ! temperature difference for snow age
     &,  fT      ! temperature weight for snow age
     &,  ydaysh  ! southern hemisphere day of year for pond fraction
 
      integer (kind=int_kind) ::
     &   i       ! longitude index
     &,  j       ! latitude index
     &,  ij      ! horizontal index, combines i and j loops
     &,  k       ! level index
     &,  icells2 ! number of cells in delta-Eddington calculation
 
      integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &   indxi2  ! compressed indices for delta-Eddington cells
     &,  indxj2  

!-----------------------------------------------------------------------
 
      ! compute fraction of nir down direct to total:
      do j = 1, ny_block
      do i = 1, nx_block
         fnidr(i,j) = c0
         if( swidr(i,j) + swidf(i,j) > puny ) then
            fnidr(i,j) = swidr(i,j)/(swidr(i,j)+swidf(i,j))
         endif
      enddo                ! i
      enddo                ! j
 
      ! zero storage albedos and fluxes for accumulation over 
      ! surface types
      fi(:,:) = c0
      fp(:,:) = c0
      hin(:,:) = c0
      hsn(:,:) = c0
      hpn(:,:) = c0
      alvdfn(:,:)   = c0
      alidfn(:,:)   = c0
      alvdrn(:,:)   = c0
      alidrn(:,:)   = c0
      fswsfc(:,:)   = c0
      fswint(:,:)   = c0
      fswthru(:,:)  = c0
      Iswabs(:,:,:) = c0
 
      ! compute shortwave radiation accounting for snow/ice (both snow over 
      ! ice and bare ice) and ponded ice (if any):
 
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      ! southern hemisphere pond day (6 month phase difference with northern)
      ydaysh  = yday + c365/c2
      if(ydaysh > c365 ) ydaysh = ydaysh - c365
      if(ydaysh < c365 ) ydaysh = ydaysh + c365

      ! find snow/ice properties for ice points
      icells2 = 0
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (aicen(i,j) > puny .and. coszen(i,j) > puny) then
            icells2 = icells2 + 1
            indxi2(icells2) = i
            indxj2(icells2) = j 

            hsn(i,j)  = vsnon(i,j) / aicen(i,j)
            hin(i,j)  = vicen(i,j) / aicen(i,j)
            hpn(i,j)  = c0
            fi(i,j)   = c1
            fp(i,j)   = c0
            if( hsn(i,j) > puny ) then
               ! snow over ice
               srftyp(i,j) = 1
            else
               ! bare ice
               srftyp(i,j) = 0
            endif
            ! snow thin; some melt ponds
            if( hsn(i,j) < hspond) then
               fs(i,j) = hsn(i,j)/hspond
               ! northern hemisphere 
               if(tlat(i,j) > c0) then
                fp(i,j)= min(max(fpmin+(yday-c165)/c167,c0),fpmax)
     &                   *(c1-fs(i,j))
               ! southern hemisphere 
               else 
                fp(i,j)= min(max(fpmin+(ydaysh-c165/c167),c0),fpmax)
     &                   *(c1-fs(i,j))
               endif
               fi(i,j) = c1-fp(i,j)
            endif
            ! snow age
            dTs = Timelt - Tsfcn(i,j)
            fT  = min(dTs/c3,c1)
            sage(i,j) = max(c1-hsn(i,j)/cp40,-cp50) - cp50*fT
         endif
      enddo                     ! ij

      ! calculate snow/ice
      call compute_dEdd
     &      (nx_block,ny_block,
     &       icells2, indxi2, indxj2, fnidr, coszen,
     &       swvdr,   swvdf,  swidr,  swidf,
     &       srftyp,  hsn,    sage,   hin,   hpn,
     &       avdr,    avdf,   aidr,   aidf,
     &       fsfc,    fint,   fthru,  Iabs)

      ! accumulate fluxes over snow/ice
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells2
         i = indxi2(ij)
         j = indxj2(ij)
         alvdfn(i,j)  = alvdfn(i,j)  + avdf(i,j) *fi(i,j)
         alidfn(i,j)  = alidfn(i,j)  + aidf(i,j) *fi(i,j)
         alvdrn(i,j)  = alvdrn(i,j)  + avdr(i,j) *fi(i,j)
         alidrn(i,j)  = alidrn(i,j)  + aidr(i,j) *fi(i,j)
         fswsfc(i,j)  = fswsfc(i,j)  + fsfc(i,j) *fi(i,j)
         fswint(i,j)  = fswint(i,j)  + fint(i,j) *fi(i,j)
         fswthru(i,j) = fswthru(i,j) + fthru(i,j)*fi(i,j)
      enddo                     ! ij

      do k = 1, nilyr
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells2
            i = indxi2(ij)
            j = indxj2(ij)
            Iswabs(i,j,k) = Iswabs(i,j,k) + Iabs(i,j,k)*fi(i,j)
         enddo                  ! ij
      enddo                     ! k
      
      ! find ponded ice properties for ice points
      icells2 = 0
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         hin(i,j) = c0
         hsn(i,j) = c0
         if (aicen(i,j) > puny .and. coszen(i,j) > puny) then
            hsn(i,j)  = vsnon(i,j) / aicen(i,j)
            hin(i,j)  = vicen(i,j) / aicen(i,j)
            ! snow thin or no snow: melt ponds
            if( hsn(i,j) < hspond ) then
               icells2 = icells2 + 1
               indxi2(icells2) = i
               indxj2(icells2) = j
               ! ponded ice
               srftyp(i,j)   = 2
               fs(i,j) = hsn(i,j)/hspond
               ! northern hemisphere
               if(tlat(i,j) > c0) then
                fp(i,j)= min(max(fpmin+(yday-c165)/c167,c0),fpmax)
     &                  *(c1-fs(i,j))
               ! southern hemisphere
               else 
                fp(i,j)= min(max(fpmin+(ydaysh-c165/c167),c0),fpmax)
     &                   *(c1-fs(i,j))
               endif              
               hpn(i,j)= fp(i,j)
               ! limit pond depth to less than 1/2 ice thickness
               hpn(i,j) = min(hpn(i,j), hin(i,j)/c2)
            endif               ! hsn < hspond
         endif                  ! aicen > puny, coszen > puny
      enddo                     ! ij

      ! calculate ponded ice
      call compute_dEdd
     &      (nx_block,ny_block,
     &       icells2, indxi2, indxj2, fnidr, coszen,
     &       swvdr,   swvdf,  swidr,  swidf,
     &       srftyp,  hsn,    sage,   hin,   hpn,
     &       avdr,    avdf,   aidr,   aidf,
     &       fsfc,    fint,   fthru,  Iabs)

      ! accumulate fluxes over ponded ice
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells2
         i = indxi2(ij)
         j = indxj2(ij)
         alvdfn(i,j) = alvdfn(i,j)  + avdf(i,j) *fp(i,j)
         alidfn(i,j) = alidfn(i,j)  + aidf(i,j) *fp(i,j)
         alvdrn(i,j) = alvdrn(i,j)  + avdr(i,j) *fp(i,j)
         alidrn(i,j) = alidrn(i,j)  + aidr(i,j) *fp(i,j)
         fswsfc(i,j) = fswsfc(i,j)  + fsfc(i,j) *fp(i,j)
         fswint(i,j) = fswint(i,j)  + fint(i,j) *fp(i,j)
         fswthru(i,j)= fswthru(i,j) + fthru(i,j)*fp(i,j) 
      enddo                     ! ij

      do k = 1, nilyr
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells2
            i = indxi2(ij)
            j = indxj2(ij)
            Iswabs(i,j,k) = Iswabs(i,j,k) + Iabs(i,j,k)*fp(i,j)
         enddo                  ! ij
      enddo                     ! k

      end subroutine shortwave_dEdd

!=======================================================================
!BOP
!
! !IROUTINE: compute_dEdd - evaluate iops for delta-Edd and compute solution
!
! !INTERFACE:
!
      subroutine compute_dEdd
     &      (nx_block,ny_block,
     &       icells,  indxi,  indxj,  fnidr, coszen,
     &       swvdr,   swvdf,  swidr,  swidf,
     &       srftyp,  hsn,    sage,   hin,   hpn,
     &       avdr ,   avdf,   aidr,   aidf,
     &       fsfc,    fint,   fthru,  Iabs)
!
! !DESCRIPTION:
!
! Evaluate snow/ice/pond inherent optical properties (iops) for evaluation 
! of delta-Eddington solution, and then calculate the solution.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR 
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  icells  ! number of ice/ocean grid cells for surface type
 
      integer (kind=int_kind), dimension(nx_block*ny_block), 
     &   intent(in)::
     &   indxi   ! compressed indices for ice/ocean cells for surface type
     &,  indxj
 
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   fnidr   ! fraction of direct to total down flux in nir
     &,  coszen  ! cosine solar zenith angle
     &,  swvdr   ! sw down, visible, direct  (W/m^2)
     &,  swvdf   ! sw down, visible, diffuse (W/m^2)
     &,  swidr   ! sw down, near IR, direct  (W/m^2)
     &,  swidf   ! sw down, near IR, diffuse (W/m^2)
 
      integer (kind=int_kind), dimension(nx_block,ny_block),
     &   intent(in)::
     &   srftyp  ! surface type over ice: (0=air, 1=snow, 2=pond)
 
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::
     &   hsn     ! snow thickness (m)
     &,  sage    ! snow age (-1 to +1)
     &,  hin     ! ice thickness (m)
     &,  hpn     ! ponded depth (m)
 
      ! albedos 
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(out) ::
     &   avdr    ! visible, direct   (fraction)
     &,  avdf    ! visible, diffuse  (fraction)
     &,  aidr    ! near-ir, direct   (fraction)
     &,  aidf    ! near-ir, diffuse  (fraction)
 
      ! fluxes absorbed in snow/sea ice and through to ocean
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(out) ::
     &   fsfc        ! SW absorbed at ice/snow surface (W m-2)
     &,  fint        ! SW absorbed in ice interior, below surface (W m-2)
     &,  fthru       ! SW through ice to ocean            (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr),
     &   intent(out) ::
     &   Iabs        ! SW absorbed in particular layer (W m-2)
! 
!EOP
!-----------------------------------------------------------------------
! Compute shortwave radiation (surface albedos and fluxes)
! August 2005  Bruce P. Briegleb
!
! Uses delta-Eddington adding method 
! 
! Snow/ice is divided into horizontally and vertically homogeneous layers.
! For each layer, a solution for reflectivity/transmissivity (R/T) for direct 
! and diffuse radiation is evaluated using the delta-Eddington method. This
! method assumes the scattering phase function can be expressed as a delta 
! function forward peak with the remainder a cosine function, in accordance 
! with the Eddington approximation. The resulting transfer equation can be 
! solved in closed form. Direct and diffuse incident radiation are distinguished, 
! with the diffuse assumed isotropic in each hemisphere. Once R/T are evaluated
! for each layer for direct and diffuse incident radiation, the layers are added 
! assuming that interlayer diffuse radiation is isotropic in each hemisphere. 
! Refraction at a specified layer interface is handled assuming an extra
! non-absorbing layer whose R/T are given by the fresnel formulas. Once the 
! layers are added, the up/down direct and diffuse fluxes at each layer 
! interface can be evaluated, and hence reflectivity, layer absorption 
! and transmission to the underlying ocean. 
!
! Divides solar spectrum into 3 intervals: 0.2-0.7, 0.7-1.19, and 1.19-5.0 
! micro-meters. The latter two are added (using a partition of incident 
! shortwave dependent on the 0.7-5.0 micro-meter direct to total incident
! flux) to give the final output of 0.2-0.7 visible and 0.7-5.0 infrared 
! albedos and fluxes.
!
! This routine specifies vertical layer optical properties based on input 
! surface type, while the called-routine solution_dEdd computes layer by layer 
! delta-Eddington reflectivity and transmissivity, and adds layers. Finally, 
! surface albedos and internal fluxes/flux divergences are evaluated in this 
! routine and saved for later use.
!
!
!      Description of the level and layer index conventions.
!
!  CCSM3 Sea Ice Model Layers       dEdd shortwave calculation
!                                 
!                                 (radiation layers and interfaces)
!                             layer index             interface index
!
!                                     ----------------------  0
!    /////////////////////            //////////////////////
!    Snow layer if present         0  air, snow or melt pond
!    /////////////////////            //////////////////////
!    ---------------------            ----------------------  1
!                                  1  \\\ surface layer \\\\
!          layer 1                    ----------------------  2
!                                  2
!    ---------------------            ----------------------  3
!
!          layer 2                 3
!
!    ---------------------            ----------------------  4
!
!          layer 3                 4
!
!    ---------------------            ----------------------  5
!
!          layer 4                 5
!
!    ---------------------            ----------------------  6
!
! The CCSM sea ice model has four vertical layers for the thermodynamic 
! heat tranfer; for the dEdd calculation, the same four vertical layers 
! are used plus a surface layer. This surface layer is implicit in the 
! CCSM3 snow/ice absorption calculation, but here is made explicit. The 
! surface layer is h1_4 thick (a specified constant) unless the first layer 
! is smaller than twice this; then it uses hi/8 for the thickness (hi=total 
! ice thickness). The dEdd layer thickness for layer 2 is hi/4 - h1_4, and
! for layers 3,4,5 it is hi/4.
!
! The radiation absorbed in the surface layer is used for the surface 
! melting calculation for bare ice; that which penetrates into radiation 
! layer 2 is used for sea ice model layer 1 heating. For the case of snow 
! and melt ponds, the radiation absorbed in the snow layer or the melt pond 
! is also used for surface melting. The net transmission in visible (vs) and 
! near-ir (ni) bands at the radiation interface 2 is equivalent to the 
! classical I_0 fraction of shortwave radiation (absorbed in the column) 
! that penetrates the surface layer.
!
!-----------------------------------------------------------------------
!
! !LOCAL PARAMETERS
!
      integer (kind=int_kind) ::
     &   i       ! longitude index
     &,  j       ! latitude index
     &,  k       ! level index
     &,  ij      ! horizontal index, combines i and j loops
     &,  ns      ! spectral index
 
      integer (kind=int_kind), parameter :: 
     &   klev    = nilyr + 1   ! number of radiation layers
     &,  klevp   = klev  + 1   ! number of radiation layers + 1
 
      integer (kind=int_kind), parameter :: 
     &   nspint  = 3     ! number of solar spectral intervals
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,nspint) ::
     &   wghtns          ! spectral weights
 
      real (kind=dbl_kind), parameter :: 
     &   cp67    = 0.67_dbl_kind   ! nir band weight parameter
     &,  cp33    = 0.33_dbl_kind   ! nir band weight parameter
     &,  cp78    = 0.78_dbl_kind   ! nir band weight parameter
     &,  cp22    = 0.22_dbl_kind   ! nir band weight parameter
     &,  cp01    = 0.01_dbl_kind   ! for ocean albedo
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,0:klev) ::
     &   tau     ! layer extinction optical depth
     &,  w0      ! layer single scattering albedo
     &,  g       ! layer asymmetry parameter
 
      ! following arrays are defined at model interfaces; 0 is the top of the
      ! layer above the sea ice; klevp is the sea ice/ocean interface.
      real (kind=dbl_kind), dimension (nx_block,ny_block,0:klevp) ::
     &   trndir  ! solar beam down transm from top
     &,  trntdr  ! total transmission to direct beam for layers above
     &,  trndif  ! diffuse transmission to diffuse beam for layers above
     &,  rupdir  ! ref to dir rad for layers below
     &,  rupdif  ! ref to dif rad for layers below
     &,  rdndif  ! ref to dif rad for layers above
 
      real (kind=dbl_kind) ::
     &   refk    ! interface multiple scattering k
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,0:klevp) ::
     &   fdirup  ! up   flux to direct beam at model interface 
     &,  fdirdn  ! down flux to direct beam at model interface
     &,  fdifup  ! up   flux to diffuse beam at model interface 
     &,  fdifdn  ! down flux to diffuse beam at model interface
 
      ! inherent optical property (iop) arrays and data
      real (kind=dbl_kind), dimension (nspint) ::
     &   ks_f    ! fresh snow grain extinction coefficient m^{-1}
     &,  ks_s    ! small snow grain extinction coefficient m^{-1}
     &,  ks_l    ! large snow grain extinction coefficient m^{-1}    
     &,  ws_f    ! fresh snow grain single scattering albedo    
     &,  ws_s    ! small snow grain single scattering albedo    
     &,  ws_l    ! large snow grain single scattering albedo    
     &,  gs_f    ! fresh snow grain asymmetry parameter        
     &,  gs_s    ! small snow grain asymmetry parameter
     &,  gs_l    ! large snow grain asymmetry parameter
 
      real (kind=dbl_kind), dimension (nspint,klev) ::
     &   ki_f    ! first year ice extinction coefficient m^{-1}
     &,  wi_f    ! first year ice single scattering albedo
     &,  gi_f    ! first year ice asymmetry parameter
 
      real (kind=dbl_kind), dimension (nspint) ::
     &   ki_p    ! first year ice under pond extinction coefficient m^{-1}
     &,  wi_p    ! first year ice under pond single scattering albedo
     &,  gi_p    ! first year ice under pond asymmetry parameter
     &,  kw      ! water extinction coefficient m^{-1}
     &,  ww      ! water single scattering albedo
     &,  gw      ! water asymmetry parameter
 
      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &  albodr   ! spectral ocean albedo to direct rad
     &, albodf   ! spectral ocean albedo to diffuse rad
     
      real (kind=dbl_kind) ::
     &  dz       ! sea ice layer thickness
     &, dzs      ! sea ice surface layer thickness
     &, dz1      ! sea ice layer thickness below surface layer
 
      real (kind=dbl_kind), parameter ::
     &  h1_4 = 0.50_dbl_kind   ! surface layer thickness (m)
 
      ! Extinction data for snow is taken from: Wiscombe and Warren, 1980
      ! and Warren and Wiscombe, 1980. Extinction for fresh snow is that 50
      ! micron radius ice spheres, aged 100 micron radius ice spheres, and
      ! that of melting snow 1000 micron radius ice spheres. The single 
      ! scattering albedos are those for ice with small dust and graphitic 
      ! carbon soot inclusions (latter order 1.5ppmw). The asymmetry 
      ! parameters are rounded to two significant digits.

      data ks_f /  10900._dbl_kind, 11000._dbl_kind, 11100._dbl_kind /
      data ks_s /   5450._dbl_kind,  5500._dbl_kind,  5550._dbl_kind /
      data ks_l /    545._dbl_kind,   550._dbl_kind,   555._dbl_kind /
 
      data ws_f /  0.9999_dbl_kind, 0.9998_dbl_kind,  0.997_dbl_kind /
      data ws_s /  0.9995_dbl_kind, 0.9992_dbl_kind,  0.995_dbl_kind /
      data ws_l /  0.9975_dbl_kind, 0.9955_dbl_kind,  0.900_dbl_kind /
 
      data gs_f /    0.89_dbl_kind,   0.89_dbl_kind,   0.90_dbl_kind /
      data gs_s /    0.89_dbl_kind,   0.89_dbl_kind,   0.90_dbl_kind /
      data gs_l /    0.89_dbl_kind,   0.90_dbl_kind,   0.93_dbl_kind /
 
      ! Following data from Bonnie Light, February 2005
 
      data ki_f /   1000._dbl_kind,  1005._dbl_kind,  2224._dbl_kind,
     &                75._dbl_kind,    82._dbl_kind,  1969._dbl_kind,
     &                60._dbl_kind,    67._dbl_kind,  1954._dbl_kind,
     &                60._dbl_kind,    67._dbl_kind,  1954._dbl_kind,
     &                60._dbl_kind,    67._dbl_kind,  1954._dbl_kind /
      data wi_f /   .9999_dbl_kind,  .9954_dbl_kind,  .4495_dbl_kind,
     &              .9972_dbl_kind,  .9125_dbl_kind,  .0381_dbl_kind,
     &              .9965_dbl_kind,  .8930_dbl_kind,  .0307_dbl_kind,
     &              .9965_dbl_kind,  .8930_dbl_kind,  .0307_dbl_kind,
     &              .9965_dbl_kind,  .8930_dbl_kind,  .0307_dbl_kind /
      data gi_f /     .94_dbl_kind,    .94_dbl_kind,    .94_dbl_kind,
     &                .94_dbl_kind,    .94_dbl_kind,    .94_dbl_kind,
     &                .94_dbl_kind,    .94_dbl_kind,    .94_dbl_kind,
     &                .94_dbl_kind,    .94_dbl_kind,    .94_dbl_kind,
     &                .94_dbl_kind,    .94_dbl_kind,    .94_dbl_kind /
 
      data ki_p /     15._dbl_kind,    22._dbl_kind,  1900._dbl_kind /
      data wi_p /    .986_dbl_kind,   .676_dbl_kind,   .008_dbl_kind /
      data gi_p /    .940_dbl_kind,   .940_dbl_kind,   .940_dbl_kind /
 
      data kw   /    0.16_dbl_kind,   11.3_dbl_kind,  1980._dbl_kind /
      data ww   /    0.01_dbl_kind,   0.01_dbl_kind,   0.01_dbl_kind /
      data gw   /    0.01_dbl_kind,   0.01_dbl_kind,   0.01_dbl_kind /
 
!-----------------------------------------------------------------------
 
      ! initialize albedos and fluxes to 0
      avdr(:,:)   = c0
      avdf(:,:)   = c0
      aidr(:,:)   = c0
      aidf(:,:)   = c0
      fsfc(:,:)   = c0
      fint(:,:)   = c0
      fthru(:,:)  = c0
      Iabs(:,:,:) = c0
 
      ! spectral weights
      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)
        wghtns(i,j,1) = c1
        wghtns(i,j,2) = cp67 + (cp78-cp67)*(c1-fnidr(i,j))
        wghtns(i,j,3) = cp33 + (cp22-cp33)*(c1-fnidr(i,j))
      enddo
 
      ! Begin spectral loop
      do ns = 1, nspint
 
        ! Set optical properties of air/snow/melt_pond/ice
        do ij = 1, icells
          i = indxi(ij)
          j = indxj(ij)
          ! air
          if( srftyp(i,j) == 0 ) then
            tau(i,j,0) = c0
            w0(i,j,0)  = c0
            g(i,j,0)   = c0
          ! snow
          else if( srftyp(i,j) == 1 ) then
            if( sage(i,j) < c0 ) then
              tau(i,j,0) = (abs(sage(i,j))*ks_f(ns)+
     &                     (c1-abs(sage(i,j)))*ks_s(ns))*hsn(i,j)
              w0(i,j,0)  = (abs(sage(i,j))*ws_f(ns)+
     &                     (c1-abs(sage(i,j)))*ws_s(ns))
              g(i,j,0)   = (abs(sage(i,j))*gs_f(ns)+
     &                     (c1-abs(sage(i,j)))*gs_s(ns))
            else if( sage(i,j) >= c0 ) then
              tau(i,j,0) = ((c1-sage(i,j))*ks_s(ns)+
     &                     sage(i,j)*ks_l(ns))*hsn(i,j)
              w0(i,j,0)  = ((c1-sage(i,j))*ws_s(ns)+
     &                     sage(i,j)*ws_l(ns))
              g(i,j,0)   = ((c1-sage(i,j))*gs_s(ns)+
     &                     sage(i,j)*gs_l(ns))
            endif
          ! pond
          else if( srftyp(i,j) == 2 ) then
            tau(i,j,0) = kw(ns)*hpn(i,j)
            w0(i,j,0)  = ww(ns)
            g(i,j,0)   = gw(ns)
          endif                 ! srftyp
        enddo                   ! ij

        ! sea ice layers 
        do ij = 1, icells
          i = indxi(ij)
          j = indxj(ij)
          ! bare or snow-covered sea ice layers
          if( srftyp(i,j) <= 1 ) then
            dzs    = h1_4
            dz     = hin(i,j)/real(nilyr)
            if( dzs > dz/c2 ) dzs = dz/c2
            dz1    = dz - dzs
            do k=1,klev
              if( k==1 ) then
                tau(i,j,k) = ki_f(ns,k)*dzs
              else if( k==2 ) then
                tau(i,j,k) = ki_f(ns,k)*dz1
              else
                tau(i,j,k) = ki_f(ns,k)*dz
              endif 
              w0(i,j,k)  = wi_f(ns,k)
              g(i,j,k)   = gi_f(ns,k)
            enddo  
          ! sea ice layers under ponds
          else if( srftyp(i,j) == 2 ) then
            dz     = hin(i,j)/real(nilyr)
            dzs = dz/c2
            do k=1,klev
              if( k <= 2 ) then
                tau(i,j,k) = ki_p(ns)*dzs
              else
                tau(i,j,k) = ki_p(ns)*dz
              endif 
              w0(i,j,k)  = wi_p(ns)
              g(i,j,k)   = gi_p(ns)
            enddo       ! k
          endif         ! srftyp
        enddo           ! ij
 
        ! Set reflectivities for ocean underlying sea ice
        if(ns == 1) then
          do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            albodr(i,j) = cp01
            albodf(i,j) = cp01
          enddo
        else if(ns >= 2) then
          do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            albodr(i,j) = c0
            albodf(i,j) = c0
          enddo
        endif
 
        ! Layer input properties now completely specified; compute the
        ! delta-Eddington solution reflectivities and transmissivities
        ! for each layer, starting from the top and working downwards;
        ! next, add the layers going downwards accounting for multiple
        ! scattering between layers, and finally start from the surface 
        ! and add successive layers upwards to the top:
 
        call solution_dEdd
     &      (nx_block, ny_block,
     &       icells, indxi,  indxj,  coszen, srftyp,
     &       tau,    w0,     g,      albodr, albodf,
     &       trndir, trntdr, trndif, rupdir, rupdif,
     &       rdndif)
 
        ! Compute up and down fluxes for each interface, using the 
        ! added layer properties at each interface:
        !              layers       interface
        !
        !       ---------------------  k
        !                 k
        !       --------------------- 
        do k=0,klevp
        do ij = 1, icells
          i = indxi(ij)
          j = indxj(ij)
          ! interface scattering
          refk          = c1/(c1 - rdndif(i,j,k)*rupdif(i,j,k))
          ! dir tran ref from below times interface scattering, plus diff
          ! tran and ref from below times interface scattering
          fdirup(i,j,k) = (trndir(i,j,k)*rupdir(i,j,k) +
     &                    (trntdr(i,j,k)-trndir(i,j,k)) *
     &                     rupdif(i,j,k))*refk
          ! dir tran plus total diff trans times interface scattering plus
          ! dir tran with up dir ref and down dif ref times interface scattering 
          fdirdn(i,j,k) = trndir(i,j,k) + (trntdr(i,j,k) - trndir(i,j,k)
     &                  + trndir(i,j,k)*rupdir(i,j,k)
     &                   *rdndif(i,j,k))*refk
          ! diffuse tran ref from below times interface scattering
          fdifup(i,j,k) = trndif(i,j,k)*rupdif(i,j,k)*refk
          ! diffuse tran times interface scattering
          fdifdn(i,j,k) = trndif(i,j,k)*refk
        enddo                   ! ij
        enddo                   ! k
 
        ! Calculate final surface albedos and fluxes-
        ! Note assumptions about surface absorption: regardless of what may overlie
        ! sea ice (air, snow or pond), absorption in the surface scattering layer is
        ! included along with snow or pond (if present) absorption. Since the 
        ! interface index of the bottom of the surface layer is 2, that index
        ! appears explicitly.

        if( ns == 1) then      ! visible

          do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            avdr(i,j)  = rupdir(i,j,0)
            avdf(i,j)  = rupdif(i,j,0)
 
            fsfc(i,j)  = fsfc(i,j) + 
     &        ((fdirdn(i,j,0)-fdirup(i,j,0))*swvdr(i,j) +
     &         (fdifdn(i,j,0)-fdifup(i,j,0))*swvdf(i,j)) -
     &        ((fdirdn(i,j,2)-fdirup(i,j,2))*swvdr(i,j) +
     &         (fdifdn(i,j,2)-fdifup(i,j,2))*swvdf(i,j))
 
            fint(i,j)  = fint(i,j) + 
     &        ((fdirdn(i,j,2)-fdirup(i,j,2))*swvdr(i,j) +
     &         (fdifdn(i,j,2)-fdifup(i,j,2))*swvdf(i,j)) -
     &        ((fdirdn(i,j,klevp)-fdirup(i,j,klevp))*swvdr(i,j) +
     &         (fdifdn(i,j,klevp)-fdifup(i,j,klevp))*swvdf(i,j))
 
            fthru(i,j)  = fthru(i,j) +
     &         (fdirdn(i,j,klevp)-fdirup(i,j,klevp))*swvdr(i,j) +
     &         (fdifdn(i,j,klevp)-fdifup(i,j,klevp))*swvdf(i,j)
 
            do k=2,klev
              Iabs(i,j,k-1) = Iabs(i,j,k-1) +
     &        ((fdirdn(i,j,k)-fdirup(i,j,k))*swvdr(i,j) +
     &         (fdifdn(i,j,k)-fdifup(i,j,k))*swvdf(i,j))
     &       -((fdirdn(i,j,k+1)-fdirup(i,j,k+1))*swvdr(i,j) +
     &         (fdifdn(i,j,k+1)-fdifup(i,j,k+1))*swvdf(i,j))
            enddo
          enddo                  ! ij

        else if(ns > 1) then  ! near IR

          do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            aidr(i,j)   = aidr(i,j) + rupdir(i,j,0)*wghtns(i,j,ns)
            aidf(i,j)   = aidf(i,j) + rupdif(i,j,0)*wghtns(i,j,ns)
 
            fsfc(i,j)  = fsfc(i,j) + 
     &      ( ((fdirdn(i,j,0)-fdirup(i,j,0))*swidr(i,j) +
     &         (fdifdn(i,j,0)-fdifup(i,j,0))*swidf(i,j)) -
     &        ((fdirdn(i,j,2)-fdirup(i,j,2))*swidr(i,j) +
     &         (fdifdn(i,j,2)-fdifup(i,j,2))*swidf(i,j)) )
     &         *wghtns(i,j,ns)
 
            fint(i,j)  = fint(i,j) + 
     &      ( ((fdirdn(i,j,2)-fdirup(i,j,2))*swidr(i,j) +
     &         (fdifdn(i,j,2)-fdifup(i,j,2))*swidf(i,j)) -
     &        ((fdirdn(i,j,klevp)-fdirup(i,j,klevp))*swidr(i,j) +
     &         (fdifdn(i,j,klevp)-fdifup(i,j,klevp))*swidf(i,j)) )
     &         *wghtns(i,j,ns)
 
            fthru(i,j) = fthru(i,j) + 
     &        ((fdirdn(i,j,klevp)-fdirup(i,j,klevp))*swidr(i,j) +
     &         (fdifdn(i,j,klevp)-fdifup(i,j,klevp))*swidf(i,j))
     &         *wghtns(i,j,ns)
 
            do k=2,klev
              Iabs(i,j,k-1) = Iabs(i,j,k-1) +
     &        ((fdirdn(i,j,k)-fdirup(i,j,k))*swidr(i,j) +
     &         (fdifdn(i,j,k)-fdifup(i,j,k))*swidf(i,j))
     &         *wghtns(i,j,ns)
     &       -((fdirdn(i,j,k+1)-fdirup(i,j,k+1))*swidr(i,j) +
     &         (fdifdn(i,j,k+1)-fdifup(i,j,k+1))*swidf(i,j))
     &         *wghtns(i,j,ns)
            enddo
          enddo                 ! ij

        endif                   ! ns = 1

      enddo                     ! ns

      end subroutine compute_dEdd
 
!=======================================================================
!BOP
!
! !IROUTINE: solution_dEdd - evaluate solution for delta-Edddington
!
! !INTERFACE:
!
      subroutine solution_dEdd
     &      (nx_block, ny_block,
     &       icells,  indxi,  indxj,  coszen, srftyp,
     &       tau,     w0,     g,      albodr, albodf,
     &       trndir,  trntdr, trndif, rupdir, rupdif,
     &       rdndif)
!
! !DESCRIPTION:
!
! Given input optical properties, evaluate the monochromatic 
! delta-Eddington solution.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:

      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  icells  ! number of ice/ocean grid cells for surface type
 
      integer (kind=int_kind), dimension (nx_block*ny_block),
     &   intent(in) ::
     &   indxi   ! compressed indices for ice/ocean cells for surface type
     &,  indxj
 
      ! cosine solar zenith angle; negative for sun below horizon
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::
     &   coszen  ! cosine solar zenith angle
 
      integer (kind=int_kind), dimension(nx_block,ny_block),
     &   intent(in) ::
     &   srftyp  ! surface type over ice: (0=air, 1=snow, 2=pond)
 
      real (kind=dbl_kind), dimension(nx_block,ny_block,0:nilyr+1),
     &   intent(in) ::
     &   tau     ! layer extinction optical depth
     &,  w0      ! layer single scattering albedo
     &,  g       ! layer asymmetry parameter
 
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in)::
     &  albodr   ! ocean albedo to direct rad
     &, albodf   ! ocean albedo to diffuse rad
 
      ! following arrays are defined at model interfaces; 0 is the top of the
      ! layer above the sea ice; klevp=nilyr+2 is the sea ice/ocean interface.

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:nilyr+2) ::
     &   trndir  ! solar beam down transm from top
     &,  trntdr  ! total transmission to direct beam for layers above
     &,  trndif  ! diffuse transmission to diffuse beam for layers above
     &,  rupdir  ! ref to dir rad for layers below
     &,  rupdif  ! ref to dif rad for layers below
     &,  rdndif  ! ref to dif rad for layers above
!
!EOP
!-----------------------------------------------------------------------
! Sea ice model delta-Eddington     August 2005   Bruce P. Briegleb
!
! Computes layer reflectivities and transmissivities, from the top down
! to the lowest interface using the delta-Eddington solutions for each 
! layer; adds layers from top down to lowest interface, and from the
! lowest interface (underlying ocean) up to the top of the sea ice 
! surface, or to the top of the layer above if snow or pond present.
!
! Note that the diffuse reflectivity and transmissivity are computed
! by integrating the direct over several gaussian angles. This is
! because the diffuse reflectivity expression sometimes is negative.
!
! Assumes monochromatic or spectrally uniform properties across a band
! for the input optical parameters.
!
! If total transmission to the interface above a particular layer is
! less than trmin, then no further delta-Eddington solutions are
! evaluated for layers below.
!
! The following describes how refraction is handled in the calculation.
!
! First, we assume that radiation is refracted when entering sea ice at 
! the base of the surface layer, or water (i.e. melt pond); we assume that
! radiation does not refract when entering snow, nor upon entering sea ice
! from a melt pond, nor upon entering the underlying ocean from sea ice. 
!
! To handle refraction, we define a "fresnel" layer, which physically 
! is of neglible vertical extent and is non-absorbing, which can be
! added to any sea ice layer or top of melt pond. The following diagram 
! shows where the fresnel layers are in relation to other layer/interfaces
! in the standard configuration: 
!
!   Radiation layers and interfaces
!
! layer index             interface index
!
!       ----------------------  0                   <<< Fresnel
!       //////////////////////
!    0  air, snow or melt pond    Air or snow       Melt pond
!       //////////////////////
!       ----------------------  1  
!    1  \\\ surface layer \\\\
!       ----------------------  2  <<< Fresnel  
!    2
!       ----------------------  3
!    3
!       ----------------------  4
!    4
!       ----------------------  5
!    5
!       ----------------------  6
!       ~~~~~~~~~~~~~~~~~~~~~~
!       ~~~~~~~ ocean ~~~~~~~~
!       ~~~~~~~~~~~~~~~~~~~~~~
!
! Layers of ice are indexed from 1 to klev (present 5) downwards; k=0 
! refers to layer above top sea ice layer. This 0 layer is either air, 
! snow or melt pond. A fresnel layer is positioned at the top of a melt 
! pond or to the top of the layer just underneath the sea ice surface 
! layer if no melt pond lies over it. The fresnel layer accounts for 
! refraction of direct beam and associated reflection and transmission 
! for solar radiation. It is assumed that radiation crossing the air/snow, 
! melt pond/sea ice and sea ice/ocean interfaces is not refracted.
!
! Some caution must be exercised for the fresnel layer, because any layer
! to which it is added is no longer a homogeneous layer, as are all other
! individual layers. For all other layers for example, the direct and diffuse 
! reflectivities/transmissivities (R/T) are the same for radiation above or 
! below the layer. This is the implication of homogeneous! But for the fresnel 
! layer this is not so. Thus, the R/T for this layer must be distinguished 
! for radiation above from that from radiation below. For generality, we 
! treat all layers to be added as inhomogeneous.
!
!-----------------------------------------------------------------------

! Local

      integer (kind=int_kind), parameter ::
     &   klev   = nilyr + 1   ! number of radiation layers
     &,  klevp  = klev  + 1   ! number of radiation layers + 1
     &,  kfrsnl = 2  ! radiation index of sea ice layer whose top is fresnel lyr
 
      ! Following variables are defined for each layer; 0 refers to layer above 
      ! sea ice. In general we must distinguish directions above and below in 
      ! the diffuse reflectivity and transmissivity, as layers are not assumed
      ! to be homogeneous (apart from the single layer delta-Edd solutions); 
      ! the direct is always from above.
      real (kind=dbl_kind), dimension (nx_block,ny_block,0:klev) ::
     &   rdir    ! layer reflectivity to direct rad
     &,  rdif_a  ! layer reflectivity to diffuse rad from above
     &,  rdif_b  ! layer reflectivity to diffuse rad from below
     &,  tdir    ! layer transmission to direct rad (inc. diffusely transmitted)
     &,  tdif_a  ! layer transmission to diffuse rad from above
     &,  tdif_b  ! layer transmission to diffuse rad from below
     &,  trnlay  ! solar beam transm for layer (direct beam only)
 
      real (kind=dbl_kind), parameter ::
     &   trmin = 0.001_dbl_kind   ! minimum total transmission allowed
 
      integer (kind=int_kind) :: 
     &   i       ! longitude index
     &,  j       ! latitude index
     &,  ij      ! longitude/latitude index
     &,  k       ! level index
 
      integer (kind=int_kind) ::
     &   icellstrmin  ! number of ice/ocean grid cells for which tr > trmin
 
      integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &   indxitrmin   ! indices for ice/ocean cells for which tr > trmin
     &,  indxjtrmin
 
      real (kind=dbl_kind) ::
     &   tautot   ! layer optical depth
     &,  wtot     ! layer single scattering albedo
     &,  gtot     ! layer asymmetry parameter
     &,  ftot     ! layer forward scattering fraction
     &,  ts       ! layer scaled extinction optical depth
     &,  ws       ! layer scaled single scattering albedo
     &,  gs       ! layer scaled asymmetry parameter
     &,  rintfc   ! reflection (multiple) at an interface
     &,  refkp1   ! interface multiple scattering for k+1
     &,  refkm1   ! interface multiple scattering for k-1
     &,  tdrrdir  ! direct tran times layer direct ref 
     &,  tdndif   ! total down diffuse = tot tran - direct tran
 
      ! perpendicular and parallel relative to pane of incidence and scattering
      real (kind=dbl_kind) ::
     &   R1       ! perpendicular polarization reflection amplitude
     &,  R2       ! parallel polarization reflection amplitude
     &,  T1       ! perpendicular polarization transmission amplitude
     &,  T2       ! parallel polarization transmission amplitude
     &,  Rf_dir_a ! fresnel reflection to direct radiation
     &,  Tf_dir_a ! fresnel transmission to direct radiation
     &,  Rf_dif_a ! fresnel reflection to diff radiation from above
     &,  Rf_dif_b ! fresnel reflection to diff radiation from below
     &,  Tf_dif_a ! fresnel transmission to diff radiation from above
     &,  Tf_dif_b ! fresnel transmission to diff radiation from below
 
      ! diffuse fresnel reflectivities pre-computed
      real (kind=dbl_kind), parameter :: 
     &   refindx = 1.310_dbl_kind  ! refractive index of ice (used for water also)
     &,  cp063   = 0.063_dbl_kind  ! diffuse fresnel reflectivity from above
     &,  cp455   = 0.455_dbl_kind  ! diffuse fresnel reflectivity from below
 
      real (kind=dbl_kind) ::
     &   mu0      ! cosine solar zenith angle incident
     &,  mu0n     ! cosine solar zenith angle in medium
 
      real (kind=dbl_kind) ::
     &   alpha    ! term in direct reflectivity and transmissivity
     &,  gamma    ! term in direct reflectivity and transmissivity
     &,  el       ! term in alpha,gamma,n,u
     &,  taus     ! scaled extinction optical depth
     &,  omgs     ! scaled single particle scattering albedo
     &,  asys     ! scaled asymmetry parameter
     &,  u        ! term in diffuse reflectivity and transmissivity
     &,  n        ! term in diffuse reflectivity and transmissivity
     &,  lm       ! temporary for el
     &,  mu       ! cosine solar zenith for either snow or water
     &,  ne       ! temporary for n
 
      real (kind=dbl_kind) ::
     &   w        ! dummy argument for statement function
     &,  uu       ! dummy argument for statement function
     &,  gg       ! dummy argument for statement function
     &,  e        ! dummy argument for statement function
     &,  f        ! dummy argument for statement function
     &,  t        ! dummy argument for statement function
     &,  et       ! dummy argument for statement function
 
      real (kind=dbl_kind) ::
     &   alp      ! temporary for alpha
     &,  gam      ! temporary for gamma
     &,  ue       ! temporary for u
     &,  arg      ! exponential argument
     &,  extins   ! extinction
     &,  amg      ! alp - gam
     &,  apg      ! alp + gam
 
      integer (kind=int_kind), parameter ::
     &   ngmax = 8  ! number of gaussian angles in hemisphere
 
      real (kind=dbl_kind), dimension (ngmax) ::
         ! gaussian angles (points)
     &   gauspt = (/ .9894009_dbl_kind,  .9445750_dbl_kind,
     &               .8656312_dbl_kind,  .7554044_dbl_kind,
     &               .6178762_dbl_kind,  .4580168_dbl_kind,
     &               .2816036_dbl_kind,  .0950125_dbl_kind /)   
         ! gaussian weights (normalized so sum is 1)
     &,  gauswt = (/ .0271525_dbl_kind,  .0622535_dbl_kind,
     &               .0951585_dbl_kind,  .1246290_dbl_kind,
     &               .1495960_dbl_kind,  .1691565_dbl_kind,
     &               .1826034_dbl_kind,  .1894506_dbl_kind /)
  
      integer (kind=int_kind) ::
     &   ng       ! gaussian integration index
 
      real (kind=dbl_kind) ::
     &   gwt      ! gaussian weight
     &,  swt      ! sum of weights
     &,  trn      ! layer transmission
     &,  rdr      ! rdir for gaussian integration
     &,  tdr      ! tdir for gaussian integration
     &,  smr      ! accumulator for rdif gaussian integration
     &,  smt      ! accumulator for tdif gaussian integration
 
      ! Delta-Eddington solution expressions
      alpha(w,uu,gg,e) = p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
      gamma(w,uu,gg,e) = p5*w*((c1 + c3*gg*(c1-w)*uu*uu)
     &                  / (c1-e*e*uu*uu))
      n(uu,et)         = ((uu+c1)*(uu+c1)/et ) - ((uu-c1)*(uu-c1)*et)
      u(w,gg,e)        = c1p5*(c1 - w*gg)/e
      el(w,gg)         = sqrt(c3*(c1-w)*(c1 - w*gg))
      taus(w,f,t)      = (c1 - w*f)*t
      omgs(w,f)        = (c1 - f)*w/(c1 - w*f)
      asys(gg,f)       = (gg - f)/(c1 - f)
 
!-----------------------------------------------------------------------
 
      ! initialize all output to 0, so that nighttime 
      ! values from previous computations are not used:
      do k=0,klevp
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         trndir(i,j,k) = c0
         trntdr(i,j,k) = c0
         trndif(i,j,k) = c0
         rupdir(i,j,k) = c0
         rupdif(i,j,k) = c0
         rdndif(i,j,k) = c0
      enddo 
      enddo
 
      ! for the layer over sea ice (either air, snow or water):
      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)
 
          ! initialize top interface of top layer:
          trndir(i,j,0) =   c1
          trntdr(i,j,0) =   c1
          trndif(i,j,0) =   c1
          rdndif(i,j,0) =   c0
 
          ! compute solar zenith angle in medium
          mu0  = max(coszen(i,j),p01)
          mu0n = sqrt(c1-((c1-mu0*mu0)/(refindx*refindx)))
 
          ! set layer properties if air:
          if( srftyp(i,j) == 0 ) then
 
            rdir(i,j,0)   = c0
            rdif_a(i,j,0) = c0
            rdif_b(i,j,0) = c0
            tdir(i,j,0)   = c1
            tdif_a(i,j,0) = c1
            tdif_b(i,j,0) = c1
            trnlay(i,j,0) = c1
 
          ! compute layer properties if snow or water:
          else if( srftyp(i,j) > 0 ) then
 
            ! Set appropriate direct beam angle:
            mu = mu0
            if( srftyp(i,j) == 2 ) mu = mu0n
 
            tautot  = tau(i,j,0)
            wtot    = w0(i,j,0)
            gtot    = g(i,j,0)
            ftot    = gtot*gtot
 
            ts   = taus(wtot,ftot,tautot)
            ws   = omgs(wtot,ftot)
            gs   = asys(gtot,ftot)
            lm   = el(ws,gs)
            ue   = u(ws,gs,lm)
 
            ! use approximate exponential
            arg  = lm*ts
            if( arg > argmax ) then
              extins   = exp_min
            else
              indx_exp = arg/dx_exp
              extins   = app_exp(indx_exp)
            endif
            ne = n(ue,extins)
 
            ! first calculation of rdif, tdif:
            rdif_a(i,j,0) = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne
            tdif_a(i,j,0) = c4*ue/ne
 
            ! evaluate rdir,tdir:
            ! use approximate exponential
            arg  = ts/mu
            if( arg > argmax ) then
              trnlay(i,j,0) = exp_min
            else
              indx_exp      = arg/dx_exp
              trnlay(i,j,0) = app_exp(indx_exp)
            endif
            alp = alpha(ws,mu,gs,lm)
            gam = gamma(ws,mu,gs,lm)
            apg = alp + gam
            amg = alp - gam
            rdir(i,j,0) = amg*(tdif_a(i,j,0)*trnlay(i,j,0) - c1) + 
     &                  apg*rdif_a(i,j,0)
            tdir(i,j,0) = apg*tdif_a(i,j,0) +
     &                  (amg*rdif_a(i,j,0) - (apg-c1))*trnlay(i,j,0)
 
            ! recalculate rdif,tdif using gaussian integration over rdir,tdir:
            swt = c0
            smr = c0
            smt = c0
            do ng=1,ngmax
              mu  = gauspt(ng)
              gwt = gauswt(ng)
              swt = swt + mu*gwt
              ! use approximate exponential
              arg = ts/mu
              if( arg > argmax ) then
                trn      = exp_min
              else
                indx_exp = arg/dx_exp
                trn      = app_exp(indx_exp)
              endif
              alp = alpha(ws,mu,gs,lm)
              gam = gamma(ws,mu,gs,lm)
              apg = alp + gam
              amg = alp - gam
              rdr = amg*(tdif_a(i,j,0)*trn-c1) + 
     &              apg*rdif_a(i,j,0)
              tdr = apg*tdif_a(i,j,0) +
     &              (amg*rdif_a(i,j,0)-(apg-c1))*trn
              smr = smr + mu*rdr*gwt
              smt = smt + mu*tdr*gwt
            enddo
            rdif_a(i,j,0) = smr/swt
            tdif_a(i,j,0) = smt/swt
 
            ! homogeneous layer
            rdif_b(i,j,0) = rdif_a(i,j,0)
            tdif_b(i,j,0) = tdif_a(i,j,0)
 
          endif ! srftyp for top layer
 
          ! set top interface terms for sea ice layer 1: note that rdndif
          ! is diffuse reflectivity from below to combined layers above.
          trndir(i,j,1) = trnlay(i,j,0)
          trntdr(i,j,1) = tdir(i,j,0)
          trndif(i,j,1) = tdif_a(i,j,0)
          rdndif(i,j,1) = rdif_b(i,j,0)
 
          ! compute fresnel reflection and transmission amplitudes
          ! for two polarizations: 1=perpendicular and 2=parallel to
          ! the plane containing incident, reflected and refracted rays.
          R1 = (mu0 - refindx*mu0n) / 
     &         (mu0 + refindx*mu0n)
          R2 = (refindx*mu0 - mu0n) / 
     &         (refindx*mu0 + mu0n)
          T1 = c2*mu0 / 
     &         (mu0 + refindx*mu0n)
          T2 = c2*mu0 / 
     &         (refindx*mu0 + mu0n)
 
          ! unpolarized light for direct beam
          Rf_dir_a = p5 * (R1*R1 + R2*R2)
          Tf_dir_a = p5 * (T1*T1 + T2*T2)*refindx*mu0n/mu0
 
          ! precalculated diffuse reflectivities and transmissivities
          ! for incident radiation above and below fresnel layer, using
          ! the direct albedos and accounting for complete internal
          ! reflection from below; precalculated because high order
          ! number of gaussian points (~256) is required for convergence:
 
          ! above
          Rf_dif_a = cp063
          Tf_dif_a = c1 - Rf_dif_a
          ! below
          Rf_dif_b = cp455
          Tf_dif_b = c1 - Rf_dif_b
 
          ! add fresnel layer to top of melt pond. The usual two-layer adding 
          ! formulas are used, with modifications for the direct transmission 
          ! in particular, since there is no equivalent exponential transmission 
          ! term for the fresnel layer. Note that the melt pond reflection and 
          ! transmission to direct radiation have already accounted for the 
          ! angular refraction from the fresnel layer.
          if( srftyp(i,j) == 2 ) then
            rintfc        = c1/(c1 - Rf_dif_b*rdif_a(i,j,0))
            tdir(i,j,0)   = Tf_dir_a*tdir(i,j,0) + 
     &                      Tf_dir_a*rdir(i,j,0)*
     &                      Rf_dif_b*rintfc*tdif_a(i,j,0)
            rdir(i,j,0)   = Rf_dir_a + 
     &                      Tf_dir_a*rdir(i,j,0)*rintfc*Tf_dif_b
            rdif_a(i,j,0) = Rf_dif_a + 
     &                      Tf_dif_a*rdif_a(i,j,0)*rintfc*Tf_dif_b
            rdif_b(i,j,0) = rdif_b(i,j,0) + 
     &                      tdif_b(i,j,0)*Rf_dif_b*rintfc*tdif_a(i,j,0)
            tdif_a(i,j,0) = Tf_dif_a*rintfc*tdif_a(i,j,0)
            tdif_b(i,j,0) = tdif_b(i,j,0)*rintfc*Tf_dif_b
            trndir(i,j,1) = Tf_dir_a*trnlay(i,j,0)
            trntdr(i,j,1) = tdir(i,j,0)
            trndif(i,j,1) = tdif_a(i,j,0)
            rdndif(i,j,1) = rdif_b(i,j,0)
          endif
 
      enddo ! end of layer over sea ice (either air, snow or water)
 
      ! now, continue down one layer at a time; if the total transmission to
      ! the interface just above a given layer is less than trmin, then no
      ! delta-eddington computation for that layer is done:
      do k=1,klev
 
         ! initialize current layer properties to zero; only if total
         ! transmission to the top interface of the current layer exceeds the
         ! minimum, will these values be computed below:
         do ij = 1, icells
           i = indxi(ij)
           j = indxj(ij)
 
               rdir(i,j,k)   =  c0
               rdif_a(i,j,k) =  c0
               rdif_b(i,j,k) =  c0
               tdir(i,j,k)   =  c0
               tdif_a(i,j,k) =  c0
               tdif_b(i,j,k) =  c0
               trnlay(i,j,k) =  c0
 
               ! Calculates the solar beam transmission, total transmission, and
               ! reflectivity for diffuse radiation from below at the top of the
               ! current layer:
               !
               !              layers       interface
               !         
               !       ---------------------  k-1 
               !                k-1
               !       ---------------------  k
               !                 k
               !       ---------------------  
              if ( k > 1 ) then
                trndir(i,j,k) = trndir(i,j,k-1)*trnlay(i,j,k-1)
                refkm1        = c1/(c1-rdndif(i,j,k-1)*rdif_a(i,j,k-1))
                tdrrdir       = trndir(i,j,k-1)*rdir(i,j,k-1)
                tdndif        = trntdr(i,j,k-1) - trndir(i,j,k-1)
                trntdr(i,j,k) = trndir(i,j,k-1)*tdir(i,j,k-1) + 
     &                         (tdndif + tdrrdir*rdndif(i,j,k-1))
     &                         *refkm1*tdif_a(i,j,k-1)
                rdndif(i,j,k) = rdif_b(i,j,k-1)  +
     &                         (tdif_b(i,j,k-1)*rdndif(i,j,k-1)
     &                         *refkm1*tdif_a(i,j,k-1))
                trndif(i,j,k) = trndif(i,j,k-1)*tdif_a(i,j,k-1)*refkm1
              endif
         enddo ! end of top interface calculation
 
         ! compute next layer delta-eddington solution only if total transmission
         ! of radiation to the interface just above the layer exceeds trmin.
         icellstrmin = 0
         do ij = 1, icells
           i = indxi(ij)
           j = indxj(ij)
           if (trntdr(i,j,k) > trmin ) then
             icellstrmin   = icellstrmin + 1
             indxitrmin(icellstrmin) = i
             indxjtrmin(icellstrmin) = j
           endif
         enddo
 
            ! calculation over layers with penetrating radiation
            do ij=1,icellstrmin
               i = indxitrmin(ij)
               j = indxjtrmin(ij)
 
               tautot  = tau(i,j,k)
               wtot    = w0(i,j,k)
               gtot    = g(i,j,k)
               ftot    = gtot*gtot
 
               ts   = taus(wtot,ftot,tautot)
               ws   = omgs(wtot,ftot)
               gs   = asys(gtot,ftot)
               lm   = el(ws,gs)
               ue   = u(ws,gs,lm)
               ! mu0 in medium
               mu0  = max(coszen(i,j),p01)
               mu0n = sqrt(c1-((c1-mu0*mu0)/(refindx*refindx)))
               ! do not change if above fresnel layer and in non-pond medium
               if( srftyp(i,j) .ne. 2 .and. k < kfrsnl ) mu0n = mu0
 
               ! use approximate exponential
               arg  = lm*ts
               if( arg > argmax ) then
                 extins   = exp_min
               else
                 indx_exp = arg/dx_exp
                 extins   = app_exp(indx_exp)
               endif
               ne = n(ue,extins)
 
               ! first calculation of rdif, tdif:
               rdif_a(i,j,k) = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne
               tdif_a(i,j,k) = c4*ue/ne
 
               ! evaluate rdir,tdir:
               ! use approximate exponential
               arg  = ts/mu0n
               if( arg > argmax ) then
                 trnlay(i,j,k) = exp_min
               else
                 indx_exp      = arg/dx_exp
                 trnlay(i,j,k) = app_exp(indx_exp)
               endif
               alp = alpha(ws,mu0n,gs,lm)
               gam = gamma(ws,mu0n,gs,lm)
               apg = alp + gam
               amg = alp - gam
               rdir(i,j,k) = amg*(tdif_a(i,j,k)*trnlay(i,j,k) - c1) + 
     &                     apg*rdif_a(i,j,k)
               tdir(i,j,k) = apg*tdif_a(i,j,k) +
     &                     (amg*rdif_a(i,j,k) - (apg-c1))*trnlay(i,j,k)
 
               ! recalculate rdif,tdif using gaussian integration over rdir,tdir:
               swt = c0
               smr = c0
               smt = c0
               do ng=1,ngmax
                 mu  = gauspt(ng)
                 gwt = gauswt(ng)
                 swt = swt + mu*gwt
                 ! use approximate exponential
                 arg = ts/mu
                 if( arg > argmax ) then
                   trn      = exp_min
                 else
                   indx_exp = arg/dx_exp
                   trn      = app_exp(indx_exp)
                 endif
                 alp = alpha(ws,mu,gs,lm)
                 gam = gamma(ws,mu,gs,lm)
                 apg = alp + gam
                 amg = alp - gam
                 rdr = amg*(tdif_a(i,j,k)*trn-c1) + 
     &                 apg*rdif_a(i,j,k)
                 tdr = apg*tdif_a(i,j,k) +
     &                 (amg*rdif_a(i,j,k)-(apg-c1))*trn
                 smr = smr + mu*rdr*gwt
                 smt = smt + mu*tdr*gwt
               enddo
               rdif_a(i,j,k) = smr/swt
               tdif_a(i,j,k) = smt/swt
 
               ! homogeneous layer
               rdif_b(i,j,k) = rdif_a(i,j,k)
               tdif_b(i,j,k) = tdif_a(i,j,k)
 
               ! add fresnel layer to top of desired layer if either 
               ! air or snow overlies ice; we ignore refraction in ice 
               ! if a melt pond overlies it:
               if( srftyp(i,j) < 2 .and. k == kfrsnl ) then
                 rintfc             = c1 /
     &                               (c1-Rf_dif_b*rdif_a(i,j,kfrsnl))
                 tdir(i,j,kfrsnl)   = Tf_dir_a*tdir(i,j,kfrsnl) + 
     &                                Tf_dir_a*rdir(i,j,kfrsnl) *
     &                                Rf_dif_b*rintfc*tdif_a(i,j,kfrsnl)
                 rdir(i,j,kfrsnl)   = Rf_dir_a + 
     &                                Tf_dir_a*rdir(i,j,kfrsnl) *
     &                                rintfc*Tf_dif_b
                 rdif_a(i,j,kfrsnl) = Rf_dif_a + 
     &                                Tf_dif_a*rdif_a(i,j,kfrsnl) *
     &                                rintfc*Tf_dif_b
                 rdif_b(i,j,kfrsnl) = rdif_b(i,j,kfrsnl) + 
     &                                tdif_b(i,j,kfrsnl)*Rf_dif_b*
     &                                rintfc*tdif_a(i,j,kfrsnl)
     &                                
                 tdif_a(i,j,kfrsnl) = Tf_dif_a*rintfc*tdif_a(i,j,kfrsnl)
                 tdif_b(i,j,kfrsnl) = tdif_b(i,j,kfrsnl)*rintfc*Tf_dif_b

                 ! update trnlay to include fresnel transmission
                 trnlay(i,j,kfrsnl) = Tf_dir_a*trnlay(i,j,kfrsnl)
               endif
            enddo ! end transmission large enough
 
      enddo ! end level k loop
 
      ! compute total direct beam transmission, total transmission, and
      ! reflectivity for diffuse radiation (from below) for all layers
      ! above the surface; note that we ignore refraction between ice
      ! and underlying ocean:
      !
      !       For k = klevp
      !
      !              layers       interface
      !
      !       ---------------------  k-1 
      !                k-1
      !       ---------------------  k
      !       \\\\\\\ ocean \\\\\\\
      k = klevp
      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)
           trndir(i,j,k) = trndir(i,j,k-1)*trnlay(i,j,k-1)
           refkm1        = c1/(c1 - rdndif(i,j,k-1)*rdif_a(i,j,k-1))
           tdrrdir       = trndir(i,j,k-1)*rdir(i,j,k-1)
           tdndif        = trntdr(i,j,k-1) - trndir(i,j,k-1)
           trntdr(i,j,k) = trndir(i,j,k-1)*tdir(i,j,k-1) + 
     &       (tdndif + tdrrdir*rdndif(i,j,k-1))*refkm1*tdif_a(i,j,k-1)
           rdndif(i,j,k) = rdif_b(i,j,k-1)  +
     &       (tdif_b(i,j,k-1)*rdndif(i,j,k-1)*refkm1*tdif_a(i,j,k-1))
           trndif(i,j,k) = trndif(i,j,k-1)*tdif_a(i,j,k-1)*refkm1
      enddo ! end surface layer
 
      ! compute reflectivity to direct and diffuse radiation for layers below by
      ! adding succesive layers starting from the surface and working upwards:
      !
      !              layers       interface
      !
      !       ---------------------  k
      !                 k
      !       ---------------------  k+1
      !                k+1
      !       ---------------------
      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)
        rupdir(i,j,klevp) = albodr(i,j)
        rupdif(i,j,klevp) = albodf(i,j)
      enddo
      do k=klev,0,-1
        do ij = 1, icells
          i = indxi(ij)
          j = indxj(ij)
          ! interface scattering
          refkp1        = c1/( c1 - rdif_b(i,j,k)*rupdif(i,j,k+1))
          ! dir from top layer plus exp tran ref from lower layer, interface
          ! scattered and tran thru top layer from below, plus diff tran ref
          ! from lower layer with interface scattering tran thru top from below
          rupdir(i,j,k) = rdir(i,j,k) + 
     &                  ( trnlay(i,j,k)*rupdir(i,j,k+1) +
     &                   (tdir(i,j,k)-trnlay(i,j,k))*rupdif(i,j,k+1) ) *
     &                    refkp1*tdif_b(i,j,k)
          ! dif from top layer from above, plus dif tran upwards reflected and
          ! interface scattered which tran top from below
          rupdif(i,j,k) = rdif_a(i,j,k) +
     &                    tdif_a(i,j,k)*rupdif(i,j,k+1)*
     &                    refkp1*tdif_b(i,j,k)
        enddo
      enddo
 
      end subroutine solution_dEdd
 
! End delta-Eddington shortwave method
!=======================================================================

      end module ice_shortwave

!=======================================================================
