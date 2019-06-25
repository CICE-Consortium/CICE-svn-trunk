c=======================================================================
c SCCS module : index.h   r3.1
c SCCS get    : 03 May 1996 / 16:59:36
c SCCS edit   : 03 May 1996 / 16:59:24
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.index.h
c=======================================================================
c  Sign convention:
c    positive value <=> downward flux
c
c  Unit convention:
c    T             ~ Kelvin
c    u,v           ~ m/s
c    pressure      ~ N/m^2 = Pa
c    heat flux     ~ W/m^2 
c    momentum flux ~ N/m^2
c    water flux    ~ (kg/s)/m^2
c
c-----------------------------------------------------------------------
c edited by E. C. Hunke 9 April 1998

      ! Faia(k) buffer (atm/ice fluxes computed by atm) indices:

      integer           
     &   kFaia_netsw ,  ! net shortwave     heat flux
     &   kFaia_lwdn  ,  ! downward longwave heat flux
     &   kFaia_prec     ! precipitation    water flux
         
      parameter( kFaia_netsw =  1)
      parameter( kFaia_lwdn  =  2)
      parameter( kFaia_prec  =  3)

c-----------------------------------------------------------------------

      ! Faid(k) buffer (atm/ice fluxes computed by drv) indices:

      integer      
     &   kFaid_taux  ,  ! wind stress, zonal
     &   kFaid_tauy     ! wind stress, meridional
          
      parameter( kFaid_taux  =  1)
      parameter( kFaid_tauy  =  2)
          
c-----------------------------------------------------------------------

      ! Faii(k) buffer (atm/ice fluxes computed by ice) indices:

      integer      
     &   kFaii_lat   ,  ! latent          heat flux
     &   kFaii_sen   ,  ! sensible        heat flux
     &   kFaii_lwup  ,  ! upward longwave heat flux
     &   kFaii_evap     ! evaporative water flux
          
      parameter( kFaii_lat   =  1)
      parameter( kFaii_sen   =  2)
      parameter( kFaii_lwup  =  3)
      parameter( kFaii_evap  =  4)
          
c
c=======================================================================

      ! Fiod(k) buffer (ice/ocn fluxes computed by drv) indices:

      integer
     &   kFiod_q        ! oceanic heat flux to ice

      parameter( kFiod_q = 3)

c-----------------------------------------------------------------------

      ! Fioi(k) buffer (ice/ocn fluxes computed by ice) indices:

      integer      
     &   kFioi_netsw ,  ! net shortwave heat flux penetrating ice
     &   kFioi_melth ,  ! heat  flux from melting ice
     &   kFioi_meltw ,  ! water flux from melting ice
     &   kFioi_tauxo ,  ! ice/ocean stress, zonal     , sent to ocn
     &   kFioi_tauyo    ! ice/ocean stress, meridional, sent to ocn
          
      parameter( kFioi_netsw =  1)
      parameter( kFioi_melth =  2)
      parameter( kFioi_meltw =  3)
      parameter( kFioi_tauxo =  4)
      parameter( kFioi_tauyo =  5)

c-----------------------------------------------------------------------

      ! Fioo(k) buffer (ice/ocn fluxes computed by ocn) indices:

      integer      
     &   kFioo_dhdx  ,  ! surface slope, zonal
     &   kFioo_dhdy     ! surface slope, meridional
          
      parameter( kFioo_dhdx  =  1)
      parameter( kFioo_dhdy  =  2)
c
c=======================================================================

      ! Sa(k) buffer (atm state) indices:

      integer      
     &   kSa_u       ,  ! wind, zonal
     &   kSa_v       ,  ! wind, meridional
     &   kSa_ptem    ,  ! potential temperature
     &   kSa_shum       ! specific humidity

      parameter( kSa_u       =  2)
      parameter( kSa_v       =  3)
      parameter( kSa_ptem    =  4)
      parameter( kSa_shum    =  5)

c-----------------------------------------------------------------------

      ! Si(k) buffer (ice state) indices:

      integer      
     &   kSi_t       ,  ! temperature
     &   kSi_ifrac   ,  ! fractional ice coverage
     &   kSi_avsdr   ,  ! albedo: visible      , direct
     &   kSi_anidr   ,  ! albedo: near infrared, direct
     &   kSi_avsdf   ,  ! albedo: visible      , diffuse
     &   kSi_anidf   ,  ! albedo: near infrared, diffuse
     &   kSi_u       ,  ! velocity, zonal
     &   kSi_v          ! velocity, meridional
          
      parameter( kSi_t       =  1)
      parameter( kSi_ifrac   =  2)
      parameter( kSi_avsdr   =  3)
      parameter( kSi_anidr   =  4)
      parameter( kSi_avsdf   =  5)
      parameter( kSi_anidf   =  6)
      parameter( kSi_u       =  7)
      parameter( kSi_v       =  8)

c-----------------------------------------------------------------------

      ! So(k) buffer (ocn state) indices:

      integer      
     &   kSo_u       ,  ! velocity, zonal
     &   kSo_v       ,  ! velocity, meridional
     &   kSo_t       ,  ! temperature
     &   kSo_s          ! salinity
          
      parameter( kSo_u       =  1)
      parameter( kSo_v       =  2)
      parameter( kSo_t       =  3)
      parameter( kSo_s       =  4)

c=======================================================================
      ! Fixx(k) buffer (ice input flux) indices:

      integer      
     &   kFixx_taux  ,  ! x surface stress
     &   kFixx_tauy  ,  ! y surface stress
     &   kFixx_dhdx  ,  ! surface slope, zonal
     &   kFixx_dhdy  ,  ! surface slope, meridional
     &   kFixx_netsw ,  ! net shortwave     heat flux
     &   kFixx_lwd   ,  ! longwave down     heat flux
     &   kFixx_water ,  ! precip + evap atm/ice    water flux
     &   kFixx_q        ! heat of fusion/melting potential from ocn

      parameter( kFixx_taux  =  1)
      parameter( kFixx_tauy  =  2)
      parameter( kFixx_dhdx  =  3)
      parameter( kFixx_dhdy  =  4)
      parameter( kFixx_netsw =  5)
      parameter( kFixx_lwd   =  6)
      parameter( kFixx_water =  7)
      parameter( kFixx_q     =  8)

c-----------------------------------------------------------------------
