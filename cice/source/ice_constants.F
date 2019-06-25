c $Id: $
c=======================================================================
!---! This module defines a variety of physical and numerical constants
!---! used throughout the ice model
!---!
!---! author Elizabeth C. Hunke
!---!        Fluid Dynamics Group, Los Alamos National Laboratory
!---!
!---! code originally based on constants.F in POP
c=======================================================================

      module ice_constants

      use ice_kinds_mod
      use ice_domain

      implicit none
      save

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter ::
     &   gravit    = 9.80616_dbl_kind   ! gravitational acceleration (m/s^2)
     &,  radius    = 6.37e6_dbl_kind    ! earth radius (m)
     &,  rhoi      = 917.0_dbl_kind          ! density of ice (kg/m^3)
     &,  rhow      = 1026.0_dbl_kind         ! density of seawater (kg/m^3)
     &,  rhos      = 330.0_dbl_kind          ! density of snow (kg/m^3)
     &,  rhofresh  = 1000.0_dbl_kind          ! density of fresh water (kg/m^3)
     &,  vonkar    = 0.4_dbl_kind            ! von Karman constant
     &,  cp_air    = 1005.0_dbl_kind         ! heat capacity of air (J/kg/K)
     &,  cpwv   =  1.81e3_dbl_kind       ! Specific heat of water vapor
     &,  stefan_boltzmann = 567.0e-10_dbl_kind !  W/m^2/K^4
     &,  Tffresh = 273.15_dbl_kind   ! Freezing temp of fresh ice (K)
     &,  depressT = 0.054_dbl_kind      ! Tf:brine salinity ratio (C/ppt)
c     &,  cp_ice = 2054._dbl_kind       ! heat capacity of fresh ice (J/kg/K)
     &,  cp_ice = 2106._dbl_kind ! specific heat of fresh ice (J kg-1 deg-1)
     &,  cp_sno = 0.0_dbl_kind   ! specific heat of snow     (J/kg/K)
     &,  cp_ocn = 4218._dbl_kind ! specific heat of ocn  (J/kg/K)
     &,  Lsub = 2.835e6_dbl_kind ! latent heat, sublimation freshwater (J/kg)
     &,  Lvap = 2.501e6_dbl_kind ! latent heat, vaporization freshwater (J/kg)
     &,  Lfresh = Lsub - Lvap   ! latent heat of melting of fresh ice (J kg-1)
     &,  Lfus = Lsub-Lvap ! latent heat of fusion freshwater          (J/kg)
     &,  Rvap = 461.           ! gas constant for H2O vapor (J kg-1 deg-1)
     &,  Timelt = 0.0_dbl_kind  ! melting temperature of ice top surface  (C)
     &,  Tsmelt = 0.0_dbl_kind  ! melting temperature of snow top surface (C)

         ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
     &,  kappav = 1.4_dbl_kind  ! vis extnctn coef in ice, wvlngth<700nm (1/m)
     &,  kappan = 17.6_dbl_kind ! vis extnctn coef in ice, wvlngth>700nm (1/m)

         ! (Briegleb JGR 97 11475-11485  July 1992)
     &,  emissivity = 0.95_dbl_kind     ! emissivity of snow and ice

     &,  kice  = 2.0340_dbl_kind   ! thermal conductivity of fresh ice(W/m/deg)
     &,  ksno  = 0.3100_dbl_kind   ! thermal conductivity of snow  (W/m/deg)
     &,  zref = 10._dbl_kind         ! reference height for stability (m)
     &,  iceruf = 0.0005_dbl_kind  ! ice surface roughness (m)
     &,  ice_ref_salinity = 4._dbl_kind  ! (psu)
c     &,  ocn_ref_salinity   = 34.7_dbl_kind  ! (psu)
c     &,  rho_air   = 1.2_dbl_kind            ! ambient air density (kg/m^3)

      !-----------------------------------------------------------------
      ! derived and miscellaneous parameters
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter ::
     &   rLfi = Lfus*rhoi   ! latent heat of fusion ice              (J/m**3)
     &,  rLfs = Lfus*rhos   ! specific latent heat of fusion snow    (J/m**3)
     &,  rLvi = Lvap*rhoi   ! specific latent heat of vapor*rhoice    (J/m**3)
     &,  rLvs = Lvap*rhos   ! specific latent heat of vapor*rhosno    (J/m**3)
     &,  rcpi = cp_ice*rhoi ! specific heat capacity of fresh ice     (J/m**3)
     &,  rcps = cp_sno*rhos ! specific heat capacity of snow          (J/m**3)
     &,  rcpidepressT = rcpi*depressT ! param for finding T(z) from q (J/m**3)
     &,  rLfidepressT = rLfi*depressT ! param for heat capacity   (J deg/m**3)
         ! heat capacity of sea ice, rhoi*C=rcpi+rLfidepressT*salinity/T**2

      !-----------------------------------------------------------------
      ! numbers
      !-----------------------------------------------------------------

      character (char_len) ::
     &  char_blank          ! empty character string

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::
     &  one                 ! array of ones which is often useful

      real (kind=dbl_kind), parameter ::
     &  c0   = 0.0_dbl_kind,
     &  c1   = 1.0_dbl_kind,
     &  c2   = 2.0_dbl_kind,
     &  c3   = 3.0_dbl_kind,
     &  c4   = 4.0_dbl_kind,
     &  c5   = 5.0_dbl_kind,
     &  c6   = 6.0_dbl_kind,
     &  c7   = 7.0_dbl_kind,
     &  c8   = 8.0_dbl_kind,
     &  c9   = 9.0_dbl_kind,
     &  c10  = 10.0_dbl_kind,
     &  c12  = 12.0_dbl_kind,
     &  c15  = 15.0_dbl_kind,
     &  c16  = 16.0_dbl_kind,
     &  c20  = 20.0_dbl_kind,
     &  c25  = 25.0_dbl_kind,
     &  c90  = 90.0_dbl_kind,
     &  c100 = 100.0_dbl_kind,
     &  c360 = 360.0_dbl_kind,
     &  c365 = 365.0_dbl_kind,
     &  c1000= 1000.0_dbl_kind,
     &  p001 = 0.001_dbl_kind,
     &  p01  = 0.01_dbl_kind,
     &  p1   = 0.1_dbl_kind,
     &  p2   = 0.2_dbl_kind,
     &  p5   = 0.5_dbl_kind,
     &  p25  = 0.25_dbl_kind,
     &  p166 = c1/c6,
     &  p333 = c1/c3,
     &  p666 = c2/c3,
     &  p111 = c1/c9,
     &  p055 = p111*p5,
     &  p027 = p055*p5,
     &  p222 = c2/c9,
     &  eps04  = 1.0e-4_dbl_kind,
     &  eps11  = 1.0e-11_dbl_kind,
     &  eps12  = 1.0e-12_dbl_kind,
     &  eps13  = 1.0e-13_dbl_kind,
     &  eps15  = 1.0e-15_dbl_kind,
     &  eps16  = 1.0e-16_dbl_kind,
c     &  puny = eps13
     &  puny = eps11

      real (kind=dbl_kind) :: 
     &  pi, pih, pi2         ! pi, pi/2 and 2pi

      !-----------------------------------------------------------------
      ! conversion factors
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter :: 
     &  cm_to_m       = 0.01_dbl_kind        ! cm to meters
     &, m_to_cm       = 100._dbl_kind        ! meters to cm
     &, m2_to_km2     = 1.e-6_dbl_kind       ! m^2 to km^2
     &, mps_to_cmpdy  = 8.64e6_dbl_kind ! m per s to cm per day
     &, mps_to_cmpyr  = mps_to_cmpdy*365._dbl_kind ! m per s to cm per yr

      real (kind=dbl_kind) :: 
     &  rad_to_deg               ! degree-radian conversion

c=======================================================================

      contains

c=======================================================================

      subroutine init_constants

!---!-------------------------------------------------------------------
!---! Initializes constants that are best defined at run time (e.g. pi)
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: n, i, j

      do j=jlo,jhi
       do i=ilo,ihi
        one(i,j) = c1
       enddo
      enddo

      pi  = c4*atan(c1)
      pi2 = c2*pi
      pih = p5*pi

      rad_to_deg = 180.0_dbl_kind/pi

      do n=1,char_len
        char_blank(n:n) = ' '
      end do

      end subroutine init_constants

c=======================================================================

      end module ice_constants

c=======================================================================
