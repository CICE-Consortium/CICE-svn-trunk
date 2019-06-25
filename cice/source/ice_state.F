c$ Id: $
c=======================================================================
!---! primary state variables in various configurations
!---! Note: other state variables are at the end of this...
!---!-------------------------------------------------------------------
!---!   main       single-column  single column,   aggregated    units
!---!                             single category  over cats      
!---!-------------------------------------------------------------------
!---! aicen(i,j,n)     ain(n)        ai             aice(i,j)     ---
!---! vicen(i,j,n)     vin(n)        vi or hi       vice(i,j)      m
!---! vsnon(i,j,n)     vsn(n)        vs or hs       vsno(i,j)      m
!---! eicen(i,j,k)     ein(k,n)      ei(k)          eice(i,j)      J/m^2
!---! esnon(i,j,n)     esn(n)        es             esno(i,j)      J/m^2
!---! Tsfcn(i,j,n)     Tsfn(n)       Tsf            Tsfc(i,j)      deg
!---!
!---! Area is dimensionless because aice is the fractional area 
!---! (normalized so that the sum over all categories, including open 
!---! water, is one).  That is why vice/vsno have units of m instead of 
!---! m^3, and eice/esno have units of J/m^2 instead of J.
!---!
!---! Note that variable names follow 2 simple rules:
!---!
!---! (1) If a variable has dimensions i and j, write 'ice' or 'sno' or 
!---! 'sfc'; else abbreviate as 'i' or 's' or 'sf'.
!---! (2) If a variable has dimension n, add an 'n'; else omit the 'n'.
!---!
!---! We will also use the variables hice = vice/aice, hsno = vsno/aice, 
!---! qice = eice/vice, and qsno = esno/vsno.  These variables will 
!---! follow the same naming convention.  They will not be state 
!---! variables but will be computed as needed in individual subroutines.
!---!
!---! authors C. M. Bitz, UW
!---!         Elizabeth C. Hunke, LANL
c=======================================================================

      module ice_state

      use ice_kinds_mod
      use ice_model_size
      use ice_domain

      implicit none

      !-----------------------------------------------------------------
      ! state of the ice aggregated over all categories
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (imt_local,jmt_local) ::
     &   aice      ! concentration of ice 
     &,  vice      ! volume per unit area of ice          (m)
     &,  vsno      ! volume per unit area of snow         (m)
     &,  Tsfc      ! temperature of ice/snow top surface  (C)
     &,  eice      ! energy of melt. of ice           (J/m^2)
     &,  esno      ! energy of melt. of snow layer    (J/m^2)
     &,  tmass     ! total mass of ice and snow
     &,  aice_init ! concentration of ice at beginning of dt (for diagnostics)

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (imt_local,jmt_local,ncat) ::
     &   aicen     ! concentration of ice 
     &,  vicen     ! volume per unit area of ice          (m)
     &,  vsnon     ! volume per unit area of snow         (m)
     &,  Tsfcn     ! temperature of ice/snow top surface  (C)
     &,  esnon     ! energy of melt. of snow layer    (J/m^2)

      real (kind=dbl_kind) ::
     &   aice0(imt_local,jmt_local)        ! concentration of open water
     &,  eicen(imt_local,jmt_local,ntilay) ! energy of melting for
                                           ! each ice layer (J/m^2)

      !-----------------------------------------------------------------
      ! state of the ice for each column
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (ncat) ::
     &   ain       ! concentration of ice 
     &,  vin       ! volume per unit area of ice          (m)
     &,  vsn       ! volume per unit area of snow         (m)
     &,  Tsfn      ! temperature of ice/snow top surface  (C)

      real (kind=dbl_kind) ::
     &   ai0       ! concentration of open water
     &,  ein(nilyr,ncat) ! ice enthalpy in 2D          (J/m^2)
     &,  esn(ncat) ! snow enthalpy in 2D               (J/m^2)

      !-----------------------------------------------------------------
      ! other variables closely related to the state of the ice
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter ::
     &   saltmax = 3.2_dbl_kind   ! max salinity, at ice base (ppm)
     &,  salnew = saltmax         ! salinity of new ice

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) :: 
     &   strength  ! ice strength                        (N/m)

c=======================================================================

      end module ice_state

c=======================================================================
