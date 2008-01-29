!=======================================================================
!BOP
!
! !MODULE: ice_state - primary state variables
!
! !DESCRIPTION:
!
! Primary state variables in various configurations
! Note: other state variables are at the end of this...
! The primary state variable names are:
!-------------------------------------------------------------------
! for each category   aggregated over     units
!                       categories
!-------------------------------------------------------------------
! aicen(i,j,n)         aice(i,j)           ---
! vicen(i,j,n)         vice(i,j)           m
! vsnon(i,j,n)         vsno(i,j)           m
! eicen(i,j,k)         eice(i,j)           J/m^2
! esnon(i,j,k)         esno(i,j)           J/m^2
! trcrn(i,j,it,n)      trcr(i,j,it)        
!
! Area is dimensionless because aice is the fractional area
! (normalized so that the sum over all categories, including open
! water, is 1.0).  That is why vice/vsno have units of m instead of
! m^3, and eice/esno have units of J/m^2 instead of J.
!
! Variable names follow these rules:
!
! (1) For 3D variables (indices i,j,n), write 'ice' or 'sno' or
!     'sfc' and put an 'n' at the end.
! (2) For 2D variables (indices i,j) aggregated over all categories,
!     write 'ice' or 'sno' or 'sfc' without the 'n'.
! (3) For 2D variables (indices i,j) associated with an individual
!     category, write 'i' or 's' instead of 'ice' or 'sno' and put an 'n'
!     at the end: e.g. hin, hsn.  These are not declared here
!     but in individual modules (e.g., ice_therm_vertical).
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors C. M. Bitz, UW
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free form source (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_state
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size
      use ice_blocks
!
!EOP
!
      implicit none
      save

      !-----------------------------------------------------------------
      ! state of the ice aggregated over all categories
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         aice  , & ! concentration of ice
         vice  , & ! volume per unit area of ice          (m)
         vsno  , & ! volume per unit area of snow         (m)
         eice  , & ! energy of melt. of ice           (J/m^2)
         esno      ! energy of melt. of snow layer    (J/m^2)

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ntrcr,max_blocks) :: &
         trcr      ! ice tracers
                   ! 1: surface temperature of ice/snow (C)
                   ! 2: meltpond volume                 (m)

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks):: &
         aice0     ! concentration of open water

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks) :: &
         apondn , & ! concentration of ponds
         hpondn     ! pond depth         (m)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ntrcr,ncat,max_blocks) :: &
         trcrn     ! tracers
                   ! 1: surface temperature of ice/snow (C)

      integer (kind=int_kind), dimension (ntrcr) :: &
         trcr_depend   ! = 0 for ice area tracers
                       ! = 1 for ice volume tracers
                       ! = 2 for snow volume tracers

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ntilyr,max_blocks) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ntslyr,max_blocks) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      !-----------------------------------------------------------------
      ! indices for tracers
      ! The maximum index should be no greater than ntrcr 
      ! (ice_domain_size) to prevent array out-of-bounds errors.
      !-----------------------------------------------------------------

      integer (kind=int_kind), parameter :: &
         nt_Tsfc  =  1, & ! ice/snow surface temperature
         nt_iage  =  2, & ! volume-weighted ice age
         nt_volpn =  2    ! melt pond volume - not used, for now

      !-----------------------------------------------------------------
      ! dynamic variables closely related to the state of the ice
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         shear    , & ! strain rate II component (1/s)
         strength     ! ice strength (N/m)

      !-----------------------------------------------------------------
      ! ice state at start of time step, saved for later in the step 
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         aice_init       ! initial concentration of ice, for diagnostics

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ncat,max_blocks) :: &
         aicen_init  , & ! initial ice concentration, for linear ITD
         vicen_init      ! initial ice volume (m), for linear ITD

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: bound_state - bound calls for ice state variables
!
! !INTERFACE:
!
      subroutine bound_state (aicen, trcrn, &
                              vicen, vsnon, &
                              eicen, esnon)
!
! !DESCRIPTION:
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
      use ice_boundary
      use ice_domain
      use ice_constants
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ncat,max_blocks), intent(inout) :: &
         aicen , & ! fractional ice area
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ntrcr,ncat,max_blocks), &
         intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ntilyr,max_blocks),intent(inout) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ntslyr,max_blocks),intent(inout) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)
!
!EOP
!
         call ice_HaloUpdate (aicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (trcrn,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (eicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (esnon,            halo_info, &
                              field_loc_center, field_type_scalar)

      end subroutine bound_state

!=======================================================================

      end module ice_state

!=======================================================================
