!  SVN:$Id$
!=======================================================================
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
! trcrn(i,j,it,n)      trcr(i,j,it)        
!
! Area is dimensionless because aice is the fractional area
! (normalized so that the sum over all categories, including open
! water, is 1.0).  That is why vice/vsno have units of m instead of m^3.
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
! authors C. M. Bitz, UW
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free form source (F90) by Elizabeth Hunke

      module ice_state

      use ice_kinds_mod
      use ice_domain_size, only: max_blocks, ncat, max_ntrcr, n_aero, &
          max_algae, max_doc, max_dic, max_aero, max_don, max_fe
      use ice_blocks, only: nx_block, ny_block

      implicit none
      private
      public :: bound_state
      save

      !-----------------------------------------------------------------
      ! state of the ice aggregated over all categories
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
         public :: &
         aice  , & ! concentration of ice
         vice  , & ! volume per unit area of ice          (m)
         vsno      ! volume per unit area of snow         (m)

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,max_ntrcr,max_blocks), public :: &
         trcr      ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), &
         public:: &
         aice0     ! concentration of open water

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), public, &
         dimension (nx_block,ny_block,max_ntrcr,ncat,max_blocks) :: &
         trcrn     ! tracers
                   ! 1: surface temperature of ice/snow (C)

      !-----------------------------------------------------------------
      ! indices and flags for tracers
      !-----------------------------------------------------------------

      integer (kind=int_kind), dimension (max_ntrcr), public :: &
         trcr_depend   ! = 0 for ice area tracers
                       ! = 1 for ice volume tracers
                       ! = 2 for snow volume tracers

      integer (kind=int_kind), public :: &
         ntrcr     ! number of tracers in use

      integer (kind=int_kind), public :: &
         nbtrcr, &    ! number of bgc tracers in use
         nbtrcr_sw    ! number of bgc tracers which impact shortwave
  
      integer (kind=int_kind), dimension(max_algae), public :: &  
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_C , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_chl   ! diatoms, phaeocystis, pico/small 

      integer (kind=int_kind), dimension(max_doc), public :: &  
         nt_bgc_DOC      !  dissolved organic carbon

      integer (kind=int_kind), dimension(max_don), public :: & 
         nt_bgc_DON         !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(max_dic), public :: &  
         nt_bgc_DIC         !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(max_fe), public :: & 
         nt_bgc_Fed,     & !  dissolved iron
         nt_bgc_Fep        !  particulate iron

      integer (kind=int_kind), dimension(max_aero), public :: &  
         nt_zaero       !  black carbon and other aerosols

      integer (kind=int_kind), public :: &  
         ntrace_start       ! index of first bio tracer
         integer (kind=int_kind), public :: &
         nt_Tsfc  , & ! ice/snow temperature
         nt_qice  , & ! volume-weighted ice enthalpy (in layers)
         nt_qsno  , & ! volume-weighted snow enthalpy (in layers)
         nt_sice  , & ! volume-weighted ice bulk salinity (CICE grid layers)
         nt_fbri  , & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
         nt_iage  , & ! volume-weighted ice age
         nt_FY    , & ! area-weighted first-year ice area
         nt_alvl  , & ! level ice area fraction
         nt_vlvl  , & ! level ice volume fraction
         nt_apnd  , & ! melt pond area fraction
         nt_hpnd  , & ! melt pond depth
         nt_ipnd  , & ! melt pond refrozen lid thickness
         nt_aero  , & ! starting index for aerosols in ice
         nt_bgc_Nit,   & ! nutrients  
         nt_bgc_Am,    & ! 
         nt_bgc_Sil,   & !
         nt_bgc_DMSPp, & ! trace gases (skeletal layer)
         nt_bgc_DMSPd, & ! 
         nt_bgc_DMS,   & ! 
         nt_bgc_PON,   & ! zooplankton and detritus  
         nt_zbgc_frac, & ! fraction of tracer in the mobile phase
         nt_bgc_S        ! Bulk salinity in fraction ice with dynamic salinity (Bio grid)

      logical (kind=log_kind), public :: &
         tr_iage,   & ! if .true., use age tracer
         tr_FY,     & ! if .true., use first-year area tracer
         tr_lvl,    & ! if .true., use level ice tracer
         tr_pond,   & ! if .true., use melt pond tracer
         tr_pond_cesm,& ! if .true., use cesm pond tracer
         tr_pond_lvl, & ! if .true., use level-ice pond tracer
         tr_pond_topo,& ! if .true., use explicit topography-based ponds
         tr_aero     ,& ! if .true., use aerosol tracers
         tr_brine       ! if .true., brine height differs from ice thickness

      !-----------------------------------------------------------------
      ! dynamic variables closely related to the state of the ice
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
         public :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         shear    , & ! strain rate II component (1/s)
         strength     ! ice strength (N/m)

      !-----------------------------------------------------------------
      ! ice state at start of time step, saved for later in the step 
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
         public :: &
         aice_init       ! initial concentration of ice, for diagnostics

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ncat,max_blocks), public :: &
         aicen_init  , & ! initial ice concentration, for linear ITD
         vicen_init  , & ! initial ice volume (m), for linear ITD
         vsnon_init      ! initial snow volume (m), for aerosol

!=======================================================================

      contains

!=======================================================================
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! author: William H. Lipscomb, LANL

      subroutine bound_state (aicen, trcrn, &
                              vicen, vsnon)

      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy
      use ice_domain, only: halo_info, maskhalo_bound, nblocks
      use ice_constants, only: field_loc_center, field_type_scalar, c0

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,ncat,max_blocks), intent(inout) :: &
         aicen , & ! fractional ice area
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,max_ntrcr,ncat,max_blocks), &
         intent(inout) :: &
         trcrn     ! ice tracers

      ! local variables

      integer (kind=int_kind) :: i, j, n, iblk

      integer (kind=int_kind), &
         dimension(nx_block,ny_block,max_blocks) :: halomask

      type (ice_halo) :: halo_info_aicemask

      call ice_HaloUpdate (aicen,            halo_info, &
                           field_loc_center, field_type_scalar)

      if (maskhalo_bound) then
         halomask(:,:,:) = 0

         !$OMP PARALLEL DO PRIVATE(iblk,n,i,j)
         do iblk = 1, nblocks
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (aicen(i,j,n,iblk) > c0) halomask(i,j,iblk) = 1
         enddo
         enddo
         enddo
         enddo
         !$OMP END PARALLEL DO

         call ice_HaloMask(halo_info_aicemask, halo_info, halomask)

         call ice_HaloUpdate (trcrn(:,:,1:ntrcr,:,:), halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloDestroy(halo_info_aicemask)

      else
         call ice_HaloUpdate (trcrn(:,:,1:ntrcr,:,:), halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info, &
                              field_loc_center, field_type_scalar)
      endif

      end subroutine bound_state

!=======================================================================

      end module ice_state

!=======================================================================
