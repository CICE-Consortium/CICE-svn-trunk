!=======================================================================
!
! Biogeochemistry variables
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_zbgc_shared

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size
      !use ice_exit
      use ice_blocks, only: nx_block, ny_block

      implicit none 

      private
      public :: remap_layers_bgc

      logical (kind=log_kind), public :: & 
         restart_hbrine     ! if .true., read hbrine from restart file

      character(char_len_long), public :: & 
         bgc_data_dir       ! directory for biogeochemistry data

      logical (kind=log_kind), dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         first_ice              ! .true. until ice is formed (false)
     
      real (kind=dbl_kind), public :: & 
         phi_snow  ,  &  ! porosity of snow
         flood_frac      ! fraction of ocean/meltwater that floods

      real (kind=dbl_kind), parameter, public :: & 
         min_salin      = p1          , & ! threshold for brine pocket treatment 
         his_min        = p01            , & ! minimum hbrine thickness
         thinS          = 0.05_dbl_kind      ! minimum ice thickness for salinity dynamics

      real (kind=dbl_kind), parameter, public :: & 
         k_o        = 3.e-8_dbl_kind      ! scaling factor permeability calculation (m^2)

      integer (kind=int_kind), parameter, public :: &
         exp_h    = 3                 ! power law for hierarchical model  

      real (kind=dbl_kind), & 
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         dh_top    ,& ! brine top change
         dh_bot    ,& ! brine bottom change
         dhi_top   ,& ! ice top change
         dhi_bot   ,& ! ice bottom change
         sice_rho     ! avg sea ice density  (kg/m^3)  ! ech: diagnostic only?

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_nbtrcr,max_blocks), public :: &
         flux_bio ,&    ! all bio fluxes to ocean
         ocean_bio,  &  ! contains all the ocean bgc tracer concentrations
         flux_bio_ai    ! all bio fluxes to ocean, averaged over grid cell

      real (kind=int_kind), dimension(max_nbtrcr), public :: &
         bgc_tracer_type  ! 1  dissolved tracers: mix like salinity
                          ! 0  tracers that cling: resist brine motion (algae)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         nit     , & ! ocean nitrate (mmol/m^3)          
         amm     , & ! ammonia/um (mmol/m^3)
         sil     , & ! silicate (mmol/m^3)
         dmsp    , & ! dmsp (mmol/m^3)
         dms     , & ! dms (mmmol/m^3     
         algalN      ! ocean algal nitrogen (mmol/m^3)

      character (char_len_long), public :: &        ! input data file names
         nit_file, &! nitrate input file
         sil_file   ! silicate input file

      character(char_len), public :: &          
         sil_data_type, & ! 'default', 'clim'
         nit_data_type, & ! 'default', 'clim'     
         bgc_flux_type    ! 'constant', 'Jin2006' describes type of ocean-ice piston velocity 

      logical (kind=log_kind), public :: & 
         tr_bgc_N_sk,       & ! if .true., nitrogen as algal tracer on ice
         tr_bgc_C_sk,       & ! if .true., carbon as algal tracer on ice
         tr_bgc_chl_sk,     & ! if .true., chlorophyll as algal tracer on ice
         tr_bgc_Nit_sk,     & ! if .true., nitrate as nutrient tracer on ice
         tr_bgc_Am_sk,      & ! if .true., ammonia/um as nutrient tracer on ice
         tr_bgc_Sil_sk,     & ! if .true., silicon as nutrient tracer on ice
         tr_bgc_DMSPp_sk,   & ! if .true., DMSPp as algal content tracer on ice
         tr_bgc_DMSPd_sk,   & ! if .true., DMSPd as precursor tracer on ice
         tr_bgc_DMS_sk,     & ! if .true., DMS as product tracer on ice
         restart_bgc,       & ! if .true., read bgc restart file
         restore_bgc,       & !            restore nitrate if true
         solve_skl_bgc        ! if .true., solve skeletal biochemistry portion of code

      ! zbgc tracer indices
      integer (kind=int_kind), public :: &
         nlt_bgc_N    ,   & ! algae 
         nlt_bgc_C    ,   & ! 
         nlt_bgc_chl  ,   & ! 
         nlt_bgc_NO   ,   & ! nutrients  
         nlt_bgc_NH   ,   & ! 
         nlt_bgc_Sil  ,   & !
         nlt_bgc_DMSPp ,   & ! trace gases (skeletal layer)
         nlt_bgc_DMSPd ,   & ! 
         nlt_bgc_DMS   ,   & ! 
         nlt_bgc_PON         ! zooplankton and detritus
!
!  Bio-variables from Diffuse_biology  (in ice_init)
!
      real (kind=dbl_kind), public :: & 
         initbio_frac    ! fraction of ocean tracer concentration used to initialize tracer

!    Bioparameters from algal_dynamic
     
    real (kind=dbl_kind), parameter, public :: &
         R_C2N      = 7.0_dbl_kind  , & ! algal C to N (mole/mole) Kristiansen 1991 (Barents), 9.0_dbl_kind
         R_gC2molC  = 12.01_dbl_kind, & ! mg/mmol C
         R_chl2N    = 3.0_dbl_kind  , & ! algal chlorophyll to N (mg/millimole)
         R_S2N      = 0.03_dbl_kind , & ! algal S to N (mole/mole)
         fr_resp    = 0.05_dbl_kind     ! respiration fraction

    real (kind=dbl_kind), parameter, public :: &
         sk_l       = 0.03_dbl_kind,  & ! skeletal layer thickness (m)
         phi_sk     = 0.30_dbl_kind     ! skeletal layer porosity

      !-----------------------------------------------------------------
      ! for bgc layer model
      !-----------------------------------------------------------------
      real (kind=dbl_kind), parameter, public :: &
         rhosi     = 940.0_dbl_kind   ! average sea ice density (919-974 kg/m^2, &
                                      ! Cox and Weeks, 1982)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         chl_net    , &  ! Total chla (mg chla/m^2) per grid cell
         grow_net   , &  ! Specific growth rate (/s) per grid cell
         PP_net     , &  ! Total production (mg C/m^2/s) per grid cell
         NO_net     , &  ! Total production (mg C/m^2/s) per grid cell
         hbri            ! brine height

      real (kind=dbl_kind), dimension (nblyr+2), public :: &
         bgrid          ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), public :: &
         igrid          ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), public :: &
         cgrid           ! CICE vertical coordinate   

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2,ncat,max_blocks), public :: &
         zphi        , & ! Porosity of layers    
         zTin         ! Temperature on ice layers interpolated on the bio grid (oC

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr,max_blocks), public :: &
         upNO    , & ! nitrate uptake rate (mmol/m^3/s)
         upNH        ! ammonium uptake rate (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr,ncat,max_blocks), public :: &
         growN            ! algal growth rate         (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr,max_blocks), public :: &
         growNp            ! algal growth rate         (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         darcy_V          ! darcy velocity positive up        (m/s)
    
!=======================================================================

      contains

!=======================================================================
!
! Remaps tracer fields in a given category from one set of layers to another.
! Grids can be very different and  so can  vertical spaces.  

      subroutine remap_layers_bgc ( ntrcr,                    &
                                   nlyrn,                    &
                                   it,                       &
                                   trcrn,    trtmp,          &
                                   nr0,      nblyr,      &
                                   hice, hinS, ice_grid,  &
                                   bio_grid, S_min)

      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice

      integer (kind=int_kind), intent(in) :: &
         ntrcr             , & ! number of tracers in use
         it                , & ! tracer index in top layer
         nr0               , & ! receiver category
         nlyrn             , & ! number of ice layers
         nblyr                 ! number of biology layers

      real (kind=dbl_kind), dimension (ntrcr), &
         intent(in) ::       &
         trcrn                 ! ice tracers

      real (kind=dbl_kind), dimension (nblyr+2), &
         intent(inout) ::    &
         trtmp                 ! temporary, remapped ice tracers

      real (kind=dbl_kind), dimension (nlyrn), intent(in) :: &
         ice_grid             ! CICE grid  cgrid(2:nilyr+1)

      real (kind=dbl_kind), dimension (nblyr), intent(in) :: &
         bio_grid             ! CICE grid  grid(2:nblyr+1)

      real(kind=dbl_kind), intent(in) :: &
         hice               ,& ! CICE ice thickness
         hinS               , & ! Brine height 
         S_min        ! for salinity on CICE grid        

      ! local variables

      integer (kind=int_kind) :: &
           kd, kr, kdr   , & ! more indices
           kdi           , & ! more indices
           n_nd          , & ! number of layers in donor
           n_nr , n_plus     ! number of layers in receiver

      real (kind=dbl_kind), dimension (nblyr+3+nlyrn) :: &
           trdr          , & ! combined tracer 
           trgrid            ! combined grid 

      real (kind=dbl_kind), dimension (nblyr+nilyr+3) :: &
           tracer         ,&  ! temporary,  ice tracers values
           dgrid           ,& ! temporary,  donor grid dimensional
           rgrid              !  temporary, receiver grid dimensional

         if ((hinS < c0) .OR. (hice < c0)) then
                  write(nu_diag, *)'Problem in remap_layers_bgc'
                  write(nu_diag, *) '(hinS < c0) .OR. (hice < c0)'
                  write(nu_diag, *) 'hinS,hice',hinS,hice
                  call abort_ice ('ice: remap_layers_bgc error')
         endif
         
         if (nr0 == 0) then ! cice to bio
            n_nd = nlyrn
            n_nr = nblyr
            n_plus = 2
            dgrid(1) = min( -hice+hinS, -hinS+hice, c0)            
            dgrid(nlyrn+2) = min(hinS,hice) 
	    tracer(1) = trcrn(it)
	    tracer(nlyrn+2) = trcrn(it+nlyrn-1)
            rgrid(nblyr+2) = min(hinS,hice)
            if (hice > hinS) then
              rgrid(1) = c0 

	     do kr = 1,n_nr
	      rgrid(kr+1) = bio_grid(kr)*hinS;
             enddo
             do kd = 1,n_nd
	      dgrid(kd+1) = (ice_grid(kd)-c1)*hice+hinS; 
              tracer(kd+1) = trcrn(it+kd-1)
             enddo
            else
              rgrid(1) = -hinS + hice 
	      do kr = 1,n_nr
	        rgrid(kr+1) = (bio_grid(kr)-c1)*hinS + hice;
              enddo
              do kd = 1,n_nd
	       dgrid(kd+1) = ice_grid(kd)*hice; 
               tracer(kd+1) = trcrn(it+kd-1)
              enddo
            endif
              
         else               ! bio to cice
            n_nd = nblyr
            n_nr = nlyrn
         
         ! bio to cice
           if (hice > hinS) then
            n_plus = 3
            tracer(1) = S_min
            tracer(2) = S_min
            dgrid(1) = -hice+hinS
            dgrid(2) = p5*(hinS-hice)
            dgrid(nblyr+3) = hinS
            tracer(nblyr+3) = trcrn(it+nblyr-1)
            rgrid(1) = -hice + hinS
            rgrid(nlyrn+2) = hinS 

           do kd = 1,n_nd
             dgrid(kd+2) = bio_grid(kd)*hinS
             tracer(kd+2) = trcrn(it+kd-1)
           enddo
           do kr = 1,n_nr
             rgrid(kr+1) = (ice_grid(kr)-c1)*hice+ hinS
           enddo


           else
            n_plus = 2
            tracer(1) = trcrn(it)
            tracer(nblyr+2) = trcrn(it+nblyr-1)
            dgrid(1) = hice-hinS
            dgrid(nblyr+2) = hice
            rgrid(nlyrn+2) = hice
            rgrid(1) = c0

            do kd = 1,n_nd
              dgrid(kd+1) = (bio_grid(kd)-c1)*hinS + hice
              tracer(kd+1) = trcrn(it+kd-1)
            enddo
            do kr = 1,n_nr
              rgrid(kr+1) = ice_grid(kr)*hice
            enddo

           endif

         endif
         kdr = 0  !combined indice
         kdi = 1  
        

         do kr = 1,n_nr
           do kd = kdi,n_nd+n_plus
             if (dgrid(kd) < rgrid(kr+1))then
                kdr = kdr+1
                trgrid(kdr) = dgrid(kd)
                trdr(kdr) = tracer(kd)
             elseif (dgrid(kd) > rgrid(kr+1)) then
                kdr = kdr + 1
                kdi = kd
                trgrid(kdr) = rgrid(kr+1)
                trtmp(kr) = trdr(kdr-1) + (rgrid(kr+1) - trgrid(kdr-1)) * &
                                   (tracer(kd) - trdr(kdr-1))/(dgrid(kd) - trgrid(kdr-1))
                trdr(kdr) = trtmp(kr)
!                EXIT
             else
                kdr = kdr+1
                kdi = kd+1
                trgrid(kdr) = rgrid(kr+1)
                trtmp(kr) = tracer(kd)              
                trdr(kdr) = tracer(kd)
!                EXIT
             endif
           enddo
         enddo

      end subroutine remap_layers_bgc

!=======================================================================

      end module ice_zbgc_shared

!=======================================================================
