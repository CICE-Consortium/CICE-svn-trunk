!=======================================================================
!
! Biogeochemistry variables
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_zbgc_public

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size
      !use ice_exit
      use ice_blocks, only: nx_block, ny_block

      implicit none 

      private
      public :: lateral_melt_bgc, calculate_qin_from_Sin, remap_layers_bgc_plus, &
           remap_layers_bgc_plus_xy

      logical (kind=log_kind), public :: & 
         restart_S     ,   &! if .true., read Salinity from restart file
         tr_bgc_S           ! if .true., S as product tracer on ice

      character(char_len_long), public :: & 
         bgc_data_dir       ! directory for biogeochemistry data

      logical (kind=log_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         Rayleigh_criteria, &   ! .true. means Ra_c was reached
         first_ice              ! .true. until ice is formed (false)  
                                !
      real (kind=dbl_kind), public :: & 
         Ra_c      ,  &  ! critical Rayleigh number for bottom convection
         grid_oS   ,  &  ! for bottom flux (zsalinity)
         l_skS     ,  &  ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
         lapidus_g , &   ! constant for artificial viscosity/diffusion during growth
         lapidus_m       ! constant for artificial diffusion during melt

      real (kind=dbl_kind), parameter, public :: & 
         min_salin = p1             , & ! threshold for brine pocket treatment 
         Dm        = 1.0e-9_dbl_kind, & ! molecular diffusion (m^2/s)
         thin      = 0.05_dbl_kind      ! minimum ice thickness for salinity dynamics

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

      real (kind=dbl_kind), & 
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         fsicen, &    ! category fsice(kg/m^2/s) 
         fsicen_g     ! salt flux from gravity drainage only

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbltrcr,max_blocks), public :: &
         flux_bio ,&    ! all bio fluxes (+ive to ocean) are included here 
         flux_bio_g, &   !gravity drainage contribution
         ocean_bio,  &  ! contains all the ocean bgc tracer concentrations
         flux_bio_gbm ,&  ! all bio fluxes (+ive to ocean) are included here
         flux_bio_g_gbm  !gravity drainage contribution

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         fsice   , & ! Total flux  of salt to ocean at time step for conservation
         fsice_g     ! Total gravity drainage flux
                        !(kg/m^2/s)

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
         nit_data_type    ! 'default', 'clim'

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
         restart_bgc  ,     & ! if .true., read bgc restart file
         scale_bgc            ! if .true., initialize bgc tracers proportionally with salinity

      logical (kind=log_kind), public :: &
         tr_bgc_NO,      &   !if .true. Nitrate tracer in ice 
         tr_bgc_N,       & ! if .true., nitrogen as algal tracer on ice
         tr_bgc_C,       & ! if .true., carbon as algal tracer on ice
         tr_bgc_chl,     & ! if .true., chlorophyll as algal tracer on ice
         tr_bgc_NH,      & ! if .true., ammonia/um as nutrient tracer on ice
         tr_bgc_Sil,     & ! if .true., silicon as nutrient tracer on ice
         tr_bgc_DMSPp,   & ! if .true., DMSPp as algal content tracer on ice
         tr_bgc_DMSPd,   & ! if .true., DMSPd as precursor tracer on ice
         tr_bgc_DMS,     & ! if .true., DMS as product tracer on ice
         tr_bgc_PON,     & ! if .true., PON as product tracer on ice
         restore_bgc,    & !            restore nitrate if true
         solve_bgc         ! if .true., solve chemistry portion of code

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
         grid_o    ,  &  ! for bottom flux        
         l_sk      ,  &  ! characteristic diffusive scale (zsalinity) (m)
         grid_o_t  , &   ! top grid point length scale
         initbio_frac    ! fraction of ocean tracer concentration used to initialize tracer

!    Bioparameters from algal_dynamic
     
    real (kind=dbl_kind), parameter, public :: &
         R_C2N      = 7.0_dbl_kind, & ! algal C to N (mole/mole) Kristiansen 1991 (Barents), 9.0_dbl_kind
         R_gC2molC  = 12.01_dbl_kind, & ! mg/mmol C
         R_chl2N    = 3.0_dbl_kind , &  ! algal chlorophyll to N (mg/millimole)
         fr_resp    = 0.05_dbl_kind     ! respiration fraction

      !-----------------------------------------------------------------
      ! for bgc layer model
      !-----------------------------------------------------------------
      real (kind=dbl_kind), public :: &
         rhosi         ! average sea ice density (919-974 kg/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         S_tot      , &  ! Total ice salinity in per grid cell (g/m^2)
         chl_net    , &  ! Total chla (mg chla/m^2) per grid cell
         PP_net     , &  ! Total production (mg C/m^2/s) per grid cell
         NO_net     , &  ! Total production (mg C/m^2/s) per grid cell
         hbri            ! brine height

      real (kind=dbl_kind), dimension (nblyr+2), public :: &
         bgrid          ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), public :: &
         igrid          ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), public :: &
         cgrid           !CICE vertical coordinate   

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2,ncat,max_blocks), public :: &
         zphi        , & ! Porosity of layers    
         zTin         ! Temperature on ice layers interpolated on the bio grid (oC

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1,ncat,max_blocks), public :: &
         zfswin         ! Shortwave flux into layers interpolated on bio grid  (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1,ncat,max_blocks), public :: &
         iDi        , & ! igrid Diffusivity          (m^2/s)    
         iki            ! Ice permeability            (m^2)     

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
!BOP
!
! !ROUTINE: calculate_qin_from_Sin  - calculate enthalpy from internal ice salinity
!
! !DESCRIPTION:
!
!  Compute the internal ice  enthalpy using new salinity and Tin
!
! !REVISION HISTORY:
!
! !INTERFACE:
!
      function calculate_qin_from_Sin (Tin, Tmltk) &
               result(qin)
              
!        
! !USES:
!
! !INPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         Tin                ,&  ! internal temperature
         Tmltk                  ! melting temperature at one level
!
! !OUTPUT PARAMETERS
!
     real (kind=dbl_kind) :: &
         qin                 ! melting temperature at one level
!
!EOP
!
    

         qin =-rhoi*(cp_ice*(Tmltk-Tin) + Lfresh*(c1-Tmltk/Tin) - cp_ocn*Tmltk)


      end function calculate_qin_from_Sin

!=======================================================================
!BOP
!
! !IROUTINE: remap_layers_bgc
!
! !INTERFACE:
!
      subroutine remap_layers_bgc (nx_block,ny_block,        &
                                   indxi,   indxj,           &
                                   icells,                   &
                                   ntrcr,                    &
                                   nlyrn,                    &
                                   it,                       &
                                   trcrn,    trtmp,          &
                                   nr0,      nblyr)
!
! !DESCRIPTION:
!
! Remaps tracer fields in a given category from one set of layers to another.
!
! !REVISION HISTORY:
!
! author: Elizabeth Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of ocean/ice cells
         ntrcr             , & ! number of tracers in use
         it                , & ! tracer index in top layer
         nr0               , & ! receiver category
         nlyrn             , & ! number of ice layers
         nblyr                 ! number of biology layers

      integer (kind=int_kind), dimension (icells), intent(in) :: &
         indxi             , & ! indices for i/j directions
         indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
         intent(in) ::       &
         trcrn                 ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
         intent(inout) ::    &
         trtmp                 ! temporary, remapped ice tracers
!
!EOP
!
      integer (kind=int_kind) :: &
           i, j, ij      , & ! horizontal indices
           kd, kr, kdr   , & ! more indices
           n_nd          , & ! number of layers in donor
           n_nr              ! number of layers in receiver

      real (kind=dbl_kind), dimension (nblyr*nlyrn) :: &
           trdr          , & ! trcrn over combined grid
           trr               ! trcrn on receiving grid

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (nr0 == 0) then ! cice to bio
            n_nd = nlyrn
            n_nr = nblyr
         else               ! bio to cice
            n_nd = nblyr
            n_nr = nlyrn
         endif

         if (n_nd /= n_nr) then ! remap
            do kd = 1, n_nd
               do kr = 1, n_nr
                  kdr = kr + (kd - 1) * n_nr
                  trdr(kdr) = trcrn(i,j,it+kd-1)
               enddo
            enddo
            do kr = 1, n_nr
               trr(kr) = c0
               do kd = 1, n_nd
                  kdr = kd + (kr - 1) * n_nd
                  trr(kr) = trr(kr) + trdr(kdr)
               enddo
               trtmp(i,j,it+kr-1) = trr(kr)/n_nd
            enddo
         else     ! fill trtmp with original trcrn values
            do kr = 1, n_nr
               trtmp(i,j,it+kr-1) = trcrn(i,j,it+kr-1)
            enddo
         endif
      enddo                  ! ij

      end subroutine remap_layers_bgc

!=======================================================================
!BOP
!
! !IROUTINE: remap_layers_bgc_plus
!
! !INTERFACE:
!
      subroutine remap_layers_bgc_plus (nx_block,ny_block,        &
                                   indxi,   indxj,           &
                                   icells,                   &
                                   ntrcr,                    &
                                   nlyrn,                    &
                                   it,                       &
                                   trcrn,    trtmp,          &
                                   nr0,      nblyr,      &
                                   hice, hinS, ice_grid,  &
                                   bio_grid, S_min)
!
! !DESCRIPTION:
!
! Remaps tracer fields in a given category from one set of layers to another.
! Grids can be very different and  so can  vertical spaces.  

! !REVISION HISTORY:
!
! author: Elizabeth Hunke, LANL
!
! !USES:

      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of ocean/ice cells
         ntrcr             , & ! number of tracers in use
         it                , & ! tracer index in top layer
         nr0               , & ! receiver category
         nlyrn             , & ! number of ice layers
         nblyr                 ! number of biology layers

      integer (kind=int_kind), dimension (icells), intent(in) :: &
         indxi             , & ! indices for i/j directions
         indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
         intent(in) ::       &
         trcrn                 ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
         intent(inout) ::    &
         trtmp                 ! temporary, remapped ice tracers

      real (kind=dbl_kind), dimension (nlyrn), intent(in) :: &
         ice_grid             ! CICE grid  cgrid(2:nilyr+1)

      real (kind=dbl_kind), dimension (nblyr), intent(in) :: &
         bio_grid             ! CICE grid  grid(2:nblyr+1)

      real(kind=dbl_kind), dimension (icells), intent(in) :: &
         hice               ,& ! CICE ice thickness
         hinS               , & ! Brine height 
         S_min        ! for salinity on CICE grid        

!
!EOP
!
      integer (kind=int_kind) :: &
           i, j, ij      , & ! horizontal indices
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

    
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if ((hinS(ij) < c0) .OR. (hice(ij) < c0)) then
                  write(nu_diag, *)'Problem in remap_layers_bgc_plus'
                  write(nu_diag, *) '(hinS(ij) < c0) .OR. (hice(ij) < c0)'
                  write(nu_diag, *) 'hinS(ij),hice(ij)',hinS(ij),hice(ij)
                  call abort_ice ('ice: remap_layers_bgc_plus error')
         endif

         
         if (nr0 == 0) then ! cice to bio
            n_nd = nlyrn
            n_nr = nblyr
            n_plus = 2
            dgrid(1) = min( -hice(ij)+hinS(ij), -hinS(ij)+hice(ij), c0)            
            dgrid(nlyrn+2) = min(hinS(ij),hice(ij)) 
	    tracer(1) = trcrn(i,j,it)
	    tracer(nlyrn+2) = trcrn(i,j,it+nlyrn-1)
            rgrid(nblyr+2) = min(hinS(ij),hice(ij))
            if (hice(ij) > hinS(ij)) then
              rgrid(1) = c0 

	     do kr = 1,n_nr
	      rgrid(kr+1) = bio_grid(kr)*hinS(ij);
             enddo
             do kd = 1,n_nd
	      dgrid(kd+1) = (ice_grid(kd)-c1)*hice(ij)+hinS(ij); 
              tracer(kd+1) = trcrn(i,j,it+kd-1)
             enddo
            else
              rgrid(1) = -hinS(ij) + hice(ij) 
	      do kr = 1,n_nr
	        rgrid(kr+1) = (bio_grid(kr)-c1)*hinS(ij) + hice(ij);
              enddo
              do kd = 1,n_nd
	       dgrid(kd+1) = ice_grid(kd)*hice(ij); 
               tracer(kd+1) = trcrn(i,j,it+kd-1)
              enddo
            endif
              
         else               ! bio to cice
            n_nd = nblyr
            n_nr = nlyrn
         
         ! bio to cice
           if (hice(ij) > hinS(ij)) then
            n_plus = 3
            tracer(1) = S_min(ij)
            tracer(2) = S_min(ij)
            dgrid(1) = -hice(ij)+hinS(ij)
            dgrid(2) = p5*(hinS(ij)-hice(ij))
            dgrid(nblyr+3) = hinS(ij)
            tracer(nblyr+3) = trcrn(i,j,it+nblyr-1)
            rgrid(1) = -hice(ij) + hinS(ij)
            rgrid(nlyrn+2) = hinS(ij) 

           do kd = 1,n_nd
             dgrid(kd+2) = bio_grid(kd)*hinS(ij)
             tracer(kd+2) = trcrn(i,j,it+kd-1)
           enddo
           do kr = 1,n_nr
             rgrid(kr+1) = (ice_grid(kr)-c1)*hice(ij)+ hinS(ij)
           enddo


           else
            n_plus = 2
            tracer(1) = trcrn(i,j,it)
            tracer(nblyr+2) = trcrn(i,j,it+nblyr-1)
            dgrid(1) = hice(ij)-hinS(ij)
            dgrid(nblyr+2) = hice(ij)
            rgrid(nlyrn+2) = hice(ij)
            rgrid(1) = c0

            do kd = 1,n_nd
              dgrid(kd+1) = (bio_grid(kd)-c1)*hinS(ij) + hice(ij)
              tracer(kd+1) = trcrn(i,j,it+kd-1)
            enddo
            do kr = 1,n_nr
              rgrid(kr+1) = ice_grid(kr)*hice(ij)
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
                trtmp(i,j,it+kr-1) = trdr(kdr-1) + (rgrid(kr+1) - trgrid(kdr-1)) * &
                                   (tracer(kd) - trdr(kdr-1))/(dgrid(kd) - trgrid(kdr-1))                
               ! if ((trtmp(i,j,it+kr-1) < c0) .AND. (it < nt_qice)) then
               !   write(nu_diag, *) 'i,j,it,kr, kd,:',i,j,it,kr, kd
               !   write(nu_diag, *) 'kdi,kdr,nr0:',kdi, kdr,nr0
               !   write(nu_diag, *) 'trcrn(i,j,it),S_min(ij),hice(ij),hinS(ij),nt_qice:',trcrn(i,j,it),S_min(ij),hice(ij),hinS(ij),nt_qice
               !   write(nu_diag, *) 'trdr(kdr-1),rgrid(kr+1),trgrid(kdr-1),nt_qice:',trdr(kdr-1),rgrid(kr+1),trgrid(kdr-1),nt_qice
               !   write(nu_diag, *) 'tracer(kd),dgrid(kd),trtmp(i,j,it+kr-1):',tracer(kd),dgrid(kd),trtmp(i,j,it+kr-1)
               !   call abort_ice ('ice: remap_layers_bgc_plus error')
               ! endif
                trdr(kdr) = trtmp(i,j,it+kr-1)
                EXIT
             else
                kdr = kdr+1
                kdi = kd+1
                trgrid(kdr) = rgrid(kr+1)
                trtmp(i,j,it+kr-1) = tracer(kd)              
                !if ((trtmp(i,j,it+kr-1)) < c0 .AND. (it < nt_qice)) then  
                !  write(nu_diag, *) 'kdi,kdr,nr0,n_nd,n_nr,n_plus:',kdi, kdr,nr0,n_nd,n_nr,n_plus
                !  write(nu_diag, *) 'trcrn(i,j,it),S_min(ij),hice(ij),hinS(ij),nt_qice:',trcrn(i,j,it),S_min(ij),hice(ij),hinS(ij),nt_qice
                !  write(nu_diag, *) 'tracer(kd),trtmp(i,j,it+kr-1),kd,i,j,it,kr:',tracer(kd),trtmp(i,j,it+kr-1),kd,i,j,it,kr
                !  write(nu_diag, *) 'trgrid(kdr),kdr,kdr,kdi:',trgrid(kdr),kdr,kdr,kdi
                !  call abort_ice ('ice: remap_layers_bgc_plus error')
                !endif
                trdr(kdr) = tracer(kd)
                EXIT
             endif
           enddo
         enddo
      enddo                  ! ij

      end subroutine remap_layers_bgc_plus

!=======================================================================
!BOP
!
! !IROUTINE: remap_layers_bgc_plus_xy
!
! !INTERFACE:
!
      subroutine remap_layers_bgc_plus_xy ( ntrcr,                    &
                                   nlyrn,                    &
                                   it,                       &
                                   trcrn,    trtmp,          &
                                   nr0,      nblyr,      &
                                   hice, hinS, ice_grid,  &
                                   bio_grid, S_min)
!
! !DESCRIPTION:
!
! Remaps tracer fields in a given category from one set of layers to another.
! Grids can be very different and  so can  vertical spaces.  

! !REVISION HISTORY:
!
! author: Elizabeth Hunke, LANL
!
! !USES:

      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         ntrcr             , & ! number of tracers in use
         it                , & ! tracer index in top layer
         nr0               , & ! receiver category
         nlyrn             , & ! number of ice layers
         nblyr                 ! number of biology layers


      real (kind=dbl_kind), dimension (ntrcr), &
         intent(in) ::       &
         trcrn                 ! ice tracers

      real (kind=dbl_kind), dimension (ntrcr), &
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

!
!EOP
!
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

    
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

 
         if ((hinS < c0) .OR. (hice < c0)) then
                  write(nu_diag, *)'Problem in remap_layers_bgc_plus'
                  write(nu_diag, *) '(hinS < c0) .OR. (hice < c0)'
                  write(nu_diag, *) 'hinS,hice',hinS,hice
                  call abort_ice ('ice: remap_layers_bgc_plus error')
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
                trtmp(it+kr-1) = trdr(kdr-1) + (rgrid(kr+1) - trgrid(kdr-1)) * &
                                   (tracer(kd) - trdr(kdr-1))/(dgrid(kd) - trgrid(kdr-1))                
               ! if ((trtmp(it+kr-1) < c0) .AND. (it < nt_qice)) then
               !   write(nu_diag, *) 'it,kr, kd,:',it,kr, kd
               !   write(nu_diag, *) 'kdi,kdr,nr0:',kdi, kdr,nr0
               !   write(nu_diag, *) 'trcrn(it),S_min,hice,hinS,nt_qice:',trcrn(it),S_min,hice,hinS,nt_qice
               !   write(nu_diag, *) 'trdr(kdr-1),rgrid(kr+1),trgrid(kdr-1),nt_qice:',trdr(kdr-1),rgrid(kr+1),trgrid(kdr-1),nt_qice
               !   write(nu_diag, *) 'tracer(kd),dgrid(kd),trtmp(it+kr-1):',tracer(kd),dgrid(kd),trtmp(it+kr-1)
               !   call abort_ice ('ice: remap_layers_bgc_plus error')
               ! endif
                trdr(kdr) = trtmp(it+kr-1)
                EXIT
             else
                kdr = kdr+1
                kdi = kd+1
                trgrid(kdr) = rgrid(kr+1)
                trtmp(it+kr-1) = tracer(kd)              
                !if ((trtmp(it+kr-1)) < c0 .AND. (it < nt_qice)) then  
                !  write(nu_diag, *) 'kdi,kdr,nr0,n_nd,n_nr,n_plus:',kdi, kdr,nr0,n_nd,n_nr,n_plus
                !  write(nu_diag, *) 'trcrn(it),S_min,hice,hinS,nt_qice:',trcrn(it),S_min,hice,hinS,nt_qice
                !  write(nu_diag, *) 'tracer(kd),trtmp(it+kr-1),kd,it,kr:',tracer(kd),trtmp(it+kr-1),kd,it,kr
                !  write(nu_diag, *) 'trgrid(kdr),kdr,kdr,kdi:',trgrid(kdr),kdr,kdr,kdi
                !  call abort_ice ('ice: remap_layers_bgc_plus error')
                !endif
                trdr(kdr) = tracer(kd)
                EXIT
             endif
           enddo
         enddo
     

      end subroutine remap_layers_bgc_plus_xy

!=======================================================================

! When sea ice melts laterally, flux bgc to ocean

      subroutine lateral_melt_bgc (nx_block, ny_block, &
                                   icells,   dt,       &
                                   indxi,    indxj,    &
                                   rside,    vicen,    &
                                   trcrn,    fsice,    &
                                   flux_bio, nbltrcr)

      use ice_state, only: nt_fbri, nt_bgc_S, hbrine,&
                           nt_bgc_N, nt_bgc_C, nt_bgc_chl, &
                           nt_bgc_NO, nt_bgc_NH, nt_bgc_Sil, &
                           nt_bgc_DMSPp, nt_bgc_DMSPd, nt_bgc_DMS, &
                           nt_bgc_PON

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of cells with aice > puny
         nbltrcr               ! number of biology tracers

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with aice > puny

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr), &
         intent(in) :: &
         trcrn     ! tracer array

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         rside     ! fraction of ice that melts laterally

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         fsice     ! salt flux from layer Salinity (kg/m^2/s)
  
      real (kind=dbl_kind), dimension(nx_block,ny_block,nbltrcr), &
         intent(inout) :: &
         flux_bio  ! biology tracer flux from layer bgc (mmol/m^2/s)

!local
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         k           , & ! layer index
         ij              ! horizontal index, combines i and j loops

      real (kind=dbl_kind) :: &
         zspace      ! bio grid spacing

       zspace = c1/(real(nblyr,kind=dbl_kind))

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (hbrine .AND. tr_bgc_S) then
            do k = 1,nblyr
              fsice(i,j) = fsice(i,j) + rhosi*trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*p001 *zspace*trcrn(i,j,nt_bgc_S+k-1)&
                            * rside(i,j)/dt
            enddo

            if (tr_bgc_NO ) then
            do k = 1,nblyr
              flux_bio(i,j,nlt_bgc_NO) = flux_bio(i,j,nlt_bgc_NO) + trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*zspace*trcrn(i,j,nt_bgc_NO+k-1)&
                            * rside(i,j)/dt
            enddo
            endif 
            if (tr_bgc_chl ) then
            do k = 1,nblyr
              flux_bio(i,j,nlt_bgc_chl) = flux_bio(i,j,nlt_bgc_chl) + trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*zspace*trcrn(i,j,nt_bgc_chl+k-1)&
                            * rside(i,j)/dt
            enddo
            endif 
            if (tr_bgc_NH ) then
            do k = 1,nblyr
              flux_bio(i,j,nlt_bgc_NH) = flux_bio(i,j,nlt_bgc_NH) + trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*zspace*trcrn(i,j,nt_bgc_NH+k-1)&
                            * rside(i,j)/dt
            enddo
            endif
 
            if (tr_bgc_C ) then
            do k = 1,nblyr
              flux_bio(i,j,nlt_bgc_C) = flux_bio(i,j,nlt_bgc_C) + trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*zspace*trcrn(i,j,nt_bgc_C+k-1)&
                            * rside(i,j)/dt
            enddo
            endif
 
            if (tr_bgc_Sil ) then
            do k = 1,nblyr
              flux_bio(i,j,nlt_bgc_Sil) = flux_bio(i,j,nlt_bgc_Sil) + trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*zspace*trcrn(i,j,nt_bgc_Sil+k-1)&
                            * rside(i,j)/dt
            enddo
            endif

 
            if (tr_bgc_DMSPp ) then
            do k = 1,nblyr
              flux_bio(i,j,nlt_bgc_DMSPp) = flux_bio(i,j,nlt_bgc_DMSPp) + trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*zspace*trcrn(i,j,nt_bgc_DMSPp+k-1)&
                            * rside(i,j)/dt
            enddo
            endif

 
            if (tr_bgc_DMSPd ) then
            do k = 1,nblyr
              flux_bio(i,j,nlt_bgc_DMSPd) = flux_bio(i,j,nlt_bgc_DMSPd) + trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*zspace*trcrn(i,j,nt_bgc_DMSPd+k-1)&
                            * rside(i,j)/dt
            enddo
            endif

            if (tr_bgc_DMS ) then
            do k = 1,nblyr
              flux_bio(i,j,nlt_bgc_DMS) = flux_bio(i,j,nlt_bgc_DMS) + trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*zspace*trcrn(i,j,nt_bgc_DMS+k-1)&
                            * rside(i,j)/dt
            enddo
            endif 

            if (tr_bgc_N ) then
            do k = 1,nblyr
              flux_bio(i,j,nlt_bgc_N) = flux_bio(i,j,nlt_bgc_N) + trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*zspace*trcrn(i,j,nt_bgc_N+k-1)&
                            * rside(i,j)/dt
            enddo
            endif

            if (tr_bgc_PON ) then
            do k = 1,nblyr
              flux_bio(i,j,nlt_bgc_PON) = flux_bio(i,j,nlt_bgc_PON) + trcrn(i,j,nt_fbri)&
                            *vicen(i,j)*zspace*trcrn(i,j,nt_bgc_PON+k-1)&
                            * rside(i,j)/dt
            enddo
            endif
 
            endif  ! hbrine and tr_bgc_S
         enddo                  ! ij

      end subroutine lateral_melt_bgc 

!=======================================================================

      end module ice_zbgc_public

!=======================================================================
