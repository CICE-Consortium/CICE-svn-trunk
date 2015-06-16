!  SVN:$Id$
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
      use ice_constants, only: p01, p1, p5, c0, c1
      use ice_domain_size, only: ncat, max_blocks, max_nbtrcr, &
                                 nblyr, nilyr, max_algae, max_doc, &
                                 max_dic, max_aero, max_don, max_nsw, &
                                 max_fe
      use ice_blocks, only: nx_block, ny_block

      implicit none 

      private
      public ::lateral_melt_bgc, calculate_qin_from_Sin, remap_zbgc, &
            zap_small_bgc  

      logical (kind=log_kind), public :: & 
         restart_S     ,   &! if .true., read Salinity from restart file 
         restart_hbrine     ! if .true., read tr_brine from restart file

      character(char_len_long), public :: & 
         bgc_data_dir       ! directory for biogeochemistry data

      logical (kind=log_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         first_ice      ! distinguishes ice that disappears (e.g. melts)
                        ! and reappears (e.g. transport) in a grid cell
                        ! during a single time step from ice that was
                        ! there the entire time step (true until ice forms)
       
      ! coupling fluxes
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_nbtrcr,max_blocks), public :: &
         flux_bio   , & ! all bio fluxes to ocean
         ocean_bio  , & ! contains all the ocean bgc tracer concentrations
         flux_bio_atm,& ! all bio fluxes to ice from atmosphere
         atm_bio,     & ! contains all atmospheric bgc/aerosol tracer concentrations
         flux_bio_ai    ! all bio fluxes (+ive to ocean) are included here
     
      ! diagnostic fluxes
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_nbtrcr,max_blocks), public :: &
         fbio_snoice, & ! fluxes from snow to ice
         fbio_atmice    ! fluxes from atm to ice

      !-----------------------------------------------------------------
      ! general biogeochemistry
      !-----------------------------------------------------------------

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,max_nsw, ncat, max_blocks), public :: &
         trcrn_sw        ! bgc tracers active in the delta-Eddington shortwave 
                         ! calculation on the shortwave grid (swgrid)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_nbtrcr, max_blocks), public :: &
         ice_bio_net  , &   ! depth integrated tracer (mmol/m^2) 
         snow_bio_net       ! depth integrated snow tracer (mmol/m^2)

      real (kind=int_kind), dimension(max_nbtrcr), public :: &
         zbgc_frac_init,&! initializes mobile fraction
         bgc_tracer_type ! described tracer in mobile or stationary phases      
                         ! < 0 is purely mobile (eg. nitrate)
                         ! > 0 has timescales for transitions between 
                         ! phases based on whether the ice is melting or growing

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_nbtrcr, max_blocks), public :: &
         ocean_bio_all, &   ! fixed order, all values even for tracers false
         atm_bio_all        ! N(1:max_algae) = 1:max_algae
                            ! Nit = max_algae + 1
                            ! DOC(1:max_doc) = max_algae + 2: max_algae + max_doc + 1
                            ! DIC(1:max_dic) = max_algae + max_doc + 2: max_algae + max_doc + 1 + max_dic
                            ! chl(1:max_algae) =  max_algae + max_doc + 2 + max_dic: &
                            !                     2*max_algae + max_doc + 1 + max_dic
                            ! Am =  2*max_algae + max_doc + 2 + max_dic
                            ! Sil=  2*max_algae + max_doc + 3 + max_dic
                            ! DMSPp=  2*max_algae + max_doc + 4 + max_dic
                            ! DMSPd=  2*max_algae + max_doc + 5 + max_dic
                            ! DMS  =  2*max_algae + max_doc + 6 + max_dic
                            ! PON  =  2*max_algae + max_doc + 7 + max_dic
                            ! DON(1:max_don)  =  2*max_algae + max_doc + 8 + max_dic:
                            !                    2*max_algae + max_doc + 7 + max_dic + max_don
                            ! Fed(1:max_fe) = 2*max_algae + max_doc + 8 + max_dic + max_don:
                            !                2*max_algae + max_doc + 7 + max_dic + max_don + max_fe
                            ! Fep(1:max_fe) = 2*max_algae + max_doc + 8 + max_dic + max_don + max_fe:
                            !                2*max_algae + max_doc + 7 + max_dic + max_don + 2*max_fe
                            ! zaero(1:max_aero) = 2*max_algae + max_doc + 8 + max_dic + max_don + 2*max_fe: 
                            !                     2*max_algae + max_doc + 7 + max_dic + max_don + 2*max_fe
                            !                     + max_aero

      integer (kind=int_kind), dimension(nx_block, ny_block,max_algae, max_blocks), public :: &
        algal_peak          ! vertical location of algal maximum, 0 if no maximum 

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         nit        , & ! ocean nitrate (mmol/m^3)          
         amm        , & ! ammonia/um (mmol/m^3)
         sil        , & ! silicate (mmol/m^3)
         dmsp       , & ! dmsp (mmol/m^3)
         dms            ! dms (mmol/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_algae, max_blocks), public :: &
         algalN         ! ocean algal nitrogen (mmol/m^3) (diatoms, phaeo, pico)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_doc, max_blocks), public :: &
         doc             ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_don, max_blocks), public :: &
         don             ! ocean don (mmol/m^3) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_dic, max_blocks), public :: &
         dic             ! ocean dic (mmol/m^3) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_fe, max_blocks), public :: &
         fed, fep        ! ocean disolved and particulate fe (nM) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_aero, max_blocks), public :: &
         zaeros          ! ocean aerosols (mmol/m^3) 

      character (char_len_long), public :: &        ! input data file names
         nit_file   , & ! nitrate input file
         sil_file       ! silicate input file

      character(char_len), public :: &          
         sil_data_type  , & ! 'default', 'clim'
         nit_data_type  , & ! 'default', 'clim'   
         fe_data_type   , & ! 'default', 'clim'     
         bgc_flux_type      ! type of ocean-ice piston velocity 
                            ! 'constant', 'Jin2006' 

      ! ocean sources/sinks
      integer (kind=int_kind), dimension(max_algae), public :: &
         nlt_bgc_N      , & ! algae 
         nlt_bgc_C      , & ! 
         nlt_bgc_chl  

      integer (kind=int_kind), dimension(max_doc), public :: &
         nlt_bgc_DOC        ! disolved organic carbon

      integer (kind=int_kind), dimension(max_don), public :: &
         nlt_bgc_DON        !

      integer (kind=int_kind), dimension(max_dic), public :: &
         nlt_bgc_DIC        ! disolved inorganic carbon

      integer (kind=int_kind), dimension(max_fe), public :: &
         nlt_bgc_Fed    , & !
         nlt_bgc_Fep        !

      integer (kind=int_kind), dimension(max_aero), public :: &
         nlt_zaero      , & ! non-reacting layer aerosols
         nlt_zaero_sw       ! points to aerosol in trcrn_sw

      integer (kind=int_kind), public :: &
         nlt_bgc_Nit   ,   & ! nutrients  
         nlt_bgc_Am    ,   & ! 
         nlt_bgc_Sil   ,   & !
         nlt_bgc_DMSPp ,   & ! trace gases (skeletal layer)
         nlt_bgc_DMSPd ,   & ! 
         nlt_bgc_DMS   ,   & ! 
         nlt_bgc_PON   ,   & ! zooplankton and detritus
         nlt_chl_sw          ! points to total chla in trcrn_sw

      integer (kind=int_kind), dimension(max_nbtrcr), public :: &
         bio_index, &        ! relates bio indices, ie.  nlt_bgc_N to nt_bgc_N or nt_bgc_N_sk
         bio_index_o         ! relates nlt_bgc_NO or nlt_bgc_nit_sk to ocean concentration index
                             ! see ocean_bio_all

      real (kind=dbl_kind), dimension(max_nbtrcr), public :: & 
         zbgc_init_frac, &   ! fraction of ocean tracer  concentration in new ice
         tau_ret,        &   ! retention timescale  (s), mobile to stationary phase
         tau_rel             ! release timescale    (s), stationary to mobile phase          

      ! bio parameters for algal_dyn

      real (kind=dbl_kind), parameter, dimension(max_algae),  public :: &
         R_C2N     = (/ 7.0_dbl_kind, 7.0_dbl_kind, 7.0_dbl_kind /),         &   ! algal C to N (mole/mole) 
                                       ! Kristiansen 1991 (Barents) 9.0
         R_Si2N    = (/ 1.4_dbl_kind, c0, c0/),                              &   ! algal C to Sil (mole/mole) 
                                       ! Kristiansen 1991 (Barents) 9.0
         R_S2N     = (/ 0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind/),       & ! algal S to N (mole/mole)
         R_chl2N   = (/ 2.1_dbl_kind,  1.1_dbl_kind,  0.84_dbl_kind/) ,      & ! algal chlorophyll to N (mg/mmol)
                 !Marchetti et al 2006,   3 umol Fe/mol C for iron limited Pseudo-nitzschia
         R_Fe2C    = (/ 3.3e-3_dbl_kind, 3.3e-3_dbl_kind, 1.0e-1_dbl_kind/), & ! algal Fe to carbon (umol/mmol)
         R_Fe2N    = (/ 2.3e-2_dbl_kind, 2.3e-2_dbl_kind, 7.0e-1_dbl_kind/), & ! algal Fe to N (umol/mmol)
         F_abs_chl = (/ 2.0_dbl_kind, 4.0_dbl_kind, 5.0_dbl_kind /)            ! to scale absorption in Dedd

      integer (kind=int_kind), dimension(max_don), public :: & 
         R_Fe2DON  = (/ 2.3e-2_dbl_kind/)   ! Fe to N of DON (nmol/umol)

      integer (kind=int_kind), dimension(max_doc), public :: &  ! increase compare to algal R_Fe2C
         R_Fe2DOC  =  (/ 1.0e-1_dbl_kind, 3.3e-2_dbl_kind, 1.0e-1_dbl_kind/)! Fe to C of DOC (nmol/umol)

      real (kind=dbl_kind), parameter, public :: &
         R_gC2molC  = 12.01_dbl_kind,  & ! mg/mmol C
         fr_resp    = 0.05_dbl_kind ,  & ! fraction of algal growth lost due to respiration
         tau_min    = 6.24e4_dbl_kind, & ! (3.12e4_dbl_kind = 6 hours (s) (or 1.25e5_dbl_kind s = 1 day) 
                                         ! rapid mobile to stationary exchanges 
         tau_max    = 6.25e5_dbl_kind   ! 6.25e5_dbl_kind = 5 days (s) 
                                         ! long time mobile to stationary exchanges 
                     ! scavenging coefficient for tracers in snow
                     ! bottom 6 are from Flanner et al., 2007
      real (kind=dbl_kind), parameter, dimension(max_nbtrcr),  public :: &
         kscavz    = (/ 0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, &
                        0.03_dbl_kind, 0.20_dbl_kind, 0.02_dbl_kind, &
                        0.02_dbl_kind, 0.01_dbl_kind, 0.01_dbl_kind/)
      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      logical (kind=log_kind), public :: & 
         tr_bgc_Nit,     & !if .true. Nitrate tracer in ice 
         tr_bgc_N,       & ! if .true., algal nitrogen tracers on ice (n_algae)
         tr_bgc_DON,     & ! if .true., DON pools are tracers on ice (n_don)
         tr_bgc_C,       & ! if .true., algal carbon tracers + DOC and DIC on ice
         tr_bgc_chl,     & ! if .true., algal chlorophyll tracers on ice
         tr_bgc_Am,      & ! if .true., ammonia/um as nutrient tracer on ice
         tr_bgc_Sil,     & ! if .true., silicon as nutrient tracer on ice
         tr_bgc_DMS,     & ! if .true., DMS as product tracer on ice
         tr_bgc_Fe,      & ! if .true., Fe as product tracer on ice
         tr_bgc_PON        ! if .true., PON as product tracer on ice

      !-----------------------------------------------------------------
      ! skeletal layer biogeochemistry
      !-----------------------------------------------------------------

      logical (kind=log_kind), public :: & 
         restart_bgc,     & ! if true, read bgc restart file
         restore_bgc,     & ! if true, restore nitrate
         skl_bgc            ! if true, solve skeletal biochemistry

      integer (kind=int_kind), public :: &
         ntrace_start      ! index for first biological tracer in trcrn

    real (kind=dbl_kind), parameter, public :: &
         sk_l       = 0.03_dbl_kind,  & ! skeletal layer thickness (m)
         phi_sk     = 0.30_dbl_kind     ! skeletal layer porosity

      !-----------------------------------------------------------------
      ! brine
      !-----------------------------------------------------------------

      integer (kind=int_kind), parameter, public :: &
         exp_h     = 3              ! power law for hierarchical model  

      real (kind=dbl_kind), parameter, public :: & 
         k_o       = 3.e-8_dbl_kind, & ! permeability scaling factor (m^2)
         rhosi     = 940.0_dbl_kind, & ! average sea ice density
                                       ! Cox and Weeks, 1982: 919-974 kg/m^2
         min_salin = p1            , & ! threshold for brine pocket treatment 
         hbr_min   = p01           , & ! minimum hbrine thickness
         thinS     = 0.05_dbl_kind     ! minimum ice thickness for brine

      real (kind=dbl_kind), public :: & 
         phi_snow   ! porosity of snow

      real (kind=dbl_kind), & 
         dimension (nx_block,ny_block,nblyr+1,ncat,max_blocks), public :: &
         Zoo        ! N losses accumulated in timestep (ie. zooplankton/bacteria)
                    ! mmol/m^3

      real (kind=dbl_kind), & 
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         dhbr_top     , & ! brine top change
         dhbr_bot         ! brine bottom change

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_blocks), public :: &
         grow_net       , & ! Specific growth rate (m/d) times vice
         PP_net         , & ! Total production (mg C/m^2/d) times aice
         hbri               ! brine height, area-averaged for comparison with hi (m)

      real (kind=dbl_kind), dimension (nblyr+2), public :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), public :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), public :: &
         cgrid            , &  ! CICE vertical coordinate   
         icgrid           , &  ! interface grid for CICE (shortwave variable)
         swgrid                ! grid for ice tracers used in dEdd scheme

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,nblyr+2,ncat,max_blocks), public :: &
         bphi           , & ! porosity of layers    
         bTiz               ! layer temperatures interpolated on bio grid (C)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         darcy_V            ! darcy velocity positive up (m/s)

      !-----------------------------------------------------------------
      ! vertical biogeochemistry  
      !-----------------------------------------------------------------

      logical (kind=log_kind), public :: &
         tr_zaero,       & ! if .true., black carbon is tracers on ice (n_zaero)
         z_tracers,      & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc,      & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae        ! if .true., algal absorptionof Shortwave is computed in the
                           ! delta-Eddington code (must use shortwave = 'dEdd')

      real (kind=dbl_kind), public :: & 
         grid_o      , & ! for bottom flux        
         l_sk        , & ! characteristic diffusive scale (zsalinity) (m)
         grid_o_t    , & ! top grid point length scale 
         initbio_frac, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav     ! multiple of ocean tracer concentration due to frazil scavenging

      real (kind=dbl_kind), parameter, public :: & 
         bphimin = 0.03_dbl_kind      ! minimum porosity for zbgc only

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         zsal_tot    , & ! Total ice salinity in per grid cell (g/m^2) 
         chl_net     , & ! Total chla (mg chla/m^2) per grid cell      
         NO_net          ! Total production (mg C/m^2/s) per grid cell  

      logical (kind=log_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         Rayleigh_criteria ! .true. means Ra_c was reached   

      real (kind=dbl_kind), public :: & 
         grid_oS     , & ! for bottom flux (zsalinity)
         l_skS           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)

      real (kind=dbl_kind), parameter, public :: & 
         viscos_dynamic = 2.2_dbl_kind   , & ! 1.8e-3_dbl_kind (pure water at 0^oC) (kg/m/s)
         min_bgc        = 0.01_dbl_kind  , & ! fraction of ocean bgc concentration in surface melt 
         Dm             = 1.0e-9_dbl_kind, & ! molecular diffusion (m^2/s)
         Ra_c           = 0.05_dbl_kind      ! critical Rayleigh number for bottom convection

      real (kind=dbl_kind), & 
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         sice_rho     ! avg sea ice density  (kg/m^3)  ! ech: diagnostic only?

      real (kind=dbl_kind), & 
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         fzsaln, &    ! category fzsal(kg/m^2/s) 
         fzsaln_g     ! salt flux from gravity drainage only

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         fzsal    , & ! Total flux  of salt to ocean at time step for conservation
         fzsal_g      ! Total gravity drainage flux

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1,ncat,max_blocks), public :: &
         zfswin       ! Shortwave flux into layers interpolated on bio grid  (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1,ncat,max_blocks), public :: &
         iDi      , & ! igrid Diffusivity (m^2/s)    
         iki          ! Ice permeability (m^2)     

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         upNO     , & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH         ! ammonium uptake rate (mmol/m^2/d) times aice

      real (kind=dbl_kind), public :: & 
         dts_b        ! zsalinity timestep
        
!=======================================================================

      contains

!=======================================================================
! 
! Compute the internal ice  enthalpy using new salinity and Tin
!

      function calculate_qin_from_Sin (Tin, Tmltk) &
               result(qin)
            
      use ice_constants, only: c1, rhoi, cp_ocn, cp_ice, Lfresh  

      real (kind=dbl_kind), intent(in) :: &
         Tin                ,&  ! internal temperature
         Tmltk                  ! melting temperature at one level

      ! local variables

     real (kind=dbl_kind) :: &
         qin                    ! melting temperature at one level   

         qin =-rhoi*(cp_ice*(Tmltk-Tin) + Lfresh*(c1-Tmltk/Tin) - cp_ocn*Tmltk)

      end function calculate_qin_from_Sin

!=======================================================================
!
! Remaps tracer fields in a given category from one set of layers to another.
! Grids can be very different and  so can  vertical spaces.  

      subroutine remap_zbgc(ntrcr,    nlyrn,    &
                            it,                 &
                            trcrn,    trtmp,    &
                            nr0,      nbyrn,    &
                            hice,     hinS,     &
                            ice_grid, bio_grid, &
                            S_min)

      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice

      integer (kind=int_kind), intent(in) :: &
         ntrcr         , & ! number of tracers in use
         it            , & ! tracer index in top layer
         nr0           , & ! receiver category
         nlyrn         , & ! number of ice layers
         nbyrn             ! number of biology layers

      real (kind=dbl_kind), dimension (ntrcr), &
         intent(in) ::       &
         trcrn             ! ice tracers

      real (kind=dbl_kind), dimension (ntrcr+2), &
         intent(inout) ::    &
         trtmp             ! temporary, remapped ice tracers

      real (kind=dbl_kind), dimension (nlyrn), intent(in) :: &
         ice_grid          ! CICE grid  cgrid(2:nilyr+1)

      real (kind=dbl_kind), dimension (nbyrn), intent(in) :: &
         bio_grid          ! CICE grid  grid(2:nbyrn+1)

      real(kind=dbl_kind), intent(in) :: &
         hice          , & ! CICE ice thickness
         hinS          , & ! brine height 
         S_min             ! for salinity on CICE grid        

      ! local variables

      integer (kind=int_kind) :: &
           kd, kr, kdr , & ! more indices
           kdi         , & ! more indices
           n_nd        , & ! number of layers in donor
           n_nr, n_plus    ! number of layers in receiver

      real (kind=dbl_kind), dimension (nbyrn+3+nlyrn) :: &
           trdr        , & ! combined tracer 
           trgrid          ! combined grid 

      real (kind=dbl_kind), dimension (nbyrn+nlyrn+3) :: &
           tracer      , & ! temporary,  ice tracers values
           dgrid       , & ! temporary,  donor grid dimensional
           rgrid           !  temporary, receiver grid dimensional

      if ((hinS < c0) .OR. (hice < c0)) then
         write(nu_diag, *)'Problem in remap_layers_bgc'
         write(nu_diag, *) '(hinS < c0) .OR. (hice < c0)'
         write(nu_diag, *) 'hinS,hice',hinS,hice
         call abort_ice ('ice: remap_layers_bgc error')
      endif
         
      if (nr0 == 0) then ! cice to bio

         n_nd            = nlyrn
         n_nr            = nbyrn
         n_plus          = 2
         dgrid (1)       = min(-hice+hinS, -hinS+hice, c0)            
         dgrid (nlyrn+2) = min(hinS, hice) 
	 tracer(1)       = trcrn(it)
	 tracer(nlyrn+2) = trcrn(it+nlyrn-1)
         rgrid (nbyrn+2) = min(hinS, hice)
         if (hice > hinS) then
            rgrid(1) = c0 
	    do kr = 1,n_nr
	       rgrid(kr+1) = bio_grid(kr)*hinS
            enddo
            do kd = 1,n_nd
	       dgrid(kd+1) = (ice_grid(kd)-c1)*hice+hinS
               tracer(kd+1) = trcrn(it+kd-1)
            enddo
         else
            rgrid(1) = -hinS + hice 
            do kr = 1,n_nr
	       rgrid(kr+1) = (bio_grid(kr)-c1)*hinS + hice
            enddo
            do kd = 1,n_nd
	       dgrid(kd+1) = ice_grid(kd)*hice
               tracer(kd+1) = trcrn(it+kd-1)
            enddo
         endif
              
      else               ! bio to cice

         n_nd = nbyrn
         n_nr = nlyrn
         if (hice > hinS) then
            n_plus          = 3
            tracer(1)       = S_min
            tracer(2)       = S_min
            dgrid (1)       = -hice+hinS
            dgrid (2)       = p5*(hinS-hice)
            dgrid (nbyrn+3) = hinS
            tracer(nbyrn+3) = trcrn(it+nbyrn-1)
            rgrid (1)       = -hice + hinS
            rgrid (nlyrn+2) = hinS 
            do kd = 1,n_nd
               dgrid(kd+2) = bio_grid(kd)*hinS
               tracer(kd+2) = trcrn(it+kd-1)
            enddo
            do kr = 1,n_nr
               rgrid(kr+1) = (ice_grid(kr)-c1)*hice+ hinS
            enddo
         else
            n_plus          = 2
            tracer(1)       = trcrn(it)
            tracer(nbyrn+2) = trcrn(it+nbyrn-1)
            dgrid (1)       = hice-hinS
            dgrid (nbyrn+2) = hice
            rgrid (nlyrn+2) = hice
            rgrid (1)       = c0
            do kd = 1,n_nd
              dgrid(kd+1) = (bio_grid(kd)-c1)*hinS + hice
              tracer(kd+1) = trcrn(it+kd-1)
            enddo
            do kr = 1,n_nr
              rgrid(kr+1) = ice_grid(kr)*hice
            enddo
         endif

      endif

      kdr = 0  !combined indices
      kdi = 1  

      do kr = 1, n_nr
         do kd = kdi, n_nd+n_plus
            if (dgrid(kd) < rgrid(kr+1)) then
               kdr = kdr+1
               trgrid(kdr) = dgrid(kd)
               trdr  (kdr) = tracer(kd)
            elseif (dgrid(kd) > rgrid(kr+1)) then
               kdr = kdr + 1
               kdi = kd
               trgrid(kdr) = rgrid(kr+1)
               trtmp (it+kr-1)  = trdr(kdr-1) &
                           + (rgrid(kr+1) - trgrid(kdr-1)) &
                           * (tracer(kd) - trdr(kdr-1)) &
                           / (dgrid(kd) - trgrid(kdr-1))
               trdr(kdr) = trtmp(it+kr-1) 
               EXIT
            else
               kdr = kdr+1
               kdi = kd+1
               trgrid(kdr) = rgrid(kr+1)
               trtmp (it+kr-1)  = tracer(kd)              
               trdr  (kdr) = tracer(kd)
               EXIT
            endif
         enddo
      enddo

      end subroutine remap_zbgc

!=======================================================================

! When sea ice melts laterally, flux bgc to ocean

      subroutine lateral_melt_bgc (nx_block, ny_block, &
                                   icells,   dt,       &
                                   indxi,    indxj,    &
                                   rside,    vicen,    &
                                   trcrn,    fzsal,    &
                                   flux_bio, nbltrcr)

      use ice_state, only: nt_fbri, nt_bgc_S
      use ice_therm_shared, only: solve_zsal
      use ice_constants, only: c1, p001
      use ice_domain_size, only: max_ntrcr

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
         fzsal     ! salt flux from layer Salinity (kg/m^2/s)
  
      real (kind=dbl_kind), dimension(nx_block,ny_block,nbltrcr), &
         intent(inout) :: &
         flux_bio  ! biology tracer flux from layer bgc (mmol/m^2/s)

      ! local variables

      integer (kind=int_kind) :: &
         i, j  , & ! horizontal indices
         k     , & ! layer index
         ij, m     ! horizontal index, combines i and j loops

      real (kind=dbl_kind) :: &
         zspace    ! bio grid spacing

       zspace = c1/(real(nblyr,kind=dbl_kind))

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells
          i = indxi(ij)
          j = indxj(ij)
            
          if (solve_zsal) then
             do k = 1,nblyr
                fzsal(i,j) = fzsal(i,j) + rhosi*trcrn(i,j,nt_fbri) &
                           * vicen(i,j)*p001*zspace*trcrn(i,j,nt_bgc_S+k-1) &
                           * rside(i,j)/dt
             enddo
          endif

          do m = 1, nbltrcr
          do k = 1,nblyr+1
             flux_bio(i,j,m) = flux_bio(i,j,m) + trcrn(i,j,nt_fbri) &
                             * vicen(i,j)*zspace*trcrn(i,j,bio_index(m)+k-1) &
                             * rside(i,j)/dt
          enddo
          enddo
       enddo                  ! ij

      end subroutine lateral_melt_bgc 

!=======================================================================

! remove tracer for very small fractional areas

      subroutine zap_small_bgc (zlevels,  dflux_bio, &
                                dt, zvol, btrcr)

      integer (kind=int_kind), intent(in) :: &
         zlevels    ! number of vertical levels in ice

      real (kind=dbl_kind), intent(in) :: &
         dt         ! time step (s)

      real (kind=dbl_kind), intent(inout) :: &
         dflux_bio  ! zapped bio tracer flux from biology (mmol/m^2/s)

      real (kind=dbl_kind), dimension (zlevels), intent(in) :: &
         btrcr  , & ! zapped bio tracer flux from biology (mmol/m^2/s)
         zvol       ! ice volume (m)

      ! local variables

      integer (kind=int_kind) :: &
         k          ! layer index

      do k = 1, zlevels
         dflux_bio = dflux_bio + btrcr(k)*zvol(k)/dt
      enddo
          
      end subroutine zap_small_bgc

!=======================================================================

      end module ice_zbgc_shared

!=======================================================================
