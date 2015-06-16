!  SVN:$Id$
!=======================================================================
!
! Biogeochemistry driver
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_zbgc

      use ice_kinds_mod
      use ice_zbgc_shared ! everything
      use ice_domain_size, only: max_nsw, ncat, max_blocks
      use ice_blocks, only: nx_block, ny_block

      implicit none 

      private
      public :: add_new_ice_bgc, init_zbgc, init_bgc, &
                init_history_bgc, biogeochemistry

!=======================================================================

      contains

!=======================================================================

! Namelist variables, set to default values; may be altered at run time
! 
! author Elizabeth C. Hunke, LANL
!        Nicole Jeffery, LANL

      subroutine init_zbgc

      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: dt
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c1, p5, c0, c5, c2, p1, rhos, rhoi, c3600
      use ice_domain_size, only: max_ntrcr, max_nbtrcr, nblyr, nilyr, nslyr, &
                           n_algae, max_algae, n_zaero, max_aero, &
                           n_doc, max_doc, n_dic, max_dic, max_don, n_don, &
                           n_fed, n_fep, max_fe
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_nml, nml_filename, get_fileunit, &
                               release_fileunit, nu_diag
      use ice_restart_shared, only: runtype
      use ice_shortwave, only: shortwave
      use ice_state, only: tr_brine, nt_fbri, ntrcr, nbtrcr, &
          trcr_depend, nbtrcr_sw, &
          nt_bgc_Nit,    nt_bgc_Am,    nt_bgc_Sil, &
          nt_bgc_DMS,    nt_bgc_PON,   nt_bgc_S, &
          nt_bgc_N,      nt_bgc_C,     nt_bgc_chl, &
          nt_bgc_DOC,    nt_bgc_DON,   nt_bgc_DIC, &
          nt_zaero  ,    nt_bgc_DMSPp, nt_bgc_DMSPd, &
          nt_bgc_Fed,    nt_bgc_Fep,   nt_zbgc_frac

      use ice_therm_shared, only: solve_zsal, ktherm

      integer (kind=int_kind) :: &
        nml_error, & ! namelist i/o error flag
        k, mm    , & ! loop index  
        ntd          ! for tracer dependency calculation 

      ! zbgc 
      real (kind=dbl_kind), parameter :: & 
         l_sk_scale = 1.e-2_dbl_kind, & ! 
         grid_o_scale = 1.e-3_dbl_kind !

      ! Transport type 
      !------------------------------------------------------------
      ! In delta Eddington, algal particles are assumed to cause no
      ! significant scattering (Brieglib and Light), only absorption
      ! in the visible spectral band (200-700 nm)
      ! Algal types: Diatoms, flagellates, Phaeocycstis
      ! DOC        : Proteins, EPS, Lipids
      !------------------------------------------------------------
      !        Tracers have mobile and  stationary phases. 
      ! ice growth allows for retention, ice melt facilitates mobility
      ! bgc_tracer_type defines the exchange timescales between these phases
      ! -1 : entirely in the mobile phase, no exchange  (this is the default)
      !  0 : retention time scale is tau_min, release time scale is tau_max
      !  1 : retention time scale is tau_max, release time scale is tau_min
      ! 0.5: retention time scale is tau_min, release time scale is tau_min
      !  2 : retention time scale is tau_max, release time scale is tau_max
      ! tau_min and tau_max are define in ice_zbgc_shared.f90
      !------------------------------------------------------------
      real (kind=dbl_kind), parameter, dimension(max_algae) :: &  
        algaltype = (/c0, c0, c0/)  ! strongly retained 
      real (kind=dbl_kind), parameter :: &  
        nitratetype = -c1, & ! purely mobile
        ammoniumtype = p5, & ! fast exchange 
        silicatetype =-c1, &
        dmspptype    = c0, & ! fast exchange (retained and released)
        dmspdtype    = -c1
      real (kind=dbl_kind), parameter, dimension(max_doc) :: &  
        doctype   = (/ c0, c0, c0/)  
      real (kind=dbl_kind), parameter, dimension(max_dic) :: &  
        dictype   = (/ -c1/) 
      real (kind=dbl_kind), parameter, dimension(max_don) :: &  
        dontype   = (/ c0/)  
      real (kind=dbl_kind), parameter, dimension(max_fe) :: &  
        fedtype    = (/ c0, c0 /)  
      real (kind=dbl_kind), parameter, dimension(max_fe) :: &  
        feptype    = (/ c0, c0/)  
      !------------------------------------------------------------
      ! Aerosol order and type should be consistent with order/type 
      ! specified in delta Eddington:  1) hydrophobic black carbon;
      ! 2) hydrophilic black carbon; 3) dust (0.05-0.5 micron);
      ! 4) dust (0.5-1.25 micron); 5) dust (1.25-2.5 micron);
      ! 6) dust (2.5-5 micron) 
      !-------------------------------------------------------------
      real (kind=dbl_kind), parameter, dimension(max_aero) :: &  
        zaerotype   = (/ 0.1_dbl_kind, 1.0_dbl_kind, 0.5_dbl_kind, &
                    1.0_dbl_kind, 1.0_dbl_kind, 1.0_dbl_kind/) 

      !-----------------------------------------------------------------
      ! namelist variables
      !-----------------------------------------------------------------

      namelist /zbgc_nml/  &
        tr_brine, restart_hbrine, tr_zaero, skl_bgc, z_tracers, &
        dEdd_algae, solve_zbgc, bgc_flux_type, &
        restore_bgc, restart_bgc, scale_bgc, solve_zsal, restart_S, &
        bgc_data_dir, sil_data_type, nit_data_type,  fe_data_type, &
        tr_bgc_Nit, tr_bgc_C, tr_bgc_chl, tr_bgc_Am, tr_bgc_Sil, &
        tr_bgc_DMS, tr_bgc_PON, tr_bgc_DON, tr_bgc_Fe, &
        grid_o, grid_o_t, l_sk, grid_oS, &   
        l_skS, phi_snow,  initbio_frac, frazil_scav

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------
      tr_brine        = .false.  ! brine height differs from ice height
      tr_zaero        = .false.  ! z aerosol tracers
      restore_bgc     = .false.  ! restore bgc if true
      solve_zsal      = .false.  ! update salinity tracer profile from solve_S_dt
      bgc_data_dir    = 'unknown_bgc_data_dir'
      sil_data_type   = 'default'
      nit_data_type   = 'default'
      fe_data_type    = 'default'
      restart_bgc     = .false.  ! biogeochemistry restart
      restart_S       = .false.  ! salinity restart
      restart_hbrine  = .false.  ! hbrine restart
      scale_bgc       = .false.  ! initial bgc tracers proportional to S  
      skl_bgc         = .false.  ! solve skeletal biochemistry 
      z_tracers       = .false.  ! solve vertically resolved tracers
      dEdd_algae      = .false.  ! dynamic algae contributes to shortwave absorption
                                 ! in delta-Eddington calculation
      solve_zbgc      = .false.  ! turn on z layer biochemistry 
      tr_bgc_PON      = .false.  !---------------------------------------------   
      tr_bgc_Nit      = .false.  ! biogeochemistry (skl or zbgc)
      tr_bgc_C        = .false.  ! if (skl_bgc = .true. then skl)
      tr_bgc_chl      = .false.  ! if (z_tracers = .true. then vertically resolved)
      tr_bgc_Sil      = .false.  ! if (z_tracers + solve_zbgc = .true. then
      tr_bgc_Am       = .false.  ! vertically resolved with  reactions  
      tr_bgc_DMS      = .false.  !------------------------------------------------
      tr_bgc_DON      = .false.  ! 
      tr_bgc_Fe       = .false.  ! 

      ! brine height parameter
      phi_snow        = p5       ! snow porosity

      !skl biology parameters
      bgc_flux_type   = 'Jin2006'! type of ocean-ice poston velocity ('constant')

      !z biology parameters  
      grid_o          = c5            ! for bottom flux        
      grid_o_t        = c5            ! for top flux        
      l_sk            = 7.0_dbl_kind  ! characteristic diffusive scale (m)   
      initbio_frac    = c1            ! fraction of ocean tracer concentration in bio tracers
      frazil_scav     = c1            ! increase in initial bio tracer from ocean scavenging 

      !z salinity  parameters
      grid_oS         = c5            ! for bottom flux         
      l_skS           = 7.0_dbl_kind  ! characteristic diffusive scale (m)  

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)

      if (my_task == master_task) then
         open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif 

         print*,'Reading zbgc_nml'
         do while (nml_error > 0)
            read(nu_nml, nml=zbgc_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         call abort_ice('ice: error reading zbgc namelist')
      endif
      call release_fileunit(nu_nml)

      !-----------------------------------------------------------------
      ! zsalinity and brine
      !-----------------------------------------------------------------
      if (solve_zsal .and. TRZS == 0) then
         write(nu_diag,*) &
            'WARNING: solve_zsal=T but 0 zsalinity  tracer  compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zsal = F'
         solve_zsal = .false.      
      elseif (solve_zsal .and. nblyr < 1)  then
         write(nu_diag,*) &
            'WARNING: solve_zsal=T but 0 zsalinity  tracer  compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zsal = F'
         solve_zsal = .false.     
      endif 
  
      if (solve_zsal) then
        tr_brine   = .true. 
        ktherm = 1
      endif

      if (tr_brine .and. TRBRI == 0 ) then
         write(nu_diag,*) &
            'WARNING: tr_brine=T but no brine height compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zsal and tr_brine = F'
         solve_zsal = .false.
         tr_brine  = .false.
      elseif (tr_brine .and. nblyr < 1 ) then
         write(nu_diag,*) &
            'WARNING: tr_brine=T but no biology layers compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zsal and tr_brine = F'
         solve_zsal = .false.
         tr_brine  = .false.
      endif 

      call broadcast_scalar(solve_zsal,         master_task)  
      call broadcast_scalar(restart_S,          master_task)  
      call broadcast_scalar(tr_brine,           master_task)
      call broadcast_scalar(restart_hbrine,     master_task) 

      if (phi_snow .le. c0) phi_snow = c1-rhos/rhoi
      call broadcast_scalar(phi_snow,           master_task)

      nt_fbri = 0
      if (tr_brine) then
          nt_fbri = ntrcr + 1   ! ice volume fraction with salt
          ntrcr = ntrcr + 1
          trcr_depend(nt_fbri) = 1  ! volume-weighted
      endif

      ntd = 0                    ! if nt_fbri /= 0 then use fbri dependency
      if (nt_fbri == 0) ntd = -1 ! otherwise make tracers depend on ice volume

      if (solve_zsal)then
          nt_bgc_S = ntrcr + 1
          ntrcr = ntrcr + nblyr
          do k = 1,nblyr
             trcr_depend(nt_bgc_S     + k - 1)  = 2+nt_fbri+ntd
          enddo
      endif 

      grid_oS = grid_oS * grid_o_scale 
      l_skS = l_skS * l_sk_scale

      call broadcast_scalar(grid_oS,            master_task)
      call broadcast_scalar(l_skS,              master_task)

      if (my_task == master_task) then
         write(nu_diag,1010) ' tr_brine                  = ', tr_brine
         if (tr_brine) then
         write(nu_diag,1010) ' restart_hbrine            = ', restart_hbrine
         write(nu_diag,1005) ' phi_snow                  = ', phi_snow
         endif
         if (solve_zsal) then
         write(nu_diag,1010) ' solve_zsal                = ', solve_zsal
         write(nu_diag,1010) ' restart_S                 = ', restart_S
         write(nu_diag,1000) ' grid_oS                   = ', grid_oS
         write(nu_diag,1005) ' l_skS                     = ', l_skS
         endif
      endif

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      if ((skl_bgc .AND. solve_zbgc) .or. (skl_bgc .AND. z_tracers)) &
              call abort_ice('error:skl_bgc &
              and solve_zbgc or z_tracers are both true')

      if (skl_bgc .AND. tr_zaero) then
         write(nu_diag,*) &
            'WARNING: skl bgc does not use vertical tracers'
         write(nu_diag,*) &
            'WARNING: setting tr_zaero = F'
         tr_zaero = .false.
      endif

      if (dEdd_algae .AND. trim(shortwave) /= 'dEdd') then 
         write(nu_diag,*) &
            'WARNING: need shortwave = dEdd if dEdd_algae = T'
         write(nu_diag,*) &
            'WARNING: setting dEdd_algae = F'
         dEdd_algae = .false.
      endif

      if (dEdd_algae .AND. .NOT. tr_bgc_N .AND. .NOT. tr_zaero) then 
         write(nu_diag,*) &
            'WARNING: need tr_bgc_N or tr_zaero if dEdd_algae = T'
         write(nu_diag,*) &
            'WARNING: setting dEdd_algae = F'
         dEdd_algae = .false.
      endif
        
      if (n_algae > max_algae) call abort_ice('error:number of algal &
            types exceeds max_algae')
      if (n_doc > max_doc) call abort_ice('error:number of doc &
            types exceeds max_doc')
      if (n_dic > max_doc) call abort_ice('error:number of dic &
            types exceeds max_dic')
      if (n_don > max_don) call abort_ice('error:number of don &
            types exceeds max_don')
      if (n_fed  > max_fe ) call abort_ice('error:number of dissolved fe &
            types exceeds max_fe ')
      if (n_fep  > max_fe ) call abort_ice('error:number of particulate fe &
            types exceeds max_fe ')
      if ((TRBGCS == 0 .and. skl_bgc) .or. (TRALG == 0 .and. skl_bgc)) then
         write(nu_diag,*) &
            'WARNING: skl_bgc=T but 0 bgc or algal tracers compiled'
         write(nu_diag,*) &
            'WARNING: setting skl_bgc = F'
         skl_bgc = .false.
      endif

      if ((TRBGCZ == 0 .and. solve_zbgc) .or. (TRALG == 0 .and. solve_zbgc)) then
         write(nu_diag,*) &
            'WARNING: solve_zbgc=T but 0 zbgc or algal tracers compiled'
         write(nu_diag,*) &
            'WARNING: setting solve_zbgc = F'
         solve_zbgc = .false.
      endif

      if (solve_zbgc .and. .not. z_tracers) z_tracers = .true.

      call broadcast_scalar(solve_zbgc,         master_task)
      call broadcast_scalar(skl_bgc,            master_task)
      call broadcast_scalar(restart_bgc,        master_task)

      if (skl_bgc .or. solve_zbgc) then
            tr_bgc_N         = .true.   ! minimum NP biogeochemistry
            tr_bgc_Nit       = .true.
      endif

      call broadcast_scalar(bgc_flux_type,      master_task)
      call broadcast_scalar(restore_bgc,        master_task)
      call broadcast_scalar(bgc_data_dir,       master_task)
      call broadcast_scalar(sil_data_type,      master_task)
      call broadcast_scalar(nit_data_type,      master_task)
      call broadcast_scalar(fe_data_type,       master_task)
      call broadcast_scalar(tr_bgc_N,           master_task)
      call broadcast_scalar(tr_bgc_C,           master_task)
      call broadcast_scalar(tr_bgc_chl,         master_task)
      call broadcast_scalar(tr_bgc_Nit,         master_task)
      call broadcast_scalar(tr_bgc_Am,          master_task)
      call broadcast_scalar(tr_bgc_Sil,         master_task)
      call broadcast_scalar(tr_bgc_DMS,         master_task) 
      call broadcast_scalar(tr_bgc_PON,         master_task) 
      call broadcast_scalar(tr_bgc_DON,         master_task) 
      call broadcast_scalar(tr_bgc_Fe,          master_task) 

      !-----------------------------------------------------------------
      ! z layer aerosols
      !-----------------------------------------------------------------
      if (tr_zaero .and. .not. z_tracers) z_tracers = .true.

      if (n_zaero > max_aero) call abort_ice('error:number of z aerosols &
            exceeds max_aero')

      grid_o = grid_o * grid_o_scale
      grid_o_t = grid_o_t * grid_o_scale
      l_sk = l_sk * l_sk_scale
         
      call broadcast_scalar(z_tracers,          master_task)
      call broadcast_scalar(tr_zaero,           master_task)
      call broadcast_scalar(grid_o,             master_task)
      call broadcast_scalar(grid_o_t,           master_task)
      call broadcast_scalar(l_sk,               master_task)
      call broadcast_scalar(scale_bgc,          master_task)
      call broadcast_scalar(initbio_frac,       master_task)
      call broadcast_scalar(frazil_scav,        master_task)

      nbtrcr = 0
      nbtrcr_sw = 0
      !--------------------------
      ! vectors of size max_algae
      !--------------------------
      nlt_bgc_N(:) = 0
      nlt_bgc_C(:) = 0
      nlt_bgc_chl(:) = 0
      !--------------------------
      ! vectors of size max_dic
      !--------------------------
      nlt_bgc_DIC(:) = 0
      !--------------------------
      ! vectors of size max_doc
      !--------------------------
      nlt_bgc_DOC(:) = 0
      !--------------------------
      ! vectors of size max_don
      !--------------------------
      nlt_bgc_DON(:) = 0
      !--------------------------
      ! vectors of size max_fe 
      !--------------------------
      nlt_bgc_Fed(:) = 0
      nlt_bgc_Fep(:) = 0
      !--------------------------
      ! vectors of size max_aero
      !--------------------------
      nlt_zaero(:) = 0
      nlt_zaero_sw(:) = 0

      nlt_bgc_Nit = 0
      nlt_bgc_Am = 0
      nlt_bgc_Sil = 0
      nlt_bgc_DMSPp = 0
      nlt_bgc_DMSPd = 0
      nlt_bgc_DMS = 0
      nlt_bgc_PON = 0
      nlt_chl_sw  = 0
      bio_index(:) = 0
      bio_index_o(:) = 0

      if (skl_bgc) then

      !-----------------------------------------------------------------
      ! assign tracer indices and dependencies
      ! bgc_tracer_type: < 0  purely mobile , >= 0 stationary 
      !------------------------------------------------------------------

      do mm = 1,n_algae
           ntrcr = ntrcr + 1             
           nt_bgc_N(mm) = ntrcr
           nbtrcr = nbtrcr + 1
           nlt_bgc_N(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)    = algaltype(mm)
           trcr_depend(nt_bgc_N(mm))  = c0
           bio_index(nlt_bgc_N(mm))   = ntrcr
           bio_index_o(nlt_bgc_N(mm)) = mm  
      enddo   !mm
      nlt_chl_sw = 1
      if (dEdd_algae) nbtrcr_sw = nilyr+nslyr+2  ! only the bottom layer 
                                                 ! will be nonzero               
      ntrcr = ntrcr + 1              
      nt_bgc_Nit = ntrcr
      nbtrcr = nbtrcr + 1
      nlt_bgc_Nit = nbtrcr
      bgc_tracer_type(nbtrcr)         = nitratetype
      trcr_depend(nt_bgc_Nit)         = c0
      bio_index(nlt_bgc_Nit)          = ntrcr
      bio_index_o(nlt_bgc_Nit)        = max_algae + 1  

      if (tr_bgc_C) then
       !
       ! Algal C is not yet distinct from algal N
       ! * Reqires  exudation and/or changing C:N ratios
       ! for implementation
       !
       ! do mm = 1,n_algae
       !    ntrcr = ntrcr + 1             
       !    nt_bgc_C(mm) = ntrcr
       !    nbtrcr = nbtrcr + 1
       !    nlt_bgc_C(mm) = nbtrcr
       !    bgc_tracer_type(nbtrcr)     = algaltype(mm) 
       !    trcr_depend(nt_bgc_C(mm))   = c0
       !    bio_index(nlt_bgc_C(mm))    = ntrcr
       !    bio_index_o(nlt_bgc_C(mm))  = max_algae + 1 + mm
       ! enddo !mm
        do mm = 1,n_doc
           ntrcr = ntrcr + 1             
           nt_bgc_DOC(mm) = ntrcr
           nbtrcr = nbtrcr + 1
           nlt_bgc_DOC(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)      = doctype(mm) 
           trcr_depend(nt_bgc_DOC(mm))  = c0
           bio_index(nlt_bgc_DOC(mm))   = ntrcr
           bio_index_o(nlt_bgc_DOC(mm)) = max_algae + 1 + mm
        enddo !mm
        do mm = 1,n_dic
           ntrcr = ntrcr + 1             
           nt_bgc_DIC(mm) = ntrcr
           nbtrcr = nbtrcr + 1
           nlt_bgc_DIC(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)      = dictype(mm)
           trcr_depend(nt_bgc_DIC(mm))  = c0
           bio_index(nlt_bgc_DIC(mm))   = ntrcr
           bio_index_o(nlt_bgc_DIC(mm)) = max_algae + max_doc+ 1 + mm
        enddo !mm

      endif !tr_bgc_C

      if (tr_bgc_chl)then
        do mm = 1,n_algae
           ntrcr = ntrcr + 1             
           nt_bgc_chl(mm) = ntrcr
           nbtrcr = nbtrcr + 1
           nlt_bgc_chl(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)      = algaltype(mm)
           trcr_depend(nt_bgc_chl(mm))  = c0
           bio_index(nlt_bgc_chl(mm))   = ntrcr
           bio_index_o(nlt_bgc_chl(mm)) = max_algae+ 1 + max_doc + max_dic + mm
         enddo !mm
      endif !tr_bgc_chl

         if (tr_bgc_Am)then
             ntrcr = ntrcr + 1
             nt_bgc_Am = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_Am= nbtrcr
             bgc_tracer_type(nbtrcr)   = ammoniumtype
             trcr_depend(nt_bgc_Am)    = c0
             bio_index(nlt_bgc_Am)     = ntrcr
             bio_index_o(nlt_bgc_Am)   = 2*max_algae + max_doc + max_dic + 2
         endif    
         if (tr_bgc_Sil)then
             ntrcr = ntrcr + 1
             nt_bgc_Sil = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_Sil = nbtrcr
             bgc_tracer_type(nbtrcr)   = silicatetype
             trcr_depend(nt_bgc_Sil)   = c0
             bio_index(nlt_bgc_Sil)    = ntrcr
             bio_index_o(nlt_bgc_Sil)  = 2*max_algae + max_doc + max_dic + 3
         endif    
         if (tr_bgc_DMS)then
             ntrcr = ntrcr + 1
             nt_bgc_DMSPp = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_DMSPp = nbtrcr
             bgc_tracer_type(nbtrcr)   = dmspptype
             trcr_depend(nt_bgc_DMSPp) = c0
             bio_index(nlt_bgc_DMSPp)  = ntrcr
             bio_index_o(nlt_bgc_DMSPp)= 2*max_algae + max_doc + max_dic + 4
             ntrcr = ntrcr + 1
             nt_bgc_DMSPd = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_DMSPd = nbtrcr
             bgc_tracer_type(nbtrcr)   = dmspdtype
             trcr_depend(nt_bgc_DMSPd) = c0
             bio_index(nlt_bgc_DMSPd)  = ntrcr
             bio_index_o(nlt_bgc_DMSPd)= 2*max_algae + max_doc + max_dic + 5
             ntrcr = ntrcr + 1
             nt_bgc_DMS = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_DMS = nbtrcr
             bgc_tracer_type(nbtrcr)   = dmspdtype
             trcr_depend(nt_bgc_DMS)   = c0
             bio_index(nlt_bgc_DMS)    = ntrcr
             bio_index_o(nlt_bgc_DMS)  = 2*max_algae + max_doc + max_dic + 6
         endif 
         if (tr_bgc_PON)then
             ntrcr = ntrcr + 1
             nt_bgc_PON = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_PON = nbtrcr
             bgc_tracer_type(nbtrcr)   = nitratetype
             trcr_depend(nt_bgc_PON)   = c0
             bio_index(nlt_bgc_PON)    = ntrcr
             bio_index_o(nlt_bgc_PON)  = 2*max_algae + max_doc + max_dic + 7
         endif     
         if (tr_bgc_DON)then
            do mm = 1,n_don
               ntrcr = ntrcr + 1             
               nt_bgc_DON(mm) = ntrcr
               nbtrcr = nbtrcr + 1
               nlt_bgc_DON(mm) = nbtrcr
               bgc_tracer_type(nbtrcr)      = dontype(mm)
               trcr_depend(nt_bgc_DON(mm))  = c0
               bio_index(nlt_bgc_DON(mm))   = ntrcr
               bio_index_o(nlt_bgc_DON(mm)) = 2*max_algae + max_doc + max_dic + 7 + mm
            enddo !mm
         endif            
         if (tr_bgc_Fe)then
            do mm = 1,n_fed
               ntrcr = ntrcr + 1             
               nt_bgc_Fed(mm) = ntrcr
               nbtrcr = nbtrcr + 1
               nlt_bgc_Fed(mm) = nbtrcr
               bgc_tracer_type(nbtrcr)     = fedtype(mm)
               trcr_depend(nt_bgc_Fed(mm))  = c0
               bio_index(nlt_bgc_Fed(mm))   = ntrcr
               bio_index_o(nlt_bgc_Fed(mm)) = 2*max_algae + max_doc + max_dic + max_don + 7 + mm
            enddo !mm
            do mm = 1,n_fep
               ntrcr = ntrcr + 1             
               nt_bgc_Fep(mm) = ntrcr
               nbtrcr = nbtrcr + 1
               nlt_bgc_Fep(mm) = nbtrcr
               bgc_tracer_type(nbtrcr)     = feptype(mm)
               trcr_depend(nt_bgc_Fep(mm))  = c0
               bio_index(nlt_bgc_Fep(mm))   = ntrcr
               bio_index_o(nlt_bgc_Fep(mm)) = 2*max_algae + max_doc + max_dic + max_don + max_fe + 7 + mm
            enddo !mm
         endif        

      if (my_task == master_task) then

         write(nu_diag,1010) ' skl_bgc                   = ', skl_bgc
         write(nu_diag,1030) ' bgc_flux_type             = ', bgc_flux_type
         write(nu_diag,1010) ' restart_bgc               = ', restart_bgc
         write(nu_diag,1010) ' restore_bgc               = ', restore_bgc
         write(nu_diag,*)    ' bgc_data_dir              = ', &
                               trim(bgc_data_dir)
         write(nu_diag,*)    ' sil_data_type             = ', &
                               trim(sil_data_type)
         write(nu_diag,*)    ' nit_data_type             = ', &
                               trim(nit_data_type)
         write(nu_diag,*)    ' fe_data_type              = ', &
                               trim(fe_data_type)
         write(nu_diag,1020) ' number of bio tracers     = ', nbtrcr
         write(nu_diag,1020) ' number of Isw tracers     = ', nbtrcr_sw
         write(nu_diag,1020) ' number of autotrophs      = ', n_algae
         write(nu_diag,1020) ' number of doc          = ', n_doc
         write(nu_diag,1020) ' number of dic          = ', n_dic
         write(nu_diag,1020) ' number of don          = ', n_don
         write(nu_diag,1020) ' number of fed          = ', n_fed
         write(nu_diag,1020) ' number of fep          = ', n_fep
         write(nu_diag,1010) ' tr_bgc_N               = ', tr_bgc_N
         write(nu_diag,1010) ' tr_bgc_C               = ', tr_bgc_C
         write(nu_diag,1010) ' tr_bgc_chl             = ', tr_bgc_chl
         write(nu_diag,1010) ' tr_bgc_Nit             = ', tr_bgc_Nit
         write(nu_diag,1010) ' tr_bgc_Am              = ', tr_bgc_Am
         write(nu_diag,1010) ' tr_bgc_Sil             = ', tr_bgc_Sil
         write(nu_diag,1010) ' tr_bgc_DMS             = ', tr_bgc_DMS
         write(nu_diag,1010) ' tr_bgc_PON             = ', tr_bgc_PON
         write(nu_diag,1010) ' tr_bgc_DON             = ', tr_bgc_DON
         write(nu_diag,1010) ' tr_bgc_Fe              = ', tr_bgc_Fe 
        
      endif   ! master_task

      !----------------------------------------------------
      ! zbgc
      !----------------------------------------------------
      elseif (z_tracers) then ! defined on nblyr+1 in ice 
                              !and 2 snow layers (snow surface + interior)

      if (tr_bgc_N) then
        do mm = 1,n_algae      
           nt_bgc_N(mm) = ntrcr + 1 
           ntrcr = ntrcr + nblyr + 3
           nbtrcr = nbtrcr + 1
           nlt_bgc_N(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)  = algaltype(mm)
           do k = 1,nblyr+1
               trcr_depend(nt_bgc_N(mm) + k - 1)  = 2+nt_fbri+ntd
           enddo
           ! snow volume
           trcr_depend(nt_bgc_N(mm) +  nblyr + 1:nt_bgc_N(mm) +  nblyr + 2) = 2 
           bio_index(nlt_bgc_N(mm)) = nt_bgc_N(mm)
           bio_index_o(nlt_bgc_N(mm)) = mm
        enddo   !mm
        if (dEdd_algae) then
           nlt_chl_sw = 1
           nbtrcr_sw =  nbtrcr_sw + nilyr+nslyr+2
        endif
      endif ! tr_bgc_N

      if (tr_bgc_Nit) then
        nt_bgc_Nit = ntrcr + 1
        ntrcr = ntrcr + nblyr+ 3
        nbtrcr = nbtrcr + 1
        nlt_bgc_Nit = nbtrcr
        bgc_tracer_type(nbtrcr) = nitratetype
        do k = 1,nblyr+1
            trcr_depend(nt_bgc_Nit + k - 1)  = 2+nt_fbri+ntd
        enddo
        trcr_depend(nt_bgc_Nit + nblyr+1:nt_bgc_Nit + nblyr+2)    = 2 
        bio_index(nlt_bgc_Nit) = nt_bgc_Nit
        bio_index_o(nlt_bgc_Nit) = max_algae + 1
      endif
         
      if (tr_bgc_C) then
       !
       ! Algal C is not yet distinct from algal N
       ! * Reqires  exudation and/or changing C:N ratios
       ! for implementation
       !
       ! do mm = 1,n_algae         
       !    nt_bgc_C(mm) = ntrcr + 1  
       !    ntrcr = ntrcr + nblyr + 3
       !    nbtrcr = nbtrcr + 1
       !    nlt_bgc_C(mm) = nbtrcr
       !    bgc_tracer_type(nbtrcr)  = algaltype(mm)
       !    do k = 1,nblyr+1
       !        trcr_depend(nt_bgc_C(mm) + k - 1)  = 2+nt_fbri+ntd
       !    enddo
       !    trcr_depend(nt_bgc_C(mm) +  nblyr + 1:nt_bgc_C(mm) +  nblyr +2) = 2 
       !    bio_index(nlt_bgc_C(mm)) = nt_bgc_C(mm)
       !    bio_index_o(nlt_bgc_C(mm)) = max_algae + 1 + mm
       ! enddo   !mm
        do mm = 1,n_doc         
           nt_bgc_DOC(mm) = ntrcr + 1 
           ntrcr = ntrcr + nblyr + 3
           nbtrcr = nbtrcr + 1
           nlt_bgc_DOC(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)  = doctype(mm)
           do k = 1,nblyr+1
               trcr_depend(nt_bgc_DOC(mm) + k - 1)  = 2+nt_fbri+ntd
           enddo
           trcr_depend(nt_bgc_DOC(mm) +  nblyr + 1:nt_bgc_DOC(mm) +  nblyr +2) = 2 
           bio_index(nlt_bgc_DOC(mm)) = nt_bgc_DOC(mm)
           bio_index_o(nlt_bgc_DOC(mm)) = max_algae + 1 + mm
        enddo   !mm
        do mm = 1,n_dic         
           nt_bgc_DIC(mm) = ntrcr + 1 
           ntrcr = ntrcr + nblyr + 3
           nbtrcr = nbtrcr + 1
           nlt_bgc_DIC(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)  = dictype(mm)
           do k = 1,nblyr+1
               trcr_depend(nt_bgc_DIC(mm) + k - 1)  = 2+nt_fbri+ntd
           enddo
           trcr_depend(nt_bgc_DIC(mm) +  nblyr + 1:nt_bgc_DIC(mm) +  nblyr +2) = 2 
           bio_index(nlt_bgc_DIC(mm)) = nt_bgc_DIC(mm)
           bio_index_o(nlt_bgc_DIC(mm)) = max_algae + max_doc + 1 + mm
        enddo   !mm
      endif !tr_bgc_C

      if (tr_bgc_chl) then
        do mm = 1,n_algae      
           nt_bgc_chl(mm) = ntrcr + 1 
           ntrcr = ntrcr + nblyr + 3 
           nbtrcr = nbtrcr + 1
           nlt_bgc_chl(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)  = algaltype(mm)
           do k = 1,nblyr+1
               trcr_depend(nt_bgc_chl(mm) + k - 1)  = 2+nt_fbri+ntd
           enddo
           trcr_depend(nt_bgc_chl(mm) +  nblyr + 1:nt_bgc_chl(mm) +  nblyr +2) = 2 
           bio_index(nlt_bgc_chl(mm)) = nt_bgc_chl(mm)
           bio_index_o(nlt_bgc_chl(mm)) = max_algae + 1 + max_doc + max_dic + mm
        enddo   !mm
      endif !tr_bgc_chl

      if (tr_bgc_Am)then
          nt_bgc_Am = ntrcr + 1
          ntrcr = ntrcr + nblyr+ 3
          nbtrcr = nbtrcr + 1
          nlt_bgc_Am = nbtrcr
          bgc_tracer_type(nbtrcr) = ammoniumtype
          do k = 1,nblyr+1
             trcr_depend(nt_bgc_Am    + k - 1)  = 2+nt_fbri+ntd
          enddo  
          trcr_depend(nt_bgc_Am + nblyr+1:nt_bgc_Am+nblyr+2)    = 2 
          bio_index(nlt_bgc_Am) = nt_bgc_Am
          bio_index_o(nlt_bgc_Am) = 2*max_algae + max_doc + max_dic + 2
      endif    
      if (tr_bgc_Sil)then
          nt_bgc_Sil = ntrcr + 1
          ntrcr = ntrcr + nblyr+3
          nbtrcr = nbtrcr + 1
          nlt_bgc_Sil = nbtrcr
          bgc_tracer_type(nbtrcr) = silicatetype
          do k = 1,nblyr+1
             trcr_depend(nt_bgc_Sil    + k - 1)  = 2+nt_fbri+ntd
          enddo
          trcr_depend(nt_bgc_Sil + nblyr+1:nt_bgc_Sil+nblyr+2)    = 2 
          bio_index(nlt_bgc_Sil) = nt_bgc_Sil
          bio_index_o(nlt_bgc_Sil) = 2*max_algae+ max_doc + max_dic + 3
      endif    
      if (tr_bgc_DMS)then   !all together
          nt_bgc_DMSPp = ntrcr + 1
          ntrcr = ntrcr + nblyr+3
          nbtrcr = nbtrcr + 1
          nlt_bgc_DMSPp = nbtrcr
          bgc_tracer_type(nbtrcr) = dmspptype
          do k = 1,nblyr+1
             trcr_depend(nt_bgc_DMSPp    + k - 1)  = 2+nt_fbri+ntd
          enddo
          trcr_depend(nt_bgc_DMSPp + nblyr+1:nt_bgc_DMSPp+nblyr+2)  = 2
          bio_index(nlt_bgc_DMSPp) = nt_bgc_DMSPp
          bio_index_o(nlt_bgc_DMSPp) = 2*max_algae+ max_doc + max_dic + 4
          nt_bgc_DMSPd = ntrcr + 1
          ntrcr = ntrcr + nblyr+3
          nbtrcr = nbtrcr + 1
          nlt_bgc_DMSPd = nbtrcr
          bgc_tracer_type(nbtrcr) = dmspdtype
          do k = 1,nblyr+1
             trcr_depend(nt_bgc_DMSPd    + k - 1)  = 2+nt_fbri+ntd
          enddo
          trcr_depend(nt_bgc_DMSPd + nblyr+1:nt_bgc_DMSPd + nblyr+2)  = 2 
          bio_index(nlt_bgc_DMSPd) = nt_bgc_DMSPd
          bio_index_o(nlt_bgc_DMSPd) = 2*max_algae+ max_doc + max_dic + 5
          nt_bgc_DMS = ntrcr + 1
          ntrcr = ntrcr + nblyr+3
          nbtrcr = nbtrcr + 1
          nlt_bgc_DMS = nbtrcr
          bgc_tracer_type(nbtrcr) = dmspdtype
          do k = 1,nblyr+1
             trcr_depend(nt_bgc_DMS    + k - 1)  = 2+nt_fbri+ntd
          enddo
          trcr_depend(nt_bgc_DMS + nblyr+1:nt_bgc_DMS + nblyr+2)  = 2 
          bio_index(nlt_bgc_DMS) = nt_bgc_DMS
          bio_index_o(nlt_bgc_DMS) = 2*max_algae+ max_doc + max_dic + 6
      endif    
      if (tr_bgc_PON)then
          nt_bgc_PON = ntrcr + 1
          ntrcr = ntrcr + nblyr+3
          nbtrcr = nbtrcr + 1
          nlt_bgc_PON = nbtrcr
          bgc_tracer_type(nbtrcr) = nitratetype
          do k = 1,nblyr+1
             trcr_depend(nt_bgc_PON    + k - 1)  = 2+nt_fbri+ntd
          enddo
          trcr_depend(nt_bgc_PON + nblyr+1:nt_bgc_PON + nblyr+2)  = 2 
          bio_index(nlt_bgc_PON) = nt_bgc_PON
          bio_index_o(nlt_bgc_PON) =  2*max_algae + max_doc + max_dic + 7
      endif
      if (tr_bgc_DON)then
        do mm = 1,n_don         
           nt_bgc_DON(mm) = ntrcr + 1  
           ntrcr = ntrcr + nblyr + 3
           nbtrcr = nbtrcr + 1
           nlt_bgc_DON(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)  = dontype(mm)
           do k = 1,nblyr+1
               trcr_depend(nt_bgc_DON(mm) + k - 1)  = 2+nt_fbri+ntd
           enddo
           trcr_depend(nt_bgc_DON(mm) +  nblyr + 1:nt_bgc_DON(mm) +  nblyr + 2) = 2 
           bio_index(nlt_bgc_DON(mm)) = nt_bgc_DON(mm)
           bio_index_o(nlt_bgc_DON(mm)) =  2*max_algae+ max_doc + max_dic + 7+ mm
        enddo   !mm
      endif   !tr_bgc_DON
      if (tr_bgc_Fe )then
        do mm = 1,n_fed          
           nt_bgc_Fed(mm) = ntrcr + 1  
           ntrcr = ntrcr + nblyr + 3
           nbtrcr = nbtrcr + 1
           nlt_bgc_Fed(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)  =  fedtype(mm)
           do k = 1,nblyr+1
               trcr_depend(nt_bgc_Fed (mm) + k - 1)  = 2+nt_fbri+ntd
           enddo
           trcr_depend(nt_bgc_Fed (mm) +  nblyr + 1:nt_bgc_Fed (mm) +  nblyr +2) = 2 
           bio_index(nlt_bgc_Fed (mm)) = nt_bgc_Fed (mm)
           bio_index_o(nlt_bgc_Fed (mm)) =  2*max_algae+ max_doc + max_dic + max_don + 7 + mm
        enddo   !mm
        do mm = 1,n_fep          
           nt_bgc_Fep(mm) = ntrcr + 1  
           ntrcr = ntrcr + nblyr + 3
           nbtrcr = nbtrcr + 1
           nlt_bgc_Fep(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)  =  feptype(mm)
           do k = 1,nblyr+1
               trcr_depend(nt_bgc_Fep (mm) + k - 1)  = 2+nt_fbri+ntd
           enddo
           trcr_depend(nt_bgc_Fep (mm) +  nblyr + 1:nt_bgc_Fep (mm) +  nblyr +2) = 2 
           bio_index(nlt_bgc_Fep (mm)) = nt_bgc_Fep (mm)
           bio_index_o(nlt_bgc_Fep (mm)) =  2*max_algae+ max_doc + max_dic + max_don + max_fe + 7 + mm
        enddo   !mm
      endif   !tr_bgc_Fe 

     !-----------------------------------------
     !  z layer aerosols
     !-----------------------------------------

      if (tr_zaero)then
        do mm = 1,n_zaero      
           nt_zaero(mm) = ntrcr + 1 
           ntrcr = ntrcr + nblyr + 3 
           nbtrcr = nbtrcr + 1
           nlt_zaero(mm) = nbtrcr
           bgc_tracer_type(nbtrcr)  = zaerotype(mm)
           do k = 1,nblyr+1
               trcr_depend(nt_zaero(mm) + k - 1)  = 2+nt_fbri+ntd
           enddo
           trcr_depend(nt_zaero(mm) +  nblyr + 1:nt_zaero(mm) +  nblyr + 2) = 2 
           bio_index(nlt_zaero(mm)) = nt_zaero(mm)
           bio_index_o(nlt_zaero(mm)) =  2*max_algae + max_doc + max_dic + max_don + 2*max_fe + 7 + mm
           if (dEdd_algae) then
                nlt_zaero_sw(mm) = nbtrcr_sw + 1
                nbtrcr_sw = nbtrcr_sw + nilyr + nslyr+2
           endif
        enddo   !mm
      endif   !tr_zaero

      nt_zbgc_frac = 0
      if (nbtrcr > 0) then
          nt_zbgc_frac = ntrcr + 1
          ntrcr = ntrcr + nbtrcr
          do k = 1,nbtrcr 
            zbgc_frac_init(k) = c1
            trcr_depend(nt_zbgc_frac+k-1) = 2+nt_fbri+ntd  ! on brine volume
            tau_ret(k) = c1
            tau_rel(k) = c1
            if (bgc_tracer_type(k) >=  c0 .and. bgc_tracer_type(k) < p5) then
              tau_ret(k) = tau_min
              tau_rel(k) = tau_max
              zbgc_frac_init(k) = c1
            elseif (bgc_tracer_type(k) >= p5 .and. bgc_tracer_type(k) < c1) then
              tau_ret(k) = tau_min
              tau_rel(k) = tau_min
              zbgc_frac_init(k) = c1
            elseif (bgc_tracer_type(k) >= c1 .and. bgc_tracer_type(k) < c2) then
              tau_ret(k) = tau_max
              tau_rel(k) = tau_min
              zbgc_frac_init(k) = c1
            elseif (bgc_tracer_type(k) >= c2 ) then
              tau_ret(k) = tau_max
              tau_rel(k) = tau_max
              zbgc_frac_init(k) = c1
            endif
          enddo
     endif

      if (my_task == master_task) then

         write(nu_diag,*)    ' sil_data_type             = ', &
                               trim(sil_data_type)
         write(nu_diag,*)    ' nit_data_type             = ', &
                               trim(nit_data_type)
         write(nu_diag,*)    ' fe_data_type              = ', &
                               trim(fe_data_type)
         write(nu_diag,*)    ' bgc_data_dir              = ', &
                               trim(bgc_data_dir)
         write(nu_diag,1010) ' restart_bgc               = ', restart_bgc
         write(nu_diag,1010) ' scale_bgc                 = ', scale_bgc
         write(nu_diag,1010) ' solve_zbgc                = ', solve_zbgc
         write(nu_diag,1020) ' number of ztracers        = ', nbtrcr
         write(nu_diag,1020) ' number of Isw tracers     = ', nbtrcr_sw
         write(nu_diag,1020) ' number of autotrophs      = ', n_algae
         write(nu_diag,1020) ' number of doc             = ', n_doc
         write(nu_diag,1020) ' number of dic             = ', n_dic
         write(nu_diag,1020) ' number of fed             = ', n_fed
         write(nu_diag,1020) ' number of fep             = ', n_fep
         write(nu_diag,1020) ' number of aerosols        = ', n_zaero
         write(nu_diag,1010) ' tr_zaero                  = ', tr_zaero
         write(nu_diag,1010) ' tr_bgc_Nit                = ', tr_bgc_Nit
         write(nu_diag,1010) ' tr_bgc_N                  = ', tr_bgc_N
         write(nu_diag,1010) ' tr_bgc_Am                 = ', tr_bgc_Am
         write(nu_diag,1010) ' tr_bgc_C                  = ', tr_bgc_C
         write(nu_diag,1010) ' tr_bgc_Sil                = ', tr_bgc_Sil
         write(nu_diag,1010) ' tr_bgc_chl                = ', tr_bgc_chl
         write(nu_diag,1010) ' tr_bgc_DMS                = ', tr_bgc_DMS
         write(nu_diag,1010) ' tr_bgc_PON                = ', tr_bgc_PON
         write(nu_diag,1010) ' tr_bgc_DON                = ', tr_bgc_DON
         write(nu_diag,1010) ' tr_bgc_Fe                 = ', tr_bgc_Fe 
         !bio parameters
         write(nu_diag,1000) ' grid_o                    = ', grid_o
         write(nu_diag,1000) ' grid_o_t                  = ', grid_o_t
         write(nu_diag,1005) ' l_sk                      = ', l_sk
         write(nu_diag,1000) ' initbio_frac              = ', initbio_frac
         write(nu_diag,1000) ' frazil_scav               = ', frazil_scav  

      endif   ! master_task
      endif  ! skl_bgc or solve_bgc

      if (nbtrcr > max_nbtrcr) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'nbtrcr > max_nbtrcr'
         write (nu_diag,*) 'nbtrcr, max_nbtrcr:',nbtrcr, max_nbtrcr
         call abort_ice ('ice: ice_zbgc error')
      endif
      if (.NOT. dEdd_algae) nbtrcr_sw = 1

      if (nbtrcr_sw > max_nsw) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'nbtrcr_sw > max_nsw'
         write (nu_diag,*) 'nbtrcr_sw, max_nsw:',nbtrcr_sw, max_nsw
         call abort_ice ('ice: ice_zbgc error')
      endif

      do k = 1, nbtrcr
         zbgc_init_frac(k) = frazil_scav
         if (bgc_tracer_type(k) < c0)  zbgc_init_frac(k) = initbio_frac
      enddo  

      if (ntrcr > max_ntrcr) then
         write(nu_diag,*) 'max_ntrcr < number of namelist tracers'
         write(nu_diag,*) 'max_ntrcr = ',max_ntrcr,' ntrcr = ',ntrcr
         call abort_ice('max_ntrcr < number of namelist tracers')
      endif                               

      if (skl_bgc .and. TRBGCS < 2) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of bgc tracers >= 2'
         write (nu_diag,*) 'number of bgc tracers compiled:',TRBGCS
         call abort_ice ('ice: ice_zbgc error')
      endif

      if (solve_zbgc .and. TRBGCZ < 2) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of zbgc tracers >= 2'
         write (nu_diag,*) 'number of bgc tracers compiled:',TRBGCZ
         call abort_ice ('ice: ice_zbgc error')
      endif

      if (tr_zaero .and. TRZAERO <  1) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of TRZAERO > 0'
         write (nu_diag,*) 'in order to solve z aerosols:',TRBGCZ
         call abort_ice ('ice: ice_zbgc error')
      endif
   
 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1005    format (a30,2x,f9.6)  ! float
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1030    format (a30,   a8)    ! character

      end subroutine init_zbgc

!=======================================================================

!  Initialize vertical profile of biogeochemistry

      subroutine init_bgc
      use ice_algae, only: read_restart_bgc, init_bgc_data
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c1, c0, c10, c5, p15, p5, p01, &
              field_type_scalar, field_loc_center, p1, c2
      use ice_domain, only: nblocks
      use ice_domain_size, only:  nblyr, nilyr,  &
             max_algae, max_don, max_doc, max_dic, max_aero, max_fe
      use ice_flux, only: sss
      use ice_calendar, only: month, dt
      use ice_shortwave, only: hs_ssl
      use ice_state, only: trcrn, aicen, vicen, vsnon, nbtrcr, &
              nt_fbri, nt_bgc_S, nbtrcr_sw, nt_sice, ntrcr, &
              nt_zbgc_frac
      use ice_therm_shared, only: solve_zsal, ktherm
 
      integer (kind=int_kind) :: &
         i, j, iblk       , & ! horizontal indices
         ixm,ixp          , & ! record numbers for neighboring months
         k                , & ! vertical index 
         n                , & ! category index 
         mm               , & ! bio tracer index
         ki               , & ! loop index
         ks, ksp

      logical (kind=log_kind) :: &
         dbug       , &      ! prints debugging output if true
         readm

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: & 
         trtmp           ! temporary, remapped tracers   
      
      real (kind=dbl_kind), dimension (nblyr+1) :: &
         zspace

      real (kind=dbl_kind) :: & 
         dvssl      , & ! volume of snow surface layer (m)
         dvint          ! volume of snow interior      (m)

      character (char_len_long) ::  &
         fieldname      !  netcdf fieldname

      dbug = .true.

      zfswin(:,:,:,:,:)     = c0   !shortwave flux on bio grid
      atm_bio_all(:,:,:,:)  = c0  
      ocean_bio_all(:,:,:,:) = c0
      ice_bio_net(:,:,:,:)  = c0   ! integrated ice tracer concentration (mmol/m^2 or mg/m^2) 
      snow_bio_net(:,:,:,:) = c0   ! integrated ice tracer concentration (mmol/m^2 or mg/m^2)
      zspace(:) = c1/real(nblyr,kind=dbl_kind)
      zspace(1) = p5*zspace(1)
      zspace(nblyr+1) = p5*zspace(nblyr+1)
      trcrn_sw(:,:,:,:,:) =c0      ! tracers active in the shortwave calculation

      !-----------------------------------------------------------------------------   
      !     Atmosphere Concentrations or defined below
      !----------------------------------------------------------------------------- 
         ! bgc tracers = c0
         atm_bio_all(:,:,:,:) = c0 ! 1:2*max_algae + max_doc + 7 + max_dic + max_don

         ! zaero  (mmol/m^3)
         atm_bio_all(:,:,2*max_algae + max_doc + 8 + max_dic + max_don + 2*max_fe,:) &
                    = 1.e-11_dbl_kind*dt/grid_o_t
         atm_bio_all(:,:,2*max_algae + max_doc + 9 + max_dic + max_don + 2*max_fe,:) &
                     = 1.e-13_dbl_kind*dt/grid_o_t
         atm_bio_all(:,:,2*max_algae + max_doc + 10 + max_dic + max_don + 2*max_fe,:) &
                     = 1.e-15_dbl_kind*dt/grid_o_t
         atm_bio_all(:,:,2*max_algae + max_doc + 11 + max_dic + max_don + 2*max_fe,:)= p15
         atm_bio_all(:,:,2*max_algae + max_doc + 12 + max_dic + max_don + 2*max_fe,:)= p15
         atm_bio_all(:,:,2*max_algae + max_doc + 13 + max_dic + max_don + 2*max_fe,:)= p15
       
      if (restart_bgc) then       
     
         call read_restart_bgc  

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,k)
      do iblk = 1, max_blocks
        do j = 1, ny_block
        do i = 1, nx_block   
               do k = 1, max_algae         
                  ocean_bio_all(i,j,k,iblk)      = algalN(i,j,k,iblk)           !N
                  ks = max_algae + 1
                  ocean_bio_all(i,j,ks + k,iblk) = R_C2N(k)*algalN(i,j,k,iblk)  !C
                  ks = max_algae + max_doc + max_dic + 1
                  ocean_bio_all(i,j,ks + k,iblk) = R_chl2N(k)*algalN(i,j,k,iblk)!chl
               enddo   
               ks = max_algae + 1
               do k = 1, max_doc
                  ocean_bio_all(i,j,ks + k,iblk) = doc(i,j,k,iblk)              !doc
               enddo  
               ks = max_algae + max_doc + 1
               do k = 1, max_dic
                  ocean_bio_all(i,j,ks + k,iblk) = dic(i,j,k,iblk)              !dic
               enddo 
               ks = 2*max_algae + max_doc + 7 + max_dic
               do k = 1, max_don
                  ocean_bio_all(i,j,ks + k,iblk) = don(i,j,k,iblk)              !don
               enddo                    
               ks = max_algae + 1
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit
               ks = 2*max_algae + max_doc + 2 + max_dic
               ocean_bio_all(i,j,ks,iblk) = amm(i,j,iblk)                       !Am
               ks = 2*max_algae + max_doc + 3 + max_dic
               ocean_bio_all(i,j,ks,iblk) = sil(i,j,iblk)                       !Sil
               ks =  2*max_algae + max_doc + 4 + max_dic
               ocean_bio_all(i,j,ks,iblk) = R_S2N(1)*algalN(i,j,1,iblk) + &    !DMSPp
                                            R_S2N(2)*algalN(i,j,2,iblk) + &
                                            R_S2N(3)*algalN(i,j,3,iblk) 
               ks =  2*max_algae + max_doc + 5 + max_dic                        
               ocean_bio_all(i,j,ks,iblk) = dmsp(i,j,iblk)                      !DMSPd
               ks =  2*max_algae + max_doc + 6 + max_dic
               ocean_bio_all(i,j,ks,iblk) = dms(i,j,iblk)                       !DMS
               ks =  2*max_algae + max_doc + 7 + max_dic
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON
               ks = 2*max_algae + max_doc + 7 + max_dic + max_don
               do k = 1, max_fe
                  ocean_bio_all(i,j,ks + k,iblk) = fed(i,j,k,iblk)               !fed
               enddo        
               ks = 2*max_algae + max_doc + 7 + max_dic + max_don + max_fe            
               do k = 1, max_fe
                  ocean_bio_all(i,j,ks + k,iblk) = fep(i,j,k,iblk)               !fep
               enddo                    
               ks = 2*max_algae + max_doc + 7 + max_dic + max_don + 2*max_fe
               do k = 1, max_aero
                  ocean_bio_all(i,j,ks+k,iblk) = c0                             !zaero
               enddo
         enddo !i
         enddo !j
        if (.not. skl_bgc) then
         do  n = 1,ncat
         do mm = 1,nbtrcr
            do k = 1, nblyr+1
               do j = 1, ny_block
               do i = 1, nx_block
                  ice_bio_net(i,j,mm,iblk) = ice_bio_net(i,j,mm,iblk) + &
                               trcrn(i,j,bio_index(mm)+k-1, n, iblk)* &
                               vicen(i,j,n,iblk)*zspace(k)*trcrn(i,j,nt_fbri,n,iblk)
               enddo   !i
               enddo   !j
             enddo     !k
           do j = 1, ny_block
           do i = 1, nx_block
               dvssl  = min(p5*vsnon(i,j,n,iblk), hs_ssl*aicen(i,j,n,iblk))   !snow surface layer
               dvint  = vsnon(i,j,n,iblk)- dvssl                              !snow interior
               snow_bio_net(i,j,mm,iblk) = snow_bio_net(i,j,mm,iblk) +  &
                           trcrn(i,j,bio_index(mm)+nblyr+1,n,iblk)*dvssl +  &
                           trcrn(i,j,bio_index(mm)+nblyr+2,n,iblk)*dvint
           enddo !i
           enddo !j
         enddo     !mm
         enddo     !n 
         endif
      enddo        !iblk
      !$OMP END PARALLEL DO
              
      else  !restart
 
      !-----------------------------------------------------------------------------   
      !   Initial Ocean Values
      !-----------------------------------------------------------------------------
      call init_bgc_data(fed(:,:,1,:),fep(:,:,1,:)) ! input dFe from file
      ki = 1
      if (trim(fe_data_type) == 'clim') ki = 2

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,k)
      do iblk = 1, max_blocks
        do j = 1, ny_block
        do i = 1, nx_block   
               amm(i,j,iblk)  = c1 !ISPOL < 1 mmol/m^3 
               dmsp(i,j,iblk) = p1  
               dms = p1    
               algalN(i,j,1,iblk) = 0.0026_dbl_kind ! ISPOL, Lannuzel 2013(pennate) 
               algalN(i,j,2,iblk) = 0.0027_dbl_kind ! ISPOL, Lannuzel 2013(Phaeocystis)
               algalN(i,j,3,iblk) = 0.0057_dbl_kind ! ISPOL, Lannuzel 2013(flagellates)  
                                    !0.024_dbl_kind ! 5% of 1 mgchl/m^3 
 
             ! possibly redefined in get_forcing_bgc  
               nit(i,j,iblk) = 15.0_dbl_kind      
               sil(i,j,iblk) = 25.0_dbl_kind      
               do k = 1, max_algae           
                  ocean_bio_all(i,j,k,iblk)      = algalN(i,j,k,iblk)           !N
                  ks = max_algae + 1
                  ocean_bio_all(i,j,ks + k,iblk) = R_C2N(k)*algalN(i,j,k,iblk)  !C
                  ks = max_algae + max_doc + max_dic + 1
                  ocean_bio_all(i,j,ks + k,iblk) = R_chl2N(k)*algalN(i,j,k,iblk)!chl
               enddo   
               ks = max_algae + 1
              
               doc(i,j,1,iblk) = 16.2_dbl_kind ! 18% saccharides
               doc(i,j,2,iblk) = 9.0_dbl_kind  ! 
               doc(i,j,3,iblk) = c1 ! 
               do k = 1, max_doc
                  ocean_bio_all(i,j,ks + k,iblk) = doc(i,j,k,iblk)              !doc
               enddo  
               ks = max_algae + max_doc + 1
               do k = 1, max_dic
                  dic(i,j,k,iblk) = c1
                  ocean_bio_all(i,j,ks + k,iblk) = dic(i,j,k,iblk)              !dic
               enddo 
               ks = 2*max_algae+ max_doc + max_dic + 7
               do k = 1, max_don
                  don(i,j,k,iblk) = 12.9_dbl_kind              
                 ! 64.3_dbl_kind ! 72% Total DOC~90 mmolC/m^3  ISPOL with N:C of 0.2
                  ocean_bio_all(i,j,ks + k,iblk) = don(i,j,k,iblk)              !don
               enddo  
               ks = max_algae + 1
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit
               ks = 2*max_algae + max_doc + 2 + max_dic
               ocean_bio_all(i,j,ks,iblk) = amm(i,j,iblk)                       !Am
               ks = 2*max_algae + max_doc + 3 + max_dic
               ocean_bio_all(i,j,ks,iblk) = sil(i,j,iblk)                       !Sil
               ks =  2*max_algae + max_doc + 4 + max_dic
               ocean_bio_all(i,j,ks,iblk) =  R_S2N(1)*algalN(i,j,1,iblk) + &    !DMSPp
                                             R_S2N(2)*algalN(i,j,2,iblk) + &
                                             R_S2N(3)*algalN(i,j,3,iblk) 
               ks =  2*max_algae + max_doc + 5 + max_dic     
               ocean_bio_all(i,j,ks,iblk) = dmsp(i,j,iblk)                      !DMSPd
               ks =  2*max_algae + max_doc + 6 + max_dic
               ocean_bio_all(i,j,ks,iblk) = dms(i,j,iblk)                       !DMS
               ks =  2*max_algae + max_doc + 7 + max_dic
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON
               ks = 2*max_algae+ max_doc + max_dic + 7 + max_don
               ocean_bio_all(i,j,ks + 1,iblk) = fed(i,j,1,iblk)                 !fed
               do k = ki, max_fe
                  fed(i,j,k,iblk) = c1 !  (nM)  Lannuzel2007 DFe, 
                                       !range 0.14-2.6 (nM) van der Merwe 2011
                  ocean_bio_all(i,j,ks + k,iblk) = fed(i,j,k,iblk)               !fed
               enddo  
               ks = 2*max_algae+ max_doc + max_dic + 7 + max_don + max_fe
               ocean_bio_all(i,j,ks+1,iblk) = fep(i,j,1,iblk)
               do k = ki, max_fe  
                  fep(i,j,k,iblk) = c2               ! (nM) van der Merwe 2011
                                                     ! (0.6 to 2.9 nM ocean)
                  ocean_bio_all(i,j,ks + k,iblk) = fep(i,j,k,iblk)               !fep
               enddo  
               ks = 2*max_algae + max_doc + 7 + max_dic + max_don + 2*max_fe
               do k = 1, max_aero
                  zaeros(i,j,k,iblk) = c0
                  ocean_bio_all(i,j,ks+k,iblk) = zaeros(i,j,k,iblk)             !zaero
               enddo
         enddo !i
         enddo !j
      enddo        !iblk

      !$OMP END PARALLEL DO
        
      !-----------------------------------------------------------------------------   
      !     Skeletal Layer Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The skeletal layer model assumes a constant 
      !  layer depth (sk_l) and porosity (phi_sk)
      !-----------------------------------------------------------------------------   
      if (skl_bgc) then
       
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, max_blocks
         do  n = 1,ncat
         do mm = 1,nbtrcr
               do j = 1, ny_block
               do i = 1, nx_block
                                   ! bulk concentration (mmol or mg per m^3, or 10^-3 mmol/m^3)
               trcrn(i,j,bio_index(mm), n, iblk) = ocean_bio_all(i,j,bio_index_o(mm),iblk)
              enddo  !i
              enddo  ! j
         enddo       ! nbtrcr
         enddo       ! n 
         enddo       ! iblk
         !$OMP END PARALLEL DO

      !-----------------------------------------------------------------------------   
      !    zbgc Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The vertical layer model uses prognosed porosity and layer depth
      !-----------------------------------------------------------------------------   

      else   !not skl_bgc
       if (scale_bgc .and. solve_zsal) then      ! bulk concentration (mmol or mg per m^3)
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, max_blocks
         do  n = 1,ncat
         do mm = 1,nbtrcr
            do k = 2, nblyr
               do j = 1, ny_block
               do i = 1, nx_block
                trcrn(i,j,bio_index(mm)+k-1,n,iblk) = &
                                   (p5*(trcrn(i,j,nt_bgc_S+k-1,n,iblk)+ trcrn(i,j,nt_bgc_S+k-2,n,iblk)) &
                                    /sss(i,j,iblk)*ocean_bio_all(i,j,bio_index_o(mm),iblk)) 
                ice_bio_net(i,j,mm,iblk) = ice_bio_net(i,j,mm,iblk) + &
                               trcrn(i,j,bio_index(mm)+k-1, n, iblk)* &
                               vicen(i,j,n,iblk)*zspace(k)*trcrn(i,j,nt_fbri,n,iblk)
               enddo        !i
               enddo        !j
            enddo  !k
            do j = 1, ny_block
            do i = 1, nx_block
                trcrn(i,j,nt_zbgc_frac-1+mm,n,iblk) = zbgc_frac_init(mm)
                trcrn(i,j,bio_index(mm),n,iblk) = (trcrn(i,j,nt_bgc_S,n,iblk)&
                                           /sss(i,j,iblk)*ocean_bio_all(i,j,bio_index_o(mm),iblk)) 
                trcrn(i,j,bio_index(mm)+nblyr,n,iblk)=  (trcrn(i,j,nt_bgc_S+nblyr-1,n,iblk)&
                                                /sss(i,j,iblk)*ocean_bio_all(i,j,bio_index_o(mm),iblk))
                trcrn(i,j,bio_index(mm)+nblyr+1:bio_index(mm)+nblyr+2,n,iblk) = c0 !snow concentration
                ice_bio_net(i,j,mm,iblk) = ice_bio_net(i,j,mm,iblk) + &
                               (trcrn(i,j,bio_index(mm), n, iblk)+ &
                                trcrn(i,j,bio_index(mm)+nblyr,n,iblk))*&
                                vicen(i,j,n,iblk)*zspace(1)*trcrn(i,j,nt_fbri,n,iblk)
             enddo        !i
             enddo        !j
         enddo           !mm
         enddo           !n 
         enddo           !iblk 
         !$OMP END PARALLEL DO
    
      elseif (scale_bgc .and. ktherm == 2) then
         trtmp(:,:,:) = c0
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, max_blocks
         do  n = 1,ncat        
          do j = 1, ny_block
          do i = 1, nx_block    
              call remap_zbgc(ntrcr,          nilyr,  &
                         nt_sice,                        &
                         trcrn(i,j,:,n,iblk), trtmp(i,j,:),&
                         0,                    nblyr+1,  &
                         c1,                   c1,       &
                         cgrid(2:nilyr+1),               &
                         igrid(1:nblyr+1),               &
                         trcrn(i,j,nt_sice,n,iblk))

           enddo  !i
           enddo  !j
           do mm = 1,nbtrcr
             do k = 1, nblyr + 1
               do j = 1, ny_block
               do i = 1, nx_block              
                 trcrn(i,j,bio_index(mm)+k-1,n,iblk) =   &
                 (trtmp(i,j,nt_sice+k-1)/sss(i,j,iblk)*ocean_bio_all(i,j,bio_index_o(mm),iblk))
                 ice_bio_net(i,j,mm,iblk) = ice_bio_net(i,j,mm,iblk) + &
                                trcrn(i,j,bio_index(mm)+k-1, n, iblk)* &
                                vicen(i,j,n,iblk)*zspace(k)*trcrn(i,j,nt_fbri,n,iblk)
                 trcrn(i,j,bio_index(mm)+nblyr+1:bio_index(mm)+nblyr+2,n,iblk) = c0 !snow concentration
               enddo        !i
               enddo        !j
            enddo  !k
          enddo           !mm
         enddo           !n 
         enddo           !iblk 

      else !  .not. scale_bgc         
     
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, max_blocks
         do  n = 1,ncat
         do mm = 1,nbtrcr
            do k = 1, nblyr+1
               do j = 1, ny_block
               do i = 1, nx_block
                  trcrn(i,j,bio_index(mm)+k-1,n,iblk) =  ocean_bio_all(i,j,bio_index_o(mm),iblk)*&
                               zbgc_init_frac(mm) 
                  ice_bio_net(i,j,mm,iblk) = ice_bio_net(i,j,mm,iblk) + &
                               trcrn(i,j,bio_index(mm)+k-1, n, iblk)* &
                               vicen(i,j,n,iblk)*zspace(k)*trcrn(i,j,nt_fbri,n,iblk)
                  trcrn(i,j,bio_index(mm)+nblyr+1:bio_index(mm)+nblyr+2,n,iblk) = c0 !snow concentration
                  trcrn(i,j,nt_zbgc_frac-1+mm,n,iblk) = zbgc_frac_init(mm)
               enddo   !i
               enddo   !j
             enddo     !k
         enddo         !mm
         enddo        !n 
         enddo         !iblk
         !$OMP END PARALLEL DO
              
       endif  ! scale_bgc
       endif  ! skl_bgc
     endif    ! restart

     end subroutine init_bgc

!=======================================================================

      subroutine biogeochemistry (dt, iblk)

      use ice_algae, only: skl_biogeochemistry, z_biogeochemistry, &
                           update_snow_bgc  
      use ice_blocks, only: block, get_block
      use ice_brine, only: preflushing_changes, compute_microS_mushy, &
                           update_hbrine, compute_microS 
      use ice_constants, only: c0, c1, puny
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: nblyr, nilyr, nslyr, n_algae, n_zaero
      use ice_flux, only: hin_old, meltbn, melttn, congeln, snoicen, &
                          sss, sst, meltsn, hmix, dsnow, dsnown,     &
                          salinz, fsnow
      use ice_grid, only: TLAT, TLON
      use ice_shortwave, only:  fswthrun, fswpenln
      use ice_state, only: aicen_init, vicen_init, aicen, vicen, vsnon, &
          trcrn, nt_fbri, tr_brine, ntrcr, nbtrcr,  nbtrcr_sw, &
          aicen_init, vicen_init, vsnon_init, &
          nt_bgc_S, nt_qice, nt_sice, nt_zbgc_frac  
      use ice_therm_shared, only: solve_zsal           
      use ice_timers, only: timer_bgc, ice_timer_start, ice_timer_stop
      use ice_zsalinity, only: solve_zsalinity, column_sum_S, &
          merge_S_fluxes 
 
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n, mm              ! thickness category index

      integer (kind=int_kind) :: &
         icells,  &         ! number of cells with aicen > puny
         fcells, &          ! number of new ice cells with aicen > puny
         nfcells            ! number of new ice cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj ,    &  ! indirect indices for cells with aicen > puny
         findxi, findxj,   &  !  for fcells
         nfindxi, nfindxj, &  !  for nfcells
         findxij, nfindxij    ! 

      real (kind=dbl_kind), dimension (nx_block*ny_block) :: &
         hin         , & ! new ice thickness
         hsn         , & ! snow thickness  (m)
         hbr_old     , & ! old brine thickness before growh/melt
         dhice       , & ! change due to sublimation/condensation (m)
         kavg        , & ! average ice permeability (m^2)
         bphi_o      , & ! surface ice porosity 
         zsal_totn  , & ! ice Salinity  (g/m^2) 
         hbrin           ! brine height

      real (kind=dbl_kind), dimension (nx_block*ny_block,nbtrcr) :: &
         zbgc_snow, &    ! snow contribution to ice surface condition 
         zbgc_atm        ! contribution from atmosphere to ice
                         !(mmol/m^3*m)

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2) :: &
      ! Defined on Bio Grid points
         bSin        , & ! salinity on the bio grid  (ppt)
         brine_sal   , & ! brine salinity (ppt)
         brine_rho       ! brine_density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1) :: &
      ! Defined on Bio Grid interfaces
         iphin       , & ! porosity 
         ibrine_sal  , & ! brine salinity  (ppt)
         ibrine_rho  , & ! brine_density (kg/m^3)
         iTin            ! Temperature on the interface grid (oC)

      real (kind=dbl_kind), dimension (nx_block*ny_block) :: & 
         sloss            ! brine flux contribution from surface runoff (g/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1,n_algae) :: &  
         grow_alg        !specific growth rate

      real (kind=dbl_kind), dimension (nx_block,ny_block,n_algae) :: &  
         grow_alg_skl        !specific growth rate

      ! for bgc sk
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: & 
         dh_bot_chl  , & ! Chlorophyll may or may not flush
         dh_top_chl  , & ! Chlorophyll may or may not flush
         darcy_V_chl     

      ! for bgc layer
      real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr) :: &
         flux_bion       !tracer flux to ocean

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1,n_algae) :: &
         upNOn      , & ! nitrate uptake rate (mmol/m^3/s)  
         upNHn          ! ammonium uptake rate (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,n_algae) :: &
         upNOn_skl  , & ! nitrate uptake rate (mmol/m^3/s)  
         upNHn_skl      ! ammonium uptake rate (mmol/m^3/s)

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &  
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &  
         istop, jstop    ! indices of grid cell where model aborts 

      if (tr_brine .or. skl_bgc) then

         call ice_timer_start(timer_bgc) ! biogeochemistry

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         ! Define ocean concentrations for tracers used in simulation
         do mm = 1,nbtrcr
         do j = jlo, jhi
         do i = ilo, ihi       
                ocean_bio(i,j,mm, iblk)    = ocean_bio_all(i,j,bio_index_o(mm),iblk)  
         enddo  !i
         enddo  !j
         enddo  !mm

         do n = 1, ncat
            hin_old(:,:,n,iblk) = c0
            do j = jlo, jhi
            do i = ilo, ihi
               if (aicen_init(i,j,n,iblk) > puny) then 
                  hin_old(i,j,n,iblk) = vicen_init(i,j,n,iblk) &
                                      / aicen_init(i,j,n,iblk)
               else
                  first_ice(i,j,n,iblk) = .true.
                  if (tr_brine) trcrn(i,j,nt_fbri,n,iblk) = c1
                  do mm = 1,nbtrcr
                    trcrn(i,j,nt_zbgc_frac-1+mm,n,iblk) = zbgc_frac_init(mm)
                  enddo
                  if (n == 1) then 
                     Rayleigh_criteria(i,j,iblk) = .false.
                  endif
                  if (solve_zsal) trcrn(i,j,nt_bgc_S:nt_bgc_S+nblyr-1,n,iblk) = c0
               endif
            enddo  !i
            enddo  !j

            icells = 0
            fcells = 0
            nfcells = 0
            indxi(:) = 0
            indxj(:) = 0
            findxi(:) = 0
            findxj(:) = 0
            findxij(:) = 0
            nfindxi(:) = 0
            nfindxj(:) = 0
            nfindxij(:) = 0

            do j = jlo, jhi
            do i = ilo, ihi
               if (aicen(i,j,n,iblk) > puny) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
                  if (first_ice(i,j,n,iblk)) then
                    fcells = fcells + 1
                    findxi(fcells) = i
                    findxj(fcells) = j
                    findxij(fcells) = icells
                  else
                    nfcells = nfcells + 1
                    nfindxi(nfcells) = i
                    nfindxj(nfcells) = j
                    nfindxij(nfcells) = icells
                  endif                    
               endif             ! aicen > puny
             enddo               ! i
             enddo               ! j
          
          flux_bion(:,:,:) = c0
          zbgc_snow(:,:) = c0 
          zbgc_atm(:,:) = c0 
          hsn(:) = c0
          hin(:) = c0
          bSin(:,:) = c0
          zsal_totn(:) = c0
          hbrin(:) = c0
          kavg(:) = c0
          bphi_o(:) = c0
          dh_top_chl(:,:) = c0
          dh_bot_chl(:,:) = c0
          darcy_V_chl(:,:) = c0
  
          if (icells > 0) then

      !-----------------------------------------------------------------
      ! brine dynamics
      !-----------------------------------------------------------------

            if (tr_brine) then 

               call preflushing_changes (nx_block,   ny_block,            &
                                icells,              n,                   &
                                indxi,               indxj,               &
                                aicen  (:,:,n,iblk),                      &
                                vicen  (:,:,n,iblk), vsnon  (:,:,n,iblk), &
                                meltbn (:,:,n,iblk), melttn (:,:,n,iblk), &
                                congeln(:,:,n,iblk), snoicen(:,:,n,iblk), &
                                hin_old(:,:,n,iblk), dhice,               & 
                                trcrn  (:,:,nt_fbri,n,iblk),              &
                                dhbr_top(:,:,n,iblk),dhbr_bot(:,:,n,iblk),&
                                hbr_old,             hin,                 &
                                hsn,                 first_ice(:,:,n,iblk))

               if (solve_zsal)  then  
                 call compute_microS (nx_block,      ny_block,            &
                                icells,              n,                   &
                                indxi,               indxj,               &   
                                trcrn(:,:,:,n,iblk), hin_old(:,:,n,iblk), &
                                hbr_old,            sss(:,:,iblk),        &
                                sst(:,:,iblk),       bTiz(:,:,:,n,iblk),  &
                                iTin,                bphi(:,:,:,n,iblk),  &
                                kavg,    bphi_o,                          &
                                Rayleigh_criteria(:,:,iblk),              &
                                first_ice(:,:,n,iblk),                    &
                                bSin,                brine_sal,           &
                                brine_rho,           iphin,               &
                                ibrine_rho,          ibrine_sal,          &
                                sice_rho(:,:,n,iblk),sloss,               &
                                salinz(:,:,1:nilyr,iblk))
               else     

               ! Requires the average ice permeability = kavg(:)
               ! and the surface ice porosity = zphi_o(:)
               ! computed in "compute_microS" or from "thermosaline_vertical"

                   call compute_microS_mushy (nx_block, ny_block,              &
                                   icells, n, indxi,    indxj,                 &   
                                   trcrn(:,:,:,n,iblk), hin_old(:,:,n,iblk),   &
                                   hbr_old,                                    &
                                   sss(:,:,iblk),       sst(:,:,iblk),         & 
                                   bTiz(:,:,:,n,iblk),  iTin,                  &
                                   bphi(:,:,:,n,iblk),  kavg,                  &
                                   bphi_o, bSin,  brine_sal,                   &
                                   brine_rho, iphin,    ibrine_rho, ibrine_sal,&
                                   sice_rho(:,:,n,iblk),iDi(:,:,:,n,iblk))
               endif !solve_zsal  

               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

               call update_hbrine (meltbn  (i,j,n,iblk), melttn(i,j,n,iblk),  &
                                   meltsn  (i,j,n,iblk), dt,                  &
                                   hin     (ij),         hsn   (ij),          &
                                   hin_old (i,j,n,iblk), hbrin (ij),          &
                                   hbr_old (ij),                              &
                                   trcrn   (i,j,nt_fbri,n,iblk),              &
                                   snoicen(i,j,n,iblk), &
                                   dhbr_top(i,j,n,iblk), dhbr_bot(i,j,n,iblk),&
                                   dh_top_chl(i,j),      dh_bot_chl(i,j),     & 
                                   kavg    (ij),         bphi_o(ij),          &
                                   darcy_V (i,j,n,iblk), darcy_V_chl(i,j),    &  
                                   bphi(i,j,2,n,iblk))
               
               hbri(i,j,iblk) = hbri(i,j,iblk) + hbrin(ij)*aicen(i,j,n,iblk)  
               enddo                     ! ij

               if (solve_zsal) then 
                  call solve_zsalinity &
                                   (nx_block,              ny_block,             &
                                   icells,      n,         dt,                   &
                                   indxi,                  indxj,                &  
                                   trcrn(:,:,nt_bgc_S:nt_bgc_S+nblyr-1,n,iblk),  &
                                   trcrn(:,:,nt_qice:nt_qice+nilyr-1,n,iblk),    &
                                   trcrn(:,:,nt_sice:nt_sice+nilyr-1,n,iblk),    &
                                   aicen(:,:,n,iblk),      vicen(:,:,n,iblk),    &  
                                   bSin,                   bTiz(:,:,:,n,iblk),   &
                                   bphi(:,:,:,n,iblk),     iphin,                &
                                   iki(:,:,:,n,iblk),      hbr_old,              &
                                   hbrin,                  hin,                  &
                                   hin_old(:,:,n,iblk),    iDi(:,:,:,n,iblk),    &
                                   darcy_V(:,:,n,iblk),    brine_sal,            & 
                                   brine_rho,              ibrine_sal,           & 
                                   ibrine_rho,                                   &
                                   Rayleigh_criteria(:,:,iblk),                  &
                                   first_ice(:,:,n,iblk),  sss(:,:,iblk),        & 
                                   sst(:,:,iblk),          dhbr_top(:,:,n,iblk), &
                                   dhbr_bot(:,:,n,iblk),   TLAT(:,:,iblk),       &
                                   TLON(:,:,iblk),         l_stop,               &
                                   istop,                  jstop,                &
                                   fzsaln(:,:,n,iblk),     fzsaln_g(:,:,n,iblk), &
                                   bphi_o)

                  call column_sum_S (nx_block,             ny_block,             &
                                    icells,                indxi,                &
                                    indxj,                 ntrcr,                &
                                    nblyr,                 vicen(:,:,n,iblk),    &
                                    trcrn(:,:,1:ntrcr,n,iblk),                   &
                                    zsal_totn)
 
                  call merge_S_fluxes (nx_block,        ny_block,               &
                                    icells,             indxi,                  &
                                    indxj,              aicen_init(:,:,n,iblk), &
                                    zsal_totn,          zsal_tot(:,:,iblk),     &
                                    dsnown(:,:,n,iblk), dsnow  (:,:,iblk),      &
                                    fzsal(:,:,iblk),    fzsaln(:,:,n,iblk),     &
                                    fzsal_g(:,:,iblk),  fzsaln_g(:,:,n,iblk))                         
               endif  !solve_zsal

        endif !tr_brine

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

        if (z_tracers) then 
                do ij = 1, icells
                     i = indxi(ij)
                     j = indxj(ij)

                     call update_snow_bgc (dt,       melttn(i,j,n,iblk),       &
                              meltsn(i,j,n,iblk),    meltbn  (i,j,n,iblk),     &
                              congeln(i,j,n,iblk),   snoicen(i,j,n,iblk),      & 
                              nbtrcr,                fsnow(i,j,iblk),          &
                              ntrcr,                 trcrn(i,j,1:ntrcr,n,iblk),&
                              bio_index,             aicen_init(i,j,n,iblk),   &
                              zbgc_snow(ij,1:nbtrcr),vicen_init(i,j,n,iblk),   &
                              vsnon_init(i,j,n,iblk),vicen(i,j,n,iblk),        &
                              vsnon(i,j,n,iblk),     aicen(i,j,n,iblk),        &
                              flux_bio_atm(i,j,1:nbtrcr,iblk), zbgc_atm(ij,1:nbtrcr)) 
                 enddo  !ij

                 call z_biogeochemistry (nx_block, ny_block,                 &
                                   icells,   n,          dt,                 &
                                   indxi,    indxj,      nbtrcr,             &
                                   n_algae,   n_zaero,                       &  
                                   fcells,  findxi,      findxj,             &   
                                   nfcells, nfindxi,     nfindxj,            &  
                                   findxij,  nfindxij,                       &  
                                   aicen (:,:,n,iblk),   vicen(:,:,n,iblk),  & 
                                   hin_old(:,:,n,iblk),                      &
                                   ocean_bio(:,:,:,iblk), flux_bion,         &
                                   bphi(:,:,:,n,iblk),   iphin,              &         
                                   trcrn(:,:,1:ntrcr,    n,iblk),            & 
                                   iDi(:,:,:,n,iblk),    sss(:,:,iblk),      &
                                   fswpenln(:,:,:,       n,iblk),            &
                                   grow_alg,             upNOn, upNHn,       &
                                   dhbr_top(:,:,         n,iblk),            &
                                   dhbr_bot(:,:,         n,iblk),            &
                                   dh_top_chl,           dh_bot_chl,         &
                                   zfswin(:,:,:,         n,iblk),            &
                                   TLAT(:,:,iblk),       TLON(:,:,iblk),     &
                                   hbrin,                hbr_old,            &
                                   darcy_V(:,:,          n,iblk),            &
                                   darcy_V_chl,  bgrid,  igrid,              &
                                   icgrid,               bphi_o,             &
                                   zbgc_snow(:,1:nbtrcr),dhice(:),           &
                                   zbgc_atm(:,1:nbtrcr), iTin,               &
                                   trcrn_sw(:,:,:,       n,iblk) ,           &
                                   max_nsw,              swgrid,             &
                                   Zoo(:,:,:,n,iblk),                        &
                                   meltbn   (:,:,        n,iblk))   

                  call merge_bgc_fluxes(nx_block, ny_block,                  &
                                         icells,    dt,                      &
                                         indxi,    indxj,           nbtrcr,  & 
                                         aicen_init(:,:,            n,iblk), &
                                         vicen     (:,:,            n,iblk), & 
                                         vsnon     (:,:,            n,iblk), & 
                                         ntrcr,    iphin,                    &
                                         trcrn     (:,:,1:ntrcr,    n,iblk), &
                                         flux_bion, flux_bio(:,:,:,   iblk), &
                                         upNOn, upNHn, upNO(:,:,    iblk),   &
                                         upNH      (:,:,            iblk),   &
                                         zbgc_snow, zbgc_atm,                &
                                         fbio_snoice(:,:,:,    iblk),        &
                                         fbio_atmice(:,:,:,    iblk),        &
                                         PP_net     (:,:,             iblk), &
                                         ice_bio_net (:,:,1:nbtrcr,   iblk), &
                                         snow_bio_net(:,:,1:nbtrcr,   iblk), &
                                         grow_alg,    grow_net  (:,:, iblk))     

           elseif (skl_bgc) then

                  call skl_biogeochemistry (nx_block, ny_block,         &
                                         icells,   dt,                  &
                                         indxi,    indxj,               &  
                                         nbtrcr,   n_algae,             &
                                         flux_bion(:,:,1:nbtrcr),       &
                                         ocean_bio(:,:,1:nbtrcr, iblk), &
                                         hmix     (:,:,          iblk), &
                                         aicen    (:,:,        n,iblk), & 
                                         meltbn   (:,:,        n,iblk), &
                                         congeln  (:,:,        n,iblk), &
                                         fswthrun (:,:,        n,iblk), &
                                         first_ice(:,:,        n,iblk), &
                                         trcrn    (:,:,1:ntrcr,n,iblk), &
                                         upNOn_skl,          upNHn_skl, &
                                         grow_alg_skl,    nbtrcr_sw,    &
                                         trcrn_sw(:,:,1:nbtrcr_sw,n,iblk), &
                                         hin)

                 call merge_bgc_fluxes_skl(nx_block, ny_block,               &
                                         icells,    ntrcr,                   &
                                         indxi,    indxj,                    &
                                         nbtrcr,   n_algae,                  &
                                         aicen_init(:,:,            n,iblk), &
                                         trcrn     (:,:,1:ntrcr,n,iblk),     &
                                         flux_bion (:,:,1:nbtrcr),           &
                                         flux_bio  (:,:,1:nbtrcr,     iblk), &
                                         PP_net    (:,:,              iblk), &
                                         upNOn_skl,               upNHn_skl, &
                                         upNO(:,:,                  iblk),   &
                                         upNH(:,:,                  iblk),   &
                                         grow_net  (:,:,            iblk),   &
                                         grow_alg_skl)
           endif  !skl_bgc

           do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               first_ice(i,j,n,iblk) = .false.
           enddo

        endif               ! icells      
        enddo               ! ncat

         call ice_timer_stop(timer_bgc) ! biogeochemistry

      endif  ! tr_brine .or. skl_bgc
 
      end subroutine biogeochemistry

!=======================================================================
!
! Aggregate flux information from all ice thickness categories
! for z layer biogeochemistry
!
      subroutine merge_bgc_fluxes (nx_block, ny_block,   &
                               icells,       dt,         &
                               indxi,        indxj,      &
                               nbtrcr,       aicen,      &    
                               vicen,        vsnon,      &
                               ntrcr,        iphin,      &
                               trcrn,                    &
                               flux_bion,    flux_bio,   &
                               upNOn,        upNHn,      &
                               upNO,         upNH,       &
                               zbgc_snown,   zbgc_atmn,  &
                               zbgc_snow,    zbgc_atm,   &
                               PP_net,       ice_bio_net,&
                               snow_bio_net, grow_alg,   &
                               grow_net)
 
      use ice_constants, only: c1, c0, p5, secday, puny
      use ice_domain_size, only: nblyr, n_algae
      use ice_state, only: tr_brine, nt_fbri, nt_bgc_N
      use ice_therm_shared, only: solve_zsal
      use ice_shortwave, only: hs_ssl
      use ice_fileunits, only: nu_diag

      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block, & ! block dimensions
          icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
          intent(in) :: &
          indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind),  intent(in):: &          
          dt              ! timestep (s)

      integer (kind=int_kind), intent(in) :: &
          ntrcr, &       ! number of tracers
          nbtrcr         ! number of biology tracer types

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
          intent(in) :: &
          trcrn         ! input tracer fields

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &          
          aicen,    & ! concentration of ice
          vicen,    & ! volume of ice (m)
          vsnon       ! volume of snow(m)

      real (kind=dbl_kind), dimension(nx_block*ny_block,nblyr+1), intent(in):: &
          iphin        ! porosity

      ! single category rates
      real (kind=dbl_kind), dimension(nx_block*ny_block,nbtrcr), intent(in):: &
          zbgc_snown, & ! bio flux from snow to ice per cat (mmol/m^3*m) 
          zbgc_atmn     ! bio flux from atm to ice per cat (mmol/m^3*m)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nbtrcr), intent(in):: &        
          flux_bion

      ! single category rates
      real (kind=dbl_kind), dimension(nx_block,ny_block,nblyr+1,n_algae), &
          intent(in):: &          
          upNOn   , & ! nitrate uptake rate per cat (mmol/m^3/s)
          upNHn       ! ammonium uptake rate per cat (mmol/m^3/s)   
 
      real (kind=dbl_kind), dimension(nx_block,ny_block,nblyr+1,n_algae), &
          intent(in):: &
          grow_alg     ! algal growth rate per cat (mmolN/m^3/s)

      ! cumulative fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block,nbtrcr), intent(inout):: &     
           flux_bio,   &
           zbgc_snow,  &  ! bio flux from snow to ice per cat (mmol/m^2/s) 
           zbgc_atm,   &  ! bio flux from atm to ice per cat (mmol/m^2/s)
           ice_bio_net, & !integrated ice tracers mmol or mg/m^2)
           snow_bio_net   !integrated snow tracers mmol or mg/m^2)

      ! cumulative variables and rates
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: & 
          PP_net , & ! net PP (mg C/m^2/d)  times aice
          grow_net   ! net specific growth (m/d) times vice

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &          
          upNO   , & ! tot nitrate uptake rate (mmol/m^2/d) times aice 
          upNH       ! tot ammonium uptake rate (mmol/m^2/d) times aice

      ! local variables

      real (kind=dbl_kind) :: &
          tmp

      integer (kind=int_kind) :: &
          ij, i, j, &   ! horizontal indices
          k, mm         ! tracer indice

      real (kind=dbl_kind), dimension (nblyr+1) :: & 
         zspace

      real (kind=dbl_kind) :: & 
         dvssl      , & ! volume of snow surface layer (m)
         dvint          ! volume of snow interior      (m)

      !-----------------------------------------------------------------
      ! Column summation
      !-----------------------------------------------------------------
      zspace(:) = c1/real(nblyr,kind=dbl_kind)
      zspace(1) = p5/real(nblyr,kind=dbl_kind)
      zspace(nblyr+1) =  p5/real(nblyr,kind=dbl_kind)

      do mm = 1, nbtrcr
      do k = 1, nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            ice_bio_net(i,j,mm) = ice_bio_net(i,j,mm) + trcrn(i,j,bio_index(mm)+k-1)*vicen(i,j)* &
                                zspace(k)*trcrn(i,j,nt_fbri)
         enddo !ij
      enddo    !k
      
      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            dvssl  = min(p5*vsnon(i,j), hs_ssl*aicen(i,j))   !snow surface layer
            dvint  = vsnon(i,j)- dvssl                       !snow interior
            snow_bio_net(i,j,mm) = snow_bio_net(i,j,mm) +  trcrn(i,j,bio_index(mm)+nblyr+1)* &
                                 dvssl + trcrn(i,j,bio_index(mm)+nblyr+2)*dvint
            flux_bio (i,j,mm)  = flux_bio(i,j,mm) + flux_bion(i,j,mm)*aicen(i,j)
            zbgc_snow (i,j,mm) = zbgc_snow(i,j,mm) + zbgc_snown(ij,mm)/dt*aicen(i,j)
            zbgc_atm (i,j,mm)  = zbgc_atm(i,j,mm)+zbgc_atmn(ij,mm)/dt*aicen(i,j)
          enddo !ij
      enddo     !mm

      if (solve_zbgc) then
      do mm = 1, n_algae
      do k = 1, nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            PP_net(i,j)   = PP_net(i,j)     + grow_alg(i,j,k,mm)*iphin(ij,k)* &
                            (c1-fr_resp)* R_C2N(mm)*R_gC2molC * vicen(i,j)* &
                            trcrn(i,j,nt_fbri)*zspace(k)*secday
            grow_net(i,j) = grow_net(i,j)  + grow_alg(i,j,k,mm)*iphin(ij,k)* &
                            vicen(i,j)*trcrn(i,j,nt_fbri)*zspace(k)*secday/&
                            (trcrn(i,j,nt_bgc_N(mm)+k-1)+puny)
            upNO (i,j)    = upNO (i,j) + upNOn(i,j,k,mm)*vicen(i,j)*zspace(k)*trcrn(i,j,nt_fbri)&
                                         *iphin(ij,k)*secday
            upNH (i,j)    = upNH (i,j) + upNHn(i,j,k,mm)*vicen(i,j)*zspace(k)*trcrn(i,j,nt_fbri)&
                                         *iphin(ij,k)*secday
         enddo                  ! ij
      enddo                     ! k
      enddo                     ! mm
      endif

      end subroutine merge_bgc_fluxes

!=======================================================================

! Aggregate flux information from all ice thickness categories
! for skeletal layer biogeochemistry
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL

      subroutine merge_bgc_fluxes_skl (nx_block, ny_block, &
                               icells,    ntrcr,           &
                               indxi,     indxj,           &
                               nbtrcr,    n_algae,         &
                               aicen,     trcrn,           &
                               flux_bion, flux_bio,        &
                               PP_net,    upNOn,           &
                               upNHn,     upNO,            &
                               upNH,      grow_net,        &
                               grow_alg)

      use ice_constants, only: c1, secday, puny
      use ice_state, only: nt_bgc_N

      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block, & ! block dimensions
          icells,   ntrcr,    & ! number of cells with aicen > puny
          nbtrcr , n_algae      ! number of bgc tracers
                                ! number of autotrophs

      integer (kind=int_kind), dimension(nx_block*ny_block), &
          intent(in) :: &
          indxi, indxj          ! compressed indices for cells with aicen > puny

      ! single category fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &          
          aicen                 ! category ice area fraction

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
          intent(in) :: &
          trcrn               ! Bulk tracer concentration (mmol N or mg/m^3)
     
      real (kind=dbl_kind), dimension(nx_block,ny_block,nbtrcr), &
          intent(in):: &
          flux_bion             ! all bio fluxes to ocean, on categories

      real (kind=dbl_kind), dimension(nx_block,ny_block,nbtrcr), &
          intent(inout):: &
          flux_bio              ! all bio fluxes to ocean, aggregated

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_algae), &
          intent(in):: & 
          grow_alg            ! algal growth rate (mmol/m^3/s) 

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_algae), &
          intent(in):: & 
          upNOn , &           ! nitrate uptake rate per cat (mmol/m^3/s)
          upNHn               ! ammonium uptake rate per cat (mmol/m^3/s)   

      ! history output
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: & 
          PP_net , &          ! Bulk net PP (mg C/m^2/s)
          grow_net            ! net specific growth (/s)
     
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout):: &          
          upNO   , & ! tot nitrate uptake rate (mmol/m^2/s) 
          upNH       ! tot ammonium uptake rate (mmol/m^2/s)
      
      ! local variables

      integer (kind=int_kind) :: &
          ij, i, j, &   ! horizontal indices
          k, mm         ! tracer indice
    
      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         do k = 1,nbtrcr
            flux_bio (i,j,k) = flux_bio(i,j,k) + flux_bion(i,j,k)*aicen(i,j)
         enddo
     enddo
     do mm = 1, n_algae
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         PP_net   (i,j) = PP_net  (i,j) &
                        + phi_sk *sk_l*grow_alg(i,j,mm)*(c1-fr_resp) &
                        * R_C2N(mm)*R_gC2molC * aicen(i,j) *secday
         grow_net (i,j) = grow_net(i,j) + phi_sk*grow_alg(i,j,mm) * aicen(i,j) * sk_l * secday/&
                        (trcrn(i,j,nt_bgc_N(mm))+puny)
         upNO (i,j)    = upNO (i,j) + upNOn(i,j,mm)*aicen(i,j) * phi_sk * sk_l * secday
         upNH (i,j)    = upNH (i,j) + upNHn(i,j,mm)*aicen(i,j) * phi_sk * sk_l * secday
      enddo                     ! ij
      enddo
      
      end subroutine merge_bgc_fluxes_skl

!=======================================================================

! Initialize bgc fields written to history files
!
! authors: Nicole Jeffery, LANL
!          Elizabeth C. Hunke, LANL

      subroutine init_history_bgc

      use ice_constants, only: c0

      upNO     (:,:,:)   = c0   
      upNH     (:,:,:)   = c0  
      zsal_tot (:,:,:)   = c0 
      PP_net   (:,:,:)   = c0
      grow_net (:,:,:)   = c0
      ice_bio_net (:,:,:,:) = c0
      snow_bio_net (:,:,:,:) = c0
      hbri     (:,:,:)   = c0
      sice_rho (:,:,:,:) = c0  
      flux_bio (:,:,:,:) = c0 
      fbio_snoice(:,:,:,:) = c0
      fbio_atmice(:,:,:,:) = c0  
      flux_bio_ai  (:,:,:,:) = c0
      trcrn_sw (:,:,:,:,:) = c0
      Zoo(:,:,:,:,:) = c0

      end subroutine init_history_bgc

!=======================================================================

! Adjust biogeochemical tracers when new frazil ice forms

      subroutine add_new_ice_bgc (nx_block,   ny_block,   dt,       &
                                  icells,     jcells,     kcells,   &
                                  indxi,      indxj,                &
                                  indxi2,     indxj2,     indxij2,  &
                                  indxi3,     indxj3,     indxij3,  &
                                  aicen_init, vicen_init, vi0_init, &
                                  aicen,      vicen,      vsnon1,   &
                                  vi0new,   &
                                  ntrcr,      trcrn,      nbtrcr,   &
                                  sss,        ocean_bio,  flux_bio, &
                                  hsurp,      &
                                  l_stop,     istop,      jstop)

      use ice_constants, only: c0, c1, puny, depressT
      use ice_domain_size, only: nilyr, nblyr, nltrcr
      use ice_itd, only: column_sum, &
                         column_conservation_check
      use ice_state, only: tr_brine, nt_fbri, nt_sice, nt_qice, nt_Tsfc
      use ice_timers, only: timer_bgc, ice_timer_start, ice_timer_stop
      use ice_therm_shared, only: calculate_Tin_from_qin, solve_zsal 

      integer (kind=int_kind), intent(in) :: &
         nx_block, & ! block dimensions
         ny_block, & ! block dimensions
         ntrcr   , & ! number of tracers in use
         icells  , & ! number of ice/ocean grid cells
         jcells  , & ! grid cell counter
         kcells      ! grid cell counter

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj, &           ! compressed i/j indices
         indxi2, indxj2, indxij2, & ! compressed i/j indices
         indxi3, indxj3, indxij3    ! compressed i/j indices

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen_init  , & ! initial concentration of ice
         vicen_init  , & ! intiial volume per unit area of ice  (m)
         aicen       , & ! concentration of ice
         vicen           ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         vsnon1          ! category 1 snow volume per unit area (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), &
         intent(inout) :: &
         trcrn           ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         sss              !sea surface salinity (ppt)

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         vi0_init    , & ! volume of new ice added to cat 1 (intial)
         vi0new          ! volume of new ice added to cat 1

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hsurp           ! thickness of new ice added to each cat

      integer (kind=int_kind), intent(in) :: &
         nbtrcr     ! number of biology tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr), &
         intent(inout) :: &
         flux_bio   ! tracer flux to ocean from biology (mmol/m^2/s) 
        
      real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr), &
         intent(in) :: &
         ocean_bio       ! ocean concentration of biological tracer

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop    ! indices of grid cell where model aborts

! local

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         n           , & ! ice category index
         k           , & ! ice layer index
         ij, m       , & ! combined i/j horizontal indices
         indxns           ! counts jcells with no snow 

      integer (kind=int_kind), dimension (nx_block*ny_block):: &
         indxnsi, indxnsj ! compressed i/j indices

      real (kind=dbl_kind), dimension (icells) :: &
         vbri1       , & ! starting volume of existing brine
         vbri_init   , & ! brine volume summed over categories
         vbri_final      ! brine volume summed over categories

      real (kind=dbl_kind), dimension(icells) :: &
         vsurp       , & ! volume of new ice added to each cat
         vtmp            ! total volume of new and old ice
        
      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) :: &
         vbrin           ! trcrn(i,j,nt_fbri,n)*vicen(i,j,n) 
       
      real (kind=dbl_kind) :: &
         Tmlts       ! melting temperature (oC)

      character (len=char_len) :: &
         fieldid         ! field identifier

      call ice_timer_start(timer_bgc) ! biogeochemistry

      !-----------------------------------------------------------------     
      ! brine
      !-----------------------------------------------------------------
      vbrin(:,:,:) = c0
      do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1,icells
         i = indxi(ij)
         j = indxj(ij)
         vbrin(i,j,n) = vicen_init(i,j,n)
         if (tr_brine) vbrin(i,j,n) =  trcrn(i,j,nt_fbri,n)*vicen_init(i,j,n)
      enddo
      enddo
     
      call column_sum (nx_block, ny_block,       &
                       icells,   indxi,   indxj, &
                       ncat,                     &
                       vbrin,    vbri_init)

      if (nbtrcr > 0) then
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         vbri_init(ij) = vbri_init(ij) + vi0_init(ij) !bgc
         do k = 1, nbtrcr  
            flux_bio(i,j,k) = flux_bio(i,j,k) &
                            - vi0_init(ij)/dt*ocean_bio(i,j,k)*zbgc_init_frac(k)
         enddo
      enddo
      else
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         vbri_init(ij) = vbri_init(ij) + vi0_init(ij) !bgc
      enddo
      endif
      !-----------------------------------------------------------------
      ! kcells:  
      ! Distribute bgc in new ice volume among all ice categories by 
      ! increasing ice thickness, leaving ice area unchanged.
      !-----------------------------------------------------------------

      do n = 1,ncat
 
         ! Diffuse_bio handles concentration changes from ice growth/melt
         ! ice area does not change for kcells
         ! add salt to the bottom 

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
        do ij = 1, kcells
           i = indxi3(ij)
           j = indxj3(ij)
           m = indxij3(ij)

           vtmp(m) = vbrin(i,j,n)
           vsurp(m) = hsurp(m) * aicen_init(i,j,n) 
           vbrin(i,j,n) = vbrin(i,j,n) + vsurp(m)
           if (tr_brine) then
              trcrn(i,j,nt_fbri,n) = c1
              if (vicen(i,j,n) > c0)  trcrn(i,j,nt_fbri,n) = vbrin(i,j,n)/vicen(i,j,n)
           endif
        enddo

        if ((nltrcr > 0) .AND. kcells > 0 ) then  
            call adjust_tracer_profile(nx_block, ny_block,            &
                                       indxi3,   indxj3,   indxij3,   &                        
                                       icells,   kcells,              &
                                       nbtrcr,  ntrcr,                &
                                       aicen_init(:,:,n),             &
                                       vbrin(:,:,n),                  &
                                       vicen(:,:,n),                  &
                                       trcrn(:,:,:,n),                &
                                       vtmp,     vsurp,    sss,       &
                                       nilyr,    nblyr,    solve_zsal,& 
                                       bgrid,    cgrid,    ocean_bio)
        endif         
      enddo              ! n

      !-----------------------------------------------------------------
      ! jcells:  
      ! Combine bgc in new ice grown in open water with category 1 ice.
      !-----------------------------------------------------------------
       
      indxnsi(:) = 0
      indxnsj(:) = 0
      indxns = 0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, jcells
         i = indxi2(ij)
         j = indxj2(ij)
         m = indxij2(ij)
         vbri1(m)     = vbrin(i,j,1) 
         vbrin(i,j,1) = vbrin(i,j,1) + vi0new(m)
         if (tr_brine) then
            trcrn(i,j,nt_fbri,1) = c1
            if (vicen(i,j,1) > c0) trcrn(i,j,nt_fbri,1) = vbrin(i,j,1)/vicen(i,j,1)
         endif
         if (vsnon1(i,j) .le. c0) then
            indxns = indxns + 1
            indxnsi(indxns) = i
            indxnsj(indxns) = j
         endif
      enddo
       
      ! Diffuse_bio handles concentration changes from ice growth/melt
      ! ice area changes for jcells
      ! add salt throughout

      if (nltrcr > 0 .AND. jcells > 0 ) then 
            call adjust_tracer_profile(nx_block, ny_block, &
                                       indxi2,   indxj2,   indxij2,  &                        
                                       icells,   jcells,             &
                                       nbtrcr,  ntrcr,              &
                                       aicen(:,:,1),                 &
                                       vbrin(:,:,1),                 &
                                       vicen(:,:,1),                 &
                                       trcrn(:,:,:,1),               &
                                       vbri1,    vi0new,   sss,      &
                                       nilyr,    nblyr,    solve_zsal,& 
                                       bgrid,    cgrid,    ocean_bio)
  
      if (solve_zsal) then
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, indxns
         i = indxnsi(ij)
         j = indxnsj(ij)
         Tmlts = -trcrn(i,j,nt_sice,1)*depressT
         trcrn(i,j,nt_Tsfc,1) =  calculate_Tin_from_qin(trcrn(i,j,nt_qice,1),Tmlts)
      enddo
      endif           ! solve_zsal 
      endif           ! nltrcr > 0 and jcells > 0

      if (tr_brine) then
         call column_sum (nx_block, ny_block,       &
                          icells,   indxi,   indxj, &
                          ncat,                     &
                          vbrin,    vbri_final)

         fieldid = 'vbrin, add_new_ice_bgc'
         call column_conservation_check (nx_block,  ny_block,      &
                                         icells,   indxi,   indxj, &
                                         fieldid,                  &
                                         vbri_init, vbri_final,    &
                                         puny,      l_stop,        &
                                         istop,     jstop)
         if (l_stop) return
       endif

      call ice_timer_stop(timer_bgc) ! biogeochemistry

      end subroutine add_new_ice_bgc

!=======================================================================
!
! Add new ice tracers to the ice bottom and adjust the vertical profile 
!
! author: Nicole Jeffery, LANL

      subroutine adjust_tracer_profile (nx_block,     ny_block,        &
                                   indxi,   indxj,    indxij,          &
                                   icells, kcells,    nbtrcr,          &
                                   ntrcr,  aicen,     vbrin,           &
                                   vicen,  trcrn,     vtmp,            &
                                   vsurp,  sss,       nilyr,           & 
                                   nblyr,  solve_zsal, bgrid,          &
                                   cgrid,  ocean_bio)

      use ice_constants, only: c1, c0, salt_loss
      use ice_state, only: nt_sice, nt_bgc_S 

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of ocean/ice cells
         kcells            , & ! cells with frazil accumulation
         ntrcr             , & ! number of tracers in use
         nilyr             , & ! number of ice layers
         nbtrcr            , & ! number of biology tracers
         nblyr                 ! number of biology layers

      integer (kind=int_kind), dimension (icells), intent(in) :: &
         indxi, indxj    , & ! indices for i/j directions
         indxij

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume of ice
         sss         ! ocean salinity (ppt)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr), &
         intent(in) :: &
         ocean_bio

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         vbrin       ! fbri*volume per unit area of ice  (m)
       
      logical (kind=log_kind), intent(in) :: &
         solve_zsal 

      real (kind=dbl_kind), dimension(icells), intent(in) :: &
         vsurp   , & ! volume of new ice added to each cat
         vtmp        ! total volume of new and old ice
        
      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid       ! thickness of new ice added to each cat

      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid       ! thickness of new ice added to each cat

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn       ! ice tracers
      
      ! local variables

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: &
         trtmp0, &      ! temporary, remapped tracers     !need extra 
         trtmp          ! temporary, remapped tracers     !need extra 
        
      real (kind=dbl_kind), dimension (kcells) :: &
         hin     , & ! ice height
         hinS_new, & ! brine height
         temp_S         

      integer (kind=int_kind) :: &
         i, j, m, &  ! grid cell counters
         k, ij, mm  

      trtmp0(:,:,:) = c0
      trtmp(:,:,:) = c0

      do k = 1, nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, kcells
            i = indxi(ij)
            j = indxj(ij)
            m = indxij(ij)
            if (vbrin(i,j) > c0 ) then

               if (solve_zsal .and. k < nblyr + 1) then
                   trcrn(i,j,nt_bgc_S+k-1)           =  &
                  (trcrn(i,j,nt_bgc_S+k-1)           * vtmp(m) &
                                + sss(i,j)*salt_loss * vsurp(m)) / vbrin(i,j)
                  trtmp0(i,j,nt_sice+k-1) = trcrn(i,j,nt_bgc_S+k-1)
               endif
               do mm = 1, nbtrcr
                   trcrn(i,j,bio_index(mm) + k-1) = &
                   (trcrn(i,j,bio_index(mm) + k-1)         * vtmp(m) &
                      +  ocean_bio(i,j,mm)*zbgc_init_frac(mm) * vsurp(m)) / vbrin(i,j) 
               enddo
            endif  
         enddo        ! ij
      enddo        ! k 
  
      if (solve_zsal) then
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, kcells
            i = indxi(ij)
            j = indxj(ij)
            m = indxij(ij)
            if (aicen(i,j) > c0) then
               hinS_new(ij)  = vbrin(i,j)/aicen(i,j)
               hin     (ij)  = vicen(i,j)/aicen(i,j)
            else
               hinS_new(ij)  = c0
               hin     (ij)  = c0
            endif
            temp_S   (ij) = min_salin   ! bio to cice
            call remap_zbgc(ntrcr,         nilyr,          &
                            nt_sice,                         &
                            trtmp0(i,j,1:ntrcr), trtmp(i,j,:),&
                            1,                nblyr,         &
                            hin(ij),          hinS_new(ij),  &
                            cgrid(2:nilyr+1),                &
                            bgrid(2:nblyr+1), temp_S(ij)     ) 
            do k = 1, nilyr
                trcrn(i,j,nt_sice+k-1) = trtmp(i,j,nt_sice+k-1)   
            enddo        !k
         enddo           !kcells
      endif              !solve_zsal

      end subroutine adjust_tracer_profile

!=======================================================================

      end module ice_zbgc

!=======================================================================
