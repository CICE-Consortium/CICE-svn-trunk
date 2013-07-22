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
      use ice_constants
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size
      use ice_communicate, only: my_task, master_task
      use ice_fileunits, only: nu_diag, nu_forcing
      use ice_algae
      use ice_zbgc_public
      use ice_brine
      use ice_state
      use ice_zsalinity
      use ice_grid, only: tlat, tlon
      use ice_shortwave, only: fswthruln, fswthrun
      use ice_therm_shared, only: solve_Sin
      use ice_domain, only: nblocks, blocks_ice
      use ice_calendar, only: istep1

      implicit none 

      private
      public :: add_new_ice_bgc, init_zbgc, init_bgc, &
           init_history_bgc, biogeochemistry

      ! zbgc
      real (kind=dbl_kind), parameter :: & 
         l_sk_scale = 1.e-2_dbl_kind, & ! 
         grid_o_scale = 1.e-3_dbl_kind !

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_zbgc - namelist variables for vertical biogeochemistry
!
! !INTERFACE:
!
      subroutine init_zbgc
!
! !DESCRIPTION:
!
! Namelist variables, set to default values; may be altered
! at run time
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!        Nicole Jeffery, LANL
!
! !USES:
!
      use ice_broadcast, only: broadcast_scalar
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_nml, nml_filename, get_fileunit, release_fileunit
      use ice_therm_shared, only: read_Sin, solve_Sin, ktherm
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
        nml_error, & ! namelist i/o error flag
        k        , & ! loop index
        ntd          ! for tracer dependency calculation

      !-----------------------------------------------------------------
      ! Namelist variables.
      !-----------------------------------------------------------------

      namelist /zbgc_nml/  &
        hbrine, bgc_data_dir, sil_data_type, nit_data_type, &
        restore_bgc,     read_Sin,      solve_Sin,  &
        solve_skl_bgc, solve_zbgc, &
        tr_bgc_N_sk, tr_bgc_C_sk, tr_bgc_chl_sk, &
        tr_bgc_Nit_sk, tr_bgc_Am_sk, tr_bgc_Sil_sk, &
        tr_bgc_DMSPp_sk, tr_bgc_DMSPd_sk, tr_bgc_DMS_sk, &
        restart_bgc, restart_S,  restart_hbrine, scale_bgc, &
        tr_bgc_NO, tr_bgc_N, tr_bgc_NH, tr_bgc_C, tr_bgc_chl, &
        tr_bgc_DMSPp, tr_bgc_DMSPd, &
        tr_bgc_DMS, tr_bgc_Sil, tr_bgc_PON, tr_bgc_S, &
        grid_o, grid_o_t, l_sk, grid_oS, &
        l_skS, phi_snow, initbio_frac

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------

      hbrine          = .false.  ! brine height differs from ice height
      restore_bgc     = .false.  ! restore bgc if true
      read_Sin        = .false.  ! update salinity tracer profile from file
      solve_Sin       = .false.  ! update salinity tracer profile from solve_S_dt
      solve_skl_bgc   = .false.  ! solve skeletal biochemistry in diffuse bio
      solve_zbgc      = .false.  ! solve layer biochemistry in diffuse bio
      bgc_data_dir    = 'unknown_bgc_data_dir'
      sil_data_type   = 'default'
      nit_data_type   = 'default'
      tr_bgc_N_sk     = .false. ! biogeochemistry, algae (skeletal)
      tr_bgc_C_sk     = .false. ! biogeochemistry, 
      tr_bgc_chl_sk   = .false. ! biogeochemistry,
      tr_bgc_Nit_sk   = .false. ! biogeochemistry, nutrients (skeletal)
      tr_bgc_Am_sk    = .false. ! biogeochemistry, 
      tr_bgc_Sil_sk   = .false. ! biogeochemistry,
      tr_bgc_DMSPp_sk  = .false. ! biogeochemistry, trace gases (skeletal)
      tr_bgc_DMSPd_sk  = .false. ! biogeochemistry, trace gases (skeletal)
      tr_bgc_DMS_sk    = .false. ! biogeochemistry, trace gases (skeletal) 
      restart_bgc      = .false. ! biogeochemistry restart
      restart_S        = .false. ! salinity restart
      restart_hbrine   = .false. ! hbrine restart
      scale_bgc        = .false. ! initial bgc tracers proportional to S
      tr_bgc_NO        = .false.  ! layer biogeochemistry
      tr_bgc_C         = .false.  ! layer biogeochemistry
      tr_bgc_chl       = .false.  ! layer biogeochemistry
      tr_bgc_Sil       = .false.  ! layer biogeochemistry
      tr_bgc_N         = .false.  ! layer biogeochemistry
      tr_bgc_NH        = .false.  ! layer biogeochemistry
      tr_bgc_DMSPp     = .false.  ! layer biogeochemistry
      tr_bgc_DMSPd     = .false.  ! layer biogeochemistry
      tr_bgc_DMS       = .false.  ! layer biogeochemistry
      tr_bgc_PON       = .false.  ! layer biogechemistry
      tr_bgc_S         = .false.  ! layer biogechemistry
      ! biology parameter
      grid_o     = c5            ! for bottom flux        
      grid_o_t   = c5            ! for top flux        
      l_sk       = 7.0_dbl_kind  ! characteristic diffusive scale (m)                                   
      initbio_frac = c1          ! fraction of ocean tracer concentration in bio tracers
      !  parameters for Salinity
      grid_oS     = c5             ! for bottom flux         
      l_skS       = 7.0_dbl_kind   ! characteristic diffusive scale (m)      
      phi_snow    = 0.5_dbl_kind   ! snow porosity

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)

      if (my_task == master_task) then
         open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif 

         do while (nml_error > 0)
            print*,'Reading zbgc_nml'
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
      ! zsalinity
      !-----------------------------------------------------------------

      call broadcast_scalar(hbrine,             master_task)
      call broadcast_scalar(read_Sin,           master_task)
      call broadcast_scalar(solve_Sin,          master_task)
      call broadcast_scalar(tr_bgc_S,           master_task)
      if (solve_Sin) then
        hbrine   = .true. ! echmod, for now
        tr_bgc_S = .true. ! echmod, for now
        ktherm = 1
      endif

      nt_fbri = c0
      if (hbrine) then
          nt_fbri = ntrcr + 1   ! ice volume fraction with salt
          ntrcr = ntrcr + 1
      endif
      if (hbrine) trcr_depend(nt_fbri) = 1

      if (tr_bgc_S)then
          nt_bgc_S = ntrcr + 1
          ntrcr = ntrcr + nblyr
      endif 

      ntd = 0                    ! if nt_fbri /= 0 then use fbri dependency
      if (nt_fbri == 0) ntd = -1 ! otherwise make tracers depend on ice volume
      do k = 1,nblyr
         if (tr_bgc_S)     trcr_depend(nt_bgc_S     + k - 1)  = 2+nt_fbri+ntd
      enddo

      call broadcast_scalar(restart_S,          master_task)
      call broadcast_scalar(restart_hbrine,     master_task)
      call broadcast_scalar(solve_skl_bgc,      master_task)
      call broadcast_scalar(solve_zbgc,         master_task)
      call broadcast_scalar(restart_bgc,        master_task)
      call broadcast_scalar(grid_oS,            master_task)
      call broadcast_scalar(l_skS,              master_task)
      call broadcast_scalar(phi_snow,           master_task)

      if (phi_snow .le. c0) phi_snow = c1-rhos/rhoi
      if (solve_skl_bgc .AND. solve_zbgc) call abort_ice('error:solve_skl_bgc &
              and solve_zbgc are both true')
      if (solve_skl_bgc) then
            tr_bgc_N_sk      = .true.   !minimum NP biogeochemistry
            tr_bgc_Nit_sk    = .true.
            tr_bgc_NO        = .false.  ! No layer biogeochemistry
            tr_bgc_C         = .false.  
            tr_bgc_chl       = .false.  
            tr_bgc_Sil       = .false.  
            tr_bgc_N         = .false. 
            tr_bgc_NH        = .false. 
            tr_bgc_DMSPp     = .false. 
            tr_bgc_DMSPd     = .false. 
            tr_bgc_DMS       = .false. 
            tr_bgc_PON       = .false. 
      elseif (solve_zbgc) then
            tr_bgc_N         = .true.   !minimum NP biochemistry
            tr_bgc_NO        = .true.   
            tr_bgc_N_sk      = .false.  ! No skl biochemistry
            tr_bgc_C_sk      = .false.
            tr_bgc_chl_sk    = .false. 
            tr_bgc_Nit_sk    = .false.
            tr_bgc_Am_sk     = .false.
            tr_bgc_Sil_sk    = .false.
            tr_bgc_DMSPp_sk  = .false.
            tr_bgc_DMSPd_sk  = .false.
            tr_bgc_DMS_sk    = .false.
      else                              ! only tr_bgc_NO allowed
            tr_bgc_C         = .false. 
            tr_bgc_chl       = .false. 
            tr_bgc_Sil       = .false. 
            tr_bgc_N         = .false. 
            tr_bgc_NH        = .false. 
            tr_bgc_DMSPp     = .false. 
            tr_bgc_DMSPd     = .false. 
            tr_bgc_DMS       = .false. 
            tr_bgc_PON       = .false. 
            tr_bgc_N_sk      = .false.
            tr_bgc_C_sk      = .false.
            tr_bgc_chl_sk    = .false.
            tr_bgc_Nit_sk    = .false.
            tr_bgc_Am_sk     = .false.
            tr_bgc_Sil_sk    = .false.
            tr_bgc_DMSPp_sk  = .false.
            tr_bgc_DMSPd_sk  = .false.
            tr_bgc_DMS_sk    = .false.
      endif

      grid_oS = grid_oS * grid_o_scale
      l_skS = l_skS * l_sk_scale

      if (my_task == master_task) then
         write(nu_diag,1010) ' hbrine                    = ', hbrine
         write(nu_diag,1010) ' read_Sin                  = ', read_Sin
         write(nu_diag,1010) ' solve_Sin                 = ', solve_Sin
         write(nu_diag,1010) ' solve_skl_bgc             = ', solve_skl_bgc
         write(nu_diag,1010) ' solve_zbgc                = ', solve_zbgc
         write(nu_diag,1010) ' restart_S                 = ', restart_S
         write(nu_diag,1010) ' restart_hbrine            = ', restart_hbrine
         write(nu_diag,1010) ' tr_bgc_S                  = ', tr_bgc_S
         write(nu_diag,1000) ' grid_oS                   = ', grid_oS
         write(nu_diag,1060) ' l_skS                     = ', l_skS
         write(nu_diag,1060) ' phi_snow                  = ', phi_snow
      endif

      !----------------------------------------------------
      ! zbgc
      !----------------------------------------------------
      call broadcast_scalar(grid_o,             master_task)
      call broadcast_scalar(grid_o_t,           master_task)
      call broadcast_scalar(l_sk,               master_task)

      grid_o = grid_o * grid_o_scale
      grid_o_t = grid_o_t * grid_o_scale
      l_sk = l_sk * l_sk_scale

      call broadcast_scalar(scale_bgc,          master_task)
      call broadcast_scalar(restore_bgc,        master_task)
      call broadcast_scalar(bgc_data_dir,       master_task)
      call broadcast_scalar(sil_data_type,      master_task)
      call broadcast_scalar(nit_data_type,      master_task)
      call broadcast_scalar(tr_bgc_N_sk,        master_task)
      call broadcast_scalar(tr_bgc_C_sk,        master_task)
      call broadcast_scalar(tr_bgc_chl_sk,      master_task)
      call broadcast_scalar(tr_bgc_Nit_sk,      master_task)
      call broadcast_scalar(tr_bgc_Am_sk,       master_task)
      call broadcast_scalar(tr_bgc_Sil_sk,      master_task)
      call broadcast_scalar(tr_bgc_DMSPp_sk,    master_task)
      call broadcast_scalar(tr_bgc_DMSPd_sk,    master_task)
      call broadcast_scalar(tr_bgc_DMS_sk,      master_task)
      call broadcast_scalar(tr_bgc_NO,          master_task)
      call broadcast_scalar(tr_bgc_C,           master_task)
      call broadcast_scalar(tr_bgc_chl,         master_task)
      call broadcast_scalar(tr_bgc_Sil,         master_task)
      call broadcast_scalar(tr_bgc_N,           master_task)
      call broadcast_scalar(tr_bgc_NH,          master_task)
      call broadcast_scalar(tr_bgc_DMSPp,       master_task)
      call broadcast_scalar(tr_bgc_DMSPd,       master_task)
      call broadcast_scalar(tr_bgc_DMS,         master_task)
      call broadcast_scalar(tr_bgc_PON,         master_task)
      call broadcast_scalar(initbio_frac,       master_task)

      if (.not. solve_zbgc .AND. .not. tr_bgc_NO .AND. .not. solve_skl_bgc) return
      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------
      if (my_task == master_task) then

         write(nu_diag,*)    ' sil_data_type             = ', &
                               trim(sil_data_type)
         write(nu_diag,*)    ' nit_data_type             = ', &
                               trim(nit_data_type)
         write(nu_diag,*)    ' bgc_data_dir              = ', &
                               trim(bgc_data_dir)
         write(nu_diag,1010) ' tr_bgc_N_sk               = ', tr_bgc_N_sk
         write(nu_diag,1010) ' tr_bgc_C_sk               = ', tr_bgc_C_sk
         write(nu_diag,1010) ' tr_bgc_chl_sk             = ', tr_bgc_chl_sk
         write(nu_diag,1010) ' tr_bgc_Nit_sk             = ', tr_bgc_Nit_sk
         write(nu_diag,1010) ' tr_bgc_Am_sk              = ', tr_bgc_Am_sk
         write(nu_diag,1010) ' tr_bgc_Sil_sk             = ', tr_bgc_Sil_sk
         write(nu_diag,1010) ' tr_bgc_DMSPp_sk           = ', tr_bgc_DMSPp_sk
         write(nu_diag,1010) ' tr_bgc_DMSPd_sk           = ', tr_bgc_DMSPd_sk
         write(nu_diag,1010) ' tr_bgc_DMS_sk             = ', tr_bgc_DMS_sk
         write(nu_diag,1010) ' restart_bgc               = ', restart_bgc
         write(nu_diag,1010) ' scale_bgc                 = ', scale_bgc
         write(nu_diag,1010) ' tr_bgc_NO                 = ', tr_bgc_NO
         write(nu_diag,1010) ' tr_bgc_N                  = ', tr_bgc_N
         write(nu_diag,1010) ' tr_bgc_NH                 = ', tr_bgc_NH
         write(nu_diag,1010) ' tr_bgc_C                  = ', tr_bgc_C
         write(nu_diag,1010) ' tr_bgc_Sil                = ', tr_bgc_Sil
         write(nu_diag,1010) ' tr_bgc_chl                = ', tr_bgc_chl
         write(nu_diag,1010) ' tr_bgc_DMSPp              = ', tr_bgc_DMSPp
         write(nu_diag,1010) ' tr_bgc_DMSPd              = ', tr_bgc_DMSPd
         write(nu_diag,1010) ' tr_bgc_DMS                = ', tr_bgc_DMS
         write(nu_diag,1010) ' tr_bgc_PON                = ', tr_bgc_PON
         !bio parameters
         write(nu_diag,1000) ' grid_o                    = ', grid_o
         write(nu_diag,1000) ' grid_o_t                  = ', grid_o_t
         write(nu_diag,1060) ' l_sk                      = ', l_sk
         write(nu_diag,1000) ' initbio_frac              = ', initbio_frac

      endif   ! master_task

!zbgc
!  skeletal layer biology model

      if (solve_skl_bgc) then

         nt_bgc_N_sk = ntrcr + 1
         ntrcr = ntrcr + 1
      
         if (tr_bgc_C_sk) then
             nt_bgc_C_sk = ntrcr + 1
             ntrcr = ntrcr + 1
         endif    
         if (tr_bgc_chl_sk)then
             nt_bgc_chl_sk = ntrcr + 1
             ntrcr = ntrcr + 1
         endif 
         nt_bgc_Nit_sk = ntrcr + 1
         ntrcr = ntrcr + 1
         if (tr_bgc_Am_sk)then
             nt_bgc_Am_sk = ntrcr + 1
             ntrcr = ntrcr + 1
         endif    
         if (tr_bgc_Sil_sk)then
             nt_bgc_Sil_sk = ntrcr + 1
             ntrcr = ntrcr + 1
         endif    
         if (tr_bgc_DMSPp_sk)then
             nt_bgc_DMSPp_sk = ntrcr + 1
             ntrcr = ntrcr + 1
         endif    
         if (tr_bgc_DMSPd_sk)then
             nt_bgc_DMSPd_sk = ntrcr + 1
             ntrcr = ntrcr + 1
         endif    
         if (tr_bgc_DMS_sk)then
             nt_bgc_DMS_sk = ntrcr + 1
             ntrcr = ntrcr + 1
         endif  
      else
         ntrace_start = -1

         if (tr_bgc_N) then
             nt_bgc_N = ntrcr + 1            
             ntrcr = ntrcr + nblyr
             ntrace_start = nt_bgc_N
         endif        
         if (tr_bgc_C) then
             nt_bgc_C = ntrcr + 1
             ntrcr = ntrcr + nblyr
             if (ntrace_start < 0) ntrace_start = nt_bgc_C
         endif    
         if (tr_bgc_chl)then
             nt_bgc_chl = ntrcr + 1
             ntrcr = ntrcr + nblyr
             if (ntrace_start < 0) ntrace_start = nt_bgc_chl
         endif    
         if (tr_bgc_NO)then
             nt_bgc_NO = ntrcr + 1
             ntrcr = ntrcr + nblyr
             if (ntrace_start < 0) ntrace_start = nt_bgc_NO
         endif    
         if (tr_bgc_NH)then
             nt_bgc_NH = ntrcr + 1
             ntrcr = ntrcr + nblyr
             if (ntrace_start < 0) ntrace_start = nt_bgc_NH
         endif    
         if (tr_bgc_Sil)then
             nt_bgc_Sil = ntrcr + 1
             ntrcr = ntrcr + nblyr
             if (ntrace_start < 0) ntrace_start = nt_bgc_Sil
         endif    
         if (tr_bgc_DMSPp)then
             nt_bgc_DMSPp = ntrcr + 1
             ntrcr = ntrcr + nblyr
             if (ntrace_start < 0) ntrace_start = nt_bgc_DMSPp
         endif    
         if (tr_bgc_DMSPd)then
             nt_bgc_DMSPd = ntrcr + 1
             ntrcr = ntrcr + nblyr
             if (ntrace_start < 0) ntrace_start = nt_bgc_DMSPd
         endif    
         if (tr_bgc_DMS)then
             nt_bgc_DMS = ntrcr + 1
             ntrcr = ntrcr + nblyr
             if (ntrace_start < 0) ntrace_start = nt_bgc_DMS
         endif    
         if (tr_bgc_PON)then
             nt_bgc_PON = ntrcr + 1
             ntrcr = ntrcr + nblyr
             if (ntrace_start < 0) ntrace_start = nt_bgc_PON
         endif
      endif  ! solve_skl_bgc
       

      if (my_task == master_task) then
         write(nu_diag,1020)'nt_bgc_N = ', nt_bgc_N
         write(nu_diag,1020)'nt_bgc_S = ', nt_bgc_S
         write(nu_diag,1020)'nt_bgc_N_sk = ', nt_bgc_N_sk
         write(nu_diag,1020)'nt_bgc_Nit_sk = ', nt_bgc_Nit_sk
         write(nu_diag,*)' '
         write(nu_diag,1020)'nblyr', nblyr
      endif

      ! BGC layer model (on bottom "skeletal" layer)
      if (tr_bgc_N_sk)     trcr_depend(nt_bgc_N_sk)     = 0 ! algae  (skeletal)
      if (tr_bgc_C_sk)     trcr_depend(nt_bgc_C_sk)     = 0 ! 
      if (tr_bgc_chl_sk)   trcr_depend(nt_bgc_chl_sk)   = 0 ! 
      if (tr_bgc_Nit_sk)   trcr_depend(nt_bgc_Nit_sk)   = 0 ! nutrients (skeletal)
      if (tr_bgc_Am_sk)    trcr_depend(nt_bgc_Am_sk)    = 0 ! 
      if (tr_bgc_Sil_sk)   trcr_depend(nt_bgc_Sil_sk)   = 0 ! 
      if (tr_bgc_DMSPp_sk) trcr_depend(nt_bgc_DMSPp_sk) = 0 ! trace gases
      if (tr_bgc_DMSPd_sk) trcr_depend(nt_bgc_DMSPd_sk) = 0 !
      if (tr_bgc_DMS_sk)   trcr_depend(nt_bgc_DMS_sk)   = 0 !
   
      ! BGC layer model (on bio grid)  Bulk concentration
      ! volume-weighted tracers determined by brine level
      do k = 1,nblyr
         if (tr_bgc_NO)    trcr_depend(nt_bgc_NO    + k - 1)  = 2+nt_fbri+ntd
         if (tr_bgc_N)     trcr_depend(nt_bgc_N     + k - 1)  = 2+nt_fbri+ntd
         if (tr_bgc_NH)    trcr_depend(nt_bgc_NH    + k - 1)  = 2+nt_fbri+ntd
         if (tr_bgc_chl)   trcr_depend(nt_bgc_chl   + k - 1)  = 2+nt_fbri+ntd
         if (tr_bgc_Sil)   trcr_depend(nt_bgc_Sil   + k - 1)  = 2+nt_fbri+ntd
         if (tr_bgc_C)     trcr_depend(nt_bgc_C     + k - 1)  = 2+nt_fbri+ntd
         if (tr_bgc_DMSPp) trcr_depend(nt_bgc_DMSPp + k - 1)  = 2+nt_fbri+ntd
         if (tr_bgc_DMSPd) trcr_depend(nt_bgc_DMSPd + k - 1)  = 2+nt_fbri+ntd
         if (tr_bgc_DMS)   trcr_depend(nt_bgc_DMS   + k - 1)  = 2+nt_fbri+ntd
         if (tr_bgc_PON)   trcr_depend(nt_bgc_PON   + k - 1)  = 2+nt_fbri+ntd
      enddo

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1060    format (a30,2x,2D13.2)! dbl precision

      end subroutine init_zbgc

!=======================================================================
!BOP
!
! !ROUTINE: init_bgc
!
! !DESCRIPTION:
!
!  Initialize vertical profile of biogeochemistry
!  
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_bgc
!
! !USES:
!
      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: nblocks
      use ice_domain_size, only: max_blocks
      use ice_flux, only:  hmix, sss
!      use ice_zbgc_public, only: upNO, upNH, growN, growNp, bgc_tracer_type
!      use ice_zbgc_public, only: zfswin, ocean_bio
      use ice_calendar, only: month, dt
      use ice_therm_shared, only: solve_Sin
      use ice_restart, only: runtype
      use ice_algae  
      use ice_exit, only: abort_ice
      use ice_read_write, only: ice_read, ice_open
!      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
           i, j, iblk       , & ! horizontal indices
           ij               , & ! horizontal index, combines i and j loops
           l                , & ! month index
           k,m               , & ! vertical index
           n                , & ! category index
           nbits

      real (kind=dbl_kind), dimension(nilyr+2) :: &
         Temp      

      real (kind=dbl_kind) :: &
         xin                , & !transformed coordinate point       
         zspace                 !grid spacing for CICE vertical grid

      logical (kind=log_kind) :: &
         dbug             ! prints debugging output if true

       
      real (kind=dbl_kind), dimension(7), parameter :: & 
         Ps  =  (/0.8828e3_dbl_kind, -2.0570e3_dbl_kind, &  !polynomial fit
                  1.8255e3_dbl_kind, -0.7523e3_dbl_kind, &  !for initial
                  0.1456e3_dbl_kind, -0.0299e3_dbl_kind, &  !Salinity
                  0.0223_dbl_kind/)  

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      !-----------------------------------------------------------------------------   
      !     BGC Layer Model
      !-----------------------------------------------------------------------------   

      dbug = .true.

      zfswin(:,:,:,:,:) = c0   !shortwave flux on bio grid
      ntraceb = 0
      nlt_bgc_NO = 0
      nlt_bgc_N = 0
      nlt_bgc_C = 0
      nlt_bgc_chl = 0
      nlt_bgc_NH = 0
      nlt_bgc_Sil = 0
      nlt_bgc_DMSPp = 0
      nlt_bgc_DMSPd = 0
      nlt_bgc_DMS = 0
      nlt_bgc_PON = 0
  
      !if (trim(runtype) == 'continue') restart_bgc = .true.

        if (tr_bgc_NO .or. tr_bgc_Nit_sk)   then  !initialize like S
           ntraceb = ntraceb + 1
           nlt_bgc_NO = ntraceb
           bgc_tracer_type(ntraceb) = c1
        endif
        if (tr_bgc_N .or. tr_bgc_N_sk)   then  !initialize like S
          ntraceb = ntraceb + 1
          nlt_bgc_N = ntraceb
           bgc_tracer_type(ntraceb) = c0
        endif
        if (tr_bgc_C .or. tr_bgc_C_sk)   then  !initialize like S
          ntraceb = ntraceb + 1
          nlt_bgc_C = ntraceb
           bgc_tracer_type(ntraceb) = c0
        endif
        if (tr_bgc_chl .or. tr_bgc_chl_sk)   then  !initialize like S
          ntraceb = ntraceb + 1
          nlt_bgc_chl = ntraceb
           bgc_tracer_type(ntraceb) = c0
        endif
        if (tr_bgc_NH .or. tr_bgc_Am_sk)   then  !initialize like S
          ntraceb = ntraceb + 1
          nlt_bgc_NH = ntraceb
           bgc_tracer_type(ntraceb) = c1
        endif
        if (tr_bgc_Sil .or. tr_bgc_Sil_sk)   then  !initialize like S
          ntraceb = ntraceb + 1
          nlt_bgc_Sil = ntraceb
           bgc_tracer_type(ntraceb) = c1
        endif
        if (tr_bgc_DMSPp .or. tr_bgc_DMSPp_sk)   then  !initialize like S
          ntraceb = ntraceb + 1
          nlt_bgc_DMSPp = ntraceb
           bgc_tracer_type(ntraceb) = c0
        endif
        if (tr_bgc_DMSPd .or. tr_bgc_DMSPd_sk)   then  !initialize like S
          ntraceb = ntraceb + 1
          nlt_bgc_DMSPd = ntraceb
           bgc_tracer_type(ntraceb) = c1
        endif
        if (tr_bgc_DMS .or. tr_bgc_DMS_sk)   then  !initialize like S
          ntraceb = ntraceb + 1
          nlt_bgc_DMS = ntraceb
           bgc_tracer_type(ntraceb) = c1
        endif
        if (tr_bgc_PON)   then  !initialize like S
          ntraceb = ntraceb + 1
          nlt_bgc_PON = ntraceb
           bgc_tracer_type(ntraceb) = c1
        endif
             
        if (ntraceb .NE. nbltrcr) then
           write (nu_diag,*) ' '
           write (nu_diag,*) 'ntraceb /= nbltrcr'
           write (nu_diag,*) 'ntraceb, nbltrcr:',ntraceb, nbltrcr
           call abort_ice ('ice: ice_init error')
        endif

      if (restart_bgc) then       

        call read_restart_bgc

      elseif (scale_bgc .AND. tr_bgc_S) then

     
      !-----------------------------------------------------------------------------   
      !     Ocean Values
      !-----------------------------------------------------------------------------   

         sil(:,:,:) = c10 !initial ocean silicate (mmol/m^3)
         nit(:,:,:) = c5  !initial ocean nitrate (mmol/m^3)
         amm(:,:,:) = c1  !initial ocean ammonia (mmol/m^3)
         dmsp(:,:,:)= R_S2N*0.15_dbl_kind !sulfur cycle product (mmol/m^3)
         dms(:,:,:) = c0 !sulfur cycle product (mmol/m^3)
         algalN(:,:,:) = 0.15_dbl_kind ! initial mixed layer algal concentration (mmol/m^3)

      !-----------------------------------------------------------------------------   
      !     Skeletal Layer Model
      !-----------------------------------------------------------------------------   

         if (tr_bgc_N_sk)     trcrn(:,:,nt_bgc_N_sk,:,:)     = 0.15_dbl_kind/phi_sk*sk_l
         if (tr_bgc_C_sk)     trcrn(:,:,nt_bgc_C_sk,:,:)     = R_C2N*0.15_dbl_kind/phi_sk*sk_l
         if (tr_bgc_chl_sk)   trcrn(:,:,nt_bgc_chl_sk,:,:)   = R_chl2N**0.15_dbl_kind/phi_sk*sk_l
         if (tr_bgc_Nit_sk)   trcrn(:,:,nt_bgc_Nit_sk,:,:)   = c5/phi_sk*sk_l
         if (tr_bgc_Am_sk)    trcrn(:,:,nt_bgc_Am_sk,:,:)    = c1/phi_sk*sk_l
         if (tr_bgc_Sil_sk)   trcrn(:,:,nt_bgc_Sil_sk,:,:)   = c10/phi_sk*sk_l
         if (tr_bgc_DMSPp_sk) trcrn(:,:,nt_bgc_DMSPp_sk,:,:) = R_S2N*0.15_dbl_kind/phi_sk*sk_l
         if (tr_bgc_DMSPd_sk) trcrn(:,:,nt_bgc_DMSPd_sk,:,:) = c0
         if (tr_bgc_DMS_sk)   trcrn(:,:,nt_bgc_DMS_sk,:,:)   = c0
 
      !-------------------------------------------------------------------
      ! silicate
      !-------------------------------------------------------------------
    
         nbits = 64                ! double precision data

         if (trim(sil_data_type) == 'clim' .AND. (tr_bgc_Sil .OR. tr_bgc_Sil_sk)) then

            sil_file = trim(bgc_data_dir)//'silicate_WOA2005_surface_monthly' ! gx1 only

            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'silicate initialized from:'
               write (nu_diag,*) trim(sil_file)
            endif

            if (my_task == master_task) &
              call ice_open (nu_forcing, sil_file, nbits)

            k = month
            call ice_read (nu_forcing, k, work1, 'rda8', dbug, &
                           field_loc_center, field_type_scalar)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  sil(i,j,iblk) = work1(i,j,iblk)
               enddo
               enddo
            enddo

           if (my_task == master_task) close(nu_forcing)

         elseif (trim(sil_data_type) == 'rct_clim' .AND. (tr_bgc_Sil .OR. tr_bgc_Sil_sk)) then
        
          !use WOA2005_surface (winter or spring) for a specific location
          !(Bering (60, 180), Okhotsk (55, 150E),  Chukchi (70, 170W) Labrador Sea (56, 50W), central(0,86)) 
          !           Just March: (25, 50, 30, 2.5,  20)
          !mmol/m^3  Apr, May, Jun spring range: (20, 40, 10, 2.5, 20)
          !          Jan, Feb, Mar winter range:  (20, 60, 25, 2.5, 20)
          
            sil(:,:,:) = 30.0_dbl_kind  !chukchi, march
         endif

      !-------------------------------------------------------------------
      ! nitrate
      !-------------------------------------------------------------------

         if (trim(nit_data_type) == 'clim' .AND. (tr_bgc_NO .OR. tr_bgc_Nit_sk)) then

            nit_file = trim(bgc_data_dir)//'nitrate_WOA2005_surface_monthly' ! gx1 only

            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'nitrate initialized from:'
               write (nu_diag,*) trim(nit_file)
            endif

            if (my_task == master_task) &
              call ice_open (nu_forcing, nit_file, nbits)

            l = month
            call ice_read (nu_forcing, l, work1, 'rda8', dbug, &
                           field_loc_center, field_type_scalar)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  nit(i,j,iblk) = work1(i,j,iblk)                 
               enddo
               enddo
            enddo

           if (my_task == master_task) close(nu_forcing)

         elseif (trim(nit_data_type) == 'rct_clim'  .AND. (tr_bgc_NO .OR. tr_bgc_Nit_sk)) then
            
 !use WOA2005_surface (winter or spring) for a specific location
 !(Bering (60, 180), Okhotsk (55, 150E),  Chukchi (70, 170W) Labrador Sea (56, 50W), central(0,86)) 
 !           Just March: (15, 25, 10, 15,  0.25)
 !mmol/m^3  Apr, May, Jun spring range: (5<, 25, 0.5, 10, 5)
 !          Jan, Feb, Mar winter range: (15, 25, 10, 13, 0.25)
            
            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'nitrate initialized from March, Chukchi Sea'
            endif
             
            nit(:,:,:) = c10  !chukchi, march

         elseif (trim(nit_data_type) == 'sss'  .AND. (tr_bgc_NO .OR. tr_bgc_Nit_sk)) then
            
            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'nitrate initialized from salinity'
            endif
            
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  nit(i,j,iblk) =  sss(i,j,iblk)        
               enddo
               enddo
            enddo
               
         endif
      
     if (tr_bgc_NO)   then  !initialize like S
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
               do k = 1, nblyr
                  trcrn(i,j,nt_bgc_NO+k-1,n,iblk) =  trcrn(i,j,nt_bgc_S+k-1,n,iblk)/sss(i,j,iblk)*nit(i,j,iblk)
               enddo      !k
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif                 
     if (tr_bgc_N)    then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_N+k-1,n,iblk)      =  algalN(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif
     if (tr_bgc_C)   then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_C+k-1,n,iblk)      = R_C2N*algalN(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif
     if (tr_bgc_chl)  then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_chl+k-1,n,iblk)      =  R_chl2N*algalN(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif  
     if (tr_bgc_NH )  then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_NH+k-1,n,iblk)      =   trcrn(i,j,nt_bgc_S+k-1,n,iblk)/sss(i,j,iblk)*amm(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif 
     if (tr_bgc_Sil) then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_Sil+k-1,n,iblk)      = trcrn(i,j,nt_bgc_S+k-1,n,iblk)/sss(i,j,iblk)*sil(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif
     if (tr_bgc_DMSPp) then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_DMSPp+k-1,n,iblk)      = dmsp(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif
     if (tr_bgc_DMSPd) then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_DMSPd+k-1,n,iblk)      =  trcrn(i,j,nt_bgc_S+k-1,n,iblk)/sss(i,j,iblk)*dmsp(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif 
     if (tr_bgc_DMS)  then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
              trcrn(i,j,nt_bgc_DMS+k-1,n,iblk)      =  trcrn(i,j,nt_bgc_S+k-1,n,iblk)/sss(i,j,iblk)*dms(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif

       ! PON shadows nitrate without reactions

      if (tr_bgc_PON )  then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_PON+k-1,n,iblk)      = trcrn(i,j,nt_bgc_S+k-1,n,iblk)/sss(i,j,iblk)*nit(i,j,iblk)
                 !algalN(i,j,iblk)
                 
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
      endif

      else ! not restarting and .not.(scale_bgc .AND. tr_bgc_S) 
         
      !-----------------------------------------------------------------------------   
      !     Ocean Values
      !-----------------------------------------------------------------------------   

         sil(:,:,:) = c10 !initial ocean silicate (mmol/m^3)
         nit(:,:,:) = c5  !initial ocean nitrate (mmol/m^3)
         amm(:,:,:) = c1  !initial ocean ammonia (mmol/m^3)
         dmsp(:,:,:)= R_S2N*0.15_dbl_kind !sulfur cycle product (mmol/m^3)
         dms(:,:,:) = c0 !sulfur cycle product (mmol/m^3)
         algalN(:,:,:) = 0.15_dbl_kind ! initial mixed layer algal concentration (mmol/m^3)

      !-----------------------------------------------------------------------------   
      !     Skeletal Layer Model
      !-----------------------------------------------------------------------------   

         if (tr_bgc_N_sk)     trcrn(:,:,nt_bgc_N_sk,:,:)     = 0.15_dbl_kind/phi_sk*sk_l
         if (tr_bgc_C_sk)     trcrn(:,:,nt_bgc_C_sk,:,:)     = R_C2N*0.15_dbl_kind/phi_sk*sk_l
         if (tr_bgc_chl_sk)   trcrn(:,:,nt_bgc_chl_sk,:,:)   = R_chl2N**0.15_dbl_kind/phi_sk*sk_l
         if (tr_bgc_Nit_sk)   trcrn(:,:,nt_bgc_Nit_sk,:,:)   = c5/phi_sk*sk_l
         if (tr_bgc_Am_sk)    trcrn(:,:,nt_bgc_Am_sk,:,:)    = c1/phi_sk*sk_l
         if (tr_bgc_Sil_sk)   trcrn(:,:,nt_bgc_Sil_sk,:,:)   = c10/phi_sk*sk_l
         if (tr_bgc_DMSPp_sk) trcrn(:,:,nt_bgc_DMSPp_sk,:,:) = R_S2N*0.15_dbl_kind/phi_sk*sk_l
         if (tr_bgc_DMSPd_sk) trcrn(:,:,nt_bgc_DMSPd_sk,:,:) = c0
         if (tr_bgc_DMS_sk)   trcrn(:,:,nt_bgc_DMS_sk,:,:)   = c0
 
 
      !-------------------------------------------------------------------
      ! silicate
      !-------------------------------------------------------------------

         nbits = 64                ! double precision data

         if (trim(sil_data_type) == 'clim' .AND. (tr_bgc_Sil .OR. tr_bgc_Sil_sk)) then  !only Arctic climatology

            sil_file = trim(bgc_data_dir)//'silicate_WOA2005_surface_monthly' ! gx1 only

            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'silicate initialized from:'
               write (nu_diag,*) trim(sil_file)
            endif

            if (my_task == master_task) &
               call ice_open (nu_forcing, sil_file, nbits)

            k = month
            call ice_read (nu_forcing, k, work1, 'rda8', dbug, &
                           field_loc_center, field_type_scalar)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  sil(i,j,iblk) = work1(i,j,iblk)
               enddo
               enddo
            enddo

           if (my_task == master_task) close(nu_forcing)

         elseif (trim(sil_data_type) == 'rct_clim'  .AND. (tr_bgc_Sil .OR. tr_bgc_Sil_sk)) then
        
          !use WOA2005_surface (winter or spring) for a specific location
          !(Bering (60, 180), Okhotsk (55, 150E),  Chukchi (70, 170W) Labrador Sea (56, 50W), central(0,86)) 
          !           Just March: (25, 50, 30, 2.5,  20)
          !mmol/m^3  Apr, May, Jun spring range: (20, 40, 10, 2.5, 20)
          !          Jan, Feb, Mar winter range:  (20, 60, 25, 2.5, 20)
          
            sil(:,:,:) = 30.0_dbl_kind  !chukchi, march

         endif

      !-------------------------------------------------------------------
      ! nitrate
      !-------------------------------------------------------------------

         if (trim(nit_data_type) == 'clim' .AND. (tr_bgc_NO .or. tr_bgc_Nit_sk)) then

            nit_file = trim(bgc_data_dir)//'nitrate_WOA2005_surface_monthly' ! gx1 only

            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'nitrate initialized from:'
               write (nu_diag,*) trim(nit_file)
            endif

            if (my_task == master_task) &
               call ice_open (nu_forcing, nit_file, nbits)

            l = month
             call ice_read (nu_forcing, l, work1, 'rda8', dbug, &
                           field_loc_center, field_type_scalar)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  nit(i,j,iblk) = work1(i,j,iblk)                 
               enddo
               enddo
            enddo
            
            if (my_task == master_task) close(nu_forcing)

         elseif (trim(nit_data_type) == 'rct_clim'  .AND. (tr_bgc_NO .OR. tr_bgc_Nit_sk)) then
            
            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'nitrate initialized from March, Chukchi Sea'
            endif
             
            nit(:,:,:) = c10 

         elseif (trim(nit_data_type) == 'sss'  .AND. (tr_bgc_NO .OR. tr_bgc_Nit_sk)) then
           
            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'nitrate initialized from salinity'
            endif
             
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  nit(i,j,iblk) =  sss(i,j,iblk)        
               enddo
               enddo
            enddo
        endif

     if (tr_bgc_NO)   then 
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
               do k = 1, nblyr
                  trcrn(i,j,nt_bgc_NO+k-1,n,iblk) =   nit(i,j,iblk)*initbio_frac
               enddo      !k
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif                 
     if (tr_bgc_N)    then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_N+k-1,n,iblk)      = algalN(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif
     if (tr_bgc_C)   then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_C+k-1,n,iblk)      = R_C2N*algalN(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif
     if (tr_bgc_chl)  then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_chl+k-1,n,iblk)      = R_chl2N*algalN(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif  
     if (tr_bgc_NH )  then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_NH+k-1,n,iblk)      =  amm(i,j,iblk)*initbio_frac
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif 
     if (tr_bgc_Sil) then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_Sil+k-1,n,iblk)      = sil(i,j,iblk)*initbio_frac   
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif
     if (tr_bgc_DMSPp) then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_DMSPp+k-1,n,iblk)      = dmsp(i,j,iblk)
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif
     if (tr_bgc_DMSPd) then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_DMSPd+k-1,n,iblk)      = dmsp(i,j,iblk)*initbio_frac  
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif 

     if (tr_bgc_DMS)  then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
              trcrn(i,j,nt_bgc_DMS+k-1,n,iblk)      = dms(i,j,iblk)*initbio_frac  
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
     endif

       ! PON shadows algalN without reactions

         if (tr_bgc_PON )  then
        do  n = 1,ncat
        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block
           do k = 1,nblyr
               trcrn(i,j,nt_bgc_PON+k-1,n,iblk)      = nit(i,j,iblk)**initbio_frac   !algalN(i,j,iblk) !
           enddo
           enddo        !i
           enddo        !j
        enddo           !iblk
        enddo           !n 
         endif
              
      endif  ! restart_bgc

      end subroutine init_bgc

!=======================================================================

      subroutine biogeochemistry (dt, iblk)

      use ice_algae, only: z_biogeochemistry, skl_biogeochemistry
      use ice_flux
      use ice_blocks, only: block, get_block

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n               ! thickness category index

      integer (kind=int_kind), save :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), save :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         hin         , &  ! new ice thickness
         hsn         , &  !snow thickness  (m)
         hinS_old    , &  ! old brine thickness before growh/melt
         zphi_o      , &  ! surface ice porosity 
         kavg        , &  ! average ice permeability (m^2)
         sloss       , &  ! brine flux contribution from surface runoff (g/m^2)
         dh_bot_chl  , &  ! Chlorophyll may or may not flush
         dh_top_chl  , &  ! Chlorophyll may or may not flush
         darcy_V_chl

!--------------------------------
! Defined on the Bio Grid points
!--------------------------------

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2) :: &
         Sinn      , &   ! salinity on the bio grid  (ppt)
         brine_sal , &   ! brine salinity (ppt)
         brine_rho       ! brine_density (kg/m^3)

!-----------------------------------
! Defined on the Bio Grid interfaces
!-----------------------------------

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2) :: &
         iphin      , &   ! porosity 
         ibrine_sal , &   ! brine salinity  (ppt)
         ibrine_rho       ! brine_density (kg/m^3)

      real (kind=dbl_kind) :: &
         pond            ! flux of water retained in ponds (kg/m^2/s)

      ! for bgc sk
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         Iavgn,  & ! Iavg (W/m^2)
         grow_Cn,& ! C growth
         fN_partn  ! N down flux (mmol/m^2/s)

      ! for bgc layer
      real (kind=dbl_kind), dimension (nx_block,ny_block,nbltrcr) :: &
         flux_bion, & !tracer flux to ocean
         flux_bio_gn  !tracer flux to ocean from gravity drainage (mmol/m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         S_totn, &  ! ice Salinity  (g/m^2)
         hbrin , &  ! brine height
         melt_frac  ! fraction of top accumulation due to meltwater
         

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr) :: &
         upNOn,  &  ! nitrate uptake rate (mmol/m^3/s)
         upNHn      ! ammonium uptake rate (mmol/m^3/s)
    
      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      logical (kind=log_kind), dimension (nx_block,ny_block) :: &
         flood_val          ! .true. if ocean water flooded the surface
                            ! use sss rather than min_salin as surface 
                            ! condition

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

       !-----------------------------------------------------------------
       ! Biogeochemistry
       !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do n = 1, ncat

             ! initialize
             hin_old(:,:,n,iblk) = c0
             flux_bion(:,:,:) = c0
             flux_bio_gn(:,:,:) = c0
             do j = jlo, jhi
             do i = ilo, ihi
               if (aicen_init(i,j,n,iblk) > puny ) then !.AND. aicen(i,j,n,iblk) > puny) then
                 ! first_ice(i,j,n,iblk) = .false.
                  hin_old(i,j,n,iblk) = vicen_init(i,j,n,iblk)/aicen_init(i,j,n,iblk)
               else  ! initialize
                  first_ice(i,j,n,iblk) = .true.
                  if (n == 1) then
                     Rayleigh_criteria(i,j,iblk) = .false.
                  endif
                  if (hbrine)   trcrn(i,j,nt_fbri,n,iblk) = c1
                  if (tr_bgc_S) trcrn(i,j,nt_bgc_S:nt_bgc_S+nblyr-1,n,iblk) = c0
               endif
             enddo
             enddo
          enddo  !ncat
          
         !  Define ocean tracer concentration
          if (solve_skl_bgc .OR. tr_bgc_NO) then
          do j = 1, ny_block
          do i = 1, nx_block
           if (tr_bgc_NO .OR. tr_bgc_Nit_sk) ocean_bio(i,j,nlt_bgc_NO,iblk) = nit(i,j,iblk)
           if (tr_bgc_chl .OR. tr_bgc_chl_sk)  ocean_bio(i,j,nlt_bgc_chl,iblk) = R_chl2N*algalN(i,j,iblk)
           if (tr_bgc_NH .OR. tr_bgc_Am_sk)  ocean_bio(i,j,nlt_bgc_NH,iblk) = amm(i,j,iblk)
           if (tr_bgc_C .OR. tr_bgc_C_sk)  ocean_bio(i,j,nlt_bgc_C,iblk) = R_C2N*algalN(i,j,iblk)
           if (tr_bgc_Sil .OR. tr_bgc_Sil_sk)  ocean_bio(i,j,nlt_bgc_Sil,iblk) = sil(i,j,iblk)
           if (tr_bgc_DMSPp .OR. tr_bgc_DMSPp_sk)  ocean_bio(i,j,nlt_bgc_DMSPp,iblk) = dmsp(i,j,iblk)
           if (tr_bgc_DMSPd .OR. tr_bgc_DMSPd_sk)  ocean_bio(i,j,nlt_bgc_DMSPd,iblk) = dmsp(i,j,iblk)
           if (tr_bgc_DMS .OR. tr_bgc_DMS_sk)  ocean_bio(i,j,nlt_bgc_DMS,iblk) = dms(i,j,iblk)
           if (tr_bgc_N .OR. tr_bgc_N_sk)  ocean_bio(i,j,nlt_bgc_N,iblk) = algalN(i,j,iblk)
           if (tr_bgc_PON) & !mimics N but with no reactions 
               ocean_bio(i,j,nlt_bgc_PON,iblk) =nit(i,j,iblk) ! algalN(i,j,iblk) !
          enddo
          enddo
          endif

        do n = 1, ncat

          hsn(:,:) = c0
          hin(:,:) = c0
          Sinn(:,:,:) = c0
          upNOn(:,:,:) = c0
          upNHn(:,:,:) = c0 
          S_totn(:,:) = c0
          hbrin(:,:) = c0
          kavg(:,:) = c0
          zphi_o(:,:) = c0
          melt_frac(:,:) = c1
          flood_val(:,:) = .false.
          dh_top_chl(:,:) = c0
          darcy_V_chl(:,:) = c0

       
          icells = 0
          do j = jlo, jhi
          do i = ilo, ihi
            if (aicen(i,j,n,iblk) > puny) then
                icells = icells + 1
                indxi(icells) = i
                indxj(icells) = j
              !  if (maxval(trcrn(i,j,nt_bgc_S:nt_bgc_S+nblyr-1,n,iblk)) .LE. c0) then
              !      first_ice(i,j,n,iblk) = .true.
              !      if (hbrine) trcrn(i,j,nt_fbri,n,iblk) = c1
              !  endif
            endif
          enddo               ! i
          enddo               ! j

          if (icells > 0) then
          if (hbrine) then 

               call preflushing_changes (nx_block, ny_block, &
                                   icells, n, indxi,    indxj,      &    
                                   aicen (:,:,n,iblk),  &
                                   vicen(:,:,n,iblk) , vsnon(:,:,n,iblk), &
                                   meltbn(:,:,n,iblk), melttn(:,:,n,iblk), &
                                   congeln(:,:,n,iblk),snoicen(:,:,n,iblk), &
                                   hin_old(:,:,n,iblk), & 
                                   trcrn(:,:,nt_fbri,n,iblk), &
                                   dh_top(:,:,n,iblk), dh_bot(:,:,n,iblk), &
                                   dh_bot_chl,&
                                   dhi_top(:,:,n,iblk),dhi_bot(:,:,n,iblk), &
                                   hinS_old, hin, hsn,&
                                   first_ice(:,:,n,iblk))

               if (tr_bgc_S)  then
                 call compute_microS (nx_block, ny_block,   &
                                   icells, n, indxi,    indxj, &   
                                   trcrn(:,:,:,n,iblk), hin_old(:,:,n,iblk), &
                                   hinS_old, &
                                   sss(:,:,iblk),sst(:,:,iblk),    & 
                                   zTin(:,:,:,n,iblk), &
                                   zphi(:,:,:,n,iblk), &
                                   kavg, zphi_o,  & 
                                   Rayleigh_criteria(:,:,iblk),&
                                   first_ice(:,:,n,iblk), Sinn, &
                                   brine_sal, brine_rho, iphin, &
                                   ibrine_rho, ibrine_sal,  &
                                   sice_rho(:,:,n,iblk),sloss)
               else          
          !-------------------------------------------------------------
          !  NOTE: Requires the average ice permeability = kavg(:,:)
          !  and the surface ice porosity = zphi_o(:,:)
          !  computed in "compute_microS" or from "thermosaline_vertical"
          !--------------------------------------------------------------

                   call compute_microS_mushy (nx_block, ny_block,   &
                                   icells, n, indxi,    indxj, &   
                                   trcrn(:,:,:,n,iblk), hin_old(:,:,n,iblk), &
                                   hinS_old, &
                                   sss(:,:,iblk),sst(:,:,iblk),    & 
                                   zTin(:,:,:,n,iblk), &
                                   zphi(:,:,:,n,iblk), &
                                   kavg, zphi_o, Sinn, &
                                   brine_sal, brine_rho, iphin, &
                                   ibrine_rho, ibrine_sal,  &
                                   sice_rho(:,:,n,iblk),iDi(:,:,:,n,iblk))
               endif

 
               call update_hbrine (icells, nx_block, ny_block, &
                                   indxi,    indxj,            &
                                   meltbn(:,:,n,iblk),melttn(:,:,n,iblk), &
                                   meltsn(:,:,n,iblk), dt, hin, hsn, &
                                   hin_old(:,:,n,iblk),  first_ice(:,:,n,iblk), &
                                   hbrin, hinS_old,                &
                                   trcrn(:,:,nt_fbri,n,iblk),  &
                                   dh_top(:,:,n,iblk),dh_bot(:,:,n,iblk), &
                                   dh_top_chl, dh_bot_chl,kavg, zphi_o, &
                                   darcy_V(:,:,n,iblk),darcy_V_chl,&
                                   flood_val, melt_frac)


               if (tr_bgc_S .AND. solve_Sin) then

                                  call solve_zsalinity &
                                   (nx_block, ny_block,   &
                                   icells, n, dt, indxi, indxj, &  
                                   trcrn(:,:,nt_bgc_S:nt_bgc_S+nblyr-1,n,iblk), &
                                   trcrn(:,:,nt_qice:nt_qice+nilyr-1,n,iblk), &
                                   trcrn(:,:,nt_sice:nt_sice+nilyr-1,n,iblk), &
                                   aicen(:,:,n,iblk),  vicen(:,:,n,iblk)  ,&  
                                   Sinn, zTin(:,:,:,n,iblk),             &
                                   zphi(:,:,:,n,iblk), iphin,                 &
                                   iki(:,:,:,n,iblk), hinS_old, hbrin, hin, &
                                   hin_old(:,:,n,iblk), iDi(:,:,:,n,iblk),  darcy_V(:,:,n,iblk), &
                                   brine_sal, brine_rho, ibrine_sal, ibrine_rho, &
                                   Rayleigh_criteria(:,:,iblk), &
                                   first_ice(:,:,n,iblk), sss(:,:,iblk), sst(:,:,iblk),&
                                   dh_top(:,:,n,iblk),dh_bot(:,:,n,iblk),&  
                                   TLAT(:,:,iblk),TLON(:,:,iblk),        &
                                   l_stop, istop, jstop, fsicen(:,:,n,iblk), &
                                   fsicen_g(:,:,n,iblk), zphi_o,sloss, flood_val,&
                                   melt_frac)   

                  call column_sum_S (nx_block, ny_block, &
                                    icells,                         &
                                    indxi,              indxj,    &
                                    ntrcr, nblyr, vicen(:,:,n,iblk),  &
                                    trcrn(:,:,1:ntrcr,n,iblk),        &
                                    S_totn)
 
                  call merge_S_fluxes (nx_block,           ny_block,  &
                                    icells,                           &
                                    indxi,              indxj,        &
                                    aicen_init(:,:,n,iblk),           &
                                    S_totn, S_tot(:,:,iblk),          &
                                    hbrin, hbri(:,:,iblk),            &
                                    dsnown(:,:,n,iblk),               &
                                    dsnow  (:,:,iblk),                  &
                                    fsice(:,:,iblk), fsicen(:,:,n,iblk),&
                                    fsice_g(:,:,iblk), fsicen_g(:,:,n,iblk))         
                
               else
                    
                  call merge_hbrine (nx_block,           ny_block,    &
                                    icells,                           &
                                    indxi,              indxj,        &
                                    aicen_init(:,:,n,iblk),           &
                                    hbrin, hbri(:,:,iblk)             )   

               endif  !tr_bgc_S .AND. solve_Sin

        endif !hbrine
       !-----------------------------------------------------------------
       !    BGC
       !-----------------------------------------------------------------  
             
           if (tr_bgc_NO) then  
                 call z_biogeochemistry (nx_block, ny_block,                &
                                   icells, n,    dt,                        &
                                   indxi,    indxj,                         &  
                                   aicen (:,:,n,iblk),  vicen(:,:,n,iblk) , &
                                   hin_old(:,:,n,iblk),  nit(:,:,iblk)   ,  &
                                   flux_bion, flux_bio_gn,                  &
                                   amm(:,:,iblk)   ,  sil(:,:,iblk)   ,     &
                                   dmsp(:,:,iblk)   ,  dms(:,:,iblk)   ,    &
                                   algalN(:,:,iblk)   , zphi(:,:,:,n,iblk), &
                                   iphin, trcrn(:,:,1:ntrcr,n,iblk),        & 
                                   iDi(:,:,:,n,iblk),  sss(:,:,iblk),       &
                                   fswthruln(:,:,:,n,iblk),                 &
                                   growN(:,:,:,n,iblk), upNOn, upNHn,       &
                                   dh_top(:,:,n,iblk),dh_bot(:,:,n,iblk),   &
                                   dh_top_chl,dh_bot_chl,                   &
                                   zfswin(:,:,:,n,iblk),        &
                                   first_ice(:,:,n,iblk),                     &
                                   TLAT(:,:,iblk),TLON(:,:,iblk),           &
                                   hbrin, hinS_old,darcy_V(:,:,n,iblk),     &
                                   darcy_V_chl, bgrid, igrid, cgrid,        &
                                   flood_val, zphi_o)     

                 call merge_bgc_fluxes (nx_block,           ny_block,  &
                                    icells,                         &
                                    indxi,              indxj,    &
                                    aicen_init(:,:,n,iblk),          &
                                    vicen(:,:,n,iblk),  & !echmod
                                    ntrcr, & !echmod
                                    trcrn(:,:,1:ntrcr,n,iblk), & !echmod
                                    flux_bion, flux_bio(:,:,:,iblk), &
                                    flux_bio_gn, flux_bio_g(:,:,:,iblk), &
                                    upNOn,                           &
                                    upNHn, upNO(:,:,:,iblk),       &
                                    upNH(:,:,:,iblk), &
                                    chl_net(:,:,iblk), &
                                    PP_net(:,:,iblk), &
                                    NO_net(:,:,iblk),          &
                                    growN(:,:,:,n,iblk),  &
                                    growNp(:,:,:,iblk), grow_net(:,:,iblk))     
                     
              elseif (solve_skl_bgc) then
                 call skl_biogeochemistry (nx_block, ny_block,           &
                                   icells,  n,   dt,                     &
                                   indxi,    indxj,                      &  
                                   flux_bion, ocean_bio(:,:,:,iblk),     &
                                   hmix (:,:,iblk), aicen(:,:,n,iblk),   & 
                                   meltbn(:,:,n,iblk), congeln(:,:,n,iblk),  &
                                   fswthrun(:,:,n,iblk), first_ice(:,:,n,iblk),&
                                   trcrn(:,:,1:ntrcr,n,iblk),            &
                                   Iavgn, grow_Cn, fN_partn) 

                 call merge_bgc_fluxes_skl (nx_block, ny_block,  &
                                    icells,                      &
                                    indxi,              indxj,   &
                                    aicen_init(:,:,n,iblk),      &
                                    trcrn(:,:,nt_bgc_N_sk,n,iblk),  &  
                                    flux_bion, flux_bio(:,:,:,iblk), &
                                    PP_net(:,:,iblk),            &
                                    grow_net(:,:,iblk), grow_Cn)
              endif  

             
              do ij = 1, icells
                 i = indxi(ij)
                 j = indxj(ij)
               
                 first_ice(i,j,n,iblk) = .false.

              enddo  !icells
          endif                  ! icells      
        enddo                  ! ncat

      end subroutine biogeochemistry

!=======================================================================
!
!BOP
!
! !IROUTINE: merge_fluxes - aggregate flux information over ITD
!
! !INTERFACE:
!
      subroutine merge_bgc_fluxes (nx_block, ny_block,   &
                               icells,               &
                               indxi,    indxj,      &
                               aicen,                &    
                               vicen,  & !echmod
                               ntrcr, trcrn, &  !echmod
                               flux_bion, flux_bio, &
                               flux_bio_gn, flux_bio_g, &
                               upNOn,  &
                               upNHn, upNO, upNH, &
!echmod                               chl_netn, & !echmod
                               chl_net, &
!echmod                               PP_netn,  !echmod
                               PP_net, &
!echmod                               NO_netn, & !echmod
                               NO_net,growNn, growNp, grow_net)
!
! !DESCRIPTION:
!
! Aggregate flux information from all ice thickness categories
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! !USES:
!
!      use ice_state, only: nt_fbri, nt_bgc_N, nt_bgc_NO
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block, & ! block dimensions
          icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
          intent(in) :: &
          indxi, indxj    ! compressed indices for cells with aicen > puny

      ! single category fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &          
          aicen   !, & ! concentration of ice !echmod
!echmod          chl_netn, & ! net chlorophyll  (mg chl) !echmod
!echmod          PP_netn , & ! net PP (mg C/s) !echmod
!echmod          NO_netn     ! net NO (mmol/s) !echmod

      ! single category fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block,nbltrcr), intent(in):: &          
           flux_bion, &
           flux_bio_gn

      ! cumulative fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block,nbltrcr), intent(inout):: &          
           flux_bio, &
           flux_bio_g

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: & 
          chl_net, & ! net chlorophyll  (mg chl )  in history divide by aice
          PP_net , & ! net PP (mg C/m^2/s)  in history divide by aice
          NO_net , & ! net NO (mmol/m^2 ) in history divide by aice
          grow_net   ! net specific growth (/s) in history divide by aice
           

      real (kind=dbl_kind), dimension(nx_block,ny_block,nblyr), &
          intent(inout):: &          
          upNO   , & ! tot nitrate uptake rate (mmol/m^3/s)  in history divide by aice 
          upNH   , & ! tot ammonium uptake rate (mmol/m^3/s) in history divide by aice 
          growNp     ! specific growth  (/s) in history divide by aice 

      real (kind=dbl_kind), dimension(nx_block,ny_block,nblyr), &
          intent(in):: &          
          upNOn   , & ! nitrate uptake rate per cat (mmol/m^3/s)
          upNHn   , &    ! ammonium uptake rate per cat (mmol/m^3/s)
          growNn      ! algal growth rate per cat (/s)

! !INPUT/OUTPUT PARAMETERS: column_sum !echmod
!
      integer (kind=int_kind), intent(in) :: &
         ntrcr             ! number of tracers
      
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &          
          vicen       !volume of ice

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
           intent(in) :: &
           trcrn         ! input fields

!local
      real (kind=dbl_kind) :: &
          tmp

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: & !, intent(out) :: &
          chl_netn, & ! net chlorophyll  (mg chl)
          PP_netn, &  ! net PP (mg C/m^2/s)   * aicen
          NO_netn, &  ! net NO (mmol/m^2)   * aicen
          grow_netn   ! net specific growth (m/s) * aicen 

      integer (kind=int_kind) :: &
          ij, i, j, &   ! horizontal indices
          k             ! tracer indice

      !-----------------------------------------------------------------
      ! Column summation
      !-----------------------------------------------------------------

     chl_netn(:,:) = c0
     PP_netn(:,:) = c0
     NO_netn(:,:) = c0
     grow_netn(:,:) = c0

      if (solve_zbgc) then
      do k = 1, nblyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            tmp = c1
            if (hbrine) tmp = trcrn(i,j,nt_fbri)
           
            chl_netn(i,j) = chl_netn(i,j) + trcrn(i,j,nt_bgc_N+k-1)*R_chl2N  &  
                         * vicen(i,j)*tmp/real(nblyr,kind=dbl_kind)
            PP_netn(i,j)  = PP_netn(i,j) + trcrn(i,j,nt_bgc_N+k-1)*growNn(i,j,k)*(c1-fr_resp)*&
                           R_C2N*R_gC2molC * vicen(i,j)*tmp/real(nblyr,kind=dbl_kind) 
            NO_netn(i,j) = NO_netn(i,j) + trcrn(i,j,nt_bgc_NO+k-1) &
                         * vicen(i,j)*tmp/real(nblyr,kind=dbl_kind)                     
            grow_netn(i,j) = grow_netn(i,j) + growNn(i,j,k) &
                         * aicen(i,j)*tmp/real(nblyr,kind=dbl_kind)                     
         enddo                  ! ij
      enddo                     ! n
      else

      do k = 1, nblyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
           
            tmp = c1
            if (hbrine) tmp = trcrn(i,j,nt_fbri)

            NO_netn(i,j) = NO_netn(i,j) + trcrn(i,j,nt_bgc_NO+k-1) &
                         * vicen(i,j)*tmp/real(nblyr,kind=dbl_kind)   
         enddo                  ! ij
      enddo     
      endif                     ! solve_zbgc

      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         do k = 1,nbltrcr
           flux_bio (i,j,k)  = flux_bio(i,j,k) + flux_bion(i,j,k)*aicen(i,j)
           flux_bio_g (i,j,k)= flux_bio_g(i,j,k)+flux_bio_gn(i,j,k)*aicen(i,j)
         enddo
         chl_net  (i,j)  = chl_net (i,j) + chl_netn(i,j) 
         PP_net   (i,j)  = PP_net  (i,j) + PP_netn (i,j)  
         NO_net   (i,j)  = NO_net  (i,j) + NO_netn (i,j)  
         grow_net (i,j)  = grow_net(i,j) + grow_netn(i,j) 

         do k = 1,nblyr
            upNO (i,j,k)  = upNO (i,j,k) + upNOn   (i,j,k)*aicen(i,j)
            upNH (i,j,k)  = upNH (i,j,k) + upNHn   (i,j,k)*aicen(i,j)
            growNp (i,j,k)  = growNp (i,j,k) + growNn  (i,j,k)*aicen(i,j)
         enddo                  !k

      enddo                     ! ij
      
      end subroutine merge_bgc_fluxes

!=======================================================================
!
!BOP
!
! !IROUTINE: merge_fluxes - aggregate flux information over ITD
!
! !INTERFACE:
!
      subroutine merge_bgc_fluxes_skl (nx_block, ny_block,   &
                               icells,               &
                               indxi,    indxj,      &
                               aicen,  algal_N,        &
                               flux_bion, flux_bio,  &
                               PP_net, grow_net,     &
                               grow_Cn)
!
! !DESCRIPTION:
!
! Aggregate flux information from all ice thickness categories
! for skeletal layer biogeochemistry
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block, & ! block dimensions
          icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
          intent(in) :: &
          indxi, indxj    ! compressed indices for cells with aicen > puny

      ! single category fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &          
          aicen   !

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          algal_N       ! (mmol N/m^2) in the bottom layer
     
      ! single category fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block,nbltrcr), intent(in):: &          
          flux_bion

      ! cumulative fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block,nbltrcr), intent(inout):: &          
          flux_bio

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(in):: & 
          grow_Cn   !Specific growth (/s) 

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: & 
          PP_net , & ! Bulk net PP (mg C/m^2/s) in history divide by aice
          grow_net   ! net specific growth (/s) in history divide by aice
      

! !INPUT/OUTPUT PARAMETERS: column_sum !echmod
!
 
      integer (kind=int_kind) :: &
          ij, i, j, &   ! horizontal indices
          k             ! tracer indice
    
      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
   
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         do k = 1,nbltrcr
           flux_bio (i,j,k)  = flux_bio(i,j,k) + flux_bion(i,j,k)*aicen(i,j)
         enddo
        
         PP_net   (i,j)  = PP_net  (i,j) + algal_N(i,j)*phi_sk*grow_Cn(i,j)*(c1-fr_resp)*&
                           R_C2N*R_gC2molC *aicen(i,j)
         grow_net (i,j)  = grow_net(i,j) + grow_Cn(i,j) * phi_sk*aicen(i,j) 

      enddo                     ! ij
      
      end subroutine merge_bgc_fluxes_skl

!=======================================================================
!BOP
!
! !IROUTINE: init_history_bgc - initialize bgc history fields
!
! !INTERFACE:
!
      subroutine init_history_bgc
!
! !DESCRIPTION:
!
! Initialize bgc fields written to history files.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      upNO     (:,:,:,:) = c0
      upNH     (:,:,:,:) = c0
      S_tot    (:,:,:)   = c0
      chl_net  (:,:,:)   = c0
      PP_net   (:,:,:)   = c0
      NO_net   (:,:,:)   = c0
      grow_net (:,:,:)   = c0
      hbri     (:,:,:)   = c0
      growNp   (:,:,:,:) = c0
      growN    (:,:,:,:,:) = c0
      flux_bio (:,:,:,:) = c0
      flux_bio_g(:,:,:,:) = c0
      flux_bio_gbm  (:,:,:,:) = c0
      flux_bio_g_gbm  (:,:,:,:) = c0

      end subroutine init_history_bgc

!=======================================================================

! Adjust biogeochemical tracers when new frazil ice forms

      subroutine add_new_ice_bgc (nx_block,   ny_block,             &
                                  icells,     jcells,     kcells,   &
                                  indxi,      indxj,                &
                                  indxi2,     indxj2,     indxij2,  &
                                  indxi3,     indxj3,     indxij3,  &
                                  aicen_init, vicen_init, vi0_init, &
                                  aicen,      vicen,      vi0new,   &
                                  ntrcr,      trcrn,      nbltrcr,  &
                                  sss,        ocean_bio,            &
                                  hsurp,      &
                                  l_stop,     istop,      jstop)

      use ice_itd, only: column_sum, &
                         column_conservation_check
      use ice_therm_shared, only: solve_Sin

      integer (kind=int_kind), intent(in) :: &
         nx_block, & ! block dimensions
         ny_block, & ! block dimensions
         ntrcr   , & ! number of tracers in use
         icells  , & ! number of ice/ocean grid cells
         jcells  , & ! grid cell counter
         kcells      ! grid cell counter

      integer (kind=int_kind), intent(in) :: &
         nbltrcr     ! number of biology tracers

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj, &           ! compressed i/j indices
         indxi2, indxj2, indxij2, & ! compressed i/j indices
         indxi3, indxj3, indxij3    ! compressed i/j indices

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen_init  , & ! initial concentration of ice
         vicen_init  , & ! intiial volume per unit area of ice  (m)
         aicen       , & ! concentration of ice
         vicen           ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), &
         intent(inout) :: &
         trcrn           ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         sss             ! sea surface salinity (ppt)

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         vi0_init    , & ! volume of new ice added to cat 1 (intial)
         vi0new          ! volume of new ice added to cat 1

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbltrcr), &
         intent(in) :: &
         ocean_bio       ! ocean concentration of biological tracer

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hsurp           ! thickness of new ice added to each cat

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop    ! indices of grid cell where model aborts

! local
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         n           , & ! ice category index
         k           , & ! ice layer index
         ij, m           ! combined i/j horizontal indices

      real (kind=dbl_kind), dimension (icells) :: &
         vbri1       , & ! starting volume of existing brine
         vbri_init   , & ! brine volume summed over categories
         vbri_final      ! brine volume summed over categories

      real (kind=dbl_kind), dimension(icells) :: &
         vsurp       , & ! volume of new ice added to each cat
         vtmp            ! total volume of new and old ice
        
      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) :: &
         vbrin           ! trcrn(i,j,nt_fbri,n)*vicen(i,j,n) 

      character (len=char_len) :: &
         fieldid         ! field identifier

      vbrin(:,:,:) = c0

      do n = 1, ncat
      do ij = 1,icells
         i = indxi(ij)
         j = indxj(ij)
         vbrin(i,j,n) = vicen_init(i,j,n)
         if (hbrine) vbrin(i,j,n) =  trcrn(i,j,nt_fbri,n)*vicen_init(i,j,n)
      enddo
      enddo

      call column_sum (nx_block, ny_block,       &
                       icells,   indxi,   indxj, &
                       ncat,                     &
                       vbrin,    vbri_init)
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         vbri_init(ij) = vbri_init(ij) + vi0_init(ij) !bgc
      enddo

      !-----------------------------------------------------------------
      ! kcells:  
      ! Distribute bgc in new ice volume among all ice categories by 
      ! increasing ice thickness, leaving ice area unchanged.
      !-----------------------------------------------------------------

      do n = 1,ncat
 
         ! Diffuse_bio handles concentration changes from ice growth/melt
         ! ice area does not change for kcells
         ! add salt to the bottom 

        do ij = 1, kcells
           i = indxi3(ij)
           j = indxj3(ij)
           m = indxij3(ij)

           vtmp(m) = vbrin(i,j,n)
           vsurp(m) = hsurp(m) * aicen_init(i,j,n) 
           vbrin(i,j,n) = vbrin(i,j,n) + vsurp(m)
           if (hbrine) then
              trcrn(i,j,nt_fbri,n) = c1
              if (vicen(i,j,n) > c0)  trcrn(i,j,nt_fbri,n) = vbrin(i,j,n)/vicen(i,j,n)
           endif
        enddo

         if ((tr_bgc_S .OR. tr_bgc_NO) .AND. kcells > 0 ) then
            call adjust_tracer_profile(nx_block, ny_block,           &
                                       indxi3,   indxj3,   indxij3,  &                        
                                       icells,   kcells,             &
                                       nbltrcr,  ntrcr,              &
                                       aicen_init(:,:,n),            &
                                       vbrin(:,:,n),                 &
                                       vicen(:,:,n),                 &
                                       trcrn(:,:,:,n),               &
                                       vtmp,     vsurp,    sss,      &
                                       nilyr,    nblyr,    solve_Sin,& 
                                       bgrid,    cgrid,    ocean_bio)
         endif           ! tr_bgc_S
      enddo              ! n

      !-----------------------------------------------------------------
      ! jcells:  
      ! Combine bgc in new ice grown in open water with category 1 ice.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, jcells
         i = indxi2(ij)
         j = indxj2(ij)
         m = indxij2(ij)

         vbri1(m)     = vbrin(i,j,1) 
         vbrin(i,j,1) = vbrin(i,j,1) + vi0new(m)
         if (hbrine) then
            trcrn(i,j,nt_fbri,1) = c1
            if (vicen(i,j,1) > c0) trcrn(i,j,nt_fbri,1) = vbrin(i,j,1)/vicen(i,j,1)
         endif
      enddo

      ! Diffuse_bio handles concentration changes from ice growth/melt
      ! ice area changes for jcells
      ! add salt throughout

      if ((tr_bgc_S .OR. tr_bgc_NO) .AND. jcells > 0 ) then
            call adjust_tracer_profile(nx_block, ny_block, &
                                       indxi2,   indxj2,   indxij2,  &                        
                                       icells,   jcells,             &
                                       nbltrcr,  ntrcr,              &
                                       aicen(:,:,1),                 &
                                       vbrin(:,:,1),                 &
                                       vicen(:,:,1),                 &
                                       trcrn(:,:,:,1),               &
                                       vbri1,    vi0new,   sss,      &
                                       nilyr,    nblyr,    solve_Sin,& 
                                       bgrid,    cgrid,    ocean_bio)
      endif           ! tr_bgc_S
  
        
      if (hbrine) then
         call column_sum (nx_block, ny_block,       &
                          icells,   indxi,   indxj, &
                          ncat,                     &
                          vbrin,    vbri_final)

         fieldid = 'vbrin, add_new_ice'
         call column_conservation_check (nx_block,  ny_block,      &
                                         icells,   indxi,   indxj, &
                                         fieldid,                  &
                                         vbri_init, vbri_final,    &
                                         puny,      l_stop,        &
                                         istop,     jstop)
         if (l_stop) return
      endif

      end subroutine add_new_ice_bgc

!=======================================================================
!BOP
!
! !IROUTINE: adjust_tracer_profile
!
! !INTERFACE:
!
      subroutine adjust_tracer_profile (nx_block,ny_block,        &
                                   indxi,   indxj,  indxij,          &
                                   icells, kcells,  nbltrcr,                    &
                                   ntrcr, aicen, vbrin, vicen,          &
                                   trcrn, vtmp, vsurp, sss, &
                                   nilyr,     & 
                                   nblyr, solve_Sin, bgrid, cgrid,       &
                                   ocean_bio)
!
! !DESCRIPTION:
!
! Add tracers to the bottom and adjust the vertical profile  
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
         kcells            , & ! cells with frazil accumulation
         ntrcr             , & ! number of tracers in use
         nilyr             , & ! number of ice layers
         nbltrcr           , & ! number of biology tracers
         nblyr                 ! number of biology layers

      integer (kind=int_kind), dimension (icells), intent(in) :: &
         indxi, indxj    , & ! indices for i/j directions
         indxij

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume of ice
         sss         ! ocean salinity (ppt)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbltrcr), &
         intent(in) :: &
         ocean_bio

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         vbrin       ! fbri*volume per unit area of ice  (m)
       
      logical (kind=log_kind), intent(in) :: &
         solve_Sin 

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
!
!EOP
!
      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: &
         trtmp0      ! temporary, remapped tracers     !need extra 
        
      real (kind=dbl_kind), dimension (kcells) :: &
         hin     , & ! ice height
         hinS_new, & ! brine height
         temp_S         

      integer (kind=int_kind) :: &
         i, j, m, &
         k, ij ! grid cell counters

      trtmp0(:,:,:) = c0

      do k = 1, nblyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, kcells
            i = indxi(ij)
            j = indxj(ij)
            m = indxij(ij)

            if (vbrin(i,j) > c0 ) then

               if (tr_bgc_S) then
                   trcrn(i,j,nt_bgc_S+k-1)           =  &
                  (trcrn(i,j,nt_bgc_S+k-1)           * vtmp(m) &
                                + sss(i,j)*salt_loss * vsurp(m)) / vbrin(i,j)
                  if (solve_Sin) trtmp0(i,j,nt_sice+k-1) = trcrn(i,j,nt_bgc_S+k-1)
               endif
               if (tr_bgc_NO) &
                   trcrn(i,j,nt_bgc_NO+k-1)          = &
                  (trcrn(i,j,nt_bgc_NO+k-1)          * vtmp(m) &
                      +  ocean_bio(i,j,nlt_bgc_NO)* initbio_frac   * vsurp(m)) / vbrin(i,j) 
               if (solve_zbgc) then
               if (tr_bgc_chl) &
                   trcrn(i,j,nt_bgc_chl+k-1)         = &
                  (trcrn(i,j,nt_bgc_chl+k-1)         * vtmp(m) &
                      +  ocean_bio(i,j,nlt_bgc_chl)  * vsurp(m)) / vbrin(i,j) 
               if (tr_bgc_NH) &
                   trcrn(i,j,nt_bgc_NH+k-1)          = &
                  (trcrn(i,j,nt_bgc_NH+k-1)          * vtmp(m) &
                      +  ocean_bio(i,j,nlt_bgc_NH)* initbio_frac   * vsurp(m)) / vbrin(i,j) 
               if (tr_bgc_C) &
                   trcrn(i,j,nt_bgc_C+k-1)           = &
                  (trcrn(i,j,nt_bgc_C+k-1)           * vtmp(m) &
                      +  ocean_bio(i,j,nlt_bgc_C)    * vsurp(m)) / vbrin(i,j) 
               if (tr_bgc_Sil) &
                   trcrn(i,j,nt_bgc_Sil+k-1)         = &
                  (trcrn(i,j,nt_bgc_Sil+k-1)         * vtmp(m) &
                      +  ocean_bio(i,j,nlt_bgc_Sil)* initbio_frac  * vsurp(m)) / vbrin(i,j) 
               if (tr_bgc_DMSPp) &
                   trcrn(i,j,nt_bgc_DMSPp+k-1)       = &
                  (trcrn(i,j,nt_bgc_DMSPp+k-1)       * vtmp(m) &
                      +  ocean_bio(i,j,nlt_bgc_DMSPp)* vsurp(m)) / vbrin(i,j) 
               if (tr_bgc_DMSPd) &
                   trcrn(i,j,nt_bgc_DMSPd+k-1)       = &
                  (trcrn(i,j,nt_bgc_DMSPd+k-1)       * vtmp(m) &
                      +  ocean_bio(i,j,nlt_bgc_DMSPd) * initbio_frac* vsurp(m)) / vbrin(i,j)
               if (tr_bgc_DMS) &
                   trcrn(i,j,nt_bgc_DMS+k-1)         = &
                  (trcrn(i,j,nt_bgc_DMS+k-1)         * vtmp(m) &
                      +  ocean_bio(i,j,nlt_bgc_DMS) * initbio_frac  * vsurp(m)) / vbrin(i,j)
               if (tr_bgc_N) &
                   trcrn(i,j,nt_bgc_N+k-1)           = &
                  (trcrn(i,j,nt_bgc_N+k-1)           * vtmp(m) &
                      +  ocean_bio(i,j,nlt_bgc_N)    * vsurp(m)) / vbrin(i,j) 
               if (tr_bgc_PON) & !mimics NO but with no reactions
                   trcrn(i,j,nt_bgc_PON+k-1)         = &
                  (trcrn(i,j,nt_bgc_PON+k-1)         * vtmp(m) &
                      +  ocean_bio(i,j,nlt_bgc_PON) * initbio_frac  * vsurp(m)) / vbrin(i,j)
               endif  !solve_zbgc
            endif  
         enddo        ! ij
      enddo        ! k   

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
         temp_S   (ij) = min_salin 
      enddo
           
      if (solve_Sin) &
         call remap_layers_bgc_plus (nx_block,ny_block,        &
                                     indxi,   indxj,           &
                                     kcells,                   &
                                     ntrcr,                    &
                                     nilyr,                    &
                                     nt_sice,                  &
                                     trtmp0, trcrn, &
                                     1,        nblyr,          &
                                     hin(:),  hinS_new(:),         &
                                     cgrid(2:nilyr+1),         &
                                     bgrid(2:nblyr+1),temp_S)

      end subroutine adjust_tracer_profile

!=======================================================================

      end module ice_zbgc

!=======================================================================
