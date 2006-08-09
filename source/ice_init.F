c $Id: ice_init.F,v 1.9 2004/02/09 17:58:46 lipscomb Exp $
!=======================================================================
!BOP
!
! !MODULE:   ice_init - parameter and variable initializations
!
! !DESCRIPTION:
!
! parameter and variable initializations
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2004: Block structure added by William Lipscomb
!
! !INTERFACE:
!
      module ice_init
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_domain_size
      use ice_constants
!
!EOP
!
      implicit none
      save

      character (len=char_len) ::
     &   advection   ! type of advection scheme used
                     ! 'upwind' => 1st order donor cell scheme
                     ! 'remap' => remapping scheme

      character(len=char_len) :: 
     &   ice_ic      ! method of ice cover initialization
                     ! 'default'  => latitude and sst dependent
                     ! 'none'     => no ice 
                     ! note:  restart = .true. overwrites

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: input_data - namelist variables
!
! !INTERFACE:
!
      subroutine input_data
!
! !DESCRIPTION:
!
! Namelist variables, set to default values; may be altered
! at run time
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_broadcast
      use ice_diagnostics
      use ice_fileunits
      use ice_calendar, only: year_init, istep0, histfreq, diagfreq,
     &                        dumpfreq, dumpfreq_n, npt, dt, ndyn_dt
      use ice_restart, only: restart, restart_dir, restart_file,
     &     pointer_file
      use ice_history, only: hist_avg, history_dir, history_file,
     &     incond_dir, incond_file 
      use ice_exit
      use ice_itd, only: kitd, kcatbound
      use ice_ocean, only: oceanmixed_ice
      use ice_forcing, only:
     &    ycycle,          fyear_init,    dbug,
     &    atm_data_type,   atm_data_dir,  precip_units, 
     &    sss_data_type,   sst_data_type, ocn_data_dir,
     &    oceanmixed_file, restore_sst,   trestore 
      use ice_grid, only: grid_file, kmt_file, grid_type 
      use ice_mechred, only: kstrength, krdg_partic, krdg_redist
      use ice_dyn_evp, only: ndte, kdyn, evp_damping, yield_curve
      use ice_shortwave, only: albicev, albicei, albsnowv, albsnowi,
     &                         shortwave
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag

      character (len=6) :: chartmp

#ifdef CCSM
      !-----------------------------------------------------------------
      ! Declare namelist variables used for CCSM coupled runs and not
      ! declared elsewhere in CICE
      !-----------------------------------------------------------------

      integer (kind=int_kind) :: kcolumn  ! 1 for column model

      character(len=char_len_long) :: 
     &   runid                 ! identifier for coupled run

#endif

      !-----------------------------------------------------------------
      ! Namelist variables.
      ! NOTE: Not all of these are used by both models.
      !-----------------------------------------------------------------

      namelist /ice_nml/
     &  year_init,      istep0,          dt,            npt
     &, diagfreq       
     &, print_points,   print_global,    diag_type,     diag_file
     &, histfreq,       hist_avg,        history_dir,   history_file
     &, dumpfreq,       dumpfreq_n,      restart_file
     &, restart,        restart_dir,     pointer_file,  ice_ic
     &, grid_type,      grid_file,       kmt_file
     &, kitd,           kcatbound
     &, kdyn,           ndyn_dt,         ndte,          evp_damping
     &, yield_curve,    advection
     &, kstrength,      krdg_partic,     krdg_redist,   shortwave
     &, albicev,        albicei,         albsnowv,      albsnowi
     &, fyear_init,     ycycle 
     &, atm_data_type,  atm_data_dir,    precip_units
     &, oceanmixed_ice, sss_data_type,   sst_data_type 
     &, ocn_data_dir,   oceanmixed_file, restore_sst,   trestore
     &, latpnt,         lonpnt,          dbug
#ifdef CCSM
! These variables are used in CCSM, but not CICE
     &, runid,           runtype,             kcolumn
#endif

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------

      year_init = 1          ! initial year
      istep0 = 0             ! no. of steps taken in previous integrations,
                             ! real (dumped) or imagined (to set calendar)
      dt = 3600.0_dbl_kind   ! time step, s 
      npt = 99999            ! total number of time steps (dt) 
      diagfreq = 24          ! how often diag output is written
      print_points = .false. ! if true, print point data
      print_global = .true.  ! if true, print global diagnostic data
      diag_type = 'stdout'
      diag_file = 'ice_diag.d'
      histfreq='m'           ! output frequency
      hist_avg = .true.      ! if true, write time-averages (not snapshots)
      history_dir  = ' '     ! Write to executable dir for default
      history_file = 'iceh'  ! history file name prefix
      dumpfreq='y'           ! restart frequency option
      dumpfreq_n= 1          ! restart frequency
      restart = .false.      ! if true, read restart files for initialization
      restart_dir  = ' '     ! Write to executable dir for default
      restart_file = 'iced'  ! restart file name prefix
      pointer_file = 'ice.restart_file'
      ice_ic       = 'default'       ! latitude and sst-dependent
      grid_type    = 'rectangular'   ! define rectangular grid internally
      grid_file    = 'unknown_grid_file'
      kmt_file     = 'unknown_kmt_file'

      kitd = 1           ! type of itd conversions (0 = delta, 1 = linear)
      kcatbound = 1      ! category boundary formula (0 = old, 1 = new)
      kdyn = 1           ! type of dynamics (1 = evp)
      ndyn_dt = 1        ! dynamic time steps per thermodynamic time step
      ndte = 120         ! subcycles per dynamics timestep:  ndte=dyn_dt/dte
      evp_damping = .false.  ! if true, use damping procedure in evp dynamics
      yield_curve = 'ellipse'
      kstrength = 1          ! 1 = Rothrock 75 strength, 0 = Hibler 79
      krdg_partic = 1        ! 1 = new participation, 0 = Thorndike et al 75
      krdg_redist = 1        ! 1 = new redistribution, 0 = Hibler 80
      advection  = 'remap'   ! incremental remapping transport scheme
      shortwave = 'ccsm3'
      albicev   = 0.78_dbl_kind   ! visible ice albedo for h > ahmax
      albicei   = 0.36_dbl_kind   ! near-ir ice albedo for h > ahmax
      albsnowv  = 0.98_dbl_kind   ! cold snow albedo, visible
      albsnowi  = 0.70_dbl_kind   ! cold snow albedo, near IR

      fyear_init = 1900           ! first year of forcing cycle
      ycycle = 1                  ! number of years in forcing cycle
      atm_data_type   = 'default'
      atm_data_dir    = ' '
      precip_units    = 'mks'     ! 'mm_per_month' or
                                  ! 'mm_per_sec' = 'mks' = kg/m^2 s
      oceanmixed_ice  = .false.   ! if true, use internal ocean mixed layer
      sss_data_type   = 'default'
      sst_data_type   = 'default'
      ocn_data_dir    = ' '
      oceanmixed_file = 'unknown_oceanmixed_file' ! ocean forcing data
      restore_sst     = .false.   ! restore sst if true
      trestore        = 90        ! restoring timescale, days (0 instantaneous)
      dbug      = .false.         ! true writes diagnostics for input forcing

      latpnt(1) =  90._dbl_kind   ! latitude of diagnostic point 1 (deg)
      lonpnt(1) =   0._dbl_kind   ! longitude of point 1 (deg)
      latpnt(2) = -65._dbl_kind   ! latitude of diagnostic point 2 (deg)
      lonpnt(2) = -45._dbl_kind   ! longitude of point 2 (deg)

#ifdef CCSM
      ! The following are in CCSM but not CICE.
      runid   = 'unknown'   ! default run ID
      runtype = 'unknown'   ! default runtype
      kcolumn = 0           ! 1 = column model
#endif

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif
         do while (nml_error > 0)
            read(nu_nml, nml=ice_nml,iostat=nml_error)
            if (nml_error > 0) read(nu_nml,*)  ! for Nagware compiler
         end do
         if (nml_error == 0) close(nu_nml)
      endif

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         call abort_ice('ERROR reading ice_nml')
      endif

      if (trim(diag_type) == 'file') then
         nu_diag = 48
      else
         nu_diag = 6
      endif

      if (histfreq == '1') hist_avg = .false. ! potential conflict
#ifdef CCSM
      if (kcolumn == 1) grid_type    = 'column'
#endif
      chartmp = advection(1:6)
      if (chartmp /= 'upwind') advection = 'remap'

      call broadcast_scalar(year_init,          master_task)
      call broadcast_scalar(istep0,             master_task)
      call broadcast_scalar(dt,                 master_task)
      call broadcast_scalar(npt,                master_task)
      call broadcast_scalar(diagfreq,           master_task)
      call broadcast_scalar(print_points,       master_task)
      call broadcast_scalar(print_global,       master_task)
      call broadcast_scalar(diag_type,          master_task)
      call broadcast_scalar(diag_file,          master_task)
      call broadcast_scalar(histfreq,           master_task)
      call broadcast_scalar(hist_avg,           master_task)
      call broadcast_scalar(history_dir,        master_task)
      call broadcast_scalar(history_file,       master_task)
      call broadcast_scalar(dumpfreq,           master_task)
      call broadcast_scalar(dumpfreq_n,         master_task)
      call broadcast_scalar(restart_file,       master_task)
      call broadcast_scalar(restart,            master_task)
      call broadcast_scalar(restart_dir,        master_task)
      call broadcast_scalar(pointer_file,       master_task)
      call broadcast_scalar(ice_ic,             master_task)
      call broadcast_scalar(grid_type,          master_task)
      call broadcast_scalar(grid_file,          master_task)
      call broadcast_scalar(kmt_file,           master_task)
      call broadcast_scalar(kitd,               master_task)
      call broadcast_scalar(kcatbound,          master_task)
      call broadcast_scalar(kdyn,               master_task)
      call broadcast_scalar(ndyn_dt,            master_task)
      call broadcast_scalar(ndte,               master_task)
      call broadcast_scalar(evp_damping,        master_task)
      call broadcast_scalar(yield_curve,        master_task)
      call broadcast_scalar(kstrength,          master_task)
      call broadcast_scalar(krdg_partic,        master_task)
      call broadcast_scalar(krdg_redist,        master_task)
      call broadcast_scalar(advection,          master_task)
      call broadcast_scalar(shortwave,          master_task)
      call broadcast_scalar(albicev,            master_task)
      call broadcast_scalar(albicei,            master_task)
      call broadcast_scalar(albsnowv,           master_task)
      call broadcast_scalar(albsnowi,           master_task)

      call broadcast_scalar(fyear_init,         master_task)
      call broadcast_scalar(ycycle,             master_task)
      call broadcast_scalar(atm_data_type,      master_task)
      call broadcast_scalar(atm_data_dir,       master_task)
      call broadcast_scalar(precip_units,       master_task)
      call broadcast_scalar(oceanmixed_ice,     master_task)
      call broadcast_scalar(sss_data_type,      master_task)
      call broadcast_scalar(sst_data_type,      master_task)
      call broadcast_scalar(ocn_data_dir,       master_task)
      call broadcast_scalar(oceanmixed_file,    master_task)
      call broadcast_scalar(restore_sst,        master_task)
      call broadcast_scalar(trestore,           master_task)
      call broadcast_scalar(dbug,               master_task)
      call broadcast_array (latpnt(1:2),        master_task)
      call broadcast_array (lonpnt(1:2),        master_task)
! only master_task writes to file
!      call broadcast_scalar(nu_diag),           master_task)

#ifdef CCSM
      ! The following are in CCSM but not CICE.
      call broadcast_scalar (runid,             master_task)
      call broadcast_scalar (runtype,           master_task)
      call broadcast_scalar (kcolumn,           master_task)
#endif

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------

      if (my_task == master_task) then

         if (trim(diag_type) == 'file') then
            write(6,*) 'Diagnostic output will be in file ', diag_file
            open (nu_diag, file=diag_file, status='unknown')
         endif

         write (nu_diag,*) '--------------------------------'
         write (nu_diag,*) '  CICE model diagnostic output  '
         write (nu_diag,*) '--------------------------------'
         write(nu_diag,*) ' '
         write(nu_diag,*) ' Document ice_in namelist parameters:'
         write(nu_diag,*) ' ==================================== '
         write(nu_diag,*) ' '
#ifdef CCSM
         if (trim(runid) /= 'unknown')
     &    write(nu_diag,*)    ' runid                     = ',
     &                         trim(runid)
         if (trim(runtype) /= 'unknown')
     &    write(nu_diag,1030) ' runtype                   = ',
     &                         trim(runtype)
#endif
         write(nu_diag,1020) ' year_init                 = ', year_init
         write(nu_diag,1020) ' istep0                    = ', istep0
         write(nu_diag,1000) ' dt                        = ', dt
         write(nu_diag,1020) ' npt                       = ', npt
         write(nu_diag,1020) ' diagfreq                  = ', diagfreq
         write(nu_diag,1010) ' print_global              = ',
     &                         print_global
         write(nu_diag,1010) ' print_points              = ',
     &                         print_points
         write(nu_diag,1030) ' histfreq                  = ',
     &                         trim(histfreq)
         write(nu_diag,1010) ' hist_avg                  = ', hist_avg
         if (hist_avg) then
            write (nu_diag,*) 'History data will be averaged over 1 ',
     &                         histfreq
         else
            write (nu_diag,*) 'History data will be snapshots'
         endif
         write(nu_diag,*)    ' history_dir               = ',
     &                         trim(history_dir)
         write(nu_diag,*)    ' history_file              = ',
     &                         trim(history_file)
         write(nu_diag,1030) ' dumpfreq                  = ',
     &                         trim(dumpfreq)
         write(nu_diag,1020) ' dumpfreq_n                = ', dumpfreq_n
         write(nu_diag,1010) ' restart                   = ', restart
         write(nu_diag,*)    ' restart_dir               = ',
     &                         trim(restart_dir)
         write(nu_diag,*)    ' restart_file              = ',
     &                         trim(restart_file)
         write(nu_diag,*)    ' pointer_file              = ',
     &                         trim(pointer_file)
         write(nu_diag,1030) ' ice_ic                    = ', ice_ic
         write(nu_diag,*)    ' grid_type                 = ',
     &                         trim(grid_type)
         if (trim(grid_type) /= 'rectangular' .or.
     &       trim(grid_type) /= 'column') then
            write(nu_diag,*) ' grid_file                 = ',
     &                         trim(grid_file)
            write(nu_diag,*) ' kmt_file                  = ',
     &                         trim(kmt_file)
         endif

         write(nu_diag,1020) ' kitd                      = ', kitd
         write(nu_diag,1020) ' kcatbound                 = ',
     &                         kcatbound
         write(nu_diag,1020) ' kdyn                      = ', kdyn
         write(nu_diag,1020) ' ndyn_dt                   = ', ndyn_dt
         write(nu_diag,1020) ' ndte                      = ', ndte
         write(nu_diag,1010) ' evp_damping               = ',
     &                         evp_damping
         write(nu_diag,*)    ' yield_curve               = ',
     &                         trim(yield_curve)
         write(nu_diag,1020) ' kstrength                 = ', kstrength
         write(nu_diag,1020) ' krdg_partic               = ',
     &                         krdg_partic
         write(nu_diag,1020) ' krdg_redist               = ',
     &                         krdg_redist
         write(nu_diag,1030) ' advection                 = ',
     &                         trim(advection)
         write(nu_diag,1030) ' shortwave                 = ',
     &                         trim(shortwave)
         write(nu_diag,1000) ' albicev                   = ', albicev
         write(nu_diag,1000) ' albicei                   = ', albicei
         write(nu_diag,1000) ' albsnowv                  = ', albsnowv
         write(nu_diag,1000) ' albsnowi                  = ', albsnowi

         write(nu_diag,1020) ' fyear_init                = ',
     &                         fyear_init
         write(nu_diag,1020) ' ycycle                    = ', ycycle
         write(nu_diag,*)    ' atm_data_type             = ',
     &                         trim(atm_data_type)
         if (trim(atm_data_type) /= 'default') then
            write(nu_diag,*) ' atm_data_dir              = ',
     &                         trim(atm_data_dir)
            write(nu_diag,*) ' precip_units              = ',
     &                         trim(precip_units)
         endif 

         write(nu_diag,1010) ' oceanmixed_ice            = ',
     &                         oceanmixed_ice
         write(nu_diag,*)    ' sss_data_type             = ',
     &                         trim(sss_data_type)
         write(nu_diag,*)    ' sst_data_type             = ',
     &                         trim(sst_data_type)
         if (trim(sss_data_type) /= 'default' .or.
     &       trim(sst_data_type) /= 'default') then
            write(nu_diag,*) ' ocn_data_dir              = ',
     &                         trim(ocn_data_dir)
         endif 
         if (trim(sss_data_type) == 'ncar' .or.
     &       trim(sst_data_type) == 'ncar') then
            write(nu_diag,*) ' oceanmixed_file           = ',
     &                         trim(oceanmixed_file)
         endif

#ifdef coupled
         if( oceanmixed_ice ) then
            write (nu_diag,*) 'WARNING WARNING WARNING WARNING '
            write (nu_diag,*) '*Coupled and oceanmixed flags are  *'
            write (nu_diag,*) '*BOTH ON.  Ocean data received from*'
            write (nu_diag,*) '*coupler will be altered by mixed  *'
            write (nu_diag,*) '*layer routine!                    *'
            write (nu_diag,*) ' '
         endif
#endif

#ifdef CCSM
         write(nu_diag,1010) ' kcolumn                   = ', kcolumn
#endif

         write (nu_diag,*) ' '
         write (nu_diag,'(a30,2f8.2)') 'Diagnostic point 1: lat, lon =',
     &                      latpnt(1), lonpnt(1)
         write (nu_diag,'(a30,2f8.2)') 'Diagnostic point 2: lat, lon =',
     &                      latpnt(2), lonpnt(2)

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1030    format (a30,   a8)    ! character

         write (nu_diag,*) ' '
         if (grid_type  /=  'displaced_pole' .and.
     &       grid_type  /=  'column'         .and.
     &       grid_type  /=  'rectangular'    .and.
     &       grid_type  /=  'panarctic' ) then
            call abort_ice('ice_init: unknown grid_type')
         endif

      endif                     ! my_task = master_task

      end subroutine input_data

!=======================================================================
!BOP
!
! !IROUTINE: init_state - initialize ice state variables
!
! !INTERFACE:
!
      subroutine init_state
!
! !DESCRIPTION:
!
! Initialize state for the itd model
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL
!
! !USES:
!
      use ice_blocks
      use ice_domain
      use ice_flux, only: sst, Tair
      use ice_grid
      use ice_state
      use ice_itd
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  it              ! tracer index
     &,  iblk            ! block index

      !-----------------------------------------------------------------
      ! Check number of layers in ice and snow.
      !-----------------------------------------------------------------

      if (nilyr < 2) then
         write (nu_diag,*) 'nilyr =', nilyr
         write (nu_diag,*) 'Must have at least two ice layers'
         call abort_ice('Not enough ice layers')
      endif

      if (nslyr < 1) then
         write (nu_diag,*) 'nslyr =', nslyr
         write (nu_diag,*) 'Must have at least one snow layer'
         call abort_ice('Not enough snow layers')
      endif

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Set tracer types
!lipscomb - Do this later based on tracer indices, e.g. it_age
      !-----------------------------------------------------------------

      trcr_depend(1) = 0   ! ice/snow surface temperature
!!      trcr_depend(2) = 1   ! volume-weighted ice age


      !-----------------------------------------------------------------
      ! Set state variables
      !-----------------------------------------------------------------

         call set_state_var (nx_block,            ny_block,
     &                       tmask(:,:,    iblk), ULAT (:,:,    iblk),
     &                       Tair (:,:,    iblk),
     &                       sst  (:,:,    iblk), trcr_depend, 
     &                       aicen(:,:,  :,iblk), trcrn(:,:,:,:,iblk),
     &                       vicen(:,:,  :,iblk), vsnon(:,:,  :,iblk),
     &                       eicen(:,:,  :,iblk), esnon(:,:,  :,iblk))

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

         aice(:,:,iblk) = c0
         vice(:,:,iblk) = c0
         vsno(:,:,iblk) = c0
         eice(:,:,iblk) = c0
         esno(:,:,iblk) = c0
         do it = 1, ntrcr
            trcr(:,:,it,iblk) = c0
         enddo

         call aggregate (nx_block, ny_block,
     &                   aicen(:,:,:,iblk),
     &                   trcrn(:,:,:,:,iblk),
     &                   vicen(:,:,:,iblk),
     &                   vsnon(:,:,:,iblk),
     &                   eicen(:,:,:,iblk),
     &                   esnon(:,:,:,iblk),
     &                   aice (:,:,  iblk),
     &                   trcr (:,:,:,iblk),
     &                   vice (:,:,  iblk),
     &                   vsno (:,:,  iblk),
     &                   eice (:,:,  iblk),
     &                   esno (:,:,  iblk),
     &                   aice0(:,:,  iblk),
     &                   tmask(:,:,  iblk),
     &                   trcr_depend)

         aice_init(:,:,iblk) = aice(:,:,iblk)

      enddo                     ! iblk

      !-----------------------------------------------------------------
      ! ghost cell updates
      !-----------------------------------------------------------------

      call bound_state (aicen, trcrn,
     &                  vicen, vsnon, 
     &                  eicen, esnon)

      end subroutine init_state

!=======================================================================
!BOP
!
! !IROUTINE: set_state_var - initialize single-category state variables
!
! !INTERFACE:
!
      subroutine set_state_var (nx_block, ny_block,
     &                          tmask,    ULAT,
     &                          Tair,
     &                          sst,      trcr_depend,
     &                          aicen,    trcrn,
     &                          vicen,    vsnon,
     &                          eicen,    esnon)
!
! !DESCRIPTION:
!
! Initialize state in each ice thickness category
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz
!          William H. Lipscomb, LANL
!
! !USES:
!
      use ice_therm_vertical, only: Tmlt
      use ice_itd, only: ilyr1, slyr1, hin_max
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions

      logical (kind=log_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   tmask      ! true for ice/ocean cells

      real (kind=dbl_kind), dimension (nx_block,ny_block), 
     &   intent(in) ::
     &   ULAT       ! latitude of velocity pts (radians)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   Tair        ! air temperature  (K)
     &,  sst         ! sea surface temperature (C) 

      integer (kind=int_kind), dimension (ntrcr), intent(inout) ::
     &   trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(out) ::
     &   aicen     ! concentration of ice
     &,  vicen     ! volume per unit area of ice          (m)
     &,  vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),
     &   intent(out) ::
     &   trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),
     &   intent(out) ::
     &   eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),
     &   intent(out) ::
     &   esnon     ! energy of melting for each ice layer  (J/m^2)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! ice layer index
     &,  n               ! thickness category index
     &,  it              ! tracer index
     &,  icells          ! number of cells initialized with ice

      integer (kind=int_kind), dimension(nx_block*ny_block) ::
     &   indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind) ::
     &   slope, Ti
     &,  ainit(ncat)

      real (kind=dbl_kind), parameter ::
     &   hsno_init = 0.20_dbl_kind    ! initial snow thickness (m)
     &,  edge_init_nh =  70._dbl_kind ! initial ice edge, N.Hem. (deg) 
     &,  edge_init_sh = -60._dbl_kind ! initial ice edge, S.Hem. (deg)


      ! Initialize state variables.
      ! If restarting, these values are overwritten.

      aicen(:,:,:) = c0
      vicen(:,:,:) = c0
      vsnon(:,:,:) = c0
      eicen(:,:,:) = c0
      esnon(:,:,:) = c0
      trcrn(:,:,1,:) = Tocnfrz  ! surface temperature
      if (ntrcr >= 2) then
         do it = 2, ntrcr
            trcrn(:,:,it,:) = c0
         enddo
      endif

      if (trim(ice_ic) == 'default') then

      !-----------------------------------------------------------------
      ! Place ice where ocean surface is cold.
      ! Note: If SST is not read from a file, then the ocean is assumed
      !       to be at its freezing point everywhere, and ice will
      !       extend to the prescribed edges.
      !-----------------------------------------------------------------

         ! initial category areas in cells with ice

         if (ncat == 1) then
            ainit(1) = 0.8_dbl_kind
         else
            do n = 1, ncat
               ainit(n) = c1 / real(ncat,kind=dbl_kind)
            enddo
         endif

         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j)) then
               ! place ice in high latitudes where ocean sfc is cold
               if ( (sst (i,j) <= Tocnfrz+p2) .and.
     &              (ULAT(i,j) < edge_init_sh/rad_to_deg .or. 
     &               ULAT(i,j) > edge_init_nh/rad_to_deg) ) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif            ! cold surface
            endif               ! tmask
         enddo                  ! i
         enddo                  ! j


         ! ice volume, snow volume, surface temperature, other tracers

         do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               aicen(i,j,n) = ainit(n)
               if (n < ncat) then
                  vicen(i,j,n) = p5*(hin_max(n-1) + hin_max(n))
     &                         * aicen(i,j,n) ! m
               else             ! n=ncat
                  vicen(i,j,n) = (hin_max(n-1) + c1)*aicen(i,j,n) ! m
               endif
               vsnon(i,j,n) = min(hsno_init*aicen(i,j,n),
     &                                   p2*vicen(i,j,n))

               ! surface temperature
               trcrn(i,j,1,n) = min(Tsmelt, Tair(i,j) - Tffresh) ! deg C

               !lipscomb - volume-weighted test tracer
!               trcrn(i,j,2,n) = c1

            enddo               ! ij

            ! ice energy

            do k = 1, nilyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  ! assume linear temp profile and compute enthalpy
                  slope = Tocnfrz - trcrn(i,j,1,n)
                  Ti = trcrn(i,j,1,n) + slope*(real(k,kind=dbl_kind)-p5)
     &                                        /real(nilyr,kind=dbl_kind)

                  eicen(i,j,ilyr1(n)+k-1) =
     &                 -(rhoi * (cp_ice*(Tmlt(k)-Ti)
     &                 + Lfresh*(c1-Tmlt(k)/Ti) - cp_ocn*Tmlt(k)))
     &                 * vicen(i,j,n)/real(nilyr,kind=dbl_kind)

               enddo            ! ij
            enddo               ! nilyr

            ! snow energy

            do k = 1, nslyr
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  Ti = min(c0, trcrn(i,j,1,n))
                  esnon(i,j,slyr1(n)+k-1) = -rhos*(Lfresh - cp_ice*Ti)
     &                                      *vsnon(i,j,n)
     &                                      /real(nslyr,kind=dbl_kind)
               enddo            ! ij
            enddo               ! nslyr
            
         enddo                  ! ncat

      endif                     ! ice_ic

      end subroutine set_state_var

!=======================================================================

      end module ice_init

!=======================================================================
