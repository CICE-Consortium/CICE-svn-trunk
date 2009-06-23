!=======================================================================
!
!BOP
!
! !MODULE: ice_history - ice model history files
!
! Output files: netCDF or binary data, Fortran unformatted dumps
!
! The following variables are currently hard-wired as snapshots 
!   (instantaneous rather than time-averages):
!   divu, shear, sig1, sig2, trsig, mlt_onset, frz_onset, hisnap, aisnap
!
! Options for histfreq: '1','h','d','m','y','x', where x means that
!   output stream will not be used (recommended for efficiency).  
! histfreq_n can be any nonnegative integer, where 0 means that the 
!   corresponding histfreq frequency will not be used.
! The flags (f_<field>) can be set to '1','h','d','m','y' or 'x', where
!   n means the field will not be written.  To output the same field at
!   more than one frequency, for instance monthy and daily, set 
!   f_<field> = 'md'.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors Tony Craig and Bruce Briegleb, NCAR
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2004 WHL: Block structure added 
! 2006 ECH: Accepted some CCSM code into mainstream CICE
!           Added ice_present, aicen, vicen; removed aice1...10, vice1...1.
!           Added histfreq_n and histfreq='h' options, removed histfreq='w'
!           Converted to free source form (F90)
!           Added option for binary output instead of netCDF
! 2009 D Bailey and ECH: Generalized for multiple frequency output
!
! !INTERFACE:
!
      module ice_history
!
! !USES:
!
      use ice_kinds_mod
      use ice_broadcast
      use ice_communicate, only: my_task, master_task
      use ice_blocks
      use ice_read_write
      use ice_fileunits
!
!EOP
!
      implicit none
      save
      
      logical (kind=log_kind) :: &
         hist_avg  ! if true, write averaged data instead of snapshots

      character (len=char_len) :: &
         history_format, & ! file format ('bin'=binary or 'nc'=netcdf)
         history_file  , & ! output file for history
         incond_file       ! output file for snapshot initial conditions

      character (len=char_len_long) :: &
         history_dir   , & ! directory name for history file
         incond_dir        ! directory for snapshot initial conditions

      character (len=char_len_long) :: &
         pointer_file      ! input pointer file for restarts

      !---------------------------------------------------------------
      ! Instructions for adding a field: (search for 'example')
      !     Here:
      ! (1) Add to frequency flags (f_<field>)
      ! (2) Add to namelist (here and also in ice_in)
      ! (3) Add to index list
      !     In init_hist:
      ! (4) Add define_hist_field call with vname, vdesc, vunit,
      !     and vcomment, vcellmeas, and conversion factor if necessary.
      ! (5) Add flag to broadcast list
      ! (6) Add accum_hist_field call with appropriate variable
      !---------------------------------------------------------------

      type, public :: ice_hist_field
          character (len=16) :: vname     ! variable name
          character (len=16) :: vunit     ! variable units
          character (len=16) :: vcoord    ! variable coordinates
          character (len=16) :: vcellmeas ! variable cell measures
          character (len=55) :: vdesc     ! variable description
          character (len=55) :: vcomment  ! variable description
          real (kind=dbl_kind) :: cona    ! multiplicative conversion factor
          real (kind=dbl_kind) :: conb    ! additive conversion factor
          character (len=1) :: vhistfreq  ! frequency of history output
          integer (kind=int_kind) :: vhistfreq_n ! number of vhistfreq intervals
      end type

      integer (kind=int_kind), parameter :: &
         max_avail_hist_fields = 600      ! Max number of history fields

      integer (kind=int_kind) :: &
         num_avail_hist_fields   = 0      ! Current number of defined fields

      type (ice_hist_field), dimension(max_avail_hist_fields) :: &
         avail_hist_fields

      character (len=16) :: vname_in     ! variable name
      character (len=55) :: vdesc_in     ! variable description
      character (len=55) :: vcomment_in  ! variable description

      integer (kind=int_kind), parameter :: &
         nvar = 11              , & ! number of grid fields that can be written
                                    !   excluding grid vertices
         ncat_hist = ncat           ! number of ice categories written <= ncat

      real (kind=real_kind) :: time_beg(max_nstrm), &
                               time_end(max_nstrm) ! bounds for averaging

      real (kind=dbl_kind), allocatable :: &
         aa (:,:,:,:)      ! field accumulations and averages
         
      real (kind=dbl_kind) :: &
         avgct(max_nstrm)   ! average sample counter

      logical (kind=log_kind) :: &
         igrd(nvar)        ! true if grid field is written to output file

      character (len=16), parameter :: &
         tstr  = 'TLON TLAT time', & ! vcoord for T cell quantities
         ustr  = 'ULON ULAT time', & ! vcoord for U cell quantities
         tcstr = 'area: tarea'   , & ! vcellmeas for T cell quantities
         ucstr = 'area: uarea'       ! vcellmeas for U cell quantities

      !---------------------------------------------------------------
      ! flags: write to output file if true or histfreq value
      !---------------------------------------------------------------

      logical (kind=log_kind) :: &
           f_tmask     = .true., &
           f_tarea     = .true., f_uarea      = .true., &
           f_dxt       = .true., f_dyt        = .true., &
           f_dxu       = .true., f_dyu        = .true., &
           f_HTN       = .true., f_HTE        = .true., &
           f_ANGLE     = .true., f_ANGLET     = .true., &
           f_bounds    = .true.

      character (len=max_nstrm) :: &
!          f_example   = 'md', &
           f_hi        = 'm', f_hs         = 'm', &
           f_Tsfc      = 'm', f_aice       = 'm', &
           f_uvel      = 'm', f_vvel       = 'm', &
           f_fswdn     = 'm', f_flwdn      = 'm', &
           f_snow      = 'm', f_snow_ai    = 'm', &
           f_rain      = 'm', f_rain_ai    = 'm', &
           f_sst       = 'm', f_sss        = 'm', &
           f_uocn      = 'm', f_vocn       = 'm', &
           f_frzmlt    = 'm', &
           f_fswfac    = 'm', &
           f_fswabs    = 'm', f_fswabs_ai  = 'm', &
           f_albsni    = 'm', &
           f_alvdr     = 'm', f_alidr      = 'm', &
           f_albice    = 'm', f_albsno     = 'm', &
           f_albpnd    = 'm', f_coszen     = 'm', &
           f_flat      = 'm', f_flat_ai    = 'm', &
           f_fsens     = 'm', f_fsens_ai   = 'm', &
           f_flwup     = 'm', f_flwup_ai   = 'm', &
           f_evap      = 'm', f_evap_ai    = 'm', &
           f_Tair      = 'm', &
           f_Tref      = 'm', f_Qref       = 'm', &
           f_congel    = 'm', f_frazil     = 'm', &
           f_snoice    = 'm', f_meltt      = 'm', &
           f_meltb     = 'm', f_meltl      = 'm', &
           f_fresh     = 'm', f_fresh_ai   = 'm', &
           f_fsalt     = 'm', f_fsalt_ai   = 'm', &
           f_fhocn     = 'm', f_fhocn_ai   = 'm', &
           f_fswthru   = 'm', f_fswthru_ai = 'm', &
           f_strairx   = 'm', f_strairy    = 'm', &
           f_strtltx   = 'm', f_strtlty    = 'm', &
           f_strcorx   = 'm', f_strcory    = 'm', &
           f_strocnx   = 'm', f_strocny    = 'm', &
           f_strintx   = 'm', f_strinty    = 'm', &
           f_strength  = 'm', f_opening    = 'm', &
           f_divu      = 'm', f_shear      = 'm', &
           f_sig1      = 'm', f_sig2       = 'm', &
           f_dvidtt    = 'm', f_dvidtd     = 'm', &
           f_daidtt    = 'm', f_daidtd     = 'm', &
           f_mlt_onset = 'm', f_frz_onset  = 'm', &
           f_dardg1dt  = 'm', f_dardg2dt   = 'm', &
           f_dvirdgdt  = 'm', f_iage       = 'x',&
           f_hisnap    = 'm', f_aisnap     = 'm', &
           f_aicen     = 'm', f_vicen      = 'm', &
           f_apondn    = 'x',&
           f_trsig     = 'm', f_icepresent = 'm', &
           f_fsurf_ai  = 'm', f_fcondtop_ai= 'm', &
           f_fmeltt_ai = 'm',                        &
           f_fsurfn_ai = 'm',f_fcondtopn_ai= 'm', &
           f_fmelttn_ai= 'm', f_flatn_ai   = 'm'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_nml /     &
           f_tmask    , &
           f_tarea    , f_uarea    , &
           f_dxt      , f_dyt      , &
           f_dxu      , f_dyu      , &
           f_HTN      , f_HTE      , &
           f_ANGLE    , f_ANGLET   , &
           f_bounds   , &
!          f_example  , &
           f_hi,        f_hs       , &
           f_Tsfc,      f_aice     , &
           f_uvel,      f_vvel     , &
           f_fswdn,     f_flwdn    , &
           f_snow,      f_snow_ai  , &     
           f_rain,      f_rain_ai  , &
           f_sst,       f_sss      , &
           f_uocn,      f_vocn     , &
           f_frzmlt                , &
           f_fswfac                , &
           f_fswabs,    f_fswabs_ai, &
           f_albsni                , &
           f_alvdr,     f_alidr    , &
           f_albice,    f_albsno   , &
           f_albpnd,    f_coszen   , &
           f_flat,      f_flat_ai  , &
           f_fsens,     f_fsens_ai , &
           f_flwup,     f_flwup_ai , &
           f_evap,      f_evap_ai  , &
           f_Tair                  , &
           f_Tref,      f_Qref     , &
           f_congel,    f_frazil   , &
           f_snoice,    f_meltt    , &
           f_meltb,     f_meltl    , &
           f_fresh,     f_fresh_ai , &  
           f_fsalt,     f_fsalt_ai , &
           f_fhocn,     f_fhocn_ai , &
           f_fswthru,   f_fswthru_ai,&
           f_strairx,   f_strairy  , &
           f_strtltx,   f_strtlty  , &
           f_strcorx,   f_strcory  , &
           f_strocnx,   f_strocny  , &
           f_strintx,   f_strinty  , &
           f_strength,  f_opening  , &
           f_divu,      f_shear    , &
           f_sig1,      f_sig2     , &
           f_dvidtt,    f_dvidtd   , &
           f_daidtt,    f_daidtd   , &
           f_mlt_onset, f_frz_onset, &
           f_dardg1dt,  f_dardg2dt , &
           f_dvirdgdt              , &
           f_hisnap,    f_aisnap   , &
           f_aicen,     f_vicen    , &
           f_iage,      &
           f_apondn   , &
           f_trsig,     f_icepresent,&
           f_fsurf_ai,  f_fcondtop_ai,&
           f_fmeltt_ai,              &
           f_fsurfn_ai,f_fcondtopn_ai,&
           f_fmelttn_ai,f_flatn_ai 

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), parameter :: &
           n_tmask      = 1,  &
           n_tarea      = 2,  &
           n_uarea      = 3,  &
           n_dxt        = 4,  &
           n_dyt        = 5,  &
           n_dxu        = 6,  & 
           n_dyu        = 7,  &
           n_HTN        = 8,  &
           n_HTE        = 9,  &
           n_ANGLE      = 10, &
           n_ANGLET     = 11, &

           n_lont_bnds  = 1, &
           n_latt_bnds  = 2, &
           n_lonu_bnds  = 3, &
           n_latu_bnds  = 4

      integer (kind=int_kind), dimension(max_nstrm) :: &
!          n_example    , &
           n_hi         , n_hs         , &
           n_Tsfc       , n_aice       , &
           n_uvel       , n_vvel       , &
           n_fswdn      , n_flwdn      , &
           n_snow       , n_snow_ai    , &
           n_rain       , n_rain_ai    , &
           n_sst        , n_sss        , &
           n_uocn       , n_vocn       , &
           n_frzmlt     , n_fswfac     , &
           n_fswabs     , n_fswabs_ai  , &
           n_albsni     , &
           n_alvdr      , n_alidr      , &
           n_albice     , n_albsno     , &
           n_albpnd     , n_coszen     , &
           n_flat       , n_flat_ai    , &
           n_fsens      , n_fsens_ai   , &
           n_flwup      , n_flwup_ai   , &
           n_evap       , n_evap_ai    , &
           n_Tair       , &
           n_Tref       , n_Qref       , &
           n_congel     , n_frazil     , &
           n_snoice     , n_meltt      , &
           n_meltb      , n_meltl      , &
           n_fresh      , n_fresh_ai   , &
           n_fsalt      , n_fsalt_ai   , &
           n_fhocn      , n_fhocn_ai   , &
           n_fswthru    , n_fswthru_ai , &
           n_strairx    , n_strairy    , &
           n_strtltx    , n_strtlty    , &
           n_strcorx    , n_strcory    , &
           n_strocnx    , n_strocny    , &
           n_strintx    , n_strinty    , &
           n_strength   , n_opening    , &
           n_divu       , n_shear      , &
           n_sig1       , n_sig2       , &
           n_dvidtt     , n_dvidtd     , &
           n_daidtt     , n_daidtd     , &
           n_mlt_onset  , n_frz_onset  , &
           n_dardg1dt   , n_dardg2dt   , &
           n_dvirdgdt   , &
           n_hisnap     , n_aisnap     , &
           n_trsig      , n_icepresent , &
           n_iage       , n_fsurf_ai   , &
           n_fcondtop_ai, n_fmeltt_ai

      ! Category dependent variables
      integer(kind=int_kind), dimension(ncat_hist,max_nstrm) :: &
           n_aicen       , &
           n_vicen       , &
           n_apondn      , &
           n_fsurfn_ai   , &
           n_fcondtopn_ai, &
           n_fmelttn_ai  , &
           n_flatn_ai    

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: init_hist - initialize history files
!
! !INTERFACE:
!
      subroutine init_hist (dt)
!
! !DESCRIPTION:
!
! Initialize history files
!
! !REVISION HISTORY:
!
! authors Tony Craig, NCAR
!         Elizabeth C. Hunke, LANL
!         C.M. Bitz, UW
!         Bruce P. Briegleb, NCAR
!         William H. Lipscomb, LANL
!
! !USES:
!
      use ice_constants
      use ice_calendar, only: yday, days_per_year, histfreq, &
           histfreq_n, nstreams
      use ice_flux, only: mlt_onset, frz_onset, albcnt
      use ice_restart, only: restart
      use ice_age, only: tr_iage
      use ice_meltpond, only: tr_pond
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: n, k, ns, ns1, ns2, lenf
      integer (kind=int_kind) :: hfreqn
      integer (kind=int_kind), dimension(max_nstrm) :: &
         ntmp
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag

      character (len=3) :: nchar
      character (len=40) :: stmp

      !-----------------------------------------------------------------
      ! read namelist
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
            read(nu_nml, nml=icefields_nml,iostat=nml_error)
            if (nml_error > 0) read(nu_nml,*)  ! for Nagware compiler
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice('ice: error reading icefields_nml')
      endif

      ! histfreq options ('1','h','d','m','y')
      nstreams = 0
      do ns = 1, max_nstrm
         if (histfreq(ns) == '1' .or. histfreq(ns) == 'h' .or. &
             histfreq(ns) == 'd' .or. histfreq(ns) == 'm' .or. &
             histfreq(ns) == 'y') then
                nstreams = nstreams + 1
         else if (histfreq(ns) /= 'x') then
             call abort_ice('ice: histfreq contains illegal element')
         endif
      enddo
      if (nstreams == 0) write (nu_diag,*) 'WARNING: No history output'
      do ns1 = 1, nstreams
         do ns2 = 1, nstreams
            if (histfreq(ns1) == histfreq(ns2) .and. ns1/=ns2 &
               .and. my_task == master_task) then
               call abort_ice('ice: histfreq elements must be unique')
            endif
         enddo
      enddo

      if (.not. tr_iage) f_iage = 'x'
      if (.not. tr_pond) f_apondn = 'x'

      ! these must be output at the same frequency because of 
      ! cos(zenith angle) averaging
      if (f_albsno(1:1) /= 'x') f_albsno = f_albice
      if (f_albpnd(1:1) /= 'x') f_albpnd = f_albice
      if (f_coszen(1:1) /= 'x') f_coszen = f_albice

      ! to prevent array-out-of-bounds when aggregating
      if (f_fmeltt_ai(1:1) /= 'x') f_fmelttn_ai = f_fmeltt_ai

#ifndef ncdf
      f_bounds = .false.
#endif

      call broadcast_scalar (f_tmask, master_task)
      call broadcast_scalar (f_tarea, master_task)
      call broadcast_scalar (f_uarea, master_task)
      call broadcast_scalar (f_dxt, master_task)
      call broadcast_scalar (f_dyt, master_task)
      call broadcast_scalar (f_dxu, master_task)
      call broadcast_scalar (f_dyu, master_task)
      call broadcast_scalar (f_HTN, master_task)
      call broadcast_scalar (f_HTE, master_task)
      call broadcast_scalar (f_ANGLE, master_task)
      call broadcast_scalar (f_ANGLET, master_task)
      call broadcast_scalar (f_bounds, master_task)

!     call broadcast_scalar (f_example, master_task)
      call broadcast_scalar (f_hi, master_task)
      call broadcast_scalar (f_hs, master_task)
      call broadcast_scalar (f_Tsfc, master_task)
      call broadcast_scalar (f_aice, master_task)
      call broadcast_scalar (f_uvel, master_task)
      call broadcast_scalar (f_vvel, master_task)
      call broadcast_scalar (f_fswdn, master_task)
      call broadcast_scalar (f_flwdn, master_task)
      call broadcast_scalar (f_snow, master_task)
      call broadcast_scalar (f_snow_ai, master_task)
      call broadcast_scalar (f_rain, master_task)
      call broadcast_scalar (f_rain_ai, master_task)
      call broadcast_scalar (f_sst, master_task)
      call broadcast_scalar (f_sss, master_task)
      call broadcast_scalar (f_uocn, master_task)
      call broadcast_scalar (f_vocn, master_task)
      call broadcast_scalar (f_frzmlt, master_task)
      call broadcast_scalar (f_fswfac, master_task)
      call broadcast_scalar (f_fswabs, master_task)
      call broadcast_scalar (f_fswabs_ai, master_task)
      call broadcast_scalar (f_albsni, master_task)
      call broadcast_scalar (f_alvdr, master_task)
      call broadcast_scalar (f_alidr, master_task)
      call broadcast_scalar (f_albice, master_task)
      call broadcast_scalar (f_albsno, master_task)
      call broadcast_scalar (f_albpnd, master_task)
      call broadcast_scalar (f_coszen, master_task)
      call broadcast_scalar (f_flat, master_task)
      call broadcast_scalar (f_flat_ai, master_task)
      call broadcast_scalar (f_fsens, master_task)
      call broadcast_scalar (f_fsens_ai, master_task)
      call broadcast_scalar (f_flwup, master_task)
      call broadcast_scalar (f_flwup_ai, master_task)
      call broadcast_scalar (f_evap, master_task)
      call broadcast_scalar (f_evap_ai, master_task)
      call broadcast_scalar (f_Tair, master_task)
      call broadcast_scalar (f_Tref, master_task)
      call broadcast_scalar (f_Qref, master_task)
      call broadcast_scalar (f_congel, master_task)
      call broadcast_scalar (f_frazil, master_task)
      call broadcast_scalar (f_snoice, master_task)
      call broadcast_scalar (f_meltt, master_task)
      call broadcast_scalar (f_meltb, master_task)
      call broadcast_scalar (f_meltl, master_task)
      call broadcast_scalar (f_fresh, master_task)
      call broadcast_scalar (f_fresh_ai, master_task)
      call broadcast_scalar (f_fsalt, master_task)
      call broadcast_scalar (f_fsalt_ai, master_task)
      call broadcast_scalar (f_fhocn, master_task)
      call broadcast_scalar (f_fhocn_ai, master_task)
      call broadcast_scalar (f_fswthru, master_task)
      call broadcast_scalar (f_fswthru_ai, master_task)
      call broadcast_scalar (f_strairx, master_task)
      call broadcast_scalar (f_strairy, master_task)
      call broadcast_scalar (f_strtltx, master_task)
      call broadcast_scalar (f_strtlty, master_task)
      call broadcast_scalar (f_strcorx, master_task)
      call broadcast_scalar (f_strcory, master_task)
      call broadcast_scalar (f_strocnx, master_task)
      call broadcast_scalar (f_strocny, master_task)
      call broadcast_scalar (f_strintx, master_task)
      call broadcast_scalar (f_strinty, master_task)
      call broadcast_scalar (f_strength, master_task)
      call broadcast_scalar (f_opening, master_task)
      call broadcast_scalar (f_divu, master_task)
      call broadcast_scalar (f_shear, master_task)
      call broadcast_scalar (f_sig1, master_task)
      call broadcast_scalar (f_sig2, master_task)
      call broadcast_scalar (f_dvidtt, master_task)
      call broadcast_scalar (f_dvidtd, master_task)
      call broadcast_scalar (f_daidtt, master_task)
      call broadcast_scalar (f_daidtd, master_task)
      call broadcast_scalar (f_mlt_onset, master_task)
      call broadcast_scalar (f_frz_onset, master_task)
      call broadcast_scalar (f_dardg1dt, master_task)
      call broadcast_scalar (f_dardg2dt, master_task)
      call broadcast_scalar (f_dvirdgdt, master_task)
      call broadcast_scalar (f_aisnap, master_task)
      call broadcast_scalar (f_hisnap, master_task)
      call broadcast_scalar (f_aicen, master_task)
      call broadcast_scalar (f_vicen, master_task)
      call broadcast_scalar (f_trsig, master_task)
      call broadcast_scalar (f_icepresent, master_task)
      call broadcast_scalar (f_fsurf_ai, master_task)
      call broadcast_scalar (f_fcondtop_ai, master_task)
      call broadcast_scalar (f_fmeltt_ai, master_task)
      call broadcast_scalar (f_fsurfn_ai, master_task)
      call broadcast_scalar (f_fcondtopn_ai, master_task)
      call broadcast_scalar (f_fmelttn_ai, master_task)
      call broadcast_scalar (f_flatn_ai, master_task)

      call broadcast_scalar (f_iage, master_task)
      call broadcast_scalar (f_apondn, master_task)

      do ns1 = 1, nstreams

!!!!! begin example
!      if (f_example(1:1) /= 'x') &
!         call define_hist_field(n_example,"example","m",tstr, tcstr, & 
!            "example: mean ice thickness",                           &
!            "ice volume per unit grid cell area", c1, c0,            &
!            ns1, f_example)
!!!!! end example

      if (f_hi(1:1) /= 'x') &
         call define_hist_field(n_hi,"hi","m",tstr, tcstr,        & 
            "grid cell mean ice thickness",                       &
            "ice volume per unit grid cell area", c1, c0,         &
            ns1, f_hi)

      if (f_hs(1:1) /= 'x') &
         call define_hist_field(n_hs,"hs","m",tstr, tcstr,        &
             "grid cell mean snow thickness",                     &
             "snow volume per unit grid cell area", c1, c0,       &
             ns1, f_hs)
      
      if (f_Tsfc(1:1) /= 'x') &
         call define_hist_field(n_Tsfc,"Tsfc","C",tstr, tcstr,    &
             "snow/ice surface temperature",                      &
             "averaged with Tf if no ice is present", c1, c0,     &
             ns1, f_Tsfc)
      
      if (f_aice(1:1) /= 'x') &
         call define_hist_field(n_aice,"aice","1",tstr, tcstr,    &
             "ice area  (aggregate)",                             &
             "none", c1, c0,                                      &
             ns1, f_aice)
      
      if (f_uvel(1:1) /= 'x') &
         call define_hist_field(n_uvel,"uvel","m/s",ustr, ucstr,  &
             "ice velocity (x)",                                  &
             "positive is x direction on U grid", c1, c0,         &
             ns1, f_uvel)
      
      if (f_vvel(1:1) /= 'x') &
         call define_hist_field(n_vvel,"vvel","m/s",ustr, ucstr,  &
             "ice velocity (y)",                                  &
             "positive is y direction on U grid", c1, c0,         &
             ns1, f_vvel)
      
      if (f_fswdn(1:1) /= 'x') &
         call define_hist_field(n_fswdn,"fswdn","W/m^2",tstr, tcstr, &
             "down solar flux",                                      &
             "positive downward", c1, c0,                            &
             ns1, f_fswdn)
      
      if (f_flwdn(1:1) /= 'x') &
         call define_hist_field(n_flwdn,"flwdn","W/m^2",tstr, tcstr, &
             "down longwave flux",                                   &
             "positive downward", c1, c0,                            &
             ns1, f_flwdn)
      
      if (f_snow(1:1) /= 'x') &
         call define_hist_field(n_snow,"snow","cm/day",tstr, tcstr, &
             "snowfall rate (cpl)",                                 &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_snow)
      
      if (f_snow_ai(1:1) /= 'x') &
         call define_hist_field(n_snow_ai,"snow_ai","cm/day",tstr, tcstr, &
             "snowfall rate",                                             &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_snow_ai)
      
      if (f_rain(1:1) /= 'x') &
         call define_hist_field(n_rain,"rain","cm/day",tstr, tcstr, &
             "rainfall rate (cpl)",                                 &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_rain)
      
      if (f_rain_ai(1:1) /= 'x') &
         call define_hist_field(n_rain_ai,"rain_ai","cm/day",tstr, tcstr, &
             "rainfall rate",                                             &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_rain_ai)
      
      if (f_sst(1:1) /= 'x') &
         call define_hist_field(n_sst,"sst","C",tstr, tcstr, &
             "sea surface temperature",                      &
             "none", c1, c0,                                 &
             ns1, f_sst)
      
      if (f_sss(1:1) /= 'x') &
         call define_hist_field(n_sss,"sss","ppt",tstr, tcstr, &
             "sea surface salinity",                           &
             "none", c1, c0,                                   &
             ns1, f_sss)
      
      if (f_uocn(1:1) /= 'x') &
         call define_hist_field(n_uocn,"uocn","m/s",ustr, ucstr, &
             "ocean current (x)",                                &
             "positive is x direction on U grid", c1, c0,        &
             ns1, f_uocn)
      
      if (f_vocn(1:1) /= 'x') &
         call define_hist_field(n_vocn,"vocn","m/s",ustr, ucstr, &
             "ocean current (y)",                                &
             "positive is y direction on U grid", c1, c0,        &
             ns1, f_vocn)
      
      if (f_frzmlt(1:1) /= 'x') &
         call define_hist_field(n_frzmlt,"frzmlt","W/m^2",tstr, tcstr, &
             "freeze/melt potential",                                  &
             "if >0, new ice forms; if <0, ice melts", c1, c0,         &
             ns1, f_frzmlt)
      
      if (f_fswfac(1:1) /= 'x') &
         call define_hist_field(n_fswfac,"fswfac","1",tstr, tcstr, &
             "shortwave scaling factor",                           &
             "ratio of netsw new:old", c1, c0,                     &
             ns1, f_fswfac)
      
      if (f_fswabs(1:1) /= 'x') &
         call define_hist_field(n_fswabs,"fswabs","W/m^2",tstr, tcstr, &
             "snow/ice/ocn absorbed solar flux (cpl)",                 &
             "positive downward", c1, c0,                              &
             ns1, f_fswabs)
      
      if (f_fswabs_ai(1:1) /= 'x') &
         call define_hist_field(n_fswabs_ai,"fswabs_ai","W/m^2",tstr, tcstr, &
             "snow/ice/ocn absorbed solar flux",                             &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_fswabs_ai)
      
      if (f_albsni(1:1) /= 'x') &
         call define_hist_field(n_albsni,"albsni","%",tstr, tcstr, &
             "snow/ice broad band albedo",                         &
             "scaled (divided) by aice", c100, c0,                 &
             ns1, f_albsni)
      
      if (f_alvdr(1:1) /= 'x') &
         call define_hist_field(n_alvdr,"alvdr","%",tstr, tcstr, &
             "visible direct albedo",                            &
             "scaled (divided) by aice", c100, c0,               &
             ns1, f_alvdr)
      
      if (f_alidr(1:1) /= 'x') &
         call define_hist_field(n_alidr,"alidr","%",tstr, tcstr, &
             "near IR direct albedo",                            &
             "scaled (divided) by aice", c100, c0,               &
             ns1, f_alidr)

      if (f_albice(1:1) /= 'x') &
         call define_hist_field(n_albice,"albice","%",tstr, tcstr, &
             "bare ice albedo",                                    &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albice)
      
      if (f_albsno(1:1) /= 'x') &
         call define_hist_field(n_albsno,"albsno","%",tstr, tcstr, &
             "snow albedo",                                        &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albsno)
      
      if (f_albpnd(1:1) /= 'x') &
         call define_hist_field(n_albpnd,"albpnd","%",tstr, tcstr, &
             "melt pond albedo",                                   &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albpnd)
      
      if (f_coszen(1:1) /= 'x') &
         call define_hist_field(n_coszen,"coszen","radian",tstr, tcstr, &
             "cosine of the zenith angle",                              &
             "negative below horizon", c1, c0,                          &
             ns1, f_coszen)

      if (f_flat(1:1) /= 'x') &
         call define_hist_field(n_flat,"flat","W/m^2",tstr, tcstr, &
             "latent heat flux (cpl)",                             &
             "positive downward", c1, c0,                          &
             ns1, f_flat)
      
      if (f_flat_ai(1:1) /= 'x') &
         call define_hist_field(n_flat_ai,"flat_ai","W/m^2",tstr, tcstr, &
             "latent heat flux",                                         &
             "weighted by ice area", c1, c0,                             &
             ns1, f_flat_ai)
      
      if (f_fsens(1:1) /= 'x') &
         call define_hist_field(n_fsens,"fsens","W/m^2",tstr, tcstr, &
             "sensible heat flux (cpl)",                             &
             "positive downward", c1, c0,                            &
             ns1, f_fsens)
      
      if (f_fsens_ai(1:1) /= 'x') &
         call define_hist_field(n_fsens_ai,"fsens_ai","W/m^2",tstr, tcstr, &
             "sensible heat flux",                                         &
             "weighted by ice area", c1, c0,                               &
             ns1, f_fsens_ai)
      
      if (f_flwup(1:1) /= 'x') &
         call define_hist_field(n_flwup,"flwup","W/m^2",tstr, tcstr, &
             "upward longwave flux (cpl)",                           &
             "positive downward", c1, c0,                            &
             ns1, f_flwup)
      
      if (f_flwup_ai(1:1) /= 'x') &
         call define_hist_field(n_flwup_ai,"flwup_ai","W/m^2",tstr, tcstr, &
             "upward longwave flux",                                       &
             "weighted by ice area", c1, c0,                               &
             ns1, f_flwup_ai)
      
      if (f_evap(1:1) /= 'x') &
         call define_hist_field(n_evap,"evap","cm/day",tstr, tcstr, &
             "evaporative water flux (cpl)",                        &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_evap)
      
      if (f_evap_ai(1:1) /= 'x') &
         call define_hist_field(n_evap_ai,"evap_ai","cm/day",tstr, tcstr, &
             "evaporative water flux",                                    &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_evap_ai)
      
      if (f_Tair(1:1) /= 'x') &
         call define_hist_field(n_Tair,"Tair","C",tstr, tcstr, &
             "air temperature",                                &
             "none", c1, -Tffresh,                             &
             ns1, f_Tair)
      
      if (f_Tref(1:1) /= 'x') &
         call define_hist_field(n_Tref,"Tref","C",tstr, tcstr, &
             "2m reference temperature",                       &
             "none", c1, -Tffresh,                             &
             ns1, f_Tref)
      
      if (f_Qref(1:1) /= 'x') &
         call define_hist_field(n_Qref,"Qref","g/kg",tstr, tcstr, &
             "2m reference specific humidity",                    &
             "none", kg_to_g, c0,                                 &
             ns1, f_Qref)
      
      if (f_congel(1:1) /= 'x') &
         call define_hist_field(n_congel,"congel","cm/day",tstr, tcstr, &
             "congelation ice growth",                                  &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_congel)
      
      if (f_frazil(1:1) /= 'x') &
         call define_hist_field(n_frazil,"frazil","cm/day",tstr, tcstr, &
             "frazil ice growth",                                       &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_frazil)
      
      if (f_snoice(1:1) /= 'x') &
         call define_hist_field(n_snoice,"snoice","cm/day",tstr, tcstr, &
             "snow-ice formation",                                      &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_snoice)
      
      if (f_meltt(1:1) /= 'x') &
         call define_hist_field(n_meltt,"meltt","cm/day",tstr, tcstr, &
             "top ice melt",                                          &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltt)
      
      if (f_meltb(1:1) /= 'x') &
         call define_hist_field(n_meltb,"meltb","cm/day",tstr, tcstr, &
             "basal ice melt",                                        &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltb)
      
      if (f_meltl(1:1) /= 'x') &
         call define_hist_field(n_meltl,"meltl","cm/day",tstr, tcstr, &
             "lateral ice melt",                                      &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltl)
      
      if (f_fresh(1:1) /= 'x') &
         call define_hist_field(n_fresh,"fresh","cm/day",tstr, tcstr,   &
             "freshwtr flx ice to ocn (cpl)",                           &
             "if positive, ocean gains fresh water",                    &
             mps_to_cmpdy/rhofresh, c0,                                 &
             ns1, f_fresh)
      
      if (f_fresh_ai(1:1) /= 'x') &
         call define_hist_field(n_fresh_ai,"fresh_ai","cm/day",tstr, tcstr, &
             "freshwtr flx ice to ocn",                                     &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,             &
             ns1, f_fresh_ai)
      
      if (f_fsalt(1:1) /= 'x') &
         call define_hist_field(n_fsalt,"fsalt","kg/m^2/s",tstr, tcstr, &
             "salt flux ice to ocn (cpl)",                              &
             "if positive, ocean gains salt", c1, c0,                   &
             ns1, f_fsalt)
      
      if (f_fsalt_ai(1:1) /= 'x') &
         call define_hist_field(n_fsalt_ai,"fsalt_ai","kg/m^2/s",tstr, tcstr, &
             "salt flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fsalt_ai)
      
      if (f_fhocn(1:1) /= 'x') &
         call define_hist_field(n_fhocn,"fhocn","W/m^2",tstr, tcstr, &
             "heat flux ice to ocn (cpl)",                           &
             "if positive, ocean gains heat", c1, c0,                &
             ns1, f_fhocn)
      
      if (f_fhocn_ai(1:1) /= 'x') &
         call define_hist_field(n_fhocn_ai,"fhocn_ai","W/m^2",tstr, tcstr, &
             "heat flux ice to ocean",                                     &
             "weighted by ice area", c1, c0,                               &
             ns1, f_fhocn_ai)
      
      if (f_fswthru(1:1) /= 'x') &
         call define_hist_field(n_fswthru,"fswthru","W/m^2",tstr, tcstr, &
             "SW thru ice to ocean (cpl)",                               &
             "if positive, ocean gains heat", c1, c0,                    &
             ns1, f_fswthru)
      
      if (f_fswthru_ai(1:1) /= 'x') &
         call define_hist_field(n_fswthru_ai,"fswthru_ai","W/m^2",tstr, tcstr,&
             "SW flux thru ice to ocean",                                     &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fswthru_ai)
      
      if (f_strairx(1:1) /= 'x') &
         call define_hist_field(n_strairx,"strairx","N/m^2",ustr, ucstr, &
             "atm/ice stress (x)",                                       &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strairx)
      
      if (f_strairy(1:1) /= 'x') &
         call define_hist_field(n_strairy,"strairy","N/m^2",ustr, ucstr, &
             "atm/ice stress (y)",                                       &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strairy)
      
      if (f_strtltx(1:1) /= 'x') &
         call define_hist_field(n_strtltx,"strtltx","N/m^2",ustr, ucstr, &
             "sea sfc tilt stress (x)",                                  &
             "none", c1, c0,                                             &
             ns1, f_strtltx)
      
      if (f_strtlty(1:1) /= 'x') &
         call define_hist_field(n_strtlty,"strtlty","N/m^2",ustr, ucstr, &
             "sea sfc tilt stress (y)",                                  &
             "none", c1, c0,                                             &
             ns1, f_strtlty)
      
      if (f_strcorx(1:1) /= 'x') &
         call define_hist_field(n_strcorx,"strcorx","N/m^2",ustr, ucstr, &
             "coriolis stress (x)",                                      &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strcorx)
      
      if (f_strcory(1:1) /= 'x') &
         call define_hist_field(n_strcory,"strcory","N/m^2",ustr, ucstr, &
             "coriolis stress (y)",                                      &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strcory)
      
      if (f_strocnx(1:1) /= 'x') &
         call define_hist_field(n_strocnx,"strocnx","N/m^2",ustr, ucstr, &
             "ocean/ice stress (x)",                                     &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strocnx)
      
      if (f_strocny(1:1) /= 'x') &
         call define_hist_field(n_strocny,"strocny","N/m^2",ustr, ucstr, &
             "ocean/ice stress (y)",                                     &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strocny)
      
      if (f_strintx(1:1) /= 'x') &
         call define_hist_field(n_strintx,"strintx","N/m^2",ustr, ucstr, &
             "internal ice stress (x)",                                  &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strintx)
      
      if (f_strinty(1:1) /= 'x') &
         call define_hist_field(n_strinty,"strinty","N/m^2",ustr, ucstr, &
             "internal ice stress (y)",                                  &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strinty)
      
      if (f_strength(1:1) /= 'x') &
         call define_hist_field(n_strength,"strength","N/m",tstr, tcstr, &
             "compressive ice strength",                                 &
             "none", c1, c0,                                             &
             ns1, f_strength)
      
      if (f_opening(1:1) /= 'x') &
         call define_hist_field(n_opening,"opening","%/day",tstr, tcstr, &
             "lead area opening rate",                                   &
             "none", secday*c100, c0,                                    &
             ns1, f_opening)
      
      if (f_divu(1:1) /= 'x') &
         call define_hist_field(n_divu,"divu","%/day",tstr, tcstr, &
             "strain rate (divergence)",                           &
             "none", secday*c100, c0,                              &
             ns1, f_divu)
      
      if (f_shear(1:1) /= 'x') &
         call define_hist_field(n_shear,"shear","%/day",tstr, tcstr, &
             "strain rate (shear)",                                  &
             "none", secday*c100, c0,                                &
             ns1, f_shear)
      
      if (f_sig1(1:1) /= 'x') &
         call define_hist_field(n_sig1,"sig1","1",ustr, ucstr, &
             "norm. principal stress 1",                       &
             "sig1 is instantaneous", c1, c0,                  &
             ns1, f_sig1)
      
      if (f_sig2(1:1) /= 'x') &
         call define_hist_field(n_sig2,"sig2","1",ustr, ucstr, &
             "norm. principal stress 2",                       &
             "sig2 is instantaneous", c1, c0,                  &
             ns1, f_sig2)
      
      if (f_dvidtt(1:1) /= 'x') &
         call define_hist_field(n_dvidtt,"dvidtt","cm/day",tstr, tcstr, &
             "volume tendency thermo",                                  &
             "none", mps_to_cmpdy, c0,                                  &
             ns1, f_dvidtt)
      
      if (f_dvidtd(1:1) /= 'x') &
         call define_hist_field(n_dvidtd,"dvidtd","cm/day",tstr, tcstr, &
             "volume tendency dynamics",                                &
             "none", mps_to_cmpdy, c0,                                  &
             ns1, f_dvidtd)
      
      if (f_daidtt(1:1) /= 'x') &
         call define_hist_field(n_daidtt,"daidtt","%/day",tstr, tcstr, &
             "area tendency thermo",                                   &
             "none", secday*c100, c0,                                  &
             ns1, f_daidtt)
      
      if (f_daidtd(1:1) /= 'x') &
         call define_hist_field(n_daidtd,"daidtd","%/day",tstr, tcstr, &
             "area tendency dynamics",                                 &
             "none", secday*c100, c0,                                  &
             ns1, f_daidtd)
      
      if (f_mlt_onset(1:1) /= 'x') &
         call define_hist_field(n_mlt_onset,"mlt_onset","day of year", &
             tstr, tcstr,"melt onset date",                            &
             "midyear restart gives erroneous dates", c1, c0,          &
             ns1, f_mlt_onset)
      
      if (f_frz_onset(1:1) /= 'x') &
         call define_hist_field(n_frz_onset,"frz_onset","day of year", &
             tstr, tcstr,"freeze onset date",                          &
             "midyear restart gives erroneous dates", c1, c0,          &
             ns1, f_frz_onset)

      if (f_dardg1dt(1:1) /= 'x') &
         call define_hist_field(n_dardg1dt,"dardg1dt","%/day",tstr, tcstr, &
             "ice area ridging rate",                                      &
             "none", secday*c100, c0,                                      &
             ns1, f_dardg1dt)
      
      if (f_dardg2dt(1:1) /= 'x') &
         call define_hist_field(n_dardg2dt,"dardg2dt","%/day",tstr, tcstr, &
             "ridge area formation rate",                                  &
             "none", secday*c100, c0,                                      &
             ns1, f_dardg2dt)
      
      if (f_dvirdgdt(1:1) /= 'x') &
         call define_hist_field(n_dvirdgdt,"dvirdgdt","cm/day",tstr, tcstr, &
             "ice volume ridging rate",                                     &
             "none", mps_to_cmpdy, c0,                                      &
             ns1, f_dvirdgdt)
      
      if (f_hisnap(1:1) /= 'x') &
         call define_hist_field(n_hisnap,"hisnap","m",tstr, tcstr, &
             "ice volume snapshot",                                &
             "none", c1, c0,                              &
             ns1, f_hisnap)
      
      if (f_aisnap(1:1) /= 'x') &
         call define_hist_field(n_aisnap,"aisnap","1",tstr, tcstr, &
             "ice area snapshot",                                  &
             "none", c1, c0,                              &
             ns1, f_aisnap)
      
      if (f_trsig(1:1) /= 'x') &
         call define_hist_field(n_trsig,"trsig","N/m^2",tstr, tcstr, &
             "internal stress tensor trace",                         &
             "ice strength approximation", c1, c0,                   &
             ns1, f_trsig)
      
      if (f_icepresent(1:1) /= 'x') &
         call define_hist_field(n_icepresent,"ice_present","1",tstr, tcstr, &
             "fraction of time-avg interval that ice is present",           &
             "ice extent flag", c1, c0,                                     &
             ns1, f_icepresent)
      
      if (f_fsurf_ai(1:1) /= 'x') &
         call define_hist_field(n_fsurf_ai,"fsurf_ai","W/m^2",tstr, tcstr, &
             "net surface heat flux",                                      &
             "positive downward, excludes conductive flux, weighted by ice area", &
             c1, c0, &
             ns1, f_fsurf_ai)
      
      if (f_fcondtop_ai(1:1) /= 'x') &
         call define_hist_field(n_fcondtop_ai,"fcondtop_ai","W/m^2", &
             tstr, tcstr,"top surface conductive heat flux",         &
             "positive downward, weighted by ice area", c1, c0,      &
             ns1, f_fcondtop_ai)
      
      if (f_fmeltt_ai(1:1) /= 'x') &
         call define_hist_field(n_fmeltt_ai,"fmeltt_ai","W/m^2",tstr, tcstr, &
             "net surface heat flux causing melt",                           &
             "always >= 0, weighted by ice area", c1, c0,                    &
             ns1, f_fmeltt_ai)
      
      ! Category variables

      if (f_aicen(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'aicen', trim(nchar) ! aicen
            stmp = 'ice area, category '   ! aicen
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_aicen(n,:),vname_in,"1",tstr, tcstr, &
                vdesc_in, "Ice range:", c1, c0,                       &
                ns1, f_aicen)
         enddo
      endif
      
      if (f_vicen(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'vicen', trim(nchar) ! vicen
            stmp = 'ice volume, category '   ! vicen
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_vicen(n,:),vname_in,"m",tstr, tcstr, &
                vdesc_in, "none", c1, c0,                               &
                ns1, f_vicen)
         enddo
      endif
      
      if (f_fsurfn_ai(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fsurfn_ai', trim(nchar) ! fsurfn_ai
            stmp = 'net surface heat flux, category '   ! fsurfn_ai
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_fsurfn_ai(n,:),vname_in,"W/m^2",    &
                tstr, tcstr, vdesc_in, "weighted by ice area", c1, c0, &
                ns1, f_fsurfn_ai)
         enddo
      endif
      
      if (f_fcondtopn_ai(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fcondtopn_ai', trim(nchar) ! fcondtopn_ai
            stmp = 'top sfc conductive heat flux, cat '   ! fcondtopn_ai
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_fcondtopn_ai(n,:),vname_in,"W/m^2", &
                tstr, tcstr, vdesc_in, "weighted by ice area", c1, c0, &
                ns1, f_fcondtopn_ai)
         enddo
      endif
      
      if (f_fmelttn_ai(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fmelttn_ai', trim(nchar) ! fmelttn_ai
            stmp = 'net sfc heat flux causing melt, cat '   ! fmelttn_ai
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_fmelttn_ai(n,:),vname_in,"W/m^2",   &
                tstr, tcstr, vdesc_in, "weighted by ice area", c1, c0, &
                ns1, f_fmelttn_ai)
         enddo
      endif
      
      if (f_flatn_ai(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'flatn_ai', trim(nchar) ! flatn_ai
            stmp = 'latent heat flux, category '   ! flatn_ai
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_flatn_ai(n,:),vname_in,"W/m^2",tstr, tcstr, &
                vdesc_in, "weighted by ice area", c1, c0,                      &
                ns1, f_flatn_ai)
         enddo
      endif      

      ! Tracers !!!!!!!!!!!

      ! Ice Age
      if (f_iage(1:1) /= 'x') &
         call define_hist_field(n_iage,"iage","years",tstr, tcstr, &
             "sea ice age",                                        &
             "none", c1/(secday*days_per_year), c0,                &
             ns1, f_iage)
       
      ! Melt ponds
      if (f_apondn(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'apond', trim(nchar) ! apondn
            stmp = 'melt pond concentration, category '   ! apondn
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_apondn(n,:),vname_in,"1",tstr, tcstr, &
                vdesc_in, "none", c1, c0,                                &
                ns1, f_apondn)
         enddo
      endif

      enddo ! ns1

      allocate(aa(nx_block,ny_block,num_avail_hist_fields,max_blocks))

      !-----------------------------------------------------------------
      ! fill igrd array with namelist values
      !-----------------------------------------------------------------

      igrd=.true.

      igrd(n_tmask     ) = f_tmask
      igrd(n_tarea     ) = f_tarea
      igrd(n_uarea     ) = f_uarea
      igrd(n_dxt       ) = f_dxt
      igrd(n_dyt       ) = f_dyt
      igrd(n_dxu       ) = f_dxu
      igrd(n_dyu       ) = f_dyu
      igrd(n_HTN       ) = f_HTN
      igrd(n_HTE       ) = f_HTE
      igrd(n_ANGLE     ) = f_ANGLE
      igrd(n_ANGLET    ) = f_ANGLET

      ntmp(:) = 0
      if (my_task == master_task) then
        write(nu_diag,*) ' '
        write(nu_diag,*) 'The following variables will be ', &
                         'written to the history tape: '
        write(nu_diag,101) 'description','units','variable','frequency','x'
        do n=1,num_avail_hist_fields
           if (avail_hist_fields(n)%vhistfreq_n /= 0) &
           write(nu_diag,100) avail_hist_fields(n)%vdesc, &
              avail_hist_fields(n)%vunit, avail_hist_fields(n)%vname, &
              avail_hist_fields(n)%vhistfreq,avail_hist_fields(n)%vhistfreq_n
           do ns = 1, nstreams
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                 ntmp(ns)=ntmp(ns)+1
           enddo
        enddo ! num_avail_hist_fields
        write(nu_diag,*) ' '
      endif
  100 format (1x,a40,2x,a16,2x,a12,1x,a1,2x,i6)
  101 format (2x,a19,10x,a16,9x,a12,2x,a,3x,a1)

      call broadcast_array(ntmp, master_task)
      do ns = 1, nstreams
         if (ntmp(ns)==0) histfreq_n(ns) = 0
      enddo

      !-----------------------------------------------------------------
      ! initialize the history arrays
      !-----------------------------------------------------------------
      aa(:,:,:,:) = c0
      avgct(:) = c0
      albcnt(:,:,:,:) = c0

      if (restart .and. yday >= c2) then
! restarting midyear gives erroneous onset dates
         mlt_onset = 999._dbl_kind 
         frz_onset = 999._dbl_kind 
      else
         mlt_onset = c0
         frz_onset = c0
      endif

      end subroutine init_hist

!=======================================================================
!
!BOP
!
! !IROUTINE: ice_write_hist - write average ice quantities or snapshots
!
! !INTERFACE:
!
      subroutine ice_write_hist (dt)
!
! !DESCRIPTION:
!
! write average ice quantities or snapshots
!
! !REVISION HISTORY:
!
! author:   Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_blocks
      use ice_domain
      use ice_grid, only: tmask, lmask_n, lmask_s
      use ice_calendar, only: new_year, secday, yday, write_history, &
                              write_ic, time, histfreq, nstreams
      use ice_state
      use ice_constants
      use ice_flux
      use ice_dyn_evp
      use ice_work, only: worka, workb
      use ice_timers
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!EOP
!
      integer (kind=int_kind) :: &
           i,j,k,n,nct,ns   , &
           iblk             , & ! block index
           ilo,ihi,jlo,jhi  , & ! beginning and end of physical domain
           nstrm                ! nstreams (1 if writing initial condition)

      real (kind=dbl_kind) :: &
           ravgct           , & ! 1/avgct
           ravgctz          , & ! 1/avgct
           ai               , & ! aice_init
           ain                  ! aicen_init

      type (block) :: &
         this_block           ! block information for current block

      !---------------------------------------------------------------
      ! increment step counter
      !---------------------------------------------------------------

      do ns = 1,nstreams
         if (.not. hist_avg .or. histfreq(ns) == '1') then  ! write snapshots
           do k = 1,num_avail_hist_fields
              if (avail_hist_fields(k)%vhistfreq == histfreq(ns)) &
                  aa(:,:,k,:) = c0
           enddo
           avgct(ns) = c1
         else                      ! write averages over time histfreq
           avgct(ns) = avgct(ns) + c1
           if (avgct(ns) == c1) time_beg(ns) = (time-dt)/int(secday)
         endif
      enddo

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

      do iblk = 1, nblocks
         workb(:,:) = aice_init(:,:,iblk)

!        if (f_example(1:1) /= 'x') &
!            call accum_hist_field(n_example,iblk, vice(:,:,iblk))
         if (f_hi     (1:1) /= 'x') &
             call accum_hist_field(n_hi,     iblk, vice(:,:,iblk))
         if (f_hs     (1:1) /= 'x') &
             call accum_hist_field(n_hs,     iblk, vsno(:,:,iblk))
         if (f_Tsfc   (1:1) /= 'x') &
             call accum_hist_field(n_Tsfc,   iblk, trcr(:,:,nt_Tsfc,iblk))
         if (f_aice   (1:1) /= 'x') &
             call accum_hist_field(n_aice,   iblk, aice(:,:,iblk))
         if (f_uvel   (1:1) /= 'x') &
             call accum_hist_field(n_uvel,   iblk, uvel(:,:,iblk))
         if (f_vvel   (1:1) /= 'x') &
             call accum_hist_field(n_vvel,   iblk, vvel(:,:,iblk))

         if (f_fswdn  (1:1) /= 'x') &
             call accum_hist_field(n_fswdn,  iblk, fsw(:,:,iblk))
         if (f_flwdn  (1:1) /= 'x') &
             call accum_hist_field(n_flwdn,  iblk, flw(:,:,iblk))
         if (f_snow   (1:1) /= 'x') &
             call accum_hist_field(n_snow,   iblk, fsnow(:,:,iblk))
         if (f_snow_ai(1:1) /= 'x') &
             call accum_hist_field(n_snow_ai,iblk, fsnow(:,:,iblk)*workb(:,:))
         if (f_rain   (1:1) /= 'x') &
             call accum_hist_field(n_rain,   iblk, frain(:,:,iblk))
         if (f_rain_ai(1:1) /= 'x') &
             call accum_hist_field(n_rain_ai,iblk, frain(:,:,iblk)*workb(:,:))

         if (f_sst    (1:1) /= 'x') &
             call accum_hist_field(n_sst,    iblk, sst(:,:,iblk))
         if (f_sss    (1:1) /= 'x') &
             call accum_hist_field(n_sss,    iblk, sss(:,:,iblk))
         if (f_uocn   (1:1) /= 'x') &
             call accum_hist_field(n_uocn,   iblk, uocn(:,:,iblk))
         if (f_vocn   (1:1) /= 'x') &
             call accum_hist_field(n_vocn,   iblk, vocn(:,:,iblk))
         if (f_frzmlt (1:1) /= 'x') &
             call accum_hist_field(n_frzmlt, iblk, frzmlt(:,:,iblk))

         if (f_fswfac (1:1) /= 'x') &
             call accum_hist_field(n_fswfac, iblk, fswfac(:,:,iblk))
         if (f_fswabs (1:1) /= 'x') &
             call accum_hist_field(n_fswabs, iblk, fswabs(:,:,iblk))
         if (f_fswabs_ai(1:1)/= 'x') &
             call accum_hist_field(n_fswabs_ai, iblk, fswabs(:,:,iblk)*workb(:,:))

         if (f_albsni (1:1) /= 'x') &
             call accum_hist_field(n_albsni, iblk, &
                                              awtvdr*alvdr(:,:,iblk) &
                                            + awtidr*alidr(:,:,iblk) &
                                            + awtvdf*alvdf(:,:,iblk) &
                                            + awtidf*alidf(:,:,iblk))
         if (f_alvdr  (1:1) /= 'x') &
             call accum_hist_field(n_alvdr,  iblk, alvdr(:,:,iblk))
         if (f_alidr  (1:1) /= 'x') &
             call accum_hist_field(n_alidr,  iblk, alidr(:,:,iblk))

         if (f_albice (1:1) /= 'x') &
             call accum_hist_field(n_albice, iblk, albice(:,:,iblk))
         if (f_albsno (1:1) /= 'x') &
             call accum_hist_field(n_albsno, iblk, albsno(:,:,iblk))
         if (f_albpnd (1:1) /= 'x') &
             call accum_hist_field(n_albpnd, iblk, albpnd(:,:,iblk))
         if (f_coszen (1:1) /= 'x') &
             call accum_hist_field(n_coszen, iblk, coszen(:,:,iblk))

         if (f_flat   (1:1) /= 'x') &
             call accum_hist_field(n_flat,   iblk, flat(:,:,iblk))
         if (f_flat_ai(1:1) /= 'x') &
             call accum_hist_field(n_flat_ai,iblk, flat(:,:,iblk)*workb(:,:))
         if (f_fsens  (1:1) /= 'x') &
             call accum_hist_field(n_fsens,   iblk, fsens(:,:,iblk))
         if (f_fsens_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsens_ai,iblk, fsens(:,:,iblk)*workb(:,:))
         if (f_flwup  (1:1) /= 'x') &
             call accum_hist_field(n_flwup,   iblk, flwout(:,:,iblk))
         if (f_flwup_ai(1:1)/= 'x') &
             call accum_hist_field(n_flwup_ai,iblk, flwout(:,:,iblk)*workb(:,:))
         if (f_evap   (1:1) /= 'x') &
             call accum_hist_field(n_evap,   iblk, evap(:,:,iblk))
         if (f_evap_ai(1:1) /= 'x') &
             call accum_hist_field(n_evap_ai,iblk, evap(:,:,iblk)*workb(:,:))

         if (f_Tair   (1:1) /= 'x') &
             call accum_hist_field(n_Tair,   iblk, Tair(:,:,iblk))
         if (f_Tref   (1:1) /= 'x') &
             call accum_hist_field(n_Tref,   iblk, Tref(:,:,iblk))
         if (f_Qref   (1:1) /= 'x') &
             call accum_hist_field(n_Qref,   iblk, Qref(:,:,iblk))
         if (f_congel (1:1) /= 'x') &
             call accum_hist_field(n_congel, iblk, congel(:,:,iblk))
         if (f_frazil (1:1) /= 'x') &
             call accum_hist_field(n_frazil, iblk, frazil(:,:,iblk))
         if (f_snoice (1:1) /= 'x') &
             call accum_hist_field(n_snoice, iblk, snoice(:,:,iblk))
         if (f_meltt  (1:1) /= 'x') &
             call accum_hist_field(n_meltt,  iblk, meltt(:,:,iblk))
         if (f_meltb  (1:1) /= 'x') &
             call accum_hist_field(n_meltb,  iblk, meltb(:,:,iblk))
         if (f_meltl  (1:1) /= 'x') &
             call accum_hist_field(n_meltl,  iblk, meltl(:,:,iblk))

         if (f_fresh  (1:1) /= 'x') &
             call accum_hist_field(n_fresh,   iblk, fresh(:,:,iblk))
         if (f_fresh_ai(1:1)/= 'x') &
             call accum_hist_field(n_fresh_ai,iblk, fresh_gbm(:,:,iblk))
         if (f_fsalt  (1:1) /= 'x') &
             call accum_hist_field(n_fsalt,   iblk, fsalt(:,:,iblk))
         if (f_fsalt_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsalt_ai,iblk, fsalt_gbm(:,:,iblk))
         if (f_fhocn  (1:1) /= 'x') &
             call accum_hist_field(n_fhocn,   iblk, fhocn(:,:,iblk))
         if (f_fhocn_ai(1:1)/= 'x') &
             call accum_hist_field(n_fhocn_ai,iblk, fhocn_gbm(:,:,iblk))
         if (f_fswthru(1:1) /= 'x') &
             call accum_hist_field(n_fswthru, iblk, fswthru(:,:,iblk))
         if (f_fswthru_ai(1:1)/= 'x') &
             call accum_hist_field(n_fswthru_ai,iblk, fswthru_gbm(:,:,iblk))
               
         if (f_strairx(1:1) /= 'x') &
             call accum_hist_field(n_strairx, iblk, strairx(:,:,iblk))
         if (f_strairy(1:1) /= 'x') &
             call accum_hist_field(n_strairy, iblk, strairy(:,:,iblk))
         if (f_strtltx(1:1) /= 'x') &
             call accum_hist_field(n_strtltx, iblk, strtltx(:,:,iblk))
         if (f_strtlty(1:1) /= 'x') &
             call accum_hist_field(n_strtlty, iblk, strtlty(:,:,iblk))
         if (f_strcorx(1:1) /= 'x') &
             call accum_hist_field(n_strcorx, iblk, fm(:,:,iblk)*vvel(:,:,iblk))
         if (f_strcory(1:1) /= 'x') &
             call accum_hist_field(n_strcory, iblk,-fm(:,:,iblk)*uvel(:,:,iblk))
         if (f_strocnx(1:1) /= 'x') &
             call accum_hist_field(n_strocnx, iblk, strocnx(:,:,iblk))
         if (f_strocny(1:1) /= 'x') &
             call accum_hist_field(n_strocny, iblk, strocny(:,:,iblk))
         if (f_strintx(1:1) /= 'x') &
             call accum_hist_field(n_strintx, iblk, strintx(:,:,iblk))
         if (f_strinty(1:1) /= 'x') &
             call accum_hist_field(n_strinty, iblk, strinty(:,:,iblk))
         if (f_strength(1:1)/= 'x') &
             call accum_hist_field(n_strength,iblk, strength(:,:,iblk))

! The following fields (divu, shear, sig1, and sig2) will be smeared
!  if averaged over more than a few days.
! Snapshots may be more useful (see below).

!        if (f_divu   (1:1) /= 'x') &
!             call accum_hist_field(n_divu,    iblk, divu(:,:,iblk))
!        if (f_shear  (1:1) /= 'x') &
!             call accum_hist_field(n_shear,   iblk, shear(:,:,iblk))
!        if (f_sig1   (1:1) /= 'x') &
!             call accum_hist_field(n_sig1,    iblk, sig1(:,:,iblk))
!        if (f_sig2   (1:1) /= 'x') &
!             call accum_hist_field(n_sig2,    iblk, sig2(:,:,iblk))
!        if (f_trsig  (1:1) /= 'x') &
!             call accum_hist_field(n_trsig,   iblk, trsig(:,:,iblk))

         if (f_dvidtt (1:1) /= 'x') &
             call accum_hist_field(n_dvidtt,  iblk, dvidtt(:,:,iblk))
         if (f_dvidtd (1:1) /= 'x') &
             call accum_hist_field(n_dvidtd,  iblk, dvidtd(:,:,iblk))
         if (f_daidtt (1:1) /= 'x') &
             call accum_hist_field(n_daidtt,  iblk, daidtt(:,:,iblk))
         if (f_daidtd (1:1) /= 'x') &
             call accum_hist_field(n_daidtd,  iblk, daidtd(:,:,iblk))

         if (f_opening(1:1) /= 'x') &
             call accum_hist_field(n_opening, iblk, opening(:,:,iblk))
         if (f_dardg1dt(1:1)/= 'x') &
             call accum_hist_field(n_dardg1dt,iblk, dardg1dt(:,:,iblk))
         if (f_dardg2dt(1:1)/= 'x') &
             call accum_hist_field(n_dardg2dt,iblk, dardg2dt(:,:,iblk))
         if (f_dvirdgdt(1:1)/= 'x') &
             call accum_hist_field(n_dvirdgdt,iblk, dvirdgdt(:,:,iblk))

         if (f_fsurf_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsurf_ai,iblk, fsurf(:,:,iblk)*workb(:,:))
         if (f_fcondtop_ai(1:1)/= 'x') &
             call accum_hist_field(n_fcondtop_ai, iblk, &
                                                 fcondtop(:,:,iblk)*workb(:,:))

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         if (f_icepresent(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) worka(i,j) = c1
           enddo
           enddo
           call accum_hist_field(n_icepresent, iblk, worka(:,:))
         endif

         nct = min(ncat, ncat_hist)
         do n=1,nct
            workb(:,:) = aicen_init(:,:,n,iblk)
            if (f_aicen   (1:1) /= 'x') &
               call accum_hist_field(n_aicen(n,:), iblk, aicen(:,:,n,iblk))
            if (f_vicen   (1:1) /= 'x') &
               call accum_hist_field(n_vicen(n,:), iblk, vicen(:,:,n,iblk))
            if (f_apondn  (1:1) /= 'x') &
               call accum_hist_field(n_apondn(n,:),iblk, apondn(:,:,n,iblk))
            if (f_fsurfn_ai(1:1)/= 'x') &
               call accum_hist_field(n_fsurfn_ai(n,:),iblk, &
                                     fsurfn(:,:,n,iblk)*workb(:,:))
            if (f_fcondtopn_ai(1:1)/= 'x') &
               call accum_hist_field(n_fcondtopn_ai(n,:),iblk, &
                                     fcondtopn(:,:,n,iblk)*workb(:,:))
            if (f_flatn_ai(1:1) /= 'x') &
               call accum_hist_field(n_flatn_ai(n,:),iblk, &
                                     flatn(:,:,n,iblk)*workb(:,:))
            ! Calculate surface heat flux that causes melt as this
            ! is what is calculated by the atmos in HadGEM3 so
            ! needed for checking purposes
            if (f_fmelttn_ai(1:1) /= 'x') &
               call accum_hist_field(n_fmelttn_ai(n,:),iblk, &
                  max(fsurfn(:,:,n,iblk) - fcondtopn(:,:,n,iblk),c0)*workb(:,:))
        enddo                    ! n

        ! Calculate aggregate surface melt flux by summing category values
        if (f_fmeltt_ai(1:1) /= 'x') then
         do ns = 1, nstreams
           if (n_fmeltt_ai(ns) /= 0) then
              worka(:,:) = c0
              do j = jlo, jhi
              do i = ilo, ihi
               if (tmask(i,j,iblk)) then
                 do n=1,nct
                    worka(i,j)  = worka(i,j) + aa(i,j,n_fmelttn_ai(n,ns),iblk)
                 enddo            ! n
               endif              ! tmask
              enddo                ! i
              enddo                ! j
              aa(:,:,n_fmeltt_ai(ns),iblk) = worka(:,:)
           endif
         enddo
        endif

      enddo                     ! iblk

      !---------------------------------------------------------------
      ! Write output files at prescribed intervals
      !---------------------------------------------------------------

      nstrm = nstreams
      if (write_ic) nstrm = 1

      do ns = 1, nstrm
      if (write_history(ns) .or. write_ic) then

      !---------------------------------------------------------------
      ! Mask out land points and convert units 
      !---------------------------------------------------------------

        ravgct = c1/avgct(ns)
        do iblk = 1, nblocks
           this_block = get_block(blocks_ice(iblk),iblk)         
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do k = 1, num_avail_hist_fields
              if (avail_hist_fields(k)%vhistfreq == histfreq(ns)) then 

              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    aa(i,j,k,iblk) = spval
                 else                            ! convert units
                    aa(i,j,k,iblk) = avail_hist_fields(k)%cona*aa(i,j,k,iblk) &
                                   * ravgct + avail_hist_fields(k)%conb
                 endif
              enddo             ! i
              enddo             ! j

              ! back out albedo/zenith angle dependence
              if (avail_hist_fields(k)%vname == 'albice') then
              do j = jlo, jhi
              do i = ilo, ihi
                 if (tmask(i,j,iblk)) then 
                    ravgctz = c0
                    if (albcnt(i,j,iblk,ns) > puny) &
                        ravgctz = c1/albcnt(i,j,iblk,ns)
                    if (n_albice(ns) /= 0) aa(i,j,n_albice(ns),iblk) = &
                       aa(i,j,n_albice(ns),iblk)*avgct(ns)*ravgctz
                    if (n_albsno(ns) /= 0) aa(i,j,n_albsno(ns),iblk) = &
                       aa(i,j,n_albsno(ns),iblk)*avgct(ns)*ravgctz
                    if (n_albpnd(ns) /= 0) aa(i,j,n_albpnd(ns),iblk) = &
                       aa(i,j,n_albpnd(ns),iblk)*avgct(ns)*ravgctz
                 endif
              enddo             ! i
              enddo             ! j
              endif

              endif
           enddo                ! k

      !---------------------------------------------------------------
      ! snapshots
      !---------------------------------------------------------------

          ! compute sig1 and sig2
        
           call principal_stress (nx_block,  ny_block,  &
                                  stressp_1 (:,:,iblk), &
                                  stressm_1 (:,:,iblk), &
                                  stress12_1(:,:,iblk), &
                                  prs_sig   (:,:,iblk), &
                                  sig1      (:,:,iblk), &
                                  sig2      (:,:,iblk))
 
           do j = jlo, jhi
           do i = ilo, ihi
              if (.not. tmask(i,j,iblk)) then ! mask out land points
                 if (n_divu     (ns) /= 0) aa(i,j,n_divu(ns),     iblk) = spval
                 if (n_shear    (ns) /= 0) aa(i,j,n_shear(ns),    iblk) = spval
                 if (n_sig1     (ns) /= 0) aa(i,j,n_sig1(ns),     iblk) = spval
                 if (n_sig2     (ns) /= 0) aa(i,j,n_sig2(ns),     iblk) = spval
                 if (n_mlt_onset(ns) /= 0) aa(i,j,n_mlt_onset(ns),iblk) = spval
                 if (n_frz_onset(ns) /= 0) aa(i,j,n_frz_onset(ns),iblk) = spval
                 if (n_hisnap   (ns) /= 0) aa(i,j,n_hisnap(ns),   iblk) = spval
                 if (n_aisnap   (ns) /= 0) aa(i,j,n_aisnap(ns),   iblk) = spval
                 if (n_trsig    (ns) /= 0) aa(i,j,n_trsig(ns),    iblk) = spval
                 if (n_iage     (ns) /= 0) aa(i,j,n_iage(ns),     iblk) = spval
              else
                 if (n_divu     (ns) /= 0) aa(i,j,n_divu(ns),iblk)      = &
                       divu (i,j,iblk)*avail_hist_fields(n_divu(ns))%cona
                 if (n_shear    (ns) /= 0) aa(i,j,n_shear(ns),iblk)     = &
                       shear(i,j,iblk)*avail_hist_fields(n_shear(ns))%cona
                 if (n_sig1     (ns) /= 0) aa(i,j,n_sig1(ns),iblk)      = &
                       sig1 (i,j,iblk)*avail_hist_fields(n_sig1(ns))%cona
                 if (n_sig2     (ns) /= 0) aa(i,j,n_sig2(ns),iblk)      = &
                       sig2 (i,j,iblk)*avail_hist_fields(n_sig2(ns))%cona
                 if (n_mlt_onset(ns) /= 0) aa(i,j,n_mlt_onset(ns),iblk) = &
                       mlt_onset(i,j,iblk)
                 if (n_frz_onset(ns) /= 0) aa(i,j,n_frz_onset(ns),iblk) = &
                       frz_onset(i,j,iblk)
                 if (n_hisnap   (ns) /= 0) aa(i,j,n_hisnap(ns),iblk)    = &
                       vice(i,j,iblk)
                 if (n_aisnap   (ns) /= 0) aa(i,j,n_aisnap(ns),iblk)    = &
                       aice(i,j,iblk)
                 if (n_trsig    (ns) /= 0) aa(i,j,n_trsig(ns),iblk )    = &
                                       p25*(stressp_1(i,j,iblk) &
                                          + stressp_2(i,j,iblk) &
                                          + stressp_3(i,j,iblk) &
                                          + stressp_4(i,j,iblk))
                 if (n_iage     (ns) /= 0) aa(i,j,n_iage(ns),iblk)  = &
                       trcr(i,j,nt_iage,iblk)*avail_hist_fields(n_iage(ns))%cona
              endif
           enddo                ! i
           enddo                ! j

        enddo                   ! iblk

        time_end(ns) = time/int(secday)

      !---------------------------------------------------------------
      ! write file
      !---------------------------------------------------------------

        call ice_timer_start(timer_readwrite)  ! reading/writing

        if (history_format == 'nc') then
          call icecdf(ns)         ! netcdf output
        else
          call icebin(ns)         ! binary output
        endif

        call ice_timer_stop(timer_readwrite)  ! reading/writing

      !---------------------------------------------------------------
      ! reset to zero
      !------------------------------------------------------------
        if (write_ic) then
           aa(:,:,:,:) = c0
           avgct(:) = c0
           albcnt(:,:,:,:) = c0
           write_ic = .false.        ! write initial condition once at most
        else
           avgct(ns) = c0
           albcnt(:,:,:,ns) = c0
        endif
!        if (write_history(ns)) albcnt(:,:,:,ns) = c0

        do k=1,num_avail_hist_fields
           if (avail_hist_fields(k)%vhistfreq == histfreq(ns)) aa(:,:,k,:) = c0
        enddo

      endif  ! write_history or write_ic
      enddo  ! nstreams

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         if (new_year) then

            do j=jlo,jhi
            do i=ilo,ihi
               ! reset NH Jan 1
               if (lmask_n(i,j,iblk)) mlt_onset(i,j,iblk) = c0
               ! reset SH Jan 1 
               if (lmask_s(i,j,iblk)) frz_onset(i,j,iblk) = c0
            enddo
            enddo
         endif                  ! new_year

         if ((yday >= 181._dbl_kind) .and. &
             (yday <  181._dbl_kind+dt/secday)) then

            do j=jlo,jhi
            do i=ilo,ihi

               ! reset SH Jul 1
               if (lmask_s(i,j,iblk)) mlt_onset(i,j,iblk) = c0

               ! reset NH Jul 1
               if (lmask_n(i,j,iblk)) frz_onset(i,j,iblk) = c0
            enddo
            enddo

         endif                  ! yday
      enddo                     ! iblk

      end subroutine ice_write_hist

!=======================================================================
!
!BOP
!
! !IROUTINE: icecdf - write netCDF history file
!
! !INTERFACE:
!
      subroutine icecdf(ns)
!
! !DESCRIPTION:
!
! write netCDF history file
!
! !REVISION HISTORY:
!
! authors:   E.C.Hunke, LANL
!            Bruce P. Briegleb, NCAR
!
! !USES:
!
#ifdef ncdf

      use ice_gather_scatter
      use ice_domain_size
      use ice_constants
      use ice_grid
      use ice_calendar, only: time, sec, idate, idate0, nyr, month, &
                              mday, write_ic, histfreq, histfreq_n, &
                              year_init, new_year, new_month, new_day, &
                              dayyr, daymo, days_per_year
      use ice_work, only: work_g1, work_gr, work_gr3, work1
      use ice_restart, only: lenstr, runid
      use ice_domain, only: distrb_info
      use ice_itd, only: c_hi_range
      use ice_exit
      use netcdf
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: ns
!
!EOP
!
      integer (kind=int_kind) :: i,j,n, &
         ncid,status,imtid,jmtid,timid,varid, &
         length,nvertexid,ivertex
      integer (kind=int_kind), dimension(3) :: dimid
      integer (kind=int_kind), dimension(3) :: dimid_nverts
      real (kind=real_kind) :: ltime
      character (char_len) :: title
      character (char_len_long) :: ncfile(max_nstrm)

      integer (kind=int_kind) :: iyear, imonth, iday
      integer (kind=int_kind) :: icategory,ind,i_aice,boundid

      character (char_len) :: start_time,current_date,current_time
      character (len=16) :: c_aice
      character (len=8) :: cdate

      ! 4 coordinate variables: TLON, TLAT, ULON, ULAT
      INTEGER (kind=int_kind), PARAMETER :: ncoord = 4

      ! 4 vertices in each grid cell
      INTEGER (kind=int_kind), PARAMETER :: nverts = 4

      ! 4 variables describe T, U grid boundaries:
      ! lont_bounds, latt_bounds, lonu_bounds, latu_bounds
      INTEGER (kind=int_kind), PARAMETER :: nvar_verts = 4

      TYPE coord_attributes         ! netcdf coordinate attributes
        character (len=11)   :: short_name
        character (len=45)   :: long_name
        character (len=20)   :: units
      END TYPE coord_attributes

      TYPE req_attributes         ! req'd netcdf attributes
        type (coord_attributes) :: req
        character (len=20)   :: coordinates
      END TYPE req_attributes

      TYPE(req_attributes), dimension(nvar) :: var
      TYPE(coord_attributes), dimension(ncoord) :: coord_var
      TYPE(coord_attributes), dimension(nvar_verts) :: var_nverts
      CHARACTER (char_len), dimension(ncoord) :: coord_bounds

      if (my_task == master_task) then

        ltime=time/int(secday)

        call construct_filename(ncfile(ns),'nc',ns)

        ! add local directory path name to ncfile
        if (write_ic) then
          ncfile(ns) = trim(incond_dir)//ncfile(ns)
        else
          ncfile(ns) = trim(history_dir)//ncfile(ns)
        endif

        ! create file
        status = nf90_create(ncfile(ns), nf90_clobber, ncid)
        if (status /= nf90_noerr) call abort_ice( &
           'ice: Error creating history ncfile '//ncfile(ns))

      !-----------------------------------------------------------------
      ! define dimensions
      !-----------------------------------------------------------------

        if (hist_avg) then
          status = nf90_def_dim(ncid,'d2',2,boundid)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error defining dim d2')
        endif

        status = nf90_def_dim(ncid,'ni',nx_global,imtid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error defining dim ni')

        status = nf90_def_dim(ncid,'nj',ny_global,jmtid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error defining dim nj')

        status = nf90_def_dim(ncid,'time',NF90_UNLIMITED,timid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error defining dim time')

        status = nf90_def_dim(ncid,'nvertices',nverts,nvertexid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error defining dim nverts')

      !-----------------------------------------------------------------
      ! define coordinate variables
      !-----------------------------------------------------------------

        status = nf90_def_var(ncid,'time',nf90_float,timid,varid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error defining var time')

        status = nf90_put_att(ncid,varid,'long_name','model time')
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: time long_name')

        write(cdate,'(i8.8)') idate0
        write(title,'(a,a,a,a,a,a,a,a)') 'days since ', &
              cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
        status = nf90_put_att(ncid,varid,'units',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: time units')

        if (days_per_year == 360) then
           status = nf90_put_att(ncid,varid,'calendar','360_day')
           if (status /= nf90_noerr) call abort_ice( &
                         'ice Error: time calendar')
        else
           status = nf90_put_att(ncid,varid,'calendar','noleap')
           if (status /= nf90_noerr) call abort_ice( &
                         'ice Error: time calendar')
        endif

        if (hist_avg) then
          status = nf90_put_att(ncid,varid,'bounds','time_bounds')
          if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: time bounds')
        endif

      !-----------------------------------------------------------------
      ! Define attributes for time bounds if hist_avg is true
      !-----------------------------------------------------------------

        if (hist_avg) then
          dimid(1) = boundid
          dimid(2) = timid
          status = nf90_def_var(ncid,'time_bounds',nf90_float,dimid(1:2),varid)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error defining var time_bounds')
          status = nf90_put_att(ncid,varid,'long_name', &
                                'boundaries for time-averaging interval')
          if (status /= nf90_noerr) call abort_ice( &
                        'ice Error: time_bounds long_name')
          write(cdate,'(i8.8)') idate0
          write(title,'(a,a,a,a,a,a,a,a)') 'days since ', &
                cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
          status = nf90_put_att(ncid,varid,'units',title)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice Error: time_bounds units')
        endif

      !-----------------------------------------------------------------
      ! define information for required time-invariant variables
      !-----------------------------------------------------------------

      ind = 0
      ind = ind + 1
      coord_var(ind) = coord_attributes('TLON', &
                       'T grid center longitude', 'degrees_east')
      coord_bounds(ind) = 'lont_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('TLAT', &
                       'T grid center latitude',  'degrees_north')
      coord_bounds(ind) = 'latt_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('ULON', &
                       'U grid center longitude', 'degrees_east')
      coord_bounds(ind) = 'lonu_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('ULAT', &
                       'U grid center latitude',  'degrees_north')
      coord_bounds(ind) = 'latu_bounds'

      !-----------------------------------------------------------------
      ! define information for optional time-invariant variables
      !-----------------------------------------------------------------

      var(n_tarea)%req = coord_attributes('tarea', &
                  'area of T grid cells', 'm^2')
      var(n_tarea)%coordinates = 'TLON TLAT'
      var(n_uarea)%req = coord_attributes('uarea', &
                  'area of U grid cells', 'm^2')
      var(n_uarea)%coordinates = 'ULON ULAT'
      var(n_dxt)%req = coord_attributes('dxt', &
                  'T cell width through middle', 'm')
      var(n_dxt)%coordinates = 'TLON TLAT'
      var(n_dyt)%req = coord_attributes('dyt', &
                  'T cell height through middle', 'm')
      var(n_dyt)%coordinates = 'TLON TLAT'
      var(n_dxu)%req = coord_attributes('dxu', &
                  'U cell width through middle', 'm')
      var(n_dxu)%coordinates = 'ULON ULAT'
      var(n_dyu)%req = coord_attributes('dyu', &
                  'U cell height through middle', 'm')
      var(n_dyu)%coordinates = 'ULON ULAT'
      var(n_HTN)%req = coord_attributes('HTN', &
                  'T cell width on North side','m')
      var(n_HTN)%coordinates = 'TLON TLAT'
      var(n_HTE)%req = coord_attributes('HTE', &
                  'T cell width on East side', 'm')
      var(n_HTE)%coordinates = 'TLON TLAT'
      var(n_ANGLE)%req = coord_attributes('ANGLE', &
                  'angle grid makes with latitude line on U grid', &
                  'radians')
      var(n_ANGLE)%coordinates = 'ULON ULAT'
      var(n_ANGLET)%req = coord_attributes('ANGLET', &
                  'angle grid makes with latitude line on T grid', &
                  'radians')
      var(n_ANGLET)%coordinates = 'TLON TLAT'

      ! These fields are required for CF compliance
      ! dimensions (nx,ny,nverts)
      var_nverts(n_lont_bnds) = coord_attributes('lont_bounds', &
                  'longitude boundaries of T cells', 'degrees_east')
      var_nverts(n_latt_bnds) = coord_attributes('latt_bounds', &
                  'latitude boundaries of T cells', 'degrees_north')
      var_nverts(n_lonu_bnds) = coord_attributes('lonu_bounds', &
                  'longitude boundaries of U cells', 'degrees_east')
      var_nverts(n_latu_bnds) = coord_attributes('latu_bounds', &
                  'latitude boundaries of U cells', 'degrees_north')

      !-----------------------------------------------------------------
      ! define attributes for time-invariant variables
      !-----------------------------------------------------------------

        dimid(1) = imtid
        dimid(2) = jmtid
        dimid(3) = timid

        do i = 1, ncoord
          status = nf90_def_var(ncid, coord_var(i)%short_name, nf90_float, &
                                dimid(1:2), varid)
          if (status /= nf90_noerr) call abort_ice( &
               'Error defining short_name for '//coord_var(i)%short_name)
          status = nf90_put_att(ncid,varid,'long_name',coord_var(i)%long_name)
          if (status /= nf90_noerr) call abort_ice( &
               'Error defining long_name for '//coord_var(i)%short_name)
          status = nf90_put_att(ncid, varid, 'units', coord_var(i)%units)
          if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//coord_var(i)%short_name)
          if (coord_var(i)%short_name == 'ULAT') then
             status = nf90_put_att(ncid,varid,'comment', &
                  'Latitude of NE corner of T grid cell')
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining comment for '//coord_var(i)%short_name)
          endif
          if (f_bounds) then
              status = nf90_put_att(ncid, varid, 'bounds', coord_bounds(i))
              if (status /= nf90_noerr) call abort_ice( &
                  'Error defining bounds for '//coord_var(i)%short_name)
          endif          
        enddo

        ! Attributes for tmask defined separately, since it has no units
        if (igrd(n_tmask)) then
           status = nf90_def_var(ncid, 'tmask', nf90_float, dimid(1:2), varid)
           if (status /= nf90_noerr) call abort_ice( &
                         'ice: Error defining var tmask')
           status = nf90_put_att(ncid,varid, 'long_name', 'ocean grid mask') 
           if (status /= nf90_noerr) call abort_ice('ice Error: tmask long_name') 
           status = nf90_put_att(ncid, varid, 'coordinates', 'TLON TLAT')
           if (status /= nf90_noerr) call abort_ice('ice Error: tmask units') 
           status = nf90_put_att(ncid,varid,'comment', '0 = land, 1 = ocean')
           if (status /= nf90_noerr) call abort_ice('ice Error: tmask comment') 
        endif

        do i = 2, nvar       ! note: n_tmask=1
          if (igrd(i)) then
             status = nf90_def_var(ncid, var(i)%req%short_name, &
                                   nf90_float, dimid(1:2), varid)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining variable '//var(i)%req%short_name)
             status = nf90_put_att(ncid,varid, 'long_name', var(i)%req%long_name)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining long_name for '//var(i)%req%short_name)
             status = nf90_put_att(ncid, varid, 'units', var(i)%req%units)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//var(i)%req%short_name)
             status = nf90_put_att(ncid, varid, 'coordinates', var(i)%coordinates)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining coordinates for '//var(i)%req%short_name)
          endif
        enddo

        ! Fields with dimensions (nverts,nx,ny)
        dimid_nverts(1) = nvertexid
        dimid_nverts(2) = imtid
        dimid_nverts(3) = jmtid
        do i = 1, nvar_verts
          if (f_bounds) then
             status = nf90_def_var(ncid, var_nverts(i)%short_name, &
                                   nf90_float,dimid_nverts, varid)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining variable '//var_nverts(i)%short_name)
             status = & 
             nf90_put_att(ncid,varid, 'long_name', var_nverts(i)%long_name)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining long_name for '//var_nverts(i)%short_name)
             status = &
             nf90_put_att(ncid, varid, 'units', var_nverts(i)%units)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//var_nverts(i)%short_name)
          endif
        enddo

        do n=1,num_avail_hist_fields
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                         nf90_float, dimid, varid)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining variable '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'units', &
                        avail_hist_fields(n)%vunit)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining units for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid, 'long_name', &
                        avail_hist_fields(n)%vdesc)


                   
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining long_name for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'coordinates', &
                        avail_hist_fields(n)%vcoord)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining coordinates for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'cell_measures', &
                        avail_hist_fields(n)%vcellmeas)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining cell measures for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining mising_value for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining _FillValue for '//avail_hist_fields(n)%vname)

      !-----------------------------------------------------------------
      ! Append ice thickness range to aicen comments
      !-----------------------------------------------------------------
            c_aice = TRIM(avail_hist_fields(n)%vname)
            i_aice = lenstr(c_aice)
            if (i_aice > 4 .and. c_aice(1:5) == 'aicen') then
              read(c_aice(6:9), '(i3)') icategory
              avail_hist_fields(n)%vcomment = &
                 'Ice range: '//c_hi_range(icategory)
            endif
            status = nf90_put_att(ncid,varid,'comment', &
                        avail_hist_fields(n)%vcomment)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining comment for '//avail_hist_fields(n)%vname)
      !-----------------------------------------------------------------
      ! Add cell_methods attribute to variables if averaged
      !-----------------------------------------------------------------
            if (hist_avg) then
              if (TRIM(avail_hist_fields(n)%vname)/='sig1' &
              .or.TRIM(avail_hist_fields(n)%vname)/='sig2') then
                status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                if (status /= nf90_noerr) call abort_ice( &
                 'Error defining cell methods for '//avail_hist_fields(n)%vname)
              endif
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg         &
                .or. n==n_divu(ns)      .or. n==n_shear(ns)     &  ! snapshots
                .or. n==n_sig1(ns)      .or. n==n_sig2(ns)      &
                .or. n==n_trsig(ns)                             &
                .or. n==n_mlt_onset(ns) .or. n==n_frz_onset(ns) &
                .or. n==n_hisnap(ns)    .or. n==n_aisnap(ns)) then
               status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
            else
               status = nf90_put_att(ncid,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields

      !-----------------------------------------------------------------
      ! global attributes
      !-----------------------------------------------------------------
      ! ... the user should change these to something useful ...
      !-----------------------------------------------------------------
#ifdef CCSMCOUPLED
        status = nf90_put_att(ncid,nf90_global,'title',runid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error in global attribute title')
#else
        title  = 'sea ice model output for CICE'
        status = nf90_put_att(ncid,nf90_global,'title',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error in global attribute title')
#endif
        title = 'Diagnostic and Prognostic Variables'
        status = nf90_put_att(ncid,nf90_global,'contents',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute contents')

        title  = 'sea ice model: Community Ice Code (CICE)'
        status = nf90_put_att(ncid,nf90_global,'source',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute source')

        write(title,'(a,i3,a)') 'All years have exactly ',int(dayyr),' days'
        status = nf90_put_att(ncid,nf90_global,'comment',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute comment')

        write(title,'(a,i8.8)') 'File written on model date ',idate
        status = nf90_put_att(ncid,nf90_global,'comment2',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute date1')

        write(title,'(a,i6)') 'seconds elapsed into model date: ',sec
        status = nf90_put_att(ncid,nf90_global,'comment3',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute date2')

        title = 'CF-1.0'
        status =  &
             nf90_put_att(ncid,nf90_global,'conventions',title)
        if (status /= nf90_noerr) call abort_ice( &
             'Error in global attribute conventions')

        call date_and_time(date=current_date, time=current_time)
        write(start_time,1000) current_date(1:4), current_date(5:6), &
                               current_date(7:8), current_time(1:2), &
                               current_time(3:4), current_time(5:8)
1000    format('This dataset was created on ', &
                a,'-',a,'-',a,' at ',a,':',a,':',a)

        status = nf90_put_att(ncid,nf90_global,'history',start_time)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute history')

      !-----------------------------------------------------------------
      ! end define mode
      !-----------------------------------------------------------------

        status = nf90_enddef(ncid)
        if (status /= nf90_noerr) call abort_ice('ice: Error in nf90_enddef')

      !-----------------------------------------------------------------
      ! write time variable
      !-----------------------------------------------------------------

        status = nf90_inq_varid(ncid,'time',varid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error getting time varid')
        status = nf90_put_var(ncid,varid,ltime)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error writing time variable')

      !-----------------------------------------------------------------
      ! write time_bounds info
      !-----------------------------------------------------------------

        if (hist_avg) then
          status = nf90_inq_varid(ncid,'time_bounds',varid)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error getting time_bounds id')
          status = nf90_put_var(ncid,varid,time_beg(ns),start=(/1/))
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error writing time_beg')
          status = nf90_put_var(ncid,varid,time_end(ns),start=(/2/))
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error writing time_end')
        endif

      endif                     ! master_task

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         allocate(work_gr(nx_global,ny_global))
      else
         allocate(work_gr(1,1))   ! to save memory
         allocate(work_g1(1,1))
      endif

      work_g1(:,:) = c0

      !-----------------------------------------------------------------
      ! write coordinate variables
      !-----------------------------------------------------------------

        do i = 1,ncoord
          call broadcast_scalar(coord_var(i)%short_name,master_task)
          SELECT CASE (coord_var(i)%short_name)
            CASE ('TLON')
              call gather_global(work_g1,TLON,master_task,distrb_info)
              if (my_task == master_task) then
              ! Convert T grid longitude from -180 -> 180 to 0 to 360
                 work_gr = work_g1*rad_to_deg + c360    ! single precision
                 where (work_gr > c360) work_gr = work_gr - c360
                 where (work_gr < c0 )  work_gr = work_gr + c360
              endif
            CASE ('TLAT')
              call gather_global(work_g1,TLAT,master_task,distrb_info)
              if (my_task == master_task) work_gr = work_g1*rad_to_deg
            CASE ('ULON')
              call gather_global(work_g1,ULON,master_task,distrb_info)
              if (my_task == master_task) work_gr = work_g1*rad_to_deg
            CASE ('ULAT')
              call gather_global(work_g1,ULAT,master_task,distrb_info)
              if (my_task == master_task) work_gr = work_g1*rad_to_deg
          END SELECT
          
          if (my_task == master_task) then
             status = nf90_inq_varid(ncid, coord_var(i)%short_name, varid)
             if (status /= nf90_noerr) call abort_ice( &
                  'ice: Error getting varid for '//coord_var(i)%short_name)
             status = nf90_put_var(ncid,varid,work_gr)
             if (status /= nf90_noerr) call abort_ice( &
                           'ice: Error writing'//coord_var(i)%short_name)
          endif
        enddo

      !-----------------------------------------------------------------
      ! write grid mask, area and rotation angle
      !-----------------------------------------------------------------

      if (igrd(n_tmask)) then
      call gather_global(work_g1, hm, master_task, distrb_info)
      if (my_task == master_task) then
        work_gr=work_g1
        status = nf90_inq_varid(ncid, 'tmask', varid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error getting varid for tmask')
        status = nf90_put_var(ncid,varid,work_gr)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error writing variable tmask')
      endif
      endif

      do i = 2, nvar       ! note: n_tmask=1
        if (igrd(i)) then
        call broadcast_scalar(var(i)%req%short_name,master_task)
        SELECT CASE (var(i)%req%short_name)
          CASE ('tarea')
            call gather_global(work_g1, tarea, master_task, distrb_info)
          CASE ('uarea')
            call gather_global(work_g1, uarea, master_task, distrb_info)
          CASE ('dxu')
            call gather_global(work_g1,   dxu, master_task, distrb_info)
          CASE ('dyu')
            call gather_global(work_g1,   dyu, master_task, distrb_info)
          CASE ('dxt')
            call gather_global(work_g1,   dxt, master_task, distrb_info)
          CASE ('dyt')
            call gather_global(work_g1,   dyt, master_task, distrb_info)
          CASE ('HTN')
            call gather_global(work_g1,   HTN, master_task, distrb_info)
          CASE ('HTE')
            call gather_global(work_g1,   HTE, master_task, distrb_info)
          CASE ('ANGLE')
            call gather_global(work_g1, ANGLE, master_task, distrb_info)
          CASE ('ANGLET')
            call gather_global(work_g1, ANGLET,master_task, distrb_info)
        END SELECT

        if (my_task == master_task) then
          work_gr=work_g1
          status = nf90_inq_varid(ncid, var(i)%req%short_name, varid)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error getting varid for '//var(i)%req%short_name)
          status = nf90_put_var(ncid,varid,work_gr)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error writing variable '//var(i)%req%short_name)
        endif
        endif
      enddo

      deallocate(work_gr)

      !----------------------------------------------------------------
      ! Write coordinates of grid box vertices
      !----------------------------------------------------------------

      if (f_bounds) then
      if (my_task==master_task) then
         allocate(work_gr3(nverts,nx_global,ny_global))
      else
         allocate(work_gr3(1,1,1))   ! to save memory
      endif

      work_gr3(:,:,:) = c0
      work1   (:,:,:) = c0

      do i = 1, nvar_verts
        call broadcast_scalar(var_nverts(i)%short_name,master_task)
        SELECT CASE (var_nverts(i)%short_name)
        CASE ('lont_bounds')
        do ivertex = 1, nverts 
           work1(:,:,:) = lont_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        CASE ('latt_bounds')
        do ivertex = 1, nverts
           work1(:,:,:) = latt_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        CASE ('lonu_bounds')
        do ivertex = 1, nverts
           work1(:,:,:) = lonu_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        CASE ('latu_bounds')
        do ivertex = 1, nverts
           work1(:,:,:) = latu_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        END SELECT

        if (my_task == master_task) then
          status = nf90_inq_varid(ncid, var_nverts(i)%short_name, varid)
          if (status /= nf90_noerr) call abort_ice( &
             'ice: Error getting varid for '//var_nverts(i)%short_name)
          status = nf90_put_var(ncid,varid,work_gr3)
          if (status /= nf90_noerr) call abort_ice( &
             'ice: Error writing variable '//var_nverts(i)%short_name)
        endif
      enddo
      deallocate(work_gr3)
      endif

      !-----------------------------------------------------------------
      ! write variable data
      !-----------------------------------------------------------------

      if (my_task==master_task) then
         allocate(work_gr3(nx_global,ny_global,1))
      else
         allocate(work_gr3(1,1,1))   ! to save memory
      endif

      do n=1,num_avail_hist_fields
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
          call gather_global(work_g1, aa(:,:,n,:), &
                             master_task, distrb_info)
          if (my_task == master_task) then
            work_gr3(:,:,1) = work_g1(:,:)
            status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
            if (status /= nf90_noerr) call abort_ice( &
               'ice: Error getting varid for '//avail_hist_fields(n)%vname)
            status  = nf90_put_var(ncid,varid,work_gr3, &
                                   count=(/nx_global,ny_global,1/))
            if (status /= nf90_noerr) call abort_ice( &
               'ice: Error writing variable '//avail_hist_fields(n)%vname)
          endif
        endif
      enddo ! num_avail_hist_fields

      deallocate(work_gr3)
      deallocate(work_g1)

      !-----------------------------------------------------------------
      ! close output dataset
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         status = nf90_close(ncid)
         if (status /= nf90_noerr) call abort_ice( &
                       'ice: Error closing netCDF history file')
         write(nu_diag,*) ' '
         write(nu_diag,*) 'Finished writing ',trim(ncfile(ns))
      endif
#endif

      end subroutine icecdf

!=======================================================================
!
!BOP
!
! !IROUTINE: icebin - write binary history file
! This routine writes fewer grid variables compared with the netcdf
! version, to reduce file size.  Grid variables can be obtained from
! the original grid input files.
!
! !INTERFACE:
!
      subroutine icebin(ns)
!
! !DESCRIPTION:
!
! write binary history file
!
! !REVISION HISTORY:
!
! authors:   E.C.Hunke, LANL
!
! !USES:
!
      use ice_gather_scatter
      use ice_domain_size
      use ice_constants
      use ice_grid
      use ice_restart, only: lenstr, runid
      use ice_itd, only: c_hi_range
      use ice_calendar, only: write_ic, dayyr, histfreq
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: ns
!
!EOP
!
      integer (kind=int_kind) :: i,j,n,nrec,nbits
      character (char_len) :: title
      character (char_len_long) :: ncfile(max_nstrm), hdrfile

      integer (kind=int_kind) :: icategory,i_aice

      character (char_len) :: current_date,current_time
      character (len=16) :: c_aice
      logical (kind=log_kind) :: diag

      diag = .false.

      if (my_task == master_task) then

        call construct_filename(ncfile(ns),'da',ns)

        ! add local directory path name to ncfile
        if (write_ic) then
          ncfile(ns) = trim(incond_dir)//ncfile(ns)
        else
          ncfile(ns) = trim(history_dir)//ncfile(ns)
        endif
        hdrfile = trim(ncfile(ns))//'.hdr'

        !-----------------------------------------------------------------
        ! create history files
        !-----------------------------------------------------------------
        nbits = 32 ! single precision
        call ice_open(nu_history, ncfile(ns), nbits) ! direct access
        open(nu_hdr,file=hdrfile,form='formatted',status='unknown') ! ascii

!echmod call ice_write(nu_history, nrec, work, rda8 or ida4, diag)

        title  = 'sea ice model: Community Ice Code (CICE)'
        write (nu_hdr, 999) 'source',title,' '

        write (nu_hdr, 999) 'file name contains model date',trim(ncfile(ns)),' '
#ifdef CCSMCOUPLED
        write (nu_hdr, 999) 'runid',runid,' '
#endif
        write (nu_hdr, 999) 'calendar','noleap',' '
        write (title,'(a,i3,a)') 'All years have exactly ',int(dayyr),' days'
        write (nu_hdr, 999) 'comment',title,' '
        write (nu_hdr, 999) 'conventions','CICE',' '
        write (nu_hdr, 997) 'missing_value',spval
        write (nu_hdr, 997) '_FillValue',spval

        call date_and_time(date=current_date, time=current_time)
        write (nu_hdr,1000) current_date(1:4), current_date(5:6), &
                            current_date(7:8), current_time(1:2), &
                            current_time(3:4), current_time(5:8)
        write (nu_hdr, *  ) ' '
        write (nu_hdr, *  ) 'Grid size:'
        write (nu_hdr, 998) '  ni',nx_global
        write (nu_hdr, 998) '  nj',ny_global

        write (nu_hdr, *  ) 'Grid variables: (left column = nrec)'
        nrec = 1
        write (nu_hdr, 996) nrec,'tarea','area of T grid cells','m^2'
        write (nu_hdr, *  ) 'History variables: (left column = nrec)'
      endif  ! my_task = master_task
      call ice_write(nu_history, nrec, tarea, 'rda4', diag)

      do n=1,num_avail_hist_fields
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then

          nrec = nrec + 1
          if (my_task == master_task) then
            write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
               trim(avail_hist_fields(n)%vdesc),trim(avail_hist_fields(n)%vunit)

            ! Append ice thickness range to aicen comments
            c_aice = TRIM(avail_hist_fields(n)%vname)
            i_aice = lenstr(c_aice)
            if (i_aice > 4 .and. c_aice(1:5) == 'aicen') then
              read(c_aice(6:9), '(i3)') icategory
              avail_hist_fields(n)%vcomment = &
                 'Ice range: '//c_hi_range(icategory)
            endif
            write (nu_hdr, 995) nrec,trim(avail_hist_fields(n)%vname), &
               trim(avail_hist_fields(n)%vcomment)

            if (histfreq(ns) == '1' .or. .not. hist_avg         &
                .or. n==n_divu(ns)      .or. n==n_shear(ns)     &  ! snapshots
                .or. n==n_sig1(ns)      .or. n==n_sig2(ns)      &
                .or. n==n_trsig(ns)                             &
                .or. n==n_mlt_onset(ns) .or. n==n_frz_onset(ns) &
                .or. n==n_hisnap(ns)    .or. n==n_aisnap(ns)) then
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','instantaneous'
            else
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','averaged'
            endif
          endif

          call ice_write(nu_history, nrec, aa(:,:,n,:), 'rda4', diag)

        endif
      enddo ! num_avail_hist_fields

995     format(i3,2x,a,' comment: ',a)
996     format(i3,2x,a,': ',a,',',2x,a)
997     format(a,': ',es13.6)
998     format(a,': ',i6)
999     format(a,': ',a,2x,a)
1000    format('This dataset was created on ', &
                a,'-',a,'-',a,' at ',a,':',a,':',a)

      if (my_task == master_task) then
        close (nu_hdr)     ! header file
        close (nu_history) ! data file
        write (nu_diag,*) ' '
        write (nu_diag,*) 'Finished writing ',trim(ncfile(ns))
      endif

      end subroutine icebin

!=======================================================================

      subroutine construct_filename(ncfile,suffix,ns)

      use ice_calendar, only: time, sec, idate, nyr, month, daymo,  &
                              mday, write_ic, histfreq, histfreq_n, &
                              year_init, new_year, new_month, new_day, &
                              dayyr, dt
      use ice_restart, only: lenstr

      character (char_len_long), intent(inout) :: ncfile
      character (len=2), intent(in) :: suffix
      integer (kind=int_kind), intent(in) :: ns

      integer (kind=int_kind) :: iyear, imonth, iday, isec

      character (char_len_long) :: tmpfile

        iyear = nyr + year_init - 1 ! set year_init=1 in ice_in to get iyear=nyr
        imonth = month
        iday = mday
        isec = sec - dt

        ! construct filename
        if (write_ic) then
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
              incond_file(1:lenstr(incond_file)),'.',iyear,'-', &
              imonth,'-',iday,'-',isec,'.',suffix
        else

         if (hist_avg) then
          if (histfreq(ns).eq.'h'.or.histfreq(ns).eq.'H') then
           ! do nothing
          elseif (new_year) then
           iyear = iyear - 1
           imonth = 12
           iday = daymo(imonth)
          elseif (new_month) then
           imonth = month - 1
           iday = daymo(imonth)
          elseif (new_day) then
           iday = iday - 1
          endif
         endif

         if (histfreq(ns) == '1') then ! instantaneous, write every dt
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
            history_file(1:lenstr(history_file)),'_inst.', &
             iyear,'-',imonth,'-',iday,'-',sec,'.',suffix

         elseif (hist_avg) then    ! write averaged data

          if (histfreq(ns).eq.'d'.or.histfreq(ns).eq.'D') then     ! daily
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,a)')  &
            history_file(1:lenstr(history_file)), &
             '.',iyear,'-',imonth,'-',iday,'.',suffix
          elseif (histfreq(ns).eq.'h'.or.histfreq(ns).eq.'H') then ! hourly
           write(ncfile,'(a,a,i2.2,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
            history_file(1:lenstr(history_file)),'_',histfreq_n(ns),'h.', &
             iyear,'-',imonth,'-',iday,'-',sec,'.',suffix
          elseif (histfreq(ns).eq.'m'.or.histfreq(ns).eq.'M') then ! monthly
           write(ncfile,'(a,a,i4.4,a,i2.2,a,a)')  &
            history_file(1:lenstr(history_file)),'.', &
             iyear,'-',imonth,'.',suffix
          elseif (histfreq(ns).eq.'y'.or.histfreq(ns).eq.'Y') then ! yearly
           write(ncfile,'(a,a,i4.4,a,a)') &
            history_file(1:lenstr(history_file)),'.', iyear,'.',suffix
          endif

         else                     ! instantaneous with histfreq > dt
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
            history_file(1:lenstr(history_file)),'_inst.', &
             iyear,'-',imonth,'-',iday,'-',sec,'.',suffix
         endif
        endif

      end subroutine construct_filename

!=======================================================================

      subroutine define_hist_field(id, vname, vunit, vcoord, vcellmeas, &
                                   vdesc, vcomment, cona, conb, &
                                   ns1, vhistfreq)

!     !DESCRIPTION:
!     Initializes description of an available field and returns location
!     in the available fields array for use in later calls.
!
!     !REVISION HISTORY:
!     2009 Created by D. Bailey following POP

!     !USES:
      use ice_exit
      use ice_calendar, only: histfreq, histfreq_n, nstreams

!     !OUTPUT PARAMETERS:

      integer (int_kind), dimension(max_nstrm), intent(out) :: &
         id                ! location in avail_fields array for use in
                           ! later routines

!     !INPUT PARAMETERS

      character (len=*), intent(in) :: &
         vname      , & ! variable names
         vunit      , & ! variable units
         vcoord     , & ! variable coordinates
         vcellmeas  , & ! variables cell measures
         vdesc      , & ! variable descriptions
         vcomment       ! variable comments

      real (kind=dbl_kind), intent(in) :: &
         cona       , & ! multiplicative conversion factor
         conb           ! additive conversion factor

      character (len=*), intent(in) :: &
         vhistfreq      ! history frequency
 
      integer (kind=int_kind), intent(in) :: &
         ns1            ! stream index

      integer (kind=int_kind) :: &
         ns         , & ! loop index
         lenf           ! length of namelist string

      character (len=40) :: stmp

      lenf = len(trim(vhistfreq))
      if (ns1 == 1) id(:) = 0

      do ns = 1, nstreams
         if (vhistfreq(ns1:ns1) == histfreq(ns)) then

            num_avail_hist_fields = num_avail_hist_fields + 1

            if (num_avail_hist_fields > max_avail_hist_fields) &
               call abort_ice("Need to increase max_avail_hist_fields")

            id(ns) = num_avail_hist_fields

            stmp = vname
            if (lenf > 1 .and. ns1 > 1) &
               write(stmp,'(a,a1,a1)') trim(stmp),'_',vhistfreq(ns1:ns1)

            avail_hist_fields(id(ns))%vname = trim(stmp)
            avail_hist_fields(id(ns))%vunit = trim(vunit)
            avail_hist_fields(id(ns))%vcoord = trim(vcoord)
            avail_hist_fields(id(ns))%vcellmeas = trim(vcellmeas)
            avail_hist_fields(id(ns))%vdesc = trim(vdesc)
            avail_hist_fields(id(ns))%vcomment = trim(vcomment)
            avail_hist_fields(id(ns))%cona = cona
            avail_hist_fields(id(ns))%conb = conb
            avail_hist_fields(id(ns))%vhistfreq = vhistfreq(ns1:ns1)
            avail_hist_fields(id(ns))%vhistfreq_n = histfreq_n(ns)

           endif
        enddo

      end subroutine define_hist_field

!=======================================================================

      subroutine accum_hist_field(id, iblk, field_accum)

!     !DESCRIPTION:
!     Accumulates a history field
!
!     !REVISION HISTORY:
!     2009 Created by D. Bailey following POP

      use ice_domain
      use ice_grid, only: tmask
      use ice_calendar, only: nstreams

!     !OUTPUT PARAMETERS:

      integer (int_kind), dimension(max_nstrm), intent(in) :: &
         id                ! location in avail_fields array for use in
                           ! later routines

      integer (kind=int_kind), intent(in) :: iblk

      real (kind=dbl_kind), intent(in) :: field_accum(nx_block,ny_block)

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: i,j, ilo, ihi, jlo, jhi, ns, idns

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

       do ns = 1, nstreams
       idns = id(ns)
       if (idns > 0) then

       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
       do i = ilo, ihi
          if (tmask(i,j,iblk)) then
             aa(i,j,idns, iblk) = aa(i,j,idns, iblk) + field_accum(i,j)
          endif
       enddo
       enddo

       endif
       enddo

      end subroutine accum_hist_field

!=======================================================================

      end module ice_history

!=======================================================================
