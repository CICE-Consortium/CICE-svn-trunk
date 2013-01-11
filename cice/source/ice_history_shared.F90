!=======================================================================
!
!BOP
!
! !MODULE: ice_history_shared - common ice model history files
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
! 2010 Alison McLaren and ECH: Added 3D capability
!
! !INTERFACE:
!
      module ice_history_shared
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
          character (len=25) :: vcoord    ! variable coordinates
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
         num_avail_hist_fields_tot  = 0, & ! Current, total number of defined fields
         num_avail_hist_fields_2D   = 0, & ! Number of 2D fields
         num_avail_hist_fields_3Dz  = 0, & ! Number of 3D fields (vertical)
         num_avail_hist_fields_3Dc  = 0, & ! Number of 3D fields (categories)
         num_avail_hist_fields_3Db  = 0, & ! Number of 3D fields (vertical biology)
         num_avail_hist_fields_4Di  = 0, & ! Number of 4D fields (categories,vertical), ice
         num_avail_hist_fields_4Ds  = 0, & ! Number of 4D fields (categories,vertical), snow
         num_avail_hist_fields_4Db  = 0    ! Number of 4D fields (categories,vertical), ice-biology

      integer (kind=int_kind) :: &        ! cumulative counts
         n2D     , & ! num_avail_hist_fields_2D
         n3Dccum , & ! n2D     + num_avail_hist_fields_3Dc
         n3Dzcum , & ! n3Dccum + num_avail_hist_fields_3Dz
         n3Dbcum , & ! n3Dzcum + num_avail_hist_fields_3Db
         n4Dicum , & ! n3Dbcum + num_avail_hist_fields_4Di
         n4Dscum , & ! n4Dicum + num_avail_hist_fields_4Ds
         n4Dbcum , & ! n4Dscum + num_avail_hist_fields_4Db
         nzlyr   , & ! vertical dimension (temp variable)
         nzlyrb      ! vertical dimension of biology grid (temp variable)

      ! for now, ice and snow have same dimensions in netcdf
      ! could set nzilyr = nilyr + nslyr and write Tin+Tsn together into Tinz
      integer (kind=int_kind), parameter :: &
         nzilyr = nilyr         , & ! vertical dimension (allows alternative grids)
         nzslyr = nslyr         , &
         nzblyr = nblyr+2

      type (ice_hist_field), dimension(max_avail_hist_fields) :: &
         avail_hist_fields

      character (len=16) :: vname_in     ! variable name
      character (len=55) :: vdesc_in     ! variable description
      character (len=55) :: vcomment_in  ! variable description

      integer (kind=int_kind), parameter :: &
         nvar = 11              , & ! number of grid fields that can be written
                                    !   excluding grid vertices
         nvarz = 4              , & ! number of category/vertical grid fields written
         ncat_hist = ncat           ! number of ice categories written <= ncat

      real (kind=real_kind) :: time_beg(max_nstrm), &
                               time_end(max_nstrm) ! bounds for averaging

      real (kind=dbl_kind), allocatable :: &
         a2D (:,:,:,:)    , & ! field accumulations/averages, 2D
         a3Dz(:,:,:,:,:)  , & ! field accumulations/averages, 3D vertical
         a3Db(:,:,:,:,:)  , & ! field accumulations/averages, 3D vertical biology
         a3Dc(:,:,:,:,:)  , & ! field accumulations/averages, 3D categories
         a4Di(:,:,:,:,:,:), & ! field accumulations/averages, 4D categories,vertical, ice
         a4Ds(:,:,:,:,:,:), & ! field accumulations/averages, 4D categories,vertical, snow
         a4Db(:,:,:,:,:,:)    ! field accumulations/averages, 4D categories,vertical, bio
         
      real (kind=dbl_kind), allocatable :: &
         Tinz4d (:,:,:,:)    , & ! array for Tin
         Tsnz4d (:,:,:,:)    , & ! array for Tsn
         Sinz4d (:,:,:,:)        ! array for Sin

      real (kind=dbl_kind) :: &
         avgct(max_nstrm)   ! average sample counter

      logical (kind=log_kind) :: &
         igrd (nvar), &        ! true if grid field is written to output file
         igrdz(nvarz)          ! true if category/vertical grid field is written

      character (len=25), parameter :: &
         tcstr = 'area: tarea'          , & ! vcellmeas for T cell quantities
         ucstr = 'area: uarea'          , & ! vcellmeas for U cell quantities
         tstr2D  = 'TLON TLAT time'     , & ! vcoord for T cell quantities, 2D
         ustr2D  = 'ULON ULAT time'     , & ! vcoord for U cell quantities, 2D
         tstr3Dz = 'TLON TLAT VGRD time', & ! vcoord for T cell quantities, 3D
         ustr3Dz = 'ULON ULAT VGRD time', & ! vcoord for U cell quantities, 3D
         tstr3Dc = 'TLON TLAT NCAT time', & ! vcoord for T cell quantities, 3D
         ustr3Dc = 'ULON ULAT NCAT time', & ! vcoord for U cell quantities, 3D
         tstr3Db = 'TLON TLAT VGRDb time', & ! vcoord for T cell quantities, 3D
         ustr3Db = 'ULON ULAT VGRDb time', & ! vcoord for U cell quantities, 3D

!ferret
         tstr4Di  = 'TLON TLAT VGRDi NCAT', & ! vcoord for T cell quantities, 4D, ice
         ustr4Di  = 'ULON ULAT VGRDi NCAT', & ! vcoord for U cell quantities, 4D, ice
         tstr4Ds  = 'TLON TLAT VGRDs NCAT', & ! vcoord for T cell quantities, 4D, snow
         ustr4Ds  = 'ULON ULAT VGRDs NCAT', & ! vcoord for U cell quantities, 4D, snow
         tstr4Db  = 'TLON TLAT VGRDb NCAT', & ! vcoord for T cell quantities, 4D, bio
         ustr4Db  = 'ULON ULAT VGRDb NCAT'    ! vcoord for U cell quantities, 4D, bio
!ferret
!         tstr4Di  = 'TLON TLAT VGRDi NCAT time', & ! ferret can not handle time 
!         ustr4Di  = 'ULON ULAT VGRDi NCAT time', & ! index on 4D variables.
!         tstr4Ds  = 'TLON TLAT VGRDs NCAT time', & ! Use 'ferret' lines instead
!         ustr4Ds  = 'ULON ULAT VGRDs NCAT time', & ! (below also)
!         tstr4Db  = 'TLON TLAT VGRDb NCAT time', & 
!         ustr4Db  = 'ULON ULAT VGRDb NCAT time'    

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
           f_bounds    = .true., f_NCAT       = .true., &
           f_VGRDi     = .true., f_VGRDs      = .true., &
           f_VGRDb     = .true.

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
           f_snoice    = 'm', f_dsnow = 'm', f_meltt      = 'm', &
           f_melts     = 'm', &
           f_meltb     = 'm', f_meltl      = 'm', &
           f_fresh     = 'm', f_fresh_ai   = 'm', &
           f_fsalt     = 'm', f_fsalt_ai   = 'm', &
           f_fsice     = 'm', f_fsice_ai   = 'm', &
           f_fsice_g   = 'm', f_fsice_g_ai = 'm', &
           f_fNO     = 'm', f_fNO_ai   = 'm', &
           f_fNO_g   = 'm', f_fNO_g_ai = 'm', &
           f_fNH     = 'm', f_fNH_ai   = 'm', &
           f_fNH_g   = 'm', f_fNH_g_ai = 'm', &
           f_fN     = 'm', f_fN_ai   = 'm', &
           f_fN_g   = 'm', f_fN_g_ai = 'm', &
           f_fSil     = 'm', f_fSil_ai   = 'm', &
           f_fSil_g   = 'm', f_fSil_g_ai = 'm', &
           f_fsicen_g  = 'm',                     &
           f_fhocn     = 'm', f_fhocn_ai   = 'm', &
           f_fswthru   = 'm', f_fswthru_ai = 'm', &
           f_strairx   = 'm', f_strairy    = 'm', &
           f_strtltx   = 'm', f_strtlty    = 'm', &
           f_strcorx   = 'm', f_strcory    = 'm', &
           f_strocnx   = 'm', f_strocny    = 'm', &
           f_strintx   = 'm', f_strinty    = 'm', &
           f_strength  = 'm', &
           f_divu      = 'm', f_shear      = 'm', &
           f_sig1      = 'm', f_sig2       = 'm', &
           f_dvidtt    = 'm', f_dvidtd     = 'm', &
           f_daidtt    = 'm', f_daidtd     = 'm', &
           f_mlt_onset = 'm', f_frz_onset  = 'm', &
           f_iage      = 'm', f_FY         = 'm', &
           f_hisnap    = 'm', f_aisnap     = 'm', &
           f_aicen     = 'x', f_vicen      = 'x', &
           f_trsig     = 'm', f_icepresent = 'm', &
           f_fsurf_ai  = 'm', f_fcondtop_ai= 'm', &
           f_fmeltt_ai = 'm',                     &
           f_fsurfn_ai = 'x' ,f_fcondtopn_ai='x', &
           f_fmelttn_ai= 'x', f_flatn_ai   = 'x', &
!          f_field3dz  = 'x',                     &
           f_Tinz      = 'x', f_Sinz       = 'x', &
           f_Tsnz      = 'x', &
           f_bgc_N_sk  = 'x',    f_bgc_C_sk= 'x', &
           f_bgc_chl_sk= 'x', f_bgc_Nit_sk = 'x', &
           f_bgc_Am_sk = 'x',  f_bgc_Sil_sk= 'x', &
           f_bgc_DMSPp_sk = 'x', f_bgc_DMSPd_sk = 'x', &
           f_bgc_DMS_sk = 'x', & 
           f_bgc_Nit_ml = 'x',  f_bgc_Am_ml = 'x', & 
           f_bgc_Sil_ml = 'x', &  
           f_bgc_DMSP_ml = 'x', f_bgc_DMS_ml = 'x', & 
           f_upNO      = 'x',   f_upNH        = 'x',   & 
           f_zTin       = 'x',      &
           f_zphi       = 'x',    &
           f_iDi       = 'x',  f_iki           = 'x',    &
           f_bgc_NO       = 'x',    &
           f_bgc_N     = 'x',  f_bgc_NH       = 'x',    &
           f_bgc_C     = 'x',  f_bgc_chl      = 'x',    &
           f_bgc_DMSPp = 'x',  f_bgc_DMSPd    = 'x',    &
           f_bgc_DMS   = 'x',  f_bgc_Sil      = 'x',   &
           f_bgc_PON   = 'x',  f_bgc_S        = 'x',   &
           f_fbri = 'x',          &
           f_growN     = 'x', f_zfswin      = 'x', &
           f_Stot = 'x', f_chlnet = 'x',  &
           f_PPnet = 'x', f_NOnet = 'x', &
           f_a11       = 'x', f_a12        = 'x' , & 
           f_e11       = 'x', f_e12        = 'x' , & 
           f_e22       = 'x',			  & 
           f_s11       = 'x', f_s12        = 'x', & 
           f_s22       = 'x',		          & 
           f_yieldstress11       = 'x', 	  & 
           f_yieldstress12       = 'x',           & 
           f_yieldstress22       = 'x'

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
           f_bounds   , f_NCAT     , &
           f_VGRDi    , f_VGRDs    , &
           f_VGRDb    , &
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
           f_snoice,  f_dsnow,   f_meltt    , &
           f_melts,                  &
           f_meltb,     f_meltl    , &
           f_fresh,     f_fresh_ai , &  
           f_fsalt,     f_fsalt_ai , &  
           f_fsice,     f_fsice_ai , &   
           f_fsice_g,   f_fsice_g_ai , &  
           f_fsicen_g, &
           f_fNO,       f_fNO_ai , &
           f_fNO_g,     f_fNO_g_ai, &
           f_fNH,       f_fNH_ai, &
           f_fNH_g,     f_fNH_g_ai, &
           f_fN,        f_fN_ai, &
           f_fN_g,      f_fN_g_ai, &
           f_fSil,      f_fSil_ai, &
           f_fSil_g,    f_fSil_g_ai, &
           f_fhocn,     f_fhocn_ai , &
           f_fswthru,   f_fswthru_ai,&
           f_strairx,   f_strairy  , &
           f_strtltx,   f_strtlty  , &
           f_strcorx,   f_strcory  , &
           f_strocnx,   f_strocny  , &
           f_strintx,   f_strinty  , &
           f_strength,  &
           f_divu,      f_shear    , &
           f_sig1,      f_sig2     , &
           f_dvidtt,    f_dvidtd   , &
           f_daidtt,    f_daidtd   , &
           f_mlt_onset, f_frz_onset, &
           f_iage,      f_FY       , &
           f_hisnap,    f_aisnap   , &
           f_aicen,     f_vicen    , &
           f_trsig,     f_icepresent,&
           f_fsurf_ai,  f_fcondtop_ai,&
           f_fmeltt_ai, &
           f_fsurfn_ai,f_fcondtopn_ai,&
           f_fmelttn_ai,f_flatn_ai,  &
!          f_field3dz,  &
           f_Tinz,      f_Sinz,      &
           f_Tsnz,  &
           f_bgc_N_sk,    f_bgc_C_sk,   f_bgc_chl_sk, & 
           f_bgc_Nit_sk,  f_bgc_Am_sk,  f_bgc_Sil_sk, &
           f_bgc_DMSPp_sk, f_bgc_DMSPd_sk, f_bgc_DMS_sk, & 
           f_bgc_Nit_ml,  f_bgc_Am_ml,  f_bgc_Sil_ml, &  
           f_bgc_DMSP_ml, f_bgc_DMS_ml, &
           f_upNO,  f_upNH, &       
           f_zTin                  , &
           f_zphi      , &
           f_iDi,       f_iki       , &
           f_bgc_NO   , f_bgc_N    , & 
           f_bgc_NH  , &
           f_bgc_C    , f_bgc_chl  , &
           f_bgc_Sil  , f_bgc_DMSPp, &
           f_bgc_DMSPd , f_bgc_DMS , &
           f_bgc_PON  , f_bgc_S, &
           f_fbri, &
           f_growN,    f_zfswin, &
           f_Stot,    f_chlnet, &
           f_PPnet, f_NOnet, &
           f_a11, 	f_a12 	   , &
           f_e11, 	f_e12	   , &
           f_e22                   , &
           f_s11, 	f_s12	   , &
           f_s22                   , &
           f_yieldstress11         , &	
           f_yieldstress12	   , &
           f_yieldstress22

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

           n_NCAT       = 1, &
           n_VGRDi      = 2, &
           n_VGRDs      = 3, &
           n_VGRDb      = 4, &

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
           n_snoice, n_dsnow     , n_meltt      , &
           n_melts      , &
           n_meltb      , n_meltl      , &
           n_fresh      , n_fresh_ai   , &
           n_fsalt      , n_fsalt_ai   , &
           n_fsice      , n_fsice_ai   , &
           n_fsice_g    , n_fsice_g_ai   , &
           n_fsicen_g, & 
           n_fNO        , n_fNO_ai , &
           n_fNO_g      , n_fNO_g_ai, &
           n_fNH        , n_fNH_ai, &
           n_fNH_g      , n_fNH_g_ai, &
           n_fN         , n_fN_ai, &
           n_fN_g       , n_fN_g_ai, &
           n_fSil       , n_fSil_ai, &
           n_fSil_g     , n_fSil_g_ai, &
           n_fhocn      , n_fhocn_ai   , &
           n_fswthru    , n_fswthru_ai , &
           n_strairx    , n_strairy    , &
           n_strtltx    , n_strtlty    , &
           n_strcorx    , n_strcory    , &
           n_strocnx    , n_strocny    , &
           n_strintx    , n_strinty    , &
           n_strength   , &
           n_divu       , n_shear      , &
           n_sig1       , n_sig2       , &
           n_dvidtt     , n_dvidtd     , &
           n_daidtt     , n_daidtd     , &
           n_mlt_onset  , n_frz_onset  , &
           n_hisnap     , n_aisnap     , &
           n_trsig      , n_icepresent , &
           n_iage       , n_FY         , &
           n_fsurf_ai   , &
           n_fcondtop_ai, n_fmeltt_ai  , &   
           n_bgc_N_sk , &
           n_bgc_C_sk, &
           n_bgc_chl_sk, &
           n_bgc_Nit_sk, &
           n_bgc_Am_sk, &
           n_bgc_Sil_sk, &
           n_bgc_Nit_ml, &
           n_bgc_Am_ml, &
           n_bgc_Sil_ml, &
           n_bgc_DMSPp_sk, &
           n_bgc_DMSPd_sk, &
           n_bgc_DMS_sk , &
           n_bgc_DMSP_ml, &
           n_bgc_DMS_ml, &
           n_aicen      , n_vicen      , &
           n_fsurfn_ai   , &
           n_fcondtopn_ai, &
           n_fmelttn_ai  , &
           n_flatn_ai    , &
!          n_field3dz    , &
           n_Tinz        , n_Sinz      , &
           n_Tsnz, &
           n_upNO,  &
           n_upNH,  &
           n_zTin         , & 
           n_zphi, &
           n_iDi, &
           n_iki, &
           n_bgc_NO, &
           n_bgc_N, &
           n_bgc_NH, &
           n_bgc_C, &
           n_bgc_chl, &
           n_bgc_DMSPp, &
           n_bgc_DMSPd, &
           n_bgc_DMS, &
           n_bgc_Sil, &
           n_bgc_PON, &
           n_bgc_S,  &
           n_fbri, &
           n_growN, &
           n_zfswin, &
           n_Stot, &
           n_chlnet, &
           n_PPnet, &
           n_NOnet, &
	   n_a11	 , n_a12	, &
	   n_e11	 , n_e12 	, &
	   n_e22	 , &
	   n_s11	 , n_s12	, &
	   n_s22	 , &
	   n_yieldstress11, n_yieldstress12,  &
	   n_yieldstress22

      interface accum_hist_field ! generic interface
           module procedure accum_hist_field_2D, &
                            accum_hist_field_3D, &
                            accum_hist_field_4D
      end interface

!=======================================================================

      contains

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

         if (hist_avg .and. histfreq(ns) /= '1') then
          if (histfreq(ns) == 'h'.or.histfreq(ns) == 'H') then
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

          if (histfreq(ns) == 'd'.or.histfreq(ns) == 'D') then     ! daily
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,a)')  &
            history_file(1:lenstr(history_file)), &
             '.',iyear,'-',imonth,'-',iday,'.',suffix
          elseif (histfreq(ns) == 'h'.or.histfreq(ns) == 'H') then ! hourly
           write(ncfile,'(a,a,i2.2,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
            history_file(1:lenstr(history_file)),'_',histfreq_n(ns),'h.', &
             iyear,'-',imonth,'-',iday,'-',sec,'.',suffix
          elseif (histfreq(ns) == 'm'.or.histfreq(ns) == 'M') then ! monthly
           write(ncfile,'(a,a,i4.4,a,i2.2,a,a)')  &
            history_file(1:lenstr(history_file)),'.', &
             iyear,'-',imonth,'.',suffix
          elseif (histfreq(ns) == 'y'.or.histfreq(ns) == 'Y') then ! yearly
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

            num_avail_hist_fields_tot = num_avail_hist_fields_tot + 1

            if (vcoord(11:14) == 'time') then
               num_avail_hist_fields_2D  = num_avail_hist_fields_2D + 1
            elseif (vcoord(11:14) == 'NCAT' .and. vcoord(16:19) == 'time') then
               num_avail_hist_fields_3Dc = num_avail_hist_fields_3Dc + 1
            elseif (vcoord(11:15) == 'VGRDi' .and. vcoord(17:20) == 'time') then
               num_avail_hist_fields_3Dz = num_avail_hist_fields_3Dz + 1
            elseif (vcoord(11:15) == 'VGRDb' .and. vcoord(17:20) == 'time') then
               num_avail_hist_fields_3Db = num_avail_hist_fields_3Db + 1
            elseif (vcoord(11:15) == 'VGRDi' .and. vcoord(17:20) == 'NCAT') then
               num_avail_hist_fields_4Di = num_avail_hist_fields_4Di + 1
            elseif (vcoord(11:15) == 'VGRDs' .and. vcoord(17:20) == 'NCAT') then
               num_avail_hist_fields_4Ds = num_avail_hist_fields_4Ds + 1
            elseif (vcoord(11:15) == 'VGRDb' .and. vcoord(17:20) == 'NCAT') then
               num_avail_hist_fields_4Db = num_avail_hist_fields_4Db + 1
            endif

            if (num_avail_hist_fields_tot > max_avail_hist_fields) &
               call abort_ice("Need to increase max_avail_hist_fields")

            if (num_avail_hist_fields_tot /= &
                num_avail_hist_fields_2D  + &
                num_avail_hist_fields_3Dc + &
                num_avail_hist_fields_3Dz + &
                num_avail_hist_fields_3Db + &
                num_avail_hist_fields_4Di + &
                num_avail_hist_fields_4Ds + &
                num_avail_hist_fields_4Db)  &
               call abort_ice("num_avail_hist_fields error")

            id(ns) = num_avail_hist_fields_tot

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

      subroutine accum_hist_field_2D(id, iblk, field_accum, field)

!     !DESCRIPTION:
!     Accumulates a history field
!
!     !REVISION HISTORY:
!     2009 Created by D. Bailey following POP
!     2010 Generalized dimension of variables by N. Jeffery, E. Hunke

      use ice_domain
      use ice_grid, only: tmask
      use ice_calendar, only: nstreams

!     !OUTPUT PARAMETERS:

      integer (int_kind), dimension(max_nstrm), intent(in) :: &
         id                ! location in avail_fields array for use in
                           ! later routines
        
      integer (kind=int_kind), intent(in) :: iblk

      real (kind=dbl_kind), intent(in) :: &
         field_accum(:,:)

      real (kind=dbl_kind), intent(inout) :: &
         field(:,:,:,:)

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
             field(i,j,idns, iblk) = field(i,j,idns, iblk) + field_accum(i,j)
          endif
       enddo
       enddo

       endif
       enddo

      end subroutine accum_hist_field_2D

!=======================================================================

      subroutine accum_hist_field_3D(id, iblk, ndim, field_accum, field)

!     !DESCRIPTION:
!     Accumulates a history field
!
!     !REVISION HISTORY:
!     2009 Created by D. Bailey following POP
!     2010 Generalized dimension of variables by N. Jeffery, E. Hunke

      use ice_domain
      use ice_grid, only: tmask
      use ice_calendar, only: nstreams

!     !OUTPUT PARAMETERS:

      integer (int_kind), dimension(max_nstrm), intent(in) :: &
         id                ! location in avail_fields array for use in
                           ! later routines
        
      integer (kind=int_kind), intent(in) :: iblk

      integer (kind=int_kind), intent(in) :: &
         ndim              ! third dimension size

      real (kind=dbl_kind), intent(in) :: &
         field_accum(:,:,:)

      real (kind=dbl_kind), intent(inout) :: &
         field(:,:,:,:,:)

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: i,j,k, ilo, ihi, jlo, jhi, ns, idns

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

       do k = 1, ndim
       do j = jlo, jhi
       do i = ilo, ihi
          if (tmask(i,j,iblk)) then
             field(i,j,k,idns,iblk) = field(i,j,k,idns,iblk) + field_accum(i,j,k)
          endif
       enddo
       enddo
       enddo

       endif
       enddo

      end subroutine accum_hist_field_3D

!=======================================================================

      subroutine accum_hist_field_4D(id, iblk, ndim3, ndim4, field_accum, field)

!     !DESCRIPTION:
!     Accumulates a history field
!
!     !REVISION HISTORY:
!     2009 Created by D. Bailey following POP
!     2010 Generalized dimension of variables by N. Jeffery, E. Hunke

      use ice_domain
      use ice_grid, only: tmask
      use ice_calendar, only: nstreams

!     !OUTPUT PARAMETERS:

      integer (int_kind), dimension(max_nstrm), intent(in) :: &
         id                ! location in avail_fields array for use in
                           ! later routines
        
      integer (kind=int_kind), intent(in) :: iblk

      integer (kind=int_kind), intent(in) :: &
         ndim3  , &        ! third dimension size
         ndim4             ! fourth dimension size

      real (kind=dbl_kind), intent(in) :: &
         field_accum(:,:,:,:)

      real (kind=dbl_kind), intent(inout) :: &
         field(:,:,:,:,:,:)

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: i,j,k,n,ilo, ihi, jlo, jhi, ns, idns

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

       do k = 1, ndim4
       do n = 1, ndim3
       do j = jlo, jhi
       do i = ilo, ihi
          if (tmask(i,j,iblk)) then
             field(i,j,n,k,idns,iblk) = field(i,j,n,k,idns,iblk) + field_accum(i,j,n,k)
          endif
       enddo
       enddo
       enddo
       enddo

       endif
       enddo

      end subroutine accum_hist_field_4D

!=======================================================================

      end module ice_history_shared

!=======================================================================
