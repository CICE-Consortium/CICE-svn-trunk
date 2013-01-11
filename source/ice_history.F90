!=======================================================================
!
!BOP
!
! !MODULE: ice_history - ice model history files
!
! Driver for core history output
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
      use ice_history_shared
      use ice_history_write
      use ice_history_mechred
      use ice_history_pond
!
!EOP
!
      implicit none
      save
      
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
      use ice_dyn_evp, only: kdyn
      use ice_flux, only: mlt_onset, frz_onset, albcnt
      use ice_restart, only: restart
      use ice_zbgc_public, only: tr_bgc_N_sk, tr_bgc_C_sk, tr_bgc_chl_sk, &
                         tr_bgc_Nit_sk, tr_bgc_Am_sk, tr_bgc_Sil_sk, &
                         tr_bgc_DMSPp_sk, tr_bgc_DMSPd_sk, tr_bgc_DMS_sk, &
                         tr_bgc_NO, tr_bgc_NH, tr_bgc_N, tr_bgc_C, &
                         tr_bgc_chl, tr_bgc_Sil, tr_bgc_DMSPp, tr_bgc_DMSPd,&
                         tr_bgc_DMS, tr_bgc_PON,  &
                         nlt_bgc_N, &
                         nlt_bgc_NO, nlt_bgc_C, nlt_bgc_chl, &
                         nlt_bgc_NH, nlt_bgc_Sil, nlt_bgc_DMSPp, &
                         nlt_bgc_DMSPd, nlt_bgc_DMS, nlt_bgc_PON , tr_bgc_S
      use ice_state, only: tr_iage, tr_FY, tr_lvl, tr_pond, tr_aero, hbrine
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
      if (.not. tr_FY)   f_FY   = 'x'
      if (.not. tr_aero) then
         f_faero_atm = 'x'
         f_faero_ocn = 'x'
         f_aero      = 'x' 
         f_aeron     = 'x' ! NOTE not implemented
      endif

      if (.not. tr_bgc_N_sk)  f_bgc_N_sk = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_C_sk = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_chl_sk = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_Nit_sk = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_Am_sk = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_Sil_sk = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_DMSPp_sk = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_DMSPd_sk = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_DMS_sk = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_Nit_ml = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_Am_ml = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_Sil_ml = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_DMSP_ml = 'x'
      if (.not. tr_bgc_N_sk)  f_bgc_DMS_ml = 'x'
      if (.not. tr_bgc_NO)    f_bgc_NO  = 'x'
      if (.not. tr_bgc_NH)    f_bgc_NH  = 'x'
      if (.not. tr_bgc_N)     f_bgc_N  = 'x'
      if (.not. tr_bgc_C)     f_bgc_C  = 'x'
      if (.not. tr_bgc_chl)   f_bgc_chl  = 'x'
      if (.not. tr_bgc_Sil)   f_bgc_Sil  = 'x'
      if (.not. tr_bgc_DMSPp) f_bgc_DMSPp  = 'x'
      if (.not. tr_bgc_DMSPd) f_bgc_DMSPd  = 'x'
      if (.not. tr_bgc_DMS)   f_bgc_DMS  = 'x'
      if (.not. tr_bgc_PON)   f_bgc_PON  = 'x'
      if (.not. tr_bgc_S)     f_bgc_S  = 'x'
      if (.not. hbrine)       f_fbri = 'x'

      if (kdyn /= 2) then
           f_a11       = 'x'
           f_a12       = 'x'
           f_e11       = 'x'
           f_e12       = 'x'
           f_e22       = 'x'
           f_s11       = 'x'
           f_s12       = 'x'
           f_s22       = 'x'
           f_yieldstress11 = 'x'
           f_yieldstress12 = 'x'
           f_yieldstress22 = 'x'
      endif

      ! these must be output at the same frequency because of 
      ! cos(zenith angle) averaging
      if (f_albice(1:1) /= 'x' .and. f_albsni(1:1) /= 'x') f_albice = f_albsni
      if (f_albsno(1:1) /= 'x') f_albsno = f_albice
      if (f_albpnd(1:1) /= 'x') f_albpnd = f_albice
      if (f_coszen(1:1) /= 'x' .and. f_albice(1:1) /= 'x') f_coszen = f_albice
      if (f_coszen(1:1) /= 'x' .and. f_albsni(1:1) /= 'x') f_coszen = f_albsni

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
      call broadcast_scalar (f_NCAT, master_task)
      call broadcast_scalar (f_VGRDi, master_task)
      call broadcast_scalar (f_VGRDs, master_task)
      call broadcast_scalar (f_VGRDb, master_task)

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
      call broadcast_scalar (f_faero_atm, master_task)
      call broadcast_scalar (f_faero_ocn, master_task)
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
      call broadcast_scalar (f_dsnow, master_task)
      call broadcast_scalar (f_meltt, master_task)
      call broadcast_scalar (f_melts, master_task)
      call broadcast_scalar (f_meltb, master_task)
      call broadcast_scalar (f_meltl, master_task)
      call broadcast_scalar (f_fresh, master_task)
      call broadcast_scalar (f_fresh_ai, master_task)
      call broadcast_scalar (f_fsalt, master_task)
      call broadcast_scalar (f_fsalt_ai, master_task)
      call broadcast_scalar (f_fsice, master_task)
      call broadcast_scalar (f_fsice_ai, master_task)
      call broadcast_scalar (f_fsice_g, master_task)
      call broadcast_scalar (f_fsice_g_ai, master_task)
      call broadcast_scalar (f_fNO, master_task)
      call broadcast_scalar (f_fNO_ai, master_task)
      call broadcast_scalar (f_fNO_g, master_task)
      call broadcast_scalar (f_fNO_g_ai, master_task)
      call broadcast_scalar (f_fNH, master_task)
      call broadcast_scalar (f_fNH_ai, master_task)
      call broadcast_scalar (f_fNH_g, master_task)
      call broadcast_scalar (f_fNH_g_ai, master_task)
      call broadcast_scalar (f_fN,  master_task)
      call broadcast_scalar (f_fN_ai, master_task)
      call broadcast_scalar (f_fN_g, master_task)
      call broadcast_scalar (f_fN_g_ai, master_task)
      call broadcast_scalar (f_fSil, master_task)
      call broadcast_scalar (f_fSil_ai, master_task)
      call broadcast_scalar (f_fSil_g, master_task)
      call broadcast_scalar (f_fSil_g_ai, master_task)
      call broadcast_scalar (f_fsicen_g, master_task)
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

!      call broadcast_scalar (f_field3dz, master_task)
      call broadcast_scalar (f_Tinz, master_task)
      call broadcast_scalar (f_Sinz, master_task)
      call broadcast_scalar (f_Tsnz, master_task)

      call broadcast_scalar (f_bgc_N_sk, master_task)
      call broadcast_scalar (f_bgc_C_sk, master_task)
      call broadcast_scalar (f_bgc_chl_sk, master_task)
      call broadcast_scalar (f_bgc_Nit_sk, master_task)
      call broadcast_scalar (f_bgc_Am_sk, master_task)
      call broadcast_scalar (f_bgc_Sil_sk, master_task)
      call broadcast_scalar (f_bgc_DMSPp_sk, master_task)
      call broadcast_scalar (f_bgc_DMSPd_sk, master_task)
      call broadcast_scalar (f_bgc_DMS_sk, master_task)
      call broadcast_scalar (f_bgc_Nit_ml, master_task)
      call broadcast_scalar (f_bgc_Am_ml, master_task)
      call broadcast_scalar (f_bgc_Sil_ml, master_task)
      call broadcast_scalar (f_bgc_DMSP_ml, master_task)
      call broadcast_scalar (f_bgc_DMS_ml, master_task)

      call broadcast_scalar (f_zTin, master_task)
      call broadcast_scalar (f_zphi, master_task)
      call broadcast_scalar (f_iDi, master_task)
      call broadcast_scalar (f_iki, master_task)
      call broadcast_scalar (f_bgc_NO, master_task)
      call broadcast_scalar (f_bgc_NH, master_task)
      call broadcast_scalar (f_bgc_N, master_task)
      call broadcast_scalar (f_bgc_C, master_task)
      call broadcast_scalar (f_bgc_chl, master_task)
      call broadcast_scalar (f_bgc_Sil, master_task)
      call broadcast_scalar (f_bgc_DMSPp, master_task)
      call broadcast_scalar (f_bgc_DMSPd, master_task)
      call broadcast_scalar (f_bgc_DMS, master_task)
      call broadcast_scalar (f_bgc_PON, master_task)
      call broadcast_scalar (f_bgc_S, master_task)
      call broadcast_scalar (f_fbri, master_task)
      call broadcast_scalar (f_growN, master_task)
      call broadcast_scalar (f_zfswin, master_task)
      call broadcast_scalar (f_Stot, master_task)
      call broadcast_scalar (f_chlnet, master_task)
      call broadcast_scalar (f_PPnet, master_task)
      call broadcast_scalar (f_NOnet, master_task)
!      call broadcast_scalar (f_Tsf_ice, master_task)
      call broadcast_scalar (f_upNO, master_task)
      call broadcast_scalar (f_upNH, master_task)
      call broadcast_scalar (f_iage, master_task)
      call broadcast_scalar (f_FY, master_task)
      call broadcast_scalar (f_aero, master_task)
      call broadcast_scalar (f_aeron, master_task)

      call broadcast_scalar (f_a11, master_task)
      call broadcast_scalar (f_a12, master_task)
      call broadcast_scalar (f_e11, master_task)
      call broadcast_scalar (f_e12, master_task)
      call broadcast_scalar (f_e22, master_task)
      call broadcast_scalar (f_s11, master_task)
      call broadcast_scalar (f_s12, master_task) 
      call broadcast_scalar (f_s22, master_task)
      call broadcast_scalar (f_yieldstress11, master_task)
      call broadcast_scalar (f_yieldstress12, master_task)
      call broadcast_scalar (f_yieldstress22, master_task)


      ! 2D variables
      do ns1 = 1, nstreams

!!!!! begin example
!      if (f_example(1:1) /= 'x') &
!         call define_hist_field(n_example,"example","m",tstr2D, tcstr, & 
!            "example: mean ice thickness",                           &
!            "ice volume per unit grid cell area", c1, c0,            &
!            ns1, f_example)
!!!!! end example

      if (f_hi(1:1) /= 'x') &
         call define_hist_field(n_hi,"hi","m",tstr2D, tcstr,        & 
            "grid cell mean ice thickness",                       &
            "ice volume per unit grid cell area", c1, c0,         &
            ns1, f_hi)

      if (f_hs(1:1) /= 'x') &
         call define_hist_field(n_hs,"hs","m",tstr2D, tcstr,        &
             "grid cell mean snow thickness",                     &
             "snow volume per unit grid cell area", c1, c0,       &
             ns1, f_hs)

      
      if (f_Tsfc(1:1) /= 'x') &
         call define_hist_field(n_Tsfc,"Tsfc","C",tstr2D, tcstr,    &
             "snow/ice surface temperature",                      &
             "averaged with Tf if no ice is present", c1, c0,     &
             ns1, f_Tsfc)
      
      if (f_aice(1:1) /= 'x') &
         call define_hist_field(n_aice,"aice","1",tstr2D, tcstr,    &
             "ice area  (aggregate)",                             &
             "none", c1, c0,                                      &
             ns1, f_aice)
      
      if (f_uvel(1:1) /= 'x') &
         call define_hist_field(n_uvel,"uvel","m/s",ustr2D, ucstr,  &
             "ice velocity (x)",                                  &
             "positive is x direction on U grid", c1, c0,         &
             ns1, f_uvel)
      
      if (f_vvel(1:1) /= 'x') &
         call define_hist_field(n_vvel,"vvel","m/s",ustr2D, ucstr,  &
             "ice velocity (y)",                                  &
             "positive is y direction on U grid", c1, c0,         &
             ns1, f_vvel)
      
      if (f_fswdn(1:1) /= 'x') &
         call define_hist_field(n_fswdn,"fswdn","W/m^2",tstr2D, tcstr, &
             "down solar flux",                                      &
             "positive downward", c1, c0,                            &
             ns1, f_fswdn)
      
      if (f_flwdn(1:1) /= 'x') &
         call define_hist_field(n_flwdn,"flwdn","W/m^2",tstr2D, tcstr, &
             "down longwave flux",                                   &
             "positive downward", c1, c0,                            &
             ns1, f_flwdn)
      
      if (f_snow(1:1) /= 'x') &
         call define_hist_field(n_snow,"snow","cm/day",tstr2D, tcstr, &
             "snowfall rate (cpl)",                                 &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_snow)
      
      if (f_snow_ai(1:1) /= 'x') &
         call define_hist_field(n_snow_ai,"snow_ai","cm/day",tstr2D, tcstr, &
             "snowfall rate",                                             &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_snow_ai)
      
      if (f_rain(1:1) /= 'x') &
         call define_hist_field(n_rain,"rain","cm/day",tstr2D, tcstr, &
             "rainfall rate (cpl)",                                 &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_rain)
      
      if (f_rain_ai(1:1) /= 'x') &
         call define_hist_field(n_rain_ai,"rain_ai","cm/day",tstr2D, tcstr, &
             "rainfall rate",                                             &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_rain_ai)
      
      if (f_sst(1:1) /= 'x') &
         call define_hist_field(n_sst,"sst","C",tstr2D, tcstr, &
             "sea surface temperature",                      &
             "none", c1, c0,                                 &
             ns1, f_sst)
      
      if (f_sss(1:1) /= 'x') &
         call define_hist_field(n_sss,"sss","ppt",tstr2D, tcstr, &
             "sea surface salinity",                           &
             "none", c1, c0,                                   &
             ns1, f_sss)
      
      if (f_uocn(1:1) /= 'x') &
         call define_hist_field(n_uocn,"uocn","m/s",ustr2D, ucstr, &
             "ocean current (x)",                                &
             "positive is x direction on U grid", c1, c0,        &
             ns1, f_uocn)
      
      if (f_vocn(1:1) /= 'x') &
         call define_hist_field(n_vocn,"vocn","m/s",ustr2D, ucstr, &
             "ocean current (y)",                                &
             "positive is y direction on U grid", c1, c0,        &
             ns1, f_vocn)
      
      if (f_frzmlt(1:1) /= 'x') &
         call define_hist_field(n_frzmlt,"frzmlt","W/m^2",tstr2D, tcstr, &
             "freeze/melt potential",                                  &
             "if >0, new ice forms; if <0, ice melts", c1, c0,         &
             ns1, f_frzmlt)
      
      if (f_fswfac(1:1) /= 'x') &
         call define_hist_field(n_fswfac,"fswfac","1",tstr2D, tcstr, &
             "shortwave scaling factor",                           &
             "ratio of netsw new:old", c1, c0,                     &
             ns1, f_fswfac)
      
      if (f_fswabs(1:1) /= 'x') &
         call define_hist_field(n_fswabs,"fswabs","W/m^2",tstr2D, tcstr, &
             "snow/ice/ocn absorbed solar flux (cpl)",                 &
             "positive downward", c1, c0,                              &
             ns1, f_fswabs)
      
      if (f_fswabs_ai(1:1) /= 'x') &
         call define_hist_field(n_fswabs_ai,"fswabs_ai","W/m^2",tstr2D, tcstr, &
             "snow/ice/ocn absorbed solar flux",                             &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_fswabs_ai)
      
!      if (f_albsni(1:1) /= 'x') &
!         call define_hist_field(n_albsni,"albsni","%",tstr2D, tcstr, &
!             "snow/ice broad band albedo",                         &
!             "scaled (divided) by aice", c100, c0,                 &
!             ns1, f_albsni)
      
      if (f_albsni(1:1) /= 'x') &
         call define_hist_field(n_albsni,"albsni","%",tstr2D, tcstr, &
             "snow/ice broad band albedo",                         &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albsni)
      
      if (f_alvdr(1:1) /= 'x') &
         call define_hist_field(n_alvdr,"alvdr","%",tstr2D, tcstr, &
             "visible direct albedo",                            &
             "scaled (divided) by aice", c100, c0,               &
             ns1, f_alvdr)
      
      if (f_alidr(1:1) /= 'x') &
         call define_hist_field(n_alidr,"alidr","%",tstr2D, tcstr, &
             "near IR direct albedo",                            &
             "scaled (divided) by aice", c100, c0,               &
             ns1, f_alidr)

      if (f_albice(1:1) /= 'x') &
         call define_hist_field(n_albice,"albice","%",tstr2D, tcstr, &
             "bare ice albedo",                                    &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albice)
      
      if (f_albsno(1:1) /= 'x') &
         call define_hist_field(n_albsno,"albsno","%",tstr2D, tcstr, &
             "snow albedo",                                        &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albsno)
      
      if (f_albpnd(1:1) /= 'x') &
         call define_hist_field(n_albpnd,"albpnd","%",tstr2D, tcstr, &
             "melt pond albedo",                                   &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albpnd)
      
      if (f_coszen(1:1) /= 'x') &
         call define_hist_field(n_coszen,"coszen","radian",tstr2D, tcstr, &
             "cosine of the zenith angle",                              &
             "negative below horizon", c1, c0,                          &
             ns1, f_coszen)

      if (f_flat(1:1) /= 'x') &
         call define_hist_field(n_flat,"flat","W/m^2",tstr2D, tcstr, &
             "latent heat flux (cpl)",                             &
             "positive downward", c1, c0,                          &
             ns1, f_flat)
      
      if (f_flat_ai(1:1) /= 'x') &
         call define_hist_field(n_flat_ai,"flat_ai","W/m^2",tstr2D, tcstr, &
             "latent heat flux",                                         &
             "weighted by ice area", c1, c0,                             &
             ns1, f_flat_ai)
      
      if (f_fsens(1:1) /= 'x') &
         call define_hist_field(n_fsens,"fsens","W/m^2",tstr2D, tcstr, &
             "sensible heat flux (cpl)",                             &
             "positive downward", c1, c0,                            &
             ns1, f_fsens)
      
      if (f_fsens_ai(1:1) /= 'x') &
         call define_hist_field(n_fsens_ai,"fsens_ai","W/m^2",tstr2D, tcstr, &
             "sensible heat flux",                                         &
             "weighted by ice area", c1, c0,                               &
             ns1, f_fsens_ai)
      
      if (f_flwup(1:1) /= 'x') &
         call define_hist_field(n_flwup,"flwup","W/m^2",tstr2D, tcstr, &
             "upward longwave flux (cpl)",                           &
             "positive downward", c1, c0,                            &
             ns1, f_flwup)
      
      if (f_flwup_ai(1:1) /= 'x') &
         call define_hist_field(n_flwup_ai,"flwup_ai","W/m^2",tstr2D, tcstr, &
             "upward longwave flux",                                       &
             "weighted by ice area", c1, c0,                               &
             ns1, f_flwup_ai)
      
      if (f_evap(1:1) /= 'x') &
         call define_hist_field(n_evap,"evap","cm/day",tstr2D, tcstr, &
             "evaporative water flux (cpl)",                        &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_evap)
      
      if (f_evap_ai(1:1) /= 'x') &
         call define_hist_field(n_evap_ai,"evap_ai","cm/day",tstr2D, tcstr, &
             "evaporative water flux",                                    &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_evap_ai)
      
      if (f_Tair(1:1) /= 'x') &
         call define_hist_field(n_Tair,"Tair","C",tstr2D, tcstr, &
             "air temperature",                                &
             "none", c1, -Tffresh,                             &
             ns1, f_Tair)
      
      if (f_Tref(1:1) /= 'x') &
         call define_hist_field(n_Tref,"Tref","C",tstr2D, tcstr, &
             "2m reference temperature",                       &
             "none", c1, -Tffresh,                             &
             ns1, f_Tref)
      
      if (f_Qref(1:1) /= 'x') &
         call define_hist_field(n_Qref,"Qref","g/kg",tstr2D, tcstr, &
             "2m reference specific humidity",                    &
             "none", kg_to_g, c0,                                 &
             ns1, f_Qref)
      
      if (f_congel(1:1) /= 'x') &
         call define_hist_field(n_congel,"congel","cm/day",tstr2D, tcstr, &
             "congelation ice growth",                                  &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_congel)
      
      if (f_frazil(1:1) /= 'x') &
         call define_hist_field(n_frazil,"frazil","cm/day",tstr2D, tcstr, &
             "frazil ice growth",                                       &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_frazil)
      
      if (f_snoice(1:1) /= 'x') &
         call define_hist_field(n_snoice,"snoice","cm/day",tstr2D, tcstr, &
             "snow-ice formation",                                      &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_snoice)
           
      if (f_dsnow(1:1) /= 'x') &
         call define_hist_field(n_dsnow,"dsnow","cm/day",tstr2D, tcstr, &
             "snow formation",                                      &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_dsnow)
      
      if (f_meltt(1:1) /= 'x') &
         call define_hist_field(n_meltt,"meltt","cm/day",tstr2D, tcstr, &
             "top ice melt",                                          &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltt)
      
      if (f_melts(1:1) /= 'x') &
         call define_hist_field(n_melts,"melts","cm/day",tstr2D, tcstr, &
             "top snow melt",                                          &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_melts)
      
      if (f_meltb(1:1) /= 'x') &
         call define_hist_field(n_meltb,"meltb","cm/day",tstr2D, tcstr, &
             "basal ice melt",                                        &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltb)
      
      if (f_meltl(1:1) /= 'x') &
         call define_hist_field(n_meltl,"meltl","cm/day",tstr2D, tcstr, &
             "lateral ice melt",                                      &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltl)
      
      if (f_fresh(1:1) /= 'x') &
         call define_hist_field(n_fresh,"fresh","cm/day",tstr2D, tcstr,   &
             "freshwtr flx ice to ocn (cpl)",                           &
             "if positive, ocean gains fresh water",                    &
             mps_to_cmpdy/rhofresh, c0,                                 &
             ns1, f_fresh)
      
      if (f_fresh_ai(1:1) /= 'x') &
         call define_hist_field(n_fresh_ai,"fresh_ai","cm/day",tstr2D, tcstr, &
             "freshwtr flx ice to ocn",                                     &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,             &
             ns1, f_fresh_ai)
      
      if (f_fsalt(1:1) /= 'x') &
         call define_hist_field(n_fsalt,"fsalt","kg/m^2/s",tstr2D, tcstr, &
             "salt flux ice to ocn (cpl)",                              &
             "if positive, ocean gains salt", c1, c0,                   &
             ns1, f_fsalt)
      
      if (f_fsalt_ai(1:1) /= 'x') &
         call define_hist_field(n_fsalt_ai,"fsalt_ai","kg/m^2/s",tstr2D, tcstr, &
             "salt flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fsalt_ai)
      
      if (f_fsice(1:1) /= 'x') &
         call define_hist_field(n_fsice,"fsice","kg/m^2/s",tstr2D, tcstr, &
             "prognostic salt flux ice to ocn (cpl)",                              &
             "if positive, ocean gains salt", c1, c0,                   &
             ns1, f_fsice)
      
      if (f_fsice_ai(1:1) /= 'x') &
         call define_hist_field(n_fsice_ai,"fsice_ai","kg/m^2/s",tstr2D, tcstr, &
             "prognostic salt flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fsice_ai)
      
      if (f_fsice_g(1:1) /= 'x') &
         call define_hist_field(n_fsice_g,"fsice_g","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage salt flux ice to ocn (cpl)",                              &
             "if positive, ocean gains salt", c1, c0,                   &
             ns1, f_fsice_g)
      
      if (f_fsice_g_ai(1:1) /= 'x') &
         call define_hist_field(n_fsice_g_ai,"fsice_g_ai","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage salt flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fsice_g_ai)
      
      if (f_fNO(1:1) /= 'x') &
         call define_hist_field(n_fNO,"fNO","mmol/m^2/s",tstr2D, tcstr, &
             "nitrate flux ice to ocn (cpl)",                              &
             "if positive, ocean gains nitrate", c1, c0,                   &
             ns1, f_fNO)
      
      if (f_fNO_ai(1:1) /= 'x') &
         call define_hist_field(n_fNO_ai,"fNO_ai","mmol/m^2/s",tstr2D, tcstr, &
             "nitrate flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fNO_ai)
      
      if (f_fNO_g(1:1) /= 'x') &
         call define_hist_field(n_fNO_g,"fNO_g","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage nitrate flux ice to ocn (cpl)",                              &
             "if positive, ocean gains nitrate", c1, c0,                   &
             ns1, f_fNO_g)
      
      if (f_fNO_g_ai(1:1) /= 'x') &
         call define_hist_field(n_fNO_g_ai,"fNO_g_ai","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage nitrate flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fNO_g_ai)
            
      if (f_fNH(1:1) /= 'x') &
         call define_hist_field(n_fNH,"fNH","mmol/m^2/s",tstr2D, tcstr, &
             "ammonium flux ice to ocn (cpl)",                              &
             "if positive, ocean gains ammonium", c1, c0,                   &
             ns1, f_fNH)
      
      if (f_fNH_ai(1:1) /= 'x') &
         call define_hist_field(n_fNH_ai,"fNH_ai","mmol/m^2/s",tstr2D, tcstr, &
             "ammonium flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fNH_ai)
      
      if (f_fNH_g(1:1) /= 'x') &
         call define_hist_field(n_fNH_g,"fNH_g","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage ammonium flux ice to ocn (cpl)",                              &
             "if positive, ocean gains ammonium", c1, c0,                   &
             ns1, f_fNH_g)
      
      if (f_fNH_g_ai(1:1) /= 'x') &
         call define_hist_field(n_fNH_g_ai,"fNH_g_ai","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage ammonium flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fNH_g_ai)
                    
      if (f_fN(1:1) /= 'x') &
         call define_hist_field(n_fN,"fN","mmol/m^2/s",tstr2D, tcstr, &
             "algal N flux ice to ocn (cpl)",                              &
             "if positive, ocean gains algal N", c1, c0,                   &
             ns1, f_fN)
      
      if (f_fN_ai(1:1) /= 'x') &
         call define_hist_field(n_fN_ai,"fN_ai","mmol/m^2/s",tstr2D, tcstr, &
             "algal N flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fN_ai)
      
      if (f_fN_g(1:1) /= 'x') &
         call define_hist_field(n_fN_g,"fN_g","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage algal N flux ice to ocn (cpl)",                              &
             "if positive, ocean gains algal N", c1, c0,                   &
             ns1, f_fN_g)
      
      if (f_fN_g_ai(1:1) /= 'x') &
         call define_hist_field(n_fN_g_ai,"fN_g_ai","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage algal N flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fN_g_ai)
              
      if (f_fSil(1:1) /= 'x') &
         call define_hist_field(n_fSil,"fSil","mmol/m^2/s",tstr2D, tcstr, &
             "silicate flux ice to ocn (cpl)",                              &
             "if positive, ocean gains silicate", c1, c0,                   &
             ns1, f_fSil)
      
      if (f_fSil_ai(1:1) /= 'x') &
         call define_hist_field(n_fSil_ai,"fSil_ai","mmol/m^2/s",tstr2D, tcstr, &
             "silicate flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fSil_ai)
      
      if (f_fSil_g(1:1) /= 'x') &
         call define_hist_field(n_fSil_g,"fSil_g","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage silicate flux ice to ocn (cpl)",                              &
             "if positive, ocean gains silicate", c1, c0,                   &
             ns1, f_fSil_g)
      
      if (f_fSil_g_ai(1:1) /= 'x') &
         call define_hist_field(n_fSil_g_ai,"fSil_g_ai","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage silicate flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fSil_g_ai)
            
      if (f_Stot(1:1) /= 'x') &
         call define_hist_field(n_Stot,"S_tot","g/m^2",tstr2D, tcstr,        &
             "Total Salt content",                     &
             "In ice volume*fbri", c1, c0,       &
             ns1, f_Stot)
     
      if (f_chlnet(1:1) /= 'x') &
         call define_hist_field(n_chlnet,"chl_net","mg chl/m^2",tstr2D, tcstr,        &
             "Net Chlorophyll",                     &
             "In ice volume*fbri ", c1, c0,       &
             ns1, f_chlnet)

      if (f_PPnet(1:1) /= 'x') &
         call define_hist_field(n_PPnet,"PP_net","mg C/d/m^2",tstr2D, tcstr,        &
             "Net Primary Production",                     &
             "In ice volume*fbri ", secday, c0,       &
             ns1, f_PPnet)

      if (f_NOnet(1:1) /= 'x') &
         call define_hist_field(n_NOnet,"NO_net","mmol NO/m^2",tstr2D, tcstr,        &
             "Net Nitrate",                     &
             "In ice volume*fbri ", c1, c0,       &
             ns1, f_NOnet)

!      if (f_Tsf_ice(1:1) /= 'x') &
!         call define_hist_field(n_Tsf_ice,"Tsf_ice","C",tstr2D, tcstr,        &
!             "Ice surface temperature",                     &
!             "Averaged over categories ", c1, c0,       &
!             ns1, f_Tsf_ice)

      if (f_fhocn(1:1) /= 'x') &
         call define_hist_field(n_fhocn,"fhocn","W/m^2",tstr2D, tcstr, &
             "heat flux ice to ocn (cpl)",                           &
             "if positive, ocean gains heat", c1, c0,                &
             ns1, f_fhocn)
      
      if (f_fhocn_ai(1:1) /= 'x') &
         call define_hist_field(n_fhocn_ai,"fhocn_ai","W/m^2",tstr2D, tcstr, &
             "heat flux ice to ocean",                                     &
             "weighted by ice area", c1, c0,                               &
             ns1, f_fhocn_ai)
      
      if (f_fswthru(1:1) /= 'x') &
         call define_hist_field(n_fswthru,"fswthru","W/m^2",tstr2D, tcstr, &
             "SW thru ice to ocean (cpl)",                               &
             "if positive, ocean gains heat", c1, c0,                    &
             ns1, f_fswthru)
      
      if (f_fswthru_ai(1:1) /= 'x') &
         call define_hist_field(n_fswthru_ai,"fswthru_ai","W/m^2",tstr2D, tcstr,&
             "SW flux thru ice to ocean",                                     &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fswthru_ai)
      
      if (f_strairx(1:1) /= 'x') &
         call define_hist_field(n_strairx,"strairx","N/m^2",ustr2D, ucstr, &
             "atm/ice stress (x)",                                       &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strairx)
      
      if (f_strairy(1:1) /= 'x') &
         call define_hist_field(n_strairy,"strairy","N/m^2",ustr2D, ucstr, &
             "atm/ice stress (y)",                                       &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strairy)
      
      if (f_strtltx(1:1) /= 'x') &
         call define_hist_field(n_strtltx,"strtltx","N/m^2",ustr2D, ucstr, &
             "sea sfc tilt stress (x)",                                  &
             "none", c1, c0,                                             &
             ns1, f_strtltx)
      
      if (f_strtlty(1:1) /= 'x') &
         call define_hist_field(n_strtlty,"strtlty","N/m^2",ustr2D, ucstr, &
             "sea sfc tilt stress (y)",                                  &
             "none", c1, c0,                                             &
             ns1, f_strtlty)
      
      if (f_strcorx(1:1) /= 'x') &
         call define_hist_field(n_strcorx,"strcorx","N/m^2",ustr2D, ucstr, &
             "coriolis stress (x)",                                      &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strcorx)
      
      if (f_strcory(1:1) /= 'x') &
         call define_hist_field(n_strcory,"strcory","N/m^2",ustr2D, ucstr, &
             "coriolis stress (y)",                                      &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strcory)
      
      if (f_strocnx(1:1) /= 'x') &
         call define_hist_field(n_strocnx,"strocnx","N/m^2",ustr2D, ucstr, &
             "ocean/ice stress (x)",                                     &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strocnx)
      
      if (f_strocny(1:1) /= 'x') &
         call define_hist_field(n_strocny,"strocny","N/m^2",ustr2D, ucstr, &
             "ocean/ice stress (y)",                                     &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strocny)
      
      if (f_strintx(1:1) /= 'x') &
         call define_hist_field(n_strintx,"strintx","N/m^2",ustr2D, ucstr, &
             "internal ice stress (x)",                                  &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strintx)
      
      if (f_strinty(1:1) /= 'x') &
         call define_hist_field(n_strinty,"strinty","N/m^2",ustr2D, ucstr, &
             "internal ice stress (y)",                                  &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strinty)
      
      if (f_strength(1:1) /= 'x') &
         call define_hist_field(n_strength,"strength","N/m",tstr2D, tcstr, &
             "compressive ice strength",                                 &
             "none", c1, c0,                                             &
             ns1, f_strength)
      
      if (f_divu(1:1) /= 'x') &
         call define_hist_field(n_divu,"divu","%/day",tstr2D, tcstr, &
             "strain rate (divergence)",                           &
             "none", secday*c100, c0,                              &
             ns1, f_divu)
      
      if (f_shear(1:1) /= 'x') &
         call define_hist_field(n_shear,"shear","%/day",tstr2D, tcstr, &
             "strain rate (shear)",                                  &
             "none", secday*c100, c0,                                &
             ns1, f_shear)
      
      if (f_sig1(1:1) /= 'x') &
         call define_hist_field(n_sig1,"sig1","1",ustr2D, ucstr, &
             "norm. principal stress 1",                       &
             "sig1 is instantaneous", c1, c0,                  &
             ns1, f_sig1)
      
      if (f_sig2(1:1) /= 'x') &
         call define_hist_field(n_sig2,"sig2","1",ustr2D, ucstr, &
             "norm. principal stress 2",                       &
             "sig2 is instantaneous", c1, c0,                  &
             ns1, f_sig2)
      
      if (f_dvidtt(1:1) /= 'x') &
         call define_hist_field(n_dvidtt,"dvidtt","cm/day",tstr2D, tcstr, &
             "volume tendency thermo",                                  &
             "none", mps_to_cmpdy, c0,                                  &
             ns1, f_dvidtt)
      
      if (f_dvidtd(1:1) /= 'x') &
         call define_hist_field(n_dvidtd,"dvidtd","cm/day",tstr2D, tcstr, &
             "volume tendency dynamics",                                &
             "none", mps_to_cmpdy, c0,                                  &
             ns1, f_dvidtd)
      
      if (f_daidtt(1:1) /= 'x') &
         call define_hist_field(n_daidtt,"daidtt","%/day",tstr2D, tcstr, &
             "area tendency thermo",                                   &
             "none", secday*c100, c0,                                  &
             ns1, f_daidtt)
      
      if (f_daidtd(1:1) /= 'x') &
         call define_hist_field(n_daidtd,"daidtd","%/day",tstr2D, tcstr, &
             "area tendency dynamics",                                 &
             "none", secday*c100, c0,                                  &
             ns1, f_daidtd)
      
      if (f_mlt_onset(1:1) /= 'x') &
         call define_hist_field(n_mlt_onset,"mlt_onset","day of year", &
             tstr2D, tcstr,"melt onset date",                            &
             "midyear restart gives erroneous dates", c1, c0,          &
             ns1, f_mlt_onset)
      
      if (f_frz_onset(1:1) /= 'x') &
         call define_hist_field(n_frz_onset,"frz_onset","day of year", &
             tstr2D, tcstr,"freeze onset date",                          &
             "midyear restart gives erroneous dates", c1, c0,          &
             ns1, f_frz_onset)

      if (f_hisnap(1:1) /= 'x') &
         call define_hist_field(n_hisnap,"hisnap","m",tstr2D, tcstr, &
             "ice volume snapshot",                                &
             "none", c1, c0,                              &
             ns1, f_hisnap)
      
      if (f_aisnap(1:1) /= 'x') &
         call define_hist_field(n_aisnap,"aisnap","1",tstr2D, tcstr, &
             "ice area snapshot",                                  &
             "none", c1, c0,                              &
             ns1, f_aisnap)
      
      if (f_trsig(1:1) /= 'x') &
         call define_hist_field(n_trsig,"trsig","N/m^2",tstr2D, tcstr, &
             "internal stress tensor trace",                         &
             "ice strength approximation", c1, c0,                   &
             ns1, f_trsig)
      
      if (f_icepresent(1:1) /= 'x') &
         call define_hist_field(n_icepresent,"ice_present","1",tstr2D, tcstr, &
             "fraction of time-avg interval that ice is present",           &
             "ice extent flag", c1, c0,                                     &
             ns1, f_icepresent)
      
      if (f_fsurf_ai(1:1) /= 'x') &
         call define_hist_field(n_fsurf_ai,"fsurf_ai","W/m^2",tstr2D, tcstr, &
             "net surface heat flux",                                      &
             "positive downward, excludes conductive flux, weighted by ice area", &
             c1, c0, &
             ns1, f_fsurf_ai)
      
      if (f_fcondtop_ai(1:1) /= 'x') &
         call define_hist_field(n_fcondtop_ai,"fcondtop_ai","W/m^2", &
             tstr2D, tcstr,"top surface conductive heat flux",         &
             "positive downward, weighted by ice area", c1, c0,      &
             ns1, f_fcondtop_ai)
      
      if (f_fmeltt_ai(1:1) /= 'x') &
         call define_hist_field(n_fmeltt_ai,"fmeltt_ai","W/m^2",tstr2D, tcstr, &
             "net surface heat flux causing melt",                           &
             "always >= 0, weighted by ice area", c1, c0,                    &
             ns1, f_fmeltt_ai)
 
      if (f_a11(1:1) /= 'x') &
         call define_hist_field(n_a11,"a11"," ",tstr2D, tcstr, &
            "a11: component a11 of the structure tensor",                           &
            "none", c1, c0,            &
            ns1, f_a11)

      if (f_a12(1:1) /= 'x') &
         call define_hist_field(n_a12,"a12"," ",tstr2D, tcstr, &
            "a12: component a12 of the structure tensor",                           &
            "none", c1, c0,            &
            ns1, f_a12)

      if (f_e11(1:1) /= 'x') &
         call define_hist_field(n_e11,"e11","1/s",tstr2D, tcstr, &
            "e11: component e11 of the strain rate tensor",                           &
            "none", c1, c0,            &
            ns1, f_e11)

      if (f_e12(1:1) /= 'x') &
         call define_hist_field(n_e12,"e12","1/s",tstr2D, tcstr, &
            "e12: component e12 of the strain rate tensor",                           &
            "none", c1, c0,            &
            ns1, f_e12)

      if (f_e22(1:1) /= 'x') &
         call define_hist_field(n_e22,"e22","1/s",tstr2D, tcstr, &
            "e22: component e22 of the strain rate tensor",                           &
            "none", c1, c0,            &
            ns1, f_e22)

      if (f_s11(1:1) /= 'x') &
         call define_hist_field(n_s11,"s11","kg/s^2",tstr2D, tcstr, &
            "s11: component s11 of the stress tensor",                           &
            "none", c1, c0,            &
            ns1, f_s11)

      if (f_s12(1:1) /= 'x') &
         call define_hist_field(n_s12,"s12","kg/s^2",tstr2D, tcstr, &
            "s12: component s12 of the stress tensor",                           &
            "none", c1, c0,            &
            ns1, f_s12)

      if (f_s22(1:1) /= 'x') &
         call define_hist_field(n_s22,"s22","kg/s^2",tstr2D, tcstr, &
            "s22: component s12 of the stress tensor",                           &
            "none", c1, c0,            &
            ns1, f_s22)

      if (f_yieldstress11(1:1) /= 'x') &
         call define_hist_field(n_yieldstress11,"yieldstress11","kg/s^2",tstr2D, tcstr, &
            "yieldstress11: component 11 of the yieldstress tensor",                    &
            "none", c1, c0,            &
            ns1, f_yieldstress11)

      if (f_yieldstress12(1:1) /= 'x') &
         call define_hist_field(n_yieldstress12,"yieldstress12","kg/s^2",tstr2D, tcstr, &
            "yieldstress12: component 12 of the yieldstress tensor",                    &
            "none", c1, c0,            &
            ns1, f_yieldstress12)

      if (f_yieldstress22(1:1) /= 'x') &
         call define_hist_field(n_yieldstress22,"yieldstress22","kg/s^2",tstr2D, tcstr, &
            "yieldstress22: component 12 of the yieldstress tensor",                    &
            "none", c1, c0,            &
            ns1, f_yieldstress22)

      ! Tracers !!!!!!!!!!!

     !!! Biogeochemistry (skeletal layer tracers)
      if (f_bgc_N_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_N_sk,"algal_N","mmol/m^2",tstr2D, tcstr, &
             "ice bottom algae (nitrogen)",                                        &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns1, f_bgc_N_sk)
      if (f_bgc_C_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_C_sk,"algal_C","mmol/m^2",tstr2D, tcstr, &
             "ice bottom algae (carbon)",                                        &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns1, f_bgc_C_sk)
      if (f_bgc_chl_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_chl_sk,"algal_chl","mmol/m^2?",tstr2D, tcstr, &
             "ice bottom algae (chlorophyll)",                                        &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns1, f_bgc_chl_sk)
      if (f_bgc_Nit_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_Nit_sk,"skl_Nit","mmol/m^2",tstr2D, tcstr, &
             "skeletal nutrient (nitrate)",                                        &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns1, f_bgc_Nit_sk)
      if (f_bgc_Am_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_Am_sk,"skl_Am","mmol/m^2",tstr2D, tcstr, &
             "skeletal nutrient (ammonia/um)",                                        &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns1, f_bgc_Am_sk)
      if (f_bgc_Sil_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_Sil_sk,"skl_Sil","mmol/m^2",tstr2D, tcstr, &
             "skeletal nutrient (silicate)",                                        &
             "skelelal layer: bottom 2-3 cm", c1, c0,                &
             ns1, f_bgc_Sil_sk)
      if (f_bgc_Nit_ml(1:1) /= 'x') &
         call define_hist_field(n_bgc_Nit_ml,"ml_Nit","mmol/m^3",tstr2D, tcstr, &
             "mixed layer nutrient (nitrate)",                                        &
             "upper ocean", c1, c0,                &
             ns1, f_bgc_Nit_ml)
      if (f_bgc_Am_ml(1:1) /= 'x') &
         call define_hist_field(n_bgc_Am_ml,"ml_Am","mmol/m^3",tstr2D, tcstr, &
             "mixed layer nutrient (ammonia/um)",                                        &
             "upper ocean", c1, c0,                &
             ns1, f_bgc_Am_ml)
      if (f_bgc_Sil_ml(1:1) /= 'x') &
         call define_hist_field(n_bgc_Sil_ml,"ml_Sil","mmol/m^3",tstr2D, tcstr, &
             "mixed layer nutrient (silicate)",                                        &
             "upper ocean", c1, c0,                &
             ns1, f_bgc_Sil_ml)
      if (f_bgc_DMSPp_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_DMSPp_sk,"skl_DMSPp","mmol/m^2",tstr2D, tcstr, &
             "particulate S in algae (DMSPp)",                                        &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns1, f_bgc_DMSPp_sk)
      if (f_bgc_DMSPd_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_DMSPd_sk,"skl_DMSPd","mmol/m^2",tstr2D, tcstr, &
             "dissolved skl precursor (DSMPd)",                                        &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns1, f_bgc_DMSPd_sk)
      if (f_bgc_DMS_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_DMS_sk,"skl_DMS","mmol/m^2",tstr2D, tcstr, &
             "dissolved skl trace gas (DMS)",                                        &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns1, f_bgc_DMS_sk)
      if (f_bgc_DMSP_ml(1:1) /= 'x') &
         call define_hist_field(n_bgc_DMSP_ml,"ml_DMSP","mmol/m^3",tstr2D, tcstr, &
             "mixed layer precursor (DMSP)",                                        &
             "upper ocean", c1, c0,                &
             ns1, f_bgc_DMSP_ml)

      if (f_bgc_DMS_ml(1:1) /= 'x') &
         call define_hist_field(n_bgc_DMS_ml,"ml_DMS","mmol/m^3",tstr2D, tcstr, &
             "mixed layer trace gas (DMS)",                                        &
             "upper ocean", c1, c0,                &
             ns1, f_bgc_DMS_ml)  

      ! Tracers

      ! Ice Age
      if (f_iage(1:1) /= 'x') &
         call define_hist_field(n_iage,"iage","years",tstr2D, tcstr, &
             "sea ice age",                                        &
             "none", c1/(secday*days_per_year), c0,                &
             ns1, f_iage)

      ! First Year Ice Area
      if (f_FY(1:1) /= 'x') &
         call define_hist_field(n_FY,"FYarea"," ",tstr2D, tcstr, &
             "first-year ice area",                            &
             "weighted by ice area", c1, c0,                   &
              ns1, f_FY)

      ! Aerosols
      if (f_aero(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'aerosnossl', trim(nchar)
            call define_hist_field(n_aerosn1(n,:),vname_in,"kg/kg",   &
                tstr2D, tcstr,"snow ssl aerosol mass","none", c1, c0, &
                ns1, f_aero)
            write(vname_in,'(a,a)') 'aerosnoint', trim(nchar)
            call define_hist_field(n_aerosn2(n,:),vname_in,"kg/kg",   &
                tstr2D, tcstr,"snow int aerosol mass","none", c1, c0, &
                ns1, f_aero)
            write(vname_in,'(a,a)') 'aeroicessl', trim(nchar)
            call define_hist_field(n_aeroic1(n,:),vname_in,"kg/kg",  &
                tstr2D, tcstr,"ice ssl aerosol mass","none", c1, c0, &
                ns1, f_aero)
            write(vname_in,'(a,a)') 'aeroiceint', trim(nchar)
            call define_hist_field(n_aeroic2(n,:),vname_in,"kg/kg",  &
                tstr2D, tcstr,"ice int aerosol mass","none", c1, c0, &
                ns1, f_aero)
         enddo
      endif

      if (f_faero_atm(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'faero_atm', trim(nchar)
            call define_hist_field(n_faero_atm(n,:),vname_in,"kg/m^2 s", &
                tstr2D, tcstr,"aerosol deposition rate","none", c1, c0,  &
                ns1, f_faero_atm)
         enddo
      endif

      if (f_faero_ocn(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'faero_ocn', trim(nchar)
            call define_hist_field(n_faero_ocn(n,:),vname_in,"kg/m^2 s", &
                tstr2D, tcstr,"aerosol flux to ocean","none", c1, c0,    &
                ns1, f_faero_ocn)
         enddo
      endif

      enddo ! ns1

      ! 3D (category) variables looped separately for ordering
      do ns1 = 1, nstreams

        if (f_aicen(1:1) /= 'x') &
           call define_hist_field(n_aicen,"aicen","1",tstr3Dc, tcstr, & 
              "ice area, categories","none", c1, c0,                  &            
              ns1, f_aicen)

        if (f_vicen(1:1) /= 'x') &
           call define_hist_field(n_vicen,"vicen","m",tstr3Dc, tcstr, & 
              "ice volume, categories","none", c1, c0,                &            
              ns1, f_vicen)

        if (f_fsurfn_ai(1:1) /= 'x') &
           call define_hist_field(n_fsurfn_ai,"fsurfn_ai","W/m^2",tstr3Dc, tcstr, & 
              "net surface heat flux, categories","weighted by ice area", c1, c0, &            
              ns1, f_fsurfn_ai)
   
         if (f_fsicen_g(1:1) /= 'x') &
         call define_hist_field(n_fsicen_g,"fsicen_g","kg/m^2/s",tstr3Dc, tcstr, &
             "salt flux from gravity drainage to ocn",                              &
             "if positive, ocean gains salt", c1, c0,                   &
             ns1, f_fsicen_g)
      
        if (f_fcondtopn_ai(1:1) /= 'x') &
           call define_hist_field(n_fcondtopn_ai,"fcondtopn_ai","W/m^2",tstr3Dc, tcstr, &
              "top sfc conductive heat flux, cat","weighted by ice area", c1, c0,       &
              ns1, f_fcondtopn_ai)

        if (f_fmelttn_ai(1:1) /= 'x') &
           call define_hist_field(n_fmelttn_ai,"fmelttn_ai","W/m^2",tstr3Dc, tcstr, & 
              "net sfc heat flux causing melt, cat","weighted by ice area", c1, c0, &            
              ns1, f_fmelttn_ai)

        if (f_flatn_ai(1:1) /= 'x') &
           call define_hist_field(n_flatn_ai,"flatn_ai","W/m^2",tstr3Dc, tcstr, & 
              "latent heat flux, category","weighted by ice area", c1, c0,      &            
              ns1, f_flatn_ai)

        if (f_fbri(1:1) /= 'x') &
         call define_hist_field(n_fbri,"fbri","1",tstr3Dc, tcstr,        &
             "ice vol frac. with dynamic sal, cat",                     &
             "none", c1, c0,       &
             ns1, f_fbri)

      enddo ! ns1

      ! 3D (vertical) variables looped separately for ordering
!      do ns1 = 1, nstreams

!      if (f_field3dz(1:1) /= 'x') &
!         call define_hist_field(n_field3dz,"field3dz","1",tstr3Dz, tcstr, & 
!            "example 3dz field",                    &
!            "vertical profile", c1, c0,                  &
!            ns1, f_field3dz)

!      enddo ! ns1 

     ! !! 3D (vertical) ice biology variables !

      do ns1 = 1, nstreams
      
       if (f_upNO(1:1) /= 'x') &
         call define_hist_field(n_upNO,"upNO","mmol/m^3/d",tstr3Db, tcstr, &
             "Algal NO uptake rate",                           &
             "Positive flux is NO to N pool", secday, c0, &
             ns1, f_upNO)

       if (f_upNH(1:1) /= 'x') &
         call define_hist_field(n_upNH,"upNH","mmol/m^3/d",tstr3Db, tcstr, &
             "Algal NH uptake rate",                           &
             "Positive flux is NH to N pool", secday, c0,&
             ns1, f_upNH)
      enddo   !ns1


      ! 4D (categories, vertical) variables looped separately for ordering
      do ns1 = 1, nstreams

      if (f_Tinz(1:1) /= 'x') &
         call define_hist_field(n_Tinz,"Tinz","C",tstr4Di, tcstr, & 
            "ice internal temperatures on CICE grid",          &
            "vertical profile", c1, c0,                    &
            ns1, f_Tinz)


      if (f_Sinz(1:1) /= 'x') &
         call define_hist_field(n_Sinz,"Sinz","ppt",tstr4Di, tcstr, & 
            "ice internal bulk salinity",          &
            "vertical profile", c1, c0,                    &
            ns1, f_Sinz)

      enddo ! ns1

      do ns1 = 1, nstreams

      if (f_Tsnz(1:1) /= 'x') &
         call define_hist_field(n_Tsnz,"Tsnz","C",tstr4Ds, tcstr, & 
            "snow internal temperatures",          &
            "vertical profile", c1, c0,                    &
            ns1, f_Tsnz)

      enddo

       if (f_Tinz   (1:1) /= 'x') then
            if (allocated(Tinz4d)) deallocate(Tinz4d)
            allocate(Tinz4d(nx_block,ny_block,nzilyr,ncat_hist))
       endif

     
       if (f_Sinz   (1:1) /= 'x')  then
            if (allocated(Sinz4d)) deallocate(Sinz4d)
            allocate(Sinz4d(nx_block,ny_block,nzilyr,ncat_hist))
       endif

      
       if (f_Tsnz   (1:1) /= 'x') then
            if (allocated(Tsnz4d)) deallocate(Tsnz4d)
            allocate(Tsnz4d(nx_block,ny_block,nzslyr,ncat_hist))
       endif
       if (f_Sinz   (1:1) /= 'x') then
            if (allocated(Sinz4d)) deallocate(Sinz4d)
            allocate(Sinz4d(nx_block,ny_block,nzilyr,ncat_hist))
       endif

      
       ! (biology vertical grid)

      do ns1 = 1, nstreams
      
      ! if (f_upNO(1:1) /= 'x') &
      !   call define_hist_field(n_upNO,"upNO","mmol/m^3/d",tstr4Db, tcstr, &
      !       "Algal NO uptake rate",                           &
      !       "Positive flux is NO to N pool", secday, c0, &
      !       ns1, f_upNO)     

      ! if (f_upNH(1:1) /= 'x') &
      !   call define_hist_field(n_upNH,"upNH","mmol/m^3/d",tstr4Db, tcstr, &
      !       "Algal NH uptake rate",                           &
      !       "Positive flux is NH to N pool", secday, c0,&
      !       ns1, f_upNH)
     
       if (f_zTin(1:1) /= 'x') &
            call define_hist_field(n_zTin,"zTinb","C",tstr4Db, tcstr, &
                "ice internal temperatures on bio grid", "interpolated to bio grid", c1, c0,  &
                ns1, f_zTin)
      
      
       if (f_zphi(1:1) /= 'x') &
            call define_hist_field(n_zphi,"zphin","%",tstr4Db, tcstr, &
                "porosity", "brine volume fraction", c100, c0, &
                ns1, f_zphi)
         
       if (f_iDi(1:1) /= 'x') &
             call define_hist_field(n_iDi,"iDin","m^2/d",tstr4Db, tcstr, &
                "interface diffusivity", "on bio interface grid", secday, c0, &
                ns1, f_iDi)
      
       if (f_iki(1:1) /= 'x') &
            call define_hist_field(n_iki,"ikin","mm^2",tstr4Db, tcstr, &
                "permeability", "on bio interface grid", 1.0e6_dbl_kind, c0,&
                ns1, f_iki)

       if (f_bgc_NO(1:1) /= 'x') &
            call define_hist_field(n_bgc_NO,"bgc_NO","mmol/m^3",tstr4Db, tcstr, &
                "bulk nitrate ", "on bio grid", c1, c0,     &
                ns1, f_bgc_NO)
 
       if (f_bgc_NH(1:1) /= 'x') &
            call define_hist_field(n_bgc_NH,"bgc_NH","mmol/m^3",tstr4Db, tcstr, &
                "bulk ammonia/um ", "on bio grid", c1, c0,  &
                ns1, f_bgc_NH)

       if (f_bgc_N(1:1) /= 'x') &
            call define_hist_field(n_bgc_N,"bgc_N","mmol/m^3",tstr4Db, tcstr, &
                "bulk algal N conc. ", "on bio grid", c1, c0, &
                ns1, f_bgc_N)
     
       if (f_bgc_C(1:1) /= 'x') &
            call define_hist_field(n_bgc_C,"bgc_C","mmol/m^3",tstr4Db, tcstr, &
                "bulk algal carbon ", "on bio grid", c1, c0, &
                ns1, f_bgc_C)

       if (f_bgc_chl(1:1) /= 'x') &
            call define_hist_field(n_bgc_chl,"bgc_chl","mg/m^3",tstr4Db, tcstr, &
                "bulk algal chlorophyll ", "on bio grid", c1, c0,&
                ns1, f_bgc_chl)
      
       if (f_bgc_Sil(1:1) /= 'x') &
            call define_hist_field(n_bgc_Sil,"bgc_Sil","mmol/m^3",tstr4Db, tcstr, &
                "bulk silicate ", "on bio grid", c1, c0, &
                ns1, f_bgc_Sil)
      
       if (f_bgc_DMSPp(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMSPp,"bgc_DMSPp","mmol/m^3",tstr4Db, tcstr, &
                "bulk algal DMSP ", "on bio grid", c1, c0,&
                ns1, f_bgc_DMSPp)
      
       if (f_bgc_DMSPd(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMSPd,"bgc_DMSPd","mmol/m^3",tstr4Db, tcstr, &
                "bulk dissolved DMSP ", "on bio grid", c1, c0, &
                ns1, f_bgc_DMSPd)
  
       if (f_bgc_DMS(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMS,"bgc_DMS","mmol/m^3",tstr4Db, tcstr, &
                "bulk DMS gas ", "on bio grid", c1, c0, &
                ns1, f_bgc_DMS)
     
       if (f_bgc_PON(1:1) /= 'x') &
            call define_hist_field(n_bgc_PON,"bgc_PON","mmol/m^3",tstr4Db, tcstr, &
                "other bulk nitrogen pool ", "on bio grid valid for (2:nblyr+1)", c1, c0, &
                ns1, f_bgc_PON)
 
       if (f_bgc_S(1:1) /= 'x') &
            call define_hist_field(n_bgc_S,"bgc_S","ppt",tstr4Db, tcstr, &
                "bulk salinity", "on bio grid valid for (2:nblyr+1)", c1, c0, &
                ns1, f_bgc_S)
     
       if (f_growN(1:1) /= 'x') &
            call define_hist_field(n_growN,"growN","d^-1",tstr4Db, tcstr, &
                "Specific algal growth rate", "on bio grid valid for (2:nblyr+1)", secday , c0,  &
                ns1, f_growN)
      
       if (f_zfswin(1:1) /= 'x') &
            call define_hist_field(n_zfswin,"zfswin","W/m^2",tstr4Db, tcstr, &
                "internal ice PAR", "on bio grid", c1, c0, &
                ns1, f_zfswin)
    

      enddo  !ns1

      !-----------------------------------------------------------------
      ! other history variables
      !-----------------------------------------------------------------
      ! mechanical redistribution
      call init_hist_mechred
      ! melt ponds
      if (tr_pond) call init_hist_pond

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

      igrdz=.true.
      igrdz(n_NCAT     ) = f_NCAT
      igrdz(n_VGRDi    ) = f_VGRDi
      igrdz(n_VGRDs    ) = f_VGRDs
      igrdz(n_VGRDb    ) = f_VGRDb

      !-----------------------------------------------------------------
      ! diagnostic output
      !-----------------------------------------------------------------

      ntmp(:) = 0
      if (my_task == master_task) then
        write(nu_diag,*) ' '
        write(nu_diag,*) 'The following variables will be ', &
                         'written to the history tape: '
        write(nu_diag,101) 'description','units','variable','frequency','x'
        do n=1,num_avail_hist_fields_tot
           if (avail_hist_fields(n)%vhistfreq_n /= 0) &
           write(nu_diag,100) avail_hist_fields(n)%vdesc, &
              avail_hist_fields(n)%vunit, avail_hist_fields(n)%vname, &
              avail_hist_fields(n)%vhistfreq,avail_hist_fields(n)%vhistfreq_n
           do ns = 1, nstreams
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                 ntmp(ns)=ntmp(ns)+1
           enddo
        enddo ! num_avail_hist_fields_tot
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

      if (allocated(a2D)) deallocate(a2D)
      if (num_avail_hist_fields_2D > 0) &
      allocate(a2D(nx_block,ny_block,num_avail_hist_fields_2D,max_blocks))

      if (allocated(a3Dc)) deallocate(a3Dc)
      if (num_avail_hist_fields_3Dc > 0) &
      allocate(a3Dc(nx_block,ny_block,ncat_hist,num_avail_hist_fields_3Dc,max_blocks))

      nzlyr = max(nzilyr, nzslyr)
      if (allocated(a3Dz)) deallocate(a3Dz)
      if (num_avail_hist_fields_3Dz > 0) &
      allocate(a3Dz(nx_block,ny_block,nzlyr,num_avail_hist_fields_3Dz,max_blocks))

      nzlyrb = nzblyr
      if (allocated(a3Db)) deallocate(a3Db)
      if (num_avail_hist_fields_3Db > 0) &
      allocate(a3Db(nx_block,ny_block,nzlyrb,num_avail_hist_fields_3Db,max_blocks))

      if (allocated(a4Di)) deallocate(a4Di)
      if (num_avail_hist_fields_4Di > 0) &
      allocate(a4Di(nx_block,ny_block,nzilyr,ncat_hist,num_avail_hist_fields_4Di,max_blocks))
      if (allocated(a4Ds)) deallocate(a4Ds)
      if (num_avail_hist_fields_4Ds > 0) &
      allocate(a4Ds(nx_block,ny_block,nzslyr,ncat_hist,num_avail_hist_fields_4Ds,max_blocks))
      if (allocated(a4Db)) deallocate(a4Db)
      if (num_avail_hist_fields_4Db > 0) &
      allocate(a4Db(nx_block,ny_block,nzblyr,ncat_hist,num_avail_hist_fields_4Db,max_blocks))

      if (allocated(a2D))  a2D (:,:,:,:)     = c0
      if (allocated(a3Dc)) a3Dc(:,:,:,:,:)   = c0
      if (allocated(a3Dz)) a3Dz(:,:,:,:,:)   = c0
      if (allocated(a3Db)) a3Db(:,:,:,:,:)   = c0
      if (allocated(a4Di)) a4Di(:,:,:,:,:,:) = c0
      if (allocated(a4Ds)) a4Ds(:,:,:,:,:,:) = c0
      if (allocated(a4Db)) a4Db(:,:,:,:,:,:) = c0
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
! !IROUTINE: accum_hist - accumulate average ice quantities or snapshots
!
! !INTERFACE:
!
      subroutine accum_hist (dt)
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
      use ice_calendar, only: new_year, secday, write_history, &
                              write_ic, time, histfreq, nstreams, month, &
                              new_month
      use ice_state
      use ice_constants
      use ice_dyn_eap
      use ice_dyn_evp
      use ice_flux
      use ice_therm_shared, only: calculate_Tin_from_qin, Tmlt, ktherm
!      use ice_therm_vertical
      use ice_therm_mushy, only: temperature_mush, temperature_snow
      use ice_shortwave, only: apeffn
      use ice_timers
      use ice_zbgc_public, only:  &
                         nit, sil, amm, dmsp, dms, algalN, R_C2N, &
                         R_chl2N, nlt_bgc_N, &
                         nlt_bgc_NO, nlt_bgc_C, nlt_bgc_chl, &
                         nlt_bgc_NH, nlt_bgc_Sil, nlt_bgc_DMSPp, &
                         nlt_bgc_DMSPd, nlt_bgc_DMS, nlt_bgc_PON , &
                         S_tot, chl_net, PP_net, NO_net, &
                         zTin, zfswin, iki, iDi, zphi
      use ice_work, only: worka, workb, workz, workzn
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!EOP
!
      integer (kind=int_kind) :: &
           i,j,k,ic,n,ns,nn, &
           iblk             , & ! block index
           ilo,ihi,jlo,jhi  , & ! beginning and end of physical domain
           nstrm                ! nstreams (1 if writing initial condition)

      real (kind=dbl_kind) :: &
           ravgct           , & ! 1/avgct
           ravgctz          , & ! 1/avgct
           ai               , & ! aice_init
           ain                  ! aicen_init

      real (kind=dbl_kind) :: & 
           qn                , & ! temporary variable for enthalpy
           Tmlts                 !  temporary variable for melting temperature

      type (block) :: &
         this_block           ! block information for current block

      !---------------------------------------------------------------
      ! increment step counter
      !---------------------------------------------------------------

      n2D     = num_avail_hist_fields_2D
      n3Dccum = n2D     + num_avail_hist_fields_3Dc
      n3Dzcum = n3Dccum + num_avail_hist_fields_3Dz
      n3Dbcum = n3Dzcum + num_avail_hist_fields_3Db
      n4Dicum = n3Dbcum + num_avail_hist_fields_4Di
      n4Dscum = n4Dicum + num_avail_hist_fields_4Ds 
      n4Dbcum = n4Dscum + num_avail_hist_fields_4Db ! should equal num_avail_hist_fields_tot

      do ns = 1,nstreams
         if (.not. hist_avg .or. histfreq(ns) == '1') then  ! write snapshots
           do n = 1,n2D
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a2D(:,:,n,:) = c0
           enddo
           do n = n2D + 1, n3Dccum                  
              nn = n - n2D
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a3Dc(:,:,:,nn,:) = c0
           enddo
           do n = n3Dccum + 1, n3Dzcum
              nn = n - n3Dccum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a3Dz(:,:,:,nn,:) = c0
           enddo
           do n = n3Dzcum + 1, n3Dbcum
              nn = n - n3Dzcum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a3Db(:,:,:,nn,:) = c0
           enddo
           do n = n3Dbcum + 1, n4Dicum
              nn = n - n3Dbcum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a4Di(:,:,:,:,nn,:) = c0
           enddo
           do n = n4Dicum + 1, n4Dscum
              nn = n - n4Dicum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a4Ds(:,:,:,:,nn,:) = c0
           enddo
           do n = n4Dscum + 1, n4Dbcum
              nn = n - n4Dscum
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                  a4Db(:,:,:,:,nn,:) = c0
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
!            call accum_hist_field(n_example,iblk, vice(:,:,iblk), a2D)
         if (f_hi     (1:1) /= 'x') &
             call accum_hist_field(n_hi,     iblk, vice(:,:,iblk), a2D)
         if (f_hs     (1:1) /= 'x') &
             call accum_hist_field(n_hs,     iblk, vsno(:,:,iblk), a2D)
         if (f_Tsfc   (1:1) /= 'x') &
             call accum_hist_field(n_Tsfc,   iblk, trcr(:,:,nt_Tsfc,iblk), a2D)
         if (f_aice   (1:1) /= 'x') &
             call accum_hist_field(n_aice,   iblk, aice(:,:,iblk), a2D)
         if (f_uvel   (1:1) /= 'x') &
             call accum_hist_field(n_uvel,   iblk, uvel(:,:,iblk), a2D)
         if (f_vvel   (1:1) /= 'x') &
             call accum_hist_field(n_vvel,   iblk, vvel(:,:,iblk), a2D)

         if (f_fswdn  (1:1) /= 'x') &
             call accum_hist_field(n_fswdn,  iblk, fsw(:,:,iblk), a2D)
         if (f_flwdn  (1:1) /= 'x') &
             call accum_hist_field(n_flwdn,  iblk, flw(:,:,iblk), a2D)
         if (f_snow   (1:1) /= 'x') &
             call accum_hist_field(n_snow,   iblk, fsnow(:,:,iblk), a2D)
         if (f_snow_ai(1:1) /= 'x') &
             call accum_hist_field(n_snow_ai,iblk, fsnow(:,:,iblk)*workb(:,:), a2D)
         if (f_rain   (1:1) /= 'x') &
             call accum_hist_field(n_rain,   iblk, frain(:,:,iblk), a2D)
         if (f_rain_ai(1:1) /= 'x') &
             call accum_hist_field(n_rain_ai,iblk, frain(:,:,iblk)*workb(:,:), a2D)

         if (f_sst    (1:1) /= 'x') &
             call accum_hist_field(n_sst,    iblk, sst(:,:,iblk), a2D)
         if (f_sss    (1:1) /= 'x') &
             call accum_hist_field(n_sss,    iblk, sss(:,:,iblk), a2D)
         if (f_uocn   (1:1) /= 'x') &
             call accum_hist_field(n_uocn,   iblk, uocn(:,:,iblk), a2D)
         if (f_vocn   (1:1) /= 'x') &
             call accum_hist_field(n_vocn,   iblk, vocn(:,:,iblk), a2D)
         if (f_frzmlt (1:1) /= 'x') &
             call accum_hist_field(n_frzmlt, iblk, frzmlt(:,:,iblk), a2D)

         if (f_fswfac (1:1) /= 'x') &
             call accum_hist_field(n_fswfac, iblk, fswfac(:,:,iblk), a2D)
         if (f_fswabs (1:1) /= 'x') &
             call accum_hist_field(n_fswabs, iblk, fswabs(:,:,iblk), a2D)
         if (f_fswabs_ai(1:1)/= 'x') &
             call accum_hist_field(n_fswabs_ai, iblk, fswabs(:,:,iblk)*workb(:,:), a2D)

         if (f_albsni (1:1) /= 'x') &
             call accum_hist_field(n_albsni, iblk, &
                                  (awtvdr*alvdr(:,:,iblk) &
                                 + awtidr*alidr(:,:,iblk) &
                                 + awtvdf*alvdf(:,:,iblk) &
                                 + awtidf*alidf(:,:,iblk))*aice(:,:,iblk), a2D)
!                                              awtvdr*alvdr(:,:,iblk) &
!                                            + awtidr*alidr(:,:,iblk) &
!                                            + awtvdf*alvdf(:,:,iblk) &
!                                            + awtidf*alidf(:,:,iblk), a2D)
         if (f_alvdr  (1:1) /= 'x') &
             call accum_hist_field(n_alvdr,  iblk, alvdr(:,:,iblk), a2D)
         if (f_alidr  (1:1) /= 'x') &
             call accum_hist_field(n_alidr,  iblk, alidr(:,:,iblk), a2D)

         if (f_albice (1:1) /= 'x') &
             call accum_hist_field(n_albice, iblk, albice(:,:,iblk), a2D)
         if (f_albsno (1:1) /= 'x') &
             call accum_hist_field(n_albsno, iblk, albsno(:,:,iblk), a2D)
         if (f_albpnd (1:1) /= 'x') &
             call accum_hist_field(n_albpnd, iblk, albpnd(:,:,iblk), a2D)
         if (f_coszen (1:1) /= 'x') &
             call accum_hist_field(n_coszen, iblk, coszen(:,:,iblk), a2D)

         if (f_flat   (1:1) /= 'x') &
             call accum_hist_field(n_flat,   iblk, flat(:,:,iblk), a2D)
         if (f_flat_ai(1:1) /= 'x') &
             call accum_hist_field(n_flat_ai,iblk, flat(:,:,iblk)*workb(:,:), a2D)
         if (f_fsens  (1:1) /= 'x') &
             call accum_hist_field(n_fsens,   iblk, fsens(:,:,iblk), a2D)
         if (f_fsens_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsens_ai,iblk, fsens(:,:,iblk)*workb(:,:), a2D)
         if (f_flwup  (1:1) /= 'x') &
             call accum_hist_field(n_flwup,   iblk, flwout(:,:,iblk), a2D)
         if (f_flwup_ai(1:1)/= 'x') &
             call accum_hist_field(n_flwup_ai,iblk, flwout(:,:,iblk)*workb(:,:), a2D)
         if (f_evap   (1:1) /= 'x') &
             call accum_hist_field(n_evap,   iblk, evap(:,:,iblk), a2D)
         if (f_evap_ai(1:1) /= 'x') &
             call accum_hist_field(n_evap_ai,iblk, evap(:,:,iblk)*workb(:,:), a2D)

         if (f_Tair   (1:1) /= 'x') &
             call accum_hist_field(n_Tair,   iblk, Tair(:,:,iblk), a2D)
         if (f_Tref   (1:1) /= 'x') &
             call accum_hist_field(n_Tref,   iblk, Tref(:,:,iblk), a2D)
         if (f_Qref   (1:1) /= 'x') &
             call accum_hist_field(n_Qref,   iblk, Qref(:,:,iblk), a2D)
         if (f_congel (1:1) /= 'x') &
             call accum_hist_field(n_congel, iblk, congel(:,:,iblk), a2D)
         if (f_frazil (1:1) /= 'x') &
             call accum_hist_field(n_frazil, iblk, frazil(:,:,iblk), a2D)
         if (f_snoice (1:1) /= 'x') &
             call accum_hist_field(n_snoice, iblk, snoice(:,:,iblk), a2D)
         if (f_dsnow (1:1) /= 'x') &
             call accum_hist_field(n_dsnow, iblk, dsnow(:,:,iblk), a2D)
         if (f_meltt  (1:1) /= 'x') &
             call accum_hist_field(n_meltt,  iblk, meltt(:,:,iblk), a2D)
         if (f_melts  (1:1) /= 'x') &
              call accum_hist_field(n_melts,  iblk, melts(:,:,iblk), a2D)
         if (f_meltb  (1:1) /= 'x') &
             call accum_hist_field(n_meltb,  iblk, meltb(:,:,iblk), a2D)
         if (f_meltl  (1:1) /= 'x') &
             call accum_hist_field(n_meltl,  iblk, meltl(:,:,iblk), a2D)

         if (f_fresh  (1:1) /= 'x') &
             call accum_hist_field(n_fresh,   iblk, fresh(:,:,iblk), a2D)
         if (f_fresh_ai(1:1)/= 'x') &
             call accum_hist_field(n_fresh_ai,iblk, fresh_gbm(:,:,iblk), a2D)
         if (f_fsalt  (1:1) /= 'x') &
             call accum_hist_field(n_fsalt,   iblk, fsalt(:,:,iblk), a2D)
         if (f_fsalt_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsalt_ai,iblk, fsalt_gbm(:,:,iblk), a2D)

         if (f_fsice  (1:1) /= 'x') &
             call accum_hist_field(n_fsice,   iblk, fsice(:,:,iblk), a2D)
         if (f_fsice_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsice_ai,iblk, fsice_gbm(:,:,iblk), a2D)
         if (f_fsice_g  (1:1) /= 'x') &
             call accum_hist_field(n_fsice_g,   iblk, fsice_g(:,:,iblk), a2D)
         if (f_fsice_g_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsice_g_ai,iblk, fsice_g_gbm(:,:,iblk), a2D)

         if (f_fNO  (1:1) /= 'x') &
             call accum_hist_field(n_fNO,   iblk, flux_bio(:,:,nlt_bgc_NO,iblk), a2D)
         if (f_fNO_ai(1:1)/= 'x') &
             call accum_hist_field(n_fNO_ai,iblk, flux_bio_gbm(:,:,nlt_bgc_NO,iblk), a2D)
         if (f_fNO_g  (1:1) /= 'x') &
             call accum_hist_field(n_fNO_g,   iblk, flux_bio_g(:,:,nlt_bgc_NO,iblk), a2D)
         if (f_fNO_g_ai(1:1)/= 'x') &
             call accum_hist_field(n_fNO_g_ai,iblk, flux_bio_g_gbm(:,:,nlt_bgc_NO,iblk), a2D)

         if (f_fNH  (1:1) /= 'x') &
             call accum_hist_field(n_fNH,   iblk, flux_bio(:,:,nlt_bgc_NH,iblk), a2D)
         if (f_fNH_ai(1:1)/= 'x') &
             call accum_hist_field(n_fNH_ai,iblk, flux_bio_gbm(:,:,nlt_bgc_NH,iblk), a2D)
         if (f_fNH_g  (1:1) /= 'x') &
             call accum_hist_field(n_fNH_g,   iblk, flux_bio_g(:,:,nlt_bgc_NH,iblk), a2D)
         if (f_fNH_g_ai(1:1)/= 'x') &
             call accum_hist_field(n_fNH_g_ai,iblk, flux_bio_g_gbm(:,:,nlt_bgc_NH,iblk), a2D)

         if (f_fN  (1:1) /= 'x') &
             call accum_hist_field(n_fN,   iblk, flux_bio(:,:,nlt_bgc_N,iblk), a2D)
         if (f_fN_ai(1:1)/= 'x') &
             call accum_hist_field(n_fN_ai,iblk, flux_bio_gbm(:,:,nlt_bgc_N,iblk), a2D)
         if (f_fN_g  (1:1) /= 'x') &
             call accum_hist_field(n_fN_g,   iblk, flux_bio_g(:,:,nlt_bgc_N,iblk), a2D)
         if (f_fN_g_ai(1:1)/= 'x') &
             call accum_hist_field(n_fN_g_ai,iblk, flux_bio_g_gbm(:,:,nlt_bgc_N,iblk), a2D)


         if (f_fSil  (1:1) /= 'x') &
             call accum_hist_field(n_fSil,   iblk, flux_bio(:,:,nlt_bgc_Sil,iblk), a2D)
         if (f_fSil_ai(1:1)/= 'x') &
             call accum_hist_field(n_fSil_ai,iblk, flux_bio_gbm(:,:,nlt_bgc_Sil,iblk), a2D)
         if (f_fSil_g  (1:1) /= 'x') &
             call accum_hist_field(n_fSil_g,   iblk, flux_bio_g(:,:,nlt_bgc_Sil,iblk), a2D)
         if (f_fSil_g_ai(1:1)/= 'x') &
             call accum_hist_field(n_fSil_g_ai,iblk, flux_bio_g_gbm(:,:,nlt_bgc_Sil,iblk), a2D)

         if (f_Stot  (1:1) /= 'x') &
             call accum_hist_field(n_Stot,   iblk, S_tot(:,:,iblk), a2D)
         if (f_chlnet  (1:1) /= 'x') &
             call accum_hist_field(n_chlnet, iblk, chl_net(:,:,iblk), a2D)
         if (f_PPnet  (1:1) /= 'x') &
             call accum_hist_field(n_PPnet,   iblk, PP_net(:,:,iblk), a2D)
         if (f_NOnet  (1:1) /= 'x') &
             call accum_hist_field(n_NOnet,   iblk, NO_net(:,:,iblk), a2D)
!         if (f_Tsf_ice  (1:1) /= 'x') &
!             call accum_hist_field(n_Tsf_ice,   iblk, Tsf_ice(:,:,iblk), a2D)
         if (f_fhocn  (1:1) /= 'x') &
             call accum_hist_field(n_fhocn,   iblk, fhocn(:,:,iblk), a2D)
         if (f_fhocn_ai(1:1)/= 'x') &
             call accum_hist_field(n_fhocn_ai,iblk, fhocn_gbm(:,:,iblk), a2D)
         if (f_fswthru(1:1) /= 'x') &
             call accum_hist_field(n_fswthru, iblk, fswthru(:,:,iblk), a2D)
         if (f_fswthru_ai(1:1)/= 'x') &
             call accum_hist_field(n_fswthru_ai,iblk, fswthru_gbm(:,:,iblk), a2D)
               
         if (f_strairx(1:1) /= 'x') &
             call accum_hist_field(n_strairx, iblk, strairx(:,:,iblk), a2D)
         if (f_strairy(1:1) /= 'x') &
             call accum_hist_field(n_strairy, iblk, strairy(:,:,iblk), a2D)
         if (f_strtltx(1:1) /= 'x') &
             call accum_hist_field(n_strtltx, iblk, strtltx(:,:,iblk), a2D)
         if (f_strtlty(1:1) /= 'x') &
             call accum_hist_field(n_strtlty, iblk, strtlty(:,:,iblk), a2D)
         if (f_strcorx(1:1) /= 'x') &
             call accum_hist_field(n_strcorx, iblk, fm(:,:,iblk)*vvel(:,:,iblk), a2D)
         if (f_strcory(1:1) /= 'x') &
             call accum_hist_field(n_strcory, iblk,-fm(:,:,iblk)*uvel(:,:,iblk), a2D)
         if (f_strocnx(1:1) /= 'x') &
             call accum_hist_field(n_strocnx, iblk, strocnx(:,:,iblk), a2D)
         if (f_strocny(1:1) /= 'x') &
             call accum_hist_field(n_strocny, iblk, strocny(:,:,iblk), a2D)
         if (f_strintx(1:1) /= 'x') &
             call accum_hist_field(n_strintx, iblk, strintx(:,:,iblk), a2D)
         if (f_strinty(1:1) /= 'x') &
             call accum_hist_field(n_strinty, iblk, strinty(:,:,iblk), a2D)
         if (f_strength(1:1)/= 'x') &
             call accum_hist_field(n_strength,iblk, strength(:,:,iblk), a2D)

! The following fields (divu, shear, sig1, and sig2) will be smeared
!  if averaged over more than a few days.
! Snapshots may be more useful (see below).

!        if (f_divu   (1:1) /= 'x') &
!             call accum_hist_field(n_divu,    iblk, divu(:,:,iblk), a2D)
!        if (f_shear  (1:1) /= 'x') &
!             call accum_hist_field(n_shear,   iblk, shear(:,:,iblk), a2D)
!        if (f_sig1   (1:1) /= 'x') &
!             call accum_hist_field(n_sig1,    iblk, sig1(:,:,iblk), a2D)
!        if (f_sig2   (1:1) /= 'x') &
!             call accum_hist_field(n_sig2,    iblk, sig2(:,:,iblk), a2D)
!        if (f_trsig  (1:1) /= 'x') &
!             call accum_hist_field(n_trsig,   iblk, trsig(:,:,iblk), a2D)

         if (f_dvidtt (1:1) /= 'x') &
             call accum_hist_field(n_dvidtt,  iblk, dvidtt(:,:,iblk), a2D)
         if (f_dvidtd (1:1) /= 'x') &
             call accum_hist_field(n_dvidtd,  iblk, dvidtd(:,:,iblk), a2D)
         if (f_daidtt (1:1) /= 'x') &
             call accum_hist_field(n_daidtt,  iblk, daidtt(:,:,iblk), a2D)
         if (f_daidtd (1:1) /= 'x') &
             call accum_hist_field(n_daidtd,  iblk, daidtd(:,:,iblk), a2D)

         ! BGC
         if (f_bgc_N_sk(1:1)/= 'x') &
             call accum_hist_field(n_bgc_N_sk,iblk, &
                        trcr(:,:,nt_bgc_N_sk,iblk), a2D)
  
         if (f_bgc_C_sk(1:1)/= 'x') &
             call accum_hist_field(n_bgc_C_sk,iblk, &
                        trcr(:,:,nt_bgc_C_sk,iblk), a2D)  
         if (f_bgc_chl_sk(1:1)/= 'x') &
             call accum_hist_field(n_bgc_chl_sk,iblk, &
                        trcr(:,:,nt_bgc_chl_sk,iblk), a2D)  
         if (f_bgc_Nit_sk(1:1)/= 'x') &
             call accum_hist_field(n_bgc_Nit_sk,iblk, &
                        trcr(:,:,nt_bgc_Nit_sk,iblk), a2D)  
         if (f_bgc_Am_sk(1:1)/= 'x') &
             call accum_hist_field(n_bgc_Am_sk,iblk, &
                        trcr(:,:,nt_bgc_Am_sk,iblk), a2D)  
         if (f_bgc_Sil_sk(1:1)/= 'x') &
             call accum_hist_field(n_bgc_Sil_sk,iblk, &
                        trcr(:,:,nt_bgc_Sil_sk,iblk), a2D)  
         if (f_bgc_DMSPp_sk(1:1)/= 'x') &
             call accum_hist_field(n_bgc_DMSPp_sk,iblk, &
                        trcr(:,:,nt_bgc_DMSPp_sk,iblk), a2D)  

        if (f_bgc_DMSPd_sk(1:1)/= 'x') &
             call accum_hist_field(n_bgc_DMSPd_sk,iblk, &
                        trcr(:,:,nt_bgc_DMSPd_sk,iblk), a2D)  
         if (f_bgc_DMS_sk(1:1)/= 'x') &
             call accum_hist_field(n_bgc_DMS_sk,iblk, &
                        trcr(:,:,nt_bgc_DMS_sk,iblk), a2D)  
         if (f_bgc_Nit_ml(1:1)/= 'x') &
             call accum_hist_field(n_bgc_Nit_ml,iblk, &
                        trcr(:,:,nt_bgc_Nit_ml,iblk), a2D)  
         if (f_bgc_Am_ml(1:1)/= 'x') &
             call accum_hist_field(n_bgc_Am_ml,iblk, &
                        trcr(:,:,nt_bgc_Am_ml,iblk), a2D)  
         if (f_bgc_Sil_ml(1:1)/= 'x') &
             call accum_hist_field(n_bgc_Sil_ml,iblk, &
                        trcr(:,:,nt_bgc_Sil_ml,iblk), a2D)  
         if (f_bgc_DMSP_ml(1:1)/= 'x') &
             call accum_hist_field(n_bgc_DMSP_ml,iblk, &
                        trcr(:,:,nt_bgc_DMSP_ml,iblk), a2D)  
         if (f_bgc_DMS_ml(1:1)/= 'x') &
             call accum_hist_field(n_bgc_DMS_ml,iblk, &
                        trcr(:,:,nt_bgc_DMS_ml,iblk), a2D)  
 
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         if (f_fsurf_ai(1:1)/= 'x') &
             call accum_hist_field(n_fsurf_ai,iblk, fsurf(:,:,iblk)*workb(:,:), a2D)
         if (f_fcondtop_ai(1:1)/= 'x') &
             call accum_hist_field(n_fcondtop_ai, iblk, &
                                                 fcondtop(:,:,iblk)*workb(:,:), a2D)

         if (f_icepresent(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) worka(i,j) = c1
           enddo
           enddo
           call accum_hist_field(n_icepresent, iblk, worka(:,:), a2D)
         endif

         ! 3D category fields
         if (f_aicen   (1:1) /= 'x') &
             call accum_hist_field(n_aicen-n2D, iblk, ncat_hist, &
                                   aicen(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_vicen   (1:1) /= 'x') &
             call accum_hist_field(n_vicen-n2D, iblk, ncat_hist, &
                                   vicen(:,:,1:ncat_hist,iblk), a3Dc)

         if (f_fsurfn_ai   (1:1) /= 'x') &
             call accum_hist_field(n_fsurfn_ai-n2D, iblk, ncat_hist, &
                  fsurfn(:,:,1:ncat_hist,iblk)*aicen_init(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_fsicen_g  (1:1) /= 'x') &
             call accum_hist_field(n_fsicen_g-n2D, iblk, ncat_hist, &
                  fsicen_g(:,:,1:ncat_hist,iblk)*aicen_init(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_fcondtopn_ai   (1:1) /= 'x') &
             call accum_hist_field(n_fcondtopn_ai-n2D, iblk, ncat_hist, &
                  fcondtopn(:,:,1:ncat_hist,iblk)*aicen_init(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_flatn_ai   (1:1) /= 'x') &
             call accum_hist_field(n_flatn_ai-n2D, iblk, ncat_hist, &
                  flatn(:,:,1:ncat_hist,iblk)*aicen_init(:,:,1:ncat_hist,iblk), a3Dc)
         ! Calculate surface heat flux that causes melt (calculated by the 
         ! atmos in HadGEM3 so needed for checking purposes)
         if (f_fmelttn_ai   (1:1) /= 'x') &
             call accum_hist_field(n_fmelttn_ai-n2D, iblk, ncat_hist, &
                  max(fsurfn(:,:,1:ncat_hist,iblk) - fcondtopn(:,:,1:ncat_hist,iblk),c0) &
                      *aicen_init(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_fbri   (1:1) /= 'x') &
             call accum_hist_field(n_fbri-n2D, iblk, ncat_hist, &
                                   trcrn(:,:,nt_fbri,1:ncat_hist,iblk), a3Dc)

! example for 3D field (x,y,z)
!         if (f_field3dz   (1:1) /= 'x') &
!             call accum_hist_field(n_field3dz-n3Dccum, iblk, nzilyr, &
!                                   field3dz(:,:,1:nzilyr,iblk), a3Dz)
        
        if (f_upNO (1:1) /= 'x') then
            workz(:,:,:) = c0
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vice(i,j,iblk) > puny) then
                  workz(i,j,1) = upNO(i,j,1,iblk)
                  workz(i,j,2:nblyr+1) = upNO(i,j,1:nblyr,iblk)
                  workz(i,j,nblyr+2) = upNO(i,j,nblyr,iblk)
                  endif
                enddo !i
               enddo  !j
             call accum_hist_field(n_upNO- n3Dccum, iblk, nzblyr, &
                                   workz(:,:,1:nzblyr), a3Db)
        endif        
          if (f_upNH (1:1) /= 'x') then
            workz(:,:,:) = c0
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vice(i,j,iblk) > puny) then
                  workz(i,j,1) = upNH(i,j,1,iblk)
                  workz(i,j,2:nblyr+1) = upNH(i,j,1:nblyr,iblk)
                  workz(i,j,nblyr+2) = upNH(i,j,nblyr,iblk)
                  endif
                enddo !i
               enddo  !j
             call accum_hist_field(n_upNH- n3Dccum, iblk, nzblyr, &
                                   workz(:,:,1:nzblyr), a3Db)
        endif        
        

         ! 4D category fields
         if (f_Tinz   (1:1) /= 'x') then
            Tinz4d(:,:,:,:) = c0
            if (ktherm == 2) then
               do n = 1, ncat_hist
                  do j = jlo, jhi
                  do i = ilo, ihi
                     do k = 1, nzilyr
                        Tinz4d(i,j,k,n) = temperature_mush( &
                             trcrn(i,j,nt_qice+k-1,n,iblk), trcrn(i,j,nt_sice+k-1,n,iblk))
                     enddo
                  enddo
                  enddo
               enddo
            else
               do n = 1, ncat_hist
                  do j = jlo, jhi
                  do i = ilo, ihi
                     do k = 1, nzilyr
                        qn = trcrn(i,j,nt_qice+k-1,n,iblk)
                        Tinz4d(i,j,k,n) = calculate_Tin_from_qin(qn,Tmlt(k))

!                       Tmlts = trcrn(i,j,nt_sice+k-1,n,iblk)*(-mlt_a +  &
!                               mlt_b*SQRT(trcrn(i,j,nt_sice+k-1,n,iblk)) - &
!                               mlt_c*trcrn(i,j,nt_sice+k-1,n,iblk)) 
!                                !-trcrn(i,j,nt_sice+k-1,n,iblk)* depressT
!                       Tinz4d(i,j,k,n) =  calculate_Tin_from_qin(qn,Tmlts)
                     enddo
                  enddo
                  enddo
               enddo
            endif
            call accum_hist_field(n_Tinz-n3Dbcum, iblk, nzilyr, ncat_hist, &
                                  Tinz4d(:,:,1:nzilyr,1:ncat_hist), a4Di)
         endif
         if (f_Sinz   (1:1) /= 'x') then
            Sinz4d(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                       Sinz4d(i,j,1:nzilyr,n) = trcrn(i,j,nt_sice:nt_sice+nzilyr-1,n,iblk) 
                  endif
               enddo
               enddo
            enddo
            call accum_hist_field(n_Sinz-n3Dbcum, iblk, nzilyr, ncat_hist, &
                                  Sinz4d(:,:,1:nzilyr,1:ncat_hist), a4Di)
         endif
         
         if (f_Sinz   (1:1) /= 'x') then
            Sinz4d(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  do k = 1, nzilyr
                     Sinz4d(i,j,k,n) = trcrn(i,j,nt_sice+k-1,n,iblk)
                  enddo
               enddo
               enddo
            enddo
            call accum_hist_field(n_Sinz-n3Dzcum, iblk, nzilyr, ncat_hist, &
                                  Sinz4d(:,:,1:nzilyr,1:ncat_hist), a4Di)
         endif

         if (f_Tsnz   (1:1) /= 'x') then
            Tsnz4d(:,:,:,:) = c0
            if (ktherm == 2) then
               do n = 1, ncat_hist
                  do j = jlo, jhi
                  do i = ilo, ihi
                     do k = 1, nzslyr
                        qn = trcrn(i,j,nt_qsno+k-1,n,iblk)
                        Tsnz4d(i,j,k,n) = temperature_snow(trcrn(i,j,nt_qsno+k-1,n,iblk))
                     enddo
                  enddo
                  enddo
               enddo
            else
               do n = 1, ncat_hist
                  do j = jlo, jhi
                  do i = ilo, ihi
                     do k = 1, nzslyr
                        qn = trcrn(i,j,nt_qsno+k-1,n,iblk)
                        Tsnz4d(i,j,k,n) = (Lfresh + qn/rhos)/cp_ice
                     enddo
                  enddo
                  enddo
               enddo
            endif
            call accum_hist_field(n_Tsnz-n4Dicum, iblk, nzslyr, ncat_hist, &
                                  Tsnz4d(:,:,1:nzslyr,1:ncat_hist), a4Ds)
         endif
         
        if (f_bgc_N   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_N,n,iblk) !*aicen_init(i,j,n,iblk)
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_N:nt_bgc_N+nblyr-1,n,iblk)
                    workzn(i,j,nblyr+2,n) = algalN(i,j,iblk) !*aicen_init(i,j,n,iblk)
                  endif
                enddo !i
               enddo  !j
               
            enddo     !n
            call accum_hist_field(n_bgc_N-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif
         
        if (f_bgc_C   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_C,n,iblk) 
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_C:nt_bgc_C+nblyr-1,n,iblk) 
                    workzn(i,j,nblyr+2,n) = R_C2N*algalN(i,j,iblk)
                  
                  endif
                enddo !i
               enddo  !j
            enddo
            call accum_hist_field(n_bgc_C-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif
         

        if (f_bgc_chl   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_chl,n,iblk)
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_chl:nt_bgc_chl+nblyr-1,n,iblk) 
                    workzn(i,j,nblyr+2,n) = R_chl2N*algalN(i,j,iblk) 
                  endif
                enddo !i
               enddo  !j
               
            enddo
            call accum_hist_field(n_bgc_chl-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif
         

        if (f_bgc_NO   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_NO,n,iblk)    
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_NO:nt_bgc_NO+nblyr-1,n,iblk)   
                    workzn(i,j,nblyr+2,n) = nit(i,j,iblk)     
                  endif
                enddo !i
               enddo  !j
            enddo
            call accum_hist_field(n_bgc_NO-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif

        if (f_bgc_NH   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_NH,n,iblk) 
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_NH:nt_bgc_NH+nblyr-1,n,iblk)
                    workzn(i,j,nblyr+2,n) = amm(i,j,iblk)
                  endif
                enddo !i
               enddo  !j
               
            enddo
            call accum_hist_field(n_bgc_NH-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif

        if (f_bgc_Sil   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_Sil,n,iblk)
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_Sil:nt_bgc_Sil+nblyr-1,n,iblk) 
                    workzn(i,j,nblyr+2,n) = sil(i,j,iblk) 
                  endif
                enddo !i
               enddo  !j
            enddo
            call accum_hist_field(n_bgc_Sil-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif
  
        if (f_bgc_DMSPd   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_DMSPd,n,iblk) 
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_DMSPd:nt_bgc_DMSPd+nblyr-1,n,iblk)
                    workzn(i,j,nblyr+2,n) = dmsp(i,j,iblk)
                  endif
                enddo !i
               enddo  !j
            enddo
            call accum_hist_field(n_bgc_DMSPd-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif       

        if (f_bgc_DMSPp   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_DMSPp,n,iblk) 
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_DMSPp:nt_bgc_DMSPp+nblyr-1,n,iblk)
                    workzn(i,j,nblyr+2,n) = dmsp(i,j,iblk)
                  endif
                enddo !i
               enddo  !j
            enddo
            call accum_hist_field(n_bgc_DMSPp-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif

       if (f_bgc_DMS   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_DMS,n,iblk) 
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_DMS:nt_bgc_DMS+nblyr-1,n,iblk)
                    workzn(i,j,nblyr+2,n) = dms(i,j,iblk) 
                  endif
                enddo !i
               enddo  !j
            enddo
            call accum_hist_field(n_bgc_DMS-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif


       if (f_bgc_PON   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_PON,n,iblk) 
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_PON:nt_bgc_PON+nblyr-1,n,iblk) 
                    workzn(i,j,nblyr+2,n) = nit(i,j,iblk) 
                  endif
                enddo !i
               enddo  !j
            enddo
            call accum_hist_field(n_bgc_PON-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif

       if (f_bgc_S   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (vicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = trcrn(i,j,nt_bgc_S,n,iblk)     
                    workzn(i,j,2:nblyr+1,n) = trcrn(i,j,nt_bgc_S:nt_bgc_S+nblyr-1,n,iblk)
                    workzn(i,j,nblyr+2,n) = sss(i,j,iblk)    
                  endif
                enddo !i
               enddo  !j
            enddo
            call accum_hist_field(n_bgc_S-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif

       if (f_zTin  (1:1) /= 'x')  &
            call accum_hist_field(n_zTin-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  zTin(:,:,1:nzblyr,1:ncat_hist,iblk), a4Db)

       if (f_zfswin   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                if (aicen(i,j,n,iblk) > c0) then
                   workzn(i,j,1:nblyr+1,n) = zfswin(i,j,1:nblyr+1,n,iblk)
                   workzn(i,j,nzblyr,n)   = zfswin(i,j,nblyr+1,n,iblk)
                endif
               enddo  !j
               enddo  !i          
            enddo
            call accum_hist_field(n_zfswin-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
       endif

       if (f_iDi   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
              do k = 1,nzblyr-1
               do j = jlo, jhi
               do i = ilo, ihi
                if (aicen(i,j,n,iblk) > c0) then
                   workzn(i,j,k,n) = iDi(i,j,k,n,iblk)*(vicen(i,j,n,iblk)/aicen(i,j,n,iblk))**2
                   workzn(i,j,nzblyr,n)   = workzn(i,j,nzblyr-1,n)     
                endif
               enddo  !j
               enddo  !i
              enddo   !k        
            enddo
            call accum_hist_field(n_iDi-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
       endif
       if (f_iki   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                if (aicen(i,j,n,iblk) > c0) then
                   workzn(i,j,1:nblyr+1,n) = iki(i,j,1:nblyr+1,n,iblk)
                   workzn(i,j,nzblyr,n)   = iki(i,j,nblyr+1,n,iblk) 
                endif
               enddo  !j
               enddo  !i      
            enddo  !n
            call accum_hist_field(n_iki-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
       endif

       if (f_growN   (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                  if (aicen(i,j,n,iblk) > puny) then
                    workzn(i,j,1,n) = growN(i,j,1,n,iblk) 
                    workzn(i,j,2:nblyr+1,n) = growN(i,j,1:nblyr,n,iblk)
                    workzn(i,j,nblyr+2,n) = growN(i,j,nblyr,n,iblk) 
                  endif
              enddo!j
              enddo!i 
            enddo  !n
            call accum_hist_field(n_growN-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
       endif

       if (f_zphi  (1:1) /= 'x') then
            workzn(:,:,:,:) = c0
            do n = 1, ncat_hist
               do j = jlo, jhi
               do i = ilo, ihi
                if (aicen(i,j,n,iblk) > c0) then
                   workzn(i,j,1:nzblyr,n) = zphi(i,j,1:nzblyr,n,iblk)
                endif
               enddo  !j
               enddo  !i
            enddo     !n
            call accum_hist_field(n_zphi-n4Dscum, iblk, nzblyr, ncat_hist, &
                                  workzn(:,:,1:nzblyr,1:ncat_hist), a4Db)
         endif

        ! Calculate aggregate surface melt flux by summing category values
        if (f_fmeltt_ai(1:1) /= 'x') then
         do ns = 1, nstreams
           if (n_fmeltt_ai(ns) /= 0) then
              worka(:,:) = c0
              do j = jlo, jhi
              do i = ilo, ihi
               if (tmask(i,j,iblk)) then
                 do n=1,ncat_hist
                    worka(i,j)  = worka(i,j) + a3Dc(i,j,n,n_fmelttn_ai(ns)-n2D,iblk)
                 enddo            ! n
               endif              ! tmask
              enddo                ! i
              enddo                ! j
              a2D(:,:,n_fmeltt_ai(ns),iblk) = worka(:,:)
           endif
         enddo
        endif

        ! Aerosols
        if (f_faero_atm(1:1) /= 'x') then
           do n=1,n_aero
              call accum_hist_field(n_faero_atm(n,:),iblk, &
                                      faero_atm(:,:,n,iblk), a2D)
           enddo
        endif
        if (f_faero_ocn(1:1) /= 'x') then
           do n=1,n_aero
              call accum_hist_field(n_faero_ocn(n,:),iblk, &
                                      faero_ocn(:,:,n,iblk), a2D)
                                    
           enddo
        endif
        if (f_aero(1:1) /= 'x') then
           do n=1,n_aero
              call accum_hist_field(n_aerosn1(n,:), iblk, &
                                 trcr(:,:,nt_aero  +4*(n-1),iblk)/rhos, a2D)
              call accum_hist_field(n_aerosn2(n,:), iblk, &
                                 trcr(:,:,nt_aero+1+4*(n-1),iblk)/rhos, a2D)
              call accum_hist_field(n_aeroic1(n,:), iblk, &
                                 trcr(:,:,nt_aero+2+4*(n-1),iblk)/rhoi, a2D)
              call accum_hist_field(n_aeroic2(n,:), iblk, &
                                 trcr(:,:,nt_aero+3+4*(n-1),iblk)/rhoi, a2D)
           enddo
        endif

      !---------------------------------------------------------------
      ! accumulate other history output
      !---------------------------------------------------------------
         ! mechanical redistribution
         call accum_hist_mechred (iblk)
         ! melt ponds
         if (tr_pond) call accum_hist_pond (iblk)

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

           do n = 1, num_avail_hist_fields_2D
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then 

              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a2D(i,j,n,iblk) = spval
                 else                            ! convert units
                    a2D(i,j,n,iblk) = avail_hist_fields(n)%cona*a2D(i,j,n,iblk) &
                                   * ravgct + avail_hist_fields(n)%conb
                 endif
              enddo             ! i
              enddo             ! j

              ! back out albedo/zenith angle dependence
              if (avail_hist_fields(n)%vname(1:6) == 'albice') then
              do j = jlo, jhi
              do i = ilo, ihi
                 if (tmask(i,j,iblk)) then 
                    ravgctz = c0
                    if (albcnt(i,j,iblk,ns) > puny) &
                        ravgctz = c1/albcnt(i,j,iblk,ns)
                    if (f_albice (1:1) /= 'x' .and. n_albice(ns) /= 0) &
                       a2D(i,j,n_albice(ns),iblk) = &
                       a2D(i,j,n_albice(ns),iblk)*avgct(ns)*ravgctz
                    if (f_albsno (1:1) /= 'x' .and. n_albsno(ns) /= 0) &
                       a2D(i,j,n_albsno(ns),iblk) = &
                       a2D(i,j,n_albsno(ns),iblk)*avgct(ns)*ravgctz
                    if (f_albpnd (1:1) /= 'x' .and. n_albpnd(ns) /= 0) &
                       a2D(i,j,n_albpnd(ns),iblk) = &
                       a2D(i,j,n_albpnd(ns),iblk)*avgct(ns)*ravgctz
                 endif
              enddo             ! i
              enddo             ! j
              endif
              if (avail_hist_fields(n)%vname(1:6) == 'albsni') then
              do j = jlo, jhi
              do i = ilo, ihi
                 if (tmask(i,j,iblk)) then 
                    ravgctz = c0
                    if (albcnt(i,j,iblk,ns) > puny) &
                        ravgctz = c1/albcnt(i,j,iblk,ns)
                    if (f_albsni (1:1) /= 'x' .and. n_albsni(ns) /= 0) &
                       a2D(i,j,n_albsni(ns),iblk) = &
                       a2D(i,j,n_albsni(ns),iblk)*avgct(ns)*ravgctz
                 endif
              enddo             ! i
              enddo             ! j
              endif

              endif
           enddo                ! n

           do n = 1, num_avail_hist_fields_3Dc
              nn = n2D + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 
              do k = 1, ncat_hist
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a3Dc(i,j,k,n,iblk) = spval
                 else                            ! convert units
                    a3Dc(i,j,k,n,iblk) = avail_hist_fields(nn)%cona*a3Dc(i,j,k,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! k
              endif
           enddo                ! n

           do n = 1, num_avail_hist_fields_3Dz
              nn = n3Dccum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 
              do k = 1, nzlyr
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a3Dz(i,j,k,n,iblk) = spval
                 else                            ! convert units
                    a3Dz(i,j,k,n,iblk) = avail_hist_fields(nn)%cona*a3Dz(i,j,k,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! k
              endif
           enddo                ! n
           do n = 1, num_avail_hist_fields_3Db
              nn = n3Dzcum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 
              do k = 1, nzblyr
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a3Db(i,j,k,n,iblk) = spval
                 else                            ! convert units
                    a3Db(i,j,k,n,iblk) = avail_hist_fields(nn)%cona*a3Db(i,j,k,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! k
              endif
           enddo                ! n

           do n = 1, num_avail_hist_fields_4Di
              nn = n3Dbcum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 
              do k = 1, nzilyr
              do ic = 1, ncat_hist
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a4Di(i,j,k,ic,n,iblk) = spval
                 else                            ! convert units
                    a4Di(i,j,k,ic,n,iblk) = avail_hist_fields(nn)%cona*a4Di(i,j,k,ic,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! ic
              enddo             ! k
              endif
           enddo                ! n

           do n = 1, num_avail_hist_fields_4Ds
              nn = n4Dicum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 
              do k = 1, nzslyr
              do ic = 1, ncat_hist
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a4Ds(i,j,k,ic,n,iblk) = spval
                 else                            ! convert units
                    a4Ds(i,j,k,ic,n,iblk) = avail_hist_fields(nn)%cona*a4Ds(i,j,k,ic,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! ic
              enddo             ! k
              endif
           enddo                ! n
           do n = 1, num_avail_hist_fields_4Db
              nn = n4Dscum + n
              if (avail_hist_fields(nn)%vhistfreq == histfreq(ns)) then 
              do k = 1, nzblyr
              do ic = 1, ncat_hist
              do j = jlo, jhi
              do i = ilo, ihi
                 if (.not. tmask(i,j,iblk)) then ! mask out land points
                    a4Db(i,j,k,ic,n,iblk) = spval
                 else                            ! convert units
                    a4Db(i,j,k,ic,n,iblk) = avail_hist_fields(nn)%cona*a4Db(i,j,k,ic,n,iblk) &
                                   * ravgct + avail_hist_fields(nn)%conb
                 endif
              enddo             ! i
              enddo             ! j
              enddo             ! ic
              enddo             ! k
              endif
           enddo                ! n

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
                 if (n_divu     (ns) /= 0) a2D(i,j,n_divu(ns),     iblk) = spval
                 if (n_shear    (ns) /= 0) a2D(i,j,n_shear(ns),    iblk) = spval
                 if (n_sig1     (ns) /= 0) a2D(i,j,n_sig1(ns),     iblk) = spval
                 if (n_sig2     (ns) /= 0) a2D(i,j,n_sig2(ns),     iblk) = spval
                 if (n_mlt_onset(ns) /= 0) a2D(i,j,n_mlt_onset(ns),iblk) = spval
                 if (n_frz_onset(ns) /= 0) a2D(i,j,n_frz_onset(ns),iblk) = spval
                 if (n_hisnap   (ns) /= 0) a2D(i,j,n_hisnap(ns),   iblk) = spval
                 if (n_aisnap   (ns) /= 0) a2D(i,j,n_aisnap(ns),   iblk) = spval
                 if (n_trsig    (ns) /= 0) a2D(i,j,n_trsig(ns),    iblk) = spval
                 if (n_iage     (ns) /= 0) a2D(i,j,n_iage(ns),     iblk) = spval
                 if (n_FY       (ns) /= 0) a2D(i,j,n_FY(ns),       iblk) = spval

                 if (n_a11      (ns) /= 0) a2D(i,j,n_a11(ns),      iblk) = spval
                 if (n_a12      (ns) /= 0) a2D(i,j,n_a12(ns),      iblk) = spval
                 if (n_e11      (ns) /= 0) a2D(i,j,n_e11(ns),      iblk) = spval
                 if (n_e12      (ns) /= 0) a2D(i,j,n_e12(ns),      iblk) = spval
                 if (n_e22      (ns) /= 0) a2D(i,j,n_e22(ns),      iblk) = spval
                 if (n_s11      (ns) /= 0) a2D(i,j,n_s11(ns),      iblk) = spval
                 if (n_s12      (ns) /= 0) a2D(i,j,n_s12(ns),      iblk) = spval
                 if (n_s22      (ns) /= 0) a2D(i,j,n_s22(ns),      iblk) = spval
                 if (n_yieldstress11 (ns) /= 0) a2D(i,j,n_yieldstress11(ns),iblk) = spval
                 if (n_yieldstress12 (ns) /= 0) a2D(i,j,n_yieldstress12(ns),iblk) = spval
                 if (n_yieldstress22 (ns) /= 0) a2D(i,j,n_yieldstress22(ns),iblk) = spval
              else
                 if (n_divu     (ns) /= 0) a2D(i,j,n_divu(ns),iblk)      = &
                       divu (i,j,iblk)*avail_hist_fields(n_divu(ns))%cona
                 if (n_shear    (ns) /= 0) a2D(i,j,n_shear(ns),iblk)     = &
                       shear(i,j,iblk)*avail_hist_fields(n_shear(ns))%cona
                 if (n_sig1     (ns) /= 0) a2D(i,j,n_sig1(ns),iblk)      = &
                       sig1 (i,j,iblk)*avail_hist_fields(n_sig1(ns))%cona
                 if (n_sig2     (ns) /= 0) a2D(i,j,n_sig2(ns),iblk)      = &
                       sig2 (i,j,iblk)*avail_hist_fields(n_sig2(ns))%cona
                 if (n_mlt_onset(ns) /= 0) a2D(i,j,n_mlt_onset(ns),iblk) = &
                       mlt_onset(i,j,iblk)
                 if (n_frz_onset(ns) /= 0) a2D(i,j,n_frz_onset(ns),iblk) = &
                       frz_onset(i,j,iblk)
                 if (n_hisnap   (ns) /= 0) a2D(i,j,n_hisnap(ns),iblk)    = &
                       vice(i,j,iblk)
                 if (n_aisnap   (ns) /= 0) a2D(i,j,n_aisnap(ns),iblk)    = &
                       aice(i,j,iblk)

                 if (kdyn == 2) then  ! for EAP dynamics different time of output
                    if (n_trsig    (ns) /= 0) a2D(i,j,n_trsig(ns),iblk ) = &
                                        prs_sig(i,j,iblk)
                 else
                    if (n_trsig    (ns) /= 0) a2D(i,j,n_trsig(ns),iblk ) = &
                                       p25*(stressp_1(i,j,iblk) &
                                          + stressp_2(i,j,iblk) &
                                          + stressp_3(i,j,iblk) &
                                          + stressp_4(i,j,iblk))
                 endif

                 if (n_iage     (ns) /= 0) a2D(i,j,n_iage(ns),iblk)  = &
                       trcr(i,j,nt_iage,iblk)*avail_hist_fields(n_iage(ns))%cona
                 if (n_FY       (ns) /= 0) a2D(i,j,n_FY(ns),iblk)  = &
                       trcr(i,j,nt_FY,iblk)*avail_hist_fields(n_FY(ns))%cona

                 if (n_a11     (ns) /= 0) a2D(i,j,n_a11(ns),iblk)      = &
                       a11 (i,j,iblk)*avail_hist_fields(n_a11(ns))%cona
                 if (n_a12     (ns) /= 0) a2D(i,j,n_a12(ns),iblk)      = &
                       a12 (i,j,iblk)*avail_hist_fields(n_a12(ns))%cona
                 if (n_e11     (ns) /= 0) a2D(i,j,n_e11(ns),iblk)      = &
                       e11 (i,j,iblk)*avail_hist_fields(n_e11(ns))%cona
                 if (n_e12     (ns) /= 0) a2D(i,j,n_e12(ns),iblk)      = &
                       e12 (i,j,iblk)*avail_hist_fields(n_e12(ns))%cona
                 if (n_e22     (ns) /= 0) a2D(i,j,n_e22(ns),iblk)      = &
                       e22 (i,j,iblk)*avail_hist_fields(n_e22(ns))%cona
                 if (n_s11     (ns) /= 0) a2D(i,j,n_s11(ns),iblk)      = &
                       s11 (i,j,iblk)*avail_hist_fields(n_s11(ns))%cona
                 if (n_s12     (ns) /= 0) a2D(i,j,n_s12(ns),iblk)      = &
                       s12 (i,j,iblk)*avail_hist_fields(n_s12(ns))%cona
                 if (n_s22     (ns) /= 0) a2D(i,j,n_s22(ns),iblk)      = &
                       s22 (i,j,iblk)*avail_hist_fields(n_s22(ns))%cona
                 if (n_yieldstress11     (ns) /= 0) a2D(i,j,n_yieldstress11(ns),iblk)      = &
                       yieldstress11 (i,j,iblk)*avail_hist_fields(n_yieldstress11(ns))%cona
                 if (n_yieldstress12     (ns) /= 0) a2D(i,j,n_yieldstress12(ns),iblk)      = &
                       yieldstress12 (i,j,iblk)*avail_hist_fields(n_yieldstress12(ns))%cona
                 if (n_yieldstress22     (ns) /= 0) a2D(i,j,n_yieldstress22(ns),iblk)      = &
                       yieldstress22 (i,j,iblk)*avail_hist_fields(n_yieldstress22(ns))%cona
              endif
           enddo                ! i
           enddo                ! j

        enddo                   ! iblk

        time_end(ns) = time/int(secday)

      !---------------------------------------------------------------
      ! write file
      !---------------------------------------------------------------

      call ice_write_hist (ns)

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

         if ( (month .eq. 7) .and. new_month ) then 

            do j=jlo,jhi
            do i=ilo,ihi

               ! reset SH Jul 1
               if (lmask_s(i,j,iblk)) mlt_onset(i,j,iblk) = c0

               ! reset NH Jul 1
               if (lmask_n(i,j,iblk)) frz_onset(i,j,iblk) = c0
            enddo
            enddo

         endif                  ! 1st of July
      enddo                     ! iblk

      end subroutine accum_hist

!=======================================================================

      end module ice_history

!=======================================================================
