!=======================================================================
!
!BOP
!
! !MODULE: ice_history - ice model history files
!
! Biogeochemistry history output
!
! The following variables are currently hard-wired as snapshots 
!   (instantaneous rather than time-averages):
!   divu, shear, sig1, sig2, trsig, mlt_onset, frz_onset, hisnap, aisnap
!
! The flags (f_<field>) can be set to '1','h','d','m','y' or 'x', where
!   n means the field will not be written.  To output the same field at
!   more than one frequency, for instance monthy and daily, set 
!   f_<field> = 'md'.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_history.F90 569 2013-01-10 15:28:29Z eclare $
!
! authors Tony Craig and Bruce Briegleb, NCAR
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2012 Elizabeth Hunke split code from ice_history.F90
!
! !INTERFACE:
!
      module ice_history_bgc
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: max_nstrm, max_aero, n_aero, nblyr
      use ice_zbgc_public
!
!EOP
!
      implicit none
      private
      public :: init_hist_bgc_2D, init_hist_bgc_3Dc, init_hist_bgc_3Db, &
                init_hist_bgc_4Db, accum_hist_bgc
      save
      
      !---------------------------------------------------------------
      ! flags: write to output file if true or histfreq value
      !---------------------------------------------------------------

      character (len=max_nstrm), public :: &
           f_faero_atm    = 'm', f_faero_ocn  = 'm', &
           f_aero         = 'm', f_aeron      = 'm', &
           f_fNO          = 'm', f_fNO_ai   = 'm', &
           f_fNO_g        = 'm', f_fNO_g_ai = 'm', &
           f_fNH          = 'm', f_fNH_ai   = 'm', &
           f_fNH_g        = 'm', f_fNH_g_ai = 'm', &
           f_fN           = 'm', f_fN_ai   = 'm', &
           f_fN_g         = 'm', f_fN_g_ai = 'm', &
           f_fSil         = 'm', f_fSil_ai   = 'm', &
           f_fSil_g       = 'm', f_fSil_g_ai = 'm', &
           f_bgc_N_sk     = 'x', f_bgc_C_sk= 'x', &
           f_bgc_chl_sk   = 'x', f_bgc_Nit_sk = 'x', &
           f_bgc_Am_sk    = 'x', f_bgc_Sil_sk= 'x', &
           f_bgc_DMSPp_sk = 'x', f_bgc_DMSPd_sk = 'x', &
           f_bgc_DMS_sk   = 'x', f_bgc_Sil_ml   = 'x', & 
           f_bgc_Nit_ml   = 'x', f_bgc_Am_ml = 'x', & 
           f_bgc_DMSP_ml  = 'x', f_bgc_DMS_ml = 'x', & 
           f_upNO         = 'x', f_upNH        = 'x',   & 
           f_zTin         = 'x', f_zphi         = 'x',  &
           f_iDi          = 'x', f_iki           = 'x',    &
           f_bgc_NO       = 'x', &
           f_bgc_N        = 'x', f_bgc_NH       = 'x',    &
           f_bgc_C        = 'x', f_bgc_chl      = 'x',    &
           f_bgc_DMSPp    = 'x', f_bgc_DMSPd    = 'x',    &
           f_bgc_DMS      = 'x', f_bgc_Sil      = 'x',   &
           f_bgc_PON      = 'x', f_bgc_S        = 'x',   &
           f_fbri         = 'x', &
           f_hbri         = 'x', &
           f_growN        = 'x', f_zfswin      = 'x', &
           f_chlnet       = 'x', &
           f_PPnet        = 'x', f_NOnet = 'x', &
           f_grownet      = 'x'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_bgc_nml /     &
           f_faero_atm, f_faero_ocn, &
           f_aero,      f_aeron, &
           f_fNO,       f_fNO_ai , &
           f_fNO_g,     f_fNO_g_ai, &
           f_fNH,       f_fNH_ai, &
           f_fNH_g,     f_fNH_g_ai, &
           f_fN,        f_fN_ai, &
           f_fN_g,      f_fN_g_ai, &
           f_fSil,      f_fSil_ai, &
           f_fSil_g,    f_fSil_g_ai, &
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
           f_hbri, &
           f_growN,    f_zfswin, &
           f_chlnet, &
           f_PPnet, f_NOnet, f_grownet

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer(kind=int_kind), dimension(max_aero,max_nstrm) :: &
           n_faero_atm    , &
           n_faero_ocn    , &
           n_aerosn1      , &
           n_aerosn2      , &
           n_aeroic1      , &
           n_aeroic2

      integer(kind=int_kind), dimension(max_nstrm) :: &
           n_fNO        , n_fNO_ai , &
           n_fNO_g      , n_fNO_g_ai, &
           n_fNH        , n_fNH_ai, &
           n_fNH_g      , n_fNH_g_ai, &
           n_fN         , n_fN_ai, &
           n_fN_g       , n_fN_g_ai, &
           n_fSil       , n_fSil_ai, &
           n_fSil_g     , n_fSil_g_ai, &
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
           n_hbri, &
           n_growN, &
           n_zfswin, &
           n_chlnet, &
           n_PPnet, &
           n_NOnet, &
           n_grownet

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
      subroutine init_hist_bgc_2D
!
! !DESCRIPTION:
!
! Initialize history files
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: nstreams
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_nml, nml_filename, &
          get_fileunit, release_fileunit
      use ice_history_shared, only: tstr2D, tcstr, define_hist_field, &
          vname_in
      use ice_state, only: tr_aero, hbrine
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: n, k, ns
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag
      character (len=3) :: nchar

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
            read(nu_nml, nml=icefields_bgc_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice('ice: error reading icefields_bgc_nml')
      endif

      if (.not. tr_aero) then
         f_faero_atm = 'x'
         f_faero_ocn = 'x'
         f_aero      = 'x' 
         f_aeron     = 'x' ! NOTE not implemented
      endif
      
      if (.not. solve_skl_bgc) then
          f_bgc_N_sk = 'x'
          f_bgc_C_sk = 'x'
          f_bgc_chl_sk = 'x'
          f_bgc_Nit_sk = 'x'
          f_bgc_Am_sk = 'x'
          f_bgc_Sil_sk = 'x'
          f_bgc_DMSPp_sk = 'x'
          f_bgc_DMSPd_sk = 'x'
          f_bgc_DMS_sk = 'x'
          f_bgc_Nit_ml = 'x'
          f_bgc_Am_ml = 'x'
          f_bgc_Sil_ml = 'x'
          f_bgc_DMSP_ml = 'x'
          f_bgc_DMS_ml = 'x'
      endif  !.not. solve_skl_bgc
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
      if (.not. hbrine)  then
              f_fbri  = 'x'
              f_hbri  = 'x'
      endif

      call broadcast_scalar (f_faero_atm, master_task)
      call broadcast_scalar (f_faero_ocn, master_task)
      call broadcast_scalar (f_aero, master_task)
      call broadcast_scalar (f_aeron, master_task)

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
      call broadcast_scalar (f_hbri, master_task)
      call broadcast_scalar (f_growN, master_task)
      call broadcast_scalar (f_zfswin, master_task)
      call broadcast_scalar (f_chlnet, master_task)
      call broadcast_scalar (f_PPnet, master_task)
      call broadcast_scalar (f_NOnet, master_task)
      call broadcast_scalar (f_grownet, master_task)
      call broadcast_scalar (f_upNO, master_task)
      call broadcast_scalar (f_upNH, master_task)

      ! 2D variables
      do ns = 1, nstreams

      ! Aerosols
      if (f_aero(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'aerosnossl', trim(nchar)
            call define_hist_field(n_aerosn1(n,:),vname_in,"kg/kg",   &
                tstr2D, tcstr,"snow ssl aerosol mass","none", c1, c0, &
                ns, f_aero)
            write(vname_in,'(a,a)') 'aerosnoint', trim(nchar)
            call define_hist_field(n_aerosn2(n,:),vname_in,"kg/kg",   &
                tstr2D, tcstr,"snow int aerosol mass","none", c1, c0, &
                ns, f_aero)
            write(vname_in,'(a,a)') 'aeroicessl', trim(nchar)
            call define_hist_field(n_aeroic1(n,:),vname_in,"kg/kg",  &
                tstr2D, tcstr,"ice ssl aerosol mass","none", c1, c0, &
                ns, f_aero)
            write(vname_in,'(a,a)') 'aeroiceint', trim(nchar)
            call define_hist_field(n_aeroic2(n,:),vname_in,"kg/kg",  &
                tstr2D, tcstr,"ice int aerosol mass","none", c1, c0, &
                ns, f_aero)
         enddo
      endif

      if (f_faero_atm(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'faero_atm', trim(nchar)
            call define_hist_field(n_faero_atm(n,:),vname_in,"kg/m^2 s", &
                tstr2D, tcstr,"aerosol deposition rate","none", c1, c0,  &
                ns, f_faero_atm)
         enddo
      endif

      if (f_faero_ocn(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'faero_ocn', trim(nchar)
            call define_hist_field(n_faero_ocn(n,:),vname_in,"kg/m^2 s", &
                tstr2D, tcstr,"aerosol flux to ocean","none", c1, c0,    &
                ns, f_faero_ocn)
         enddo
      endif

      ! Biogeochemistry

      ! skeletal layer tracers
      if (solve_skl_bgc) then
      if (f_bgc_N_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_N_sk,"algal_N","mmol/m^2",tstr2D, tcstr, &
             "ice bottom algae (nitrogen)",                                      &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns, f_bgc_N_sk)
      if (f_bgc_C_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_C_sk,"algal_C","mmol/m^2",tstr2D, tcstr, &
             "ice bottom algae (carbon)",                                        &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns, f_bgc_C_sk)
      if (f_bgc_chl_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_chl_sk,"algal_chl","mmol/m^2?",tstr2D, tcstr, &
             "ice bottom algae (chlorophyll)",                                   &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns, f_bgc_chl_sk)
      if (f_bgc_Nit_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_Nit_sk,"skl_Nit","mmol/m^2",tstr2D, tcstr, &
             "skeletal nutrient (nitrate)",                                      &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns, f_bgc_Nit_sk)
      if (f_bgc_Am_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_Am_sk,"skl_Am","mmol/m^2",tstr2D, tcstr, &
             "skeletal nutrient (ammonia/um)",                                   &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns, f_bgc_Am_sk)
      if (f_bgc_Sil_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_Sil_sk,"skl_Sil","mmol/m^2",tstr2D, tcstr, &
             "skeletal nutrient (silicate)",                                     &
             "skelelal layer: bottom 2-3 cm", c1, c0,                &
             ns, f_bgc_Sil_sk)
      if (f_bgc_Nit_ml(1:1) /= 'x') &
         call define_hist_field(n_bgc_Nit_ml,"ml_Nit","mmol/m^3",tstr2D, tcstr, &
             "mixed layer nutrient (nitrate)",                                   &
             "upper ocean", c1, c0,                &
             ns, f_bgc_Nit_ml)
      if (f_bgc_Am_ml(1:1) /= 'x') &
         call define_hist_field(n_bgc_Am_ml,"ml_Am","mmol/m^3",tstr2D, tcstr, &
             "mixed layer nutrient (ammonia/um)",                                &
             "upper ocean", c1, c0,                &
             ns, f_bgc_Am_ml)
      if (f_bgc_Sil_ml(1:1) /= 'x') &
         call define_hist_field(n_bgc_Sil_ml,"ml_Sil","mmol/m^3",tstr2D, tcstr, &
             "mixed layer nutrient (silicate)",                                  &
             "upper ocean", c1, c0,                &
             ns, f_bgc_Sil_ml)
      if (f_bgc_DMSPp_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_DMSPp_sk,"skl_DMSPp","mmol/m^2",tstr2D, tcstr, &
             "particulate S in algae (DMSPp)",                                   &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns, f_bgc_DMSPp_sk)
      if (f_bgc_DMSPd_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_DMSPd_sk,"skl_DMSPd","mmol/m^2",tstr2D, tcstr, &
             "dissolved skl precursor (DSMPd)",                                  &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns, f_bgc_DMSPd_sk)
      if (f_bgc_DMS_sk(1:1) /= 'x') &
         call define_hist_field(n_bgc_DMS_sk,"skl_DMS","mmol/m^2",tstr2D, tcstr, &
             "dissolved skl trace gas (DMS)",                                    &
             "skeletal layer: bottom 2-3 cm", c1, c0,                &
             ns, f_bgc_DMS_sk)
      if (f_bgc_DMSP_ml(1:1) /= 'x') &
         call define_hist_field(n_bgc_DMSP_ml,"ml_DMSP","mmol/m^3",tstr2D, tcstr, &
             "mixed layer precursor (DMSP)",                                     &
             "upper ocean", c1, c0,                &
             ns, f_bgc_DMSP_ml)

      if (f_bgc_DMS_ml(1:1) /= 'x') &
         call define_hist_field(n_bgc_DMS_ml,"ml_DMS","mmol/m^3",tstr2D, tcstr, &
             "mixed layer trace gas (DMS)",                                      &
             "upper ocean", c1, c0,                &
             ns, f_bgc_DMS_ml) 
      endif   !solve_skl_bgc 

      ! zbgc
      if (f_fNO(1:1) /= 'x') &
         call define_hist_field(n_fNO,"fNO","mmol/m^2/s",tstr2D, tcstr, &
             "nitrate flux ice to ocn (cpl)",                              &
             "if positive, ocean gains nitrate", c1, c0,                   &
             ns, f_fNO)
      
      if (f_fNO_ai(1:1) /= 'x') &
         call define_hist_field(n_fNO_ai,"fNO_ai","mmol/m^2/s",tstr2D, tcstr, &
             "nitrate flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns, f_fNO_ai)
      
      if (f_fNO_g(1:1) /= 'x') &
         call define_hist_field(n_fNO_g,"fNO_g","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage nitrate flux ice to ocn (cpl)",                   &
             "if positive, ocean gains nitrate", c1, c0,                   &
             ns, f_fNO_g)
      
      if (f_fNO_g_ai(1:1) /= 'x') &
         call define_hist_field(n_fNO_g_ai,"fNO_g_ai","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage nitrate flux ice to ocean",                       &
             "weighted by ice area", c1, c0,                                  &
             ns, f_fNO_g_ai)
            
      if (f_fNH(1:1) /= 'x') &
         call define_hist_field(n_fNH,"fNH","mmol/m^2/s",tstr2D, tcstr, &
             "ammonium flux ice to ocn (cpl)",                              &
             "if positive, ocean gains ammonium", c1, c0,                   &
             ns, f_fNH)
      
      if (f_fNH_ai(1:1) /= 'x') &
         call define_hist_field(n_fNH_ai,"fNH_ai","mmol/m^2/s",tstr2D, tcstr, &
             "ammonium flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns, f_fNH_ai)
      
      if (f_fNH_g(1:1) /= 'x') &
         call define_hist_field(n_fNH_g,"fNH_g","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage ammonium flux ice to ocn (cpl)",                  &
             "if positive, ocean gains ammonium", c1, c0,                   &
             ns, f_fNH_g)
      
      if (f_fNH_g_ai(1:1) /= 'x') &
         call define_hist_field(n_fNH_g_ai,"fNH_g_ai","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage ammonium flux ice to ocean",                      &
             "weighted by ice area", c1, c0,                                  &
             ns, f_fNH_g_ai)
                    
      if (f_fN(1:1) /= 'x') &
         call define_hist_field(n_fN,"fN","mmol/m^2/s",tstr2D, tcstr, &
             "algal N flux ice to ocn (cpl)",                              &
             "if positive, ocean gains algal N", c1, c0,                   &
             ns, f_fN)
      
      if (f_fN_ai(1:1) /= 'x') &
         call define_hist_field(n_fN_ai,"fN_ai","mmol/m^2/s",tstr2D, tcstr, &
             "algal N flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns, f_fN_ai)
      
      if (f_fN_g(1:1) /= 'x') &
         call define_hist_field(n_fN_g,"fN_g","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage algal N flux ice to ocn (cpl)",                   &
             "if positive, ocean gains algal N", c1, c0,                   &
             ns, f_fN_g)
      
      if (f_fN_g_ai(1:1) /= 'x') &
         call define_hist_field(n_fN_g_ai,"fN_g_ai","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage algal N flux ice to ocean",                       &
             "weighted by ice area", c1, c0,                                  &
             ns, f_fN_g_ai)
              
      if (f_fSil(1:1) /= 'x') &
         call define_hist_field(n_fSil,"fSil","mmol/m^2/s",tstr2D, tcstr, &
             "silicate flux ice to ocn (cpl)",                              &
             "if positive, ocean gains silicate", c1, c0,                   &
             ns, f_fSil)
      
      if (f_fSil_ai(1:1) /= 'x') &
         call define_hist_field(n_fSil_ai,"fSil_ai","mmol/m^2/s",tstr2D, tcstr, &
             "silicate flux ice to ocean",                                        &
             "weighted by ice area", c1, c0,                                  &
             ns, f_fSil_ai)
      
      if (f_fSil_g(1:1) /= 'x') &
         call define_hist_field(n_fSil_g,"fSil_g","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage silicate flux ice to ocn (cpl)",                  &
             "if positive, ocean gains silicate", c1, c0,                   &
             ns, f_fSil_g)
      
      if (f_fSil_g_ai(1:1) /= 'x') &
         call define_hist_field(n_fSil_g_ai,"fSil_g_ai","kg/m^2/s",tstr2D, tcstr, &
             "Gravity drainage silicate flux ice to ocean",                      &
             "weighted by ice area", c1, c0,                                  &
             ns, f_fSil_g_ai)
            
      if (f_chlnet(1:1) /= 'x') &
         call define_hist_field(n_chlnet,"chl_net","mg chl/m^2",tstr2D, tcstr,        &
             "Net Chlorophyll",                     &
             "weighted by ice area ", c1, c0,       &
             ns, f_chlnet)

      if (f_NOnet(1:1) /= 'x') &
         call define_hist_field(n_NOnet,"NO_net","mmol NO/m^2",tstr2D, tcstr,        &
             "Net Nitrate",                     &
             "weighted by ice area", c1, c0,       &
             ns, f_NOnet)

     ! both skl and zbgc
       
      if (f_PPnet(1:1) /= 'x') &
         call define_hist_field(n_PPnet,"PP_net","mg C/d/m^2",tstr2D, tcstr,        &
             "Net Primary Production",                     &
             "weighted by ice area", secday, c0,       &
             ns, f_PPnet)
      if (f_grownet(1:1) /= 'x') &
         call define_hist_field(n_grownet,"grow_net","/d",tstr2D, tcstr,        &
             "Net specific growth",                     &
             "weighted by ice area", secday, c0,       &
             ns, f_grownet)
      if (f_hbri(1:1) /= 'x') &
         call define_hist_field(n_hbri,"hbrine","m",tstr2D, tcstr,        &
             "Brine height",                     &
             "distance from ice bottom to brine surface", c1, c0,       &
             ns, f_hbri)

      enddo ! nstreams
      
      end subroutine init_hist_bgc_2D

!=======================================================================
!
!BOP
!
! !IROUTINE: init_hist - initialize history files
!
! !INTERFACE:
!
      subroutine init_hist_bgc_3Dc
!
! !DESCRIPTION:
!
! Initialize history files
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: nstreams
      use ice_constants, only: c0, c1
      use ice_history_shared, only: tstr3Dc, tcstr, define_hist_field
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: ns
      
      ! 3D (category) variables must be looped separately
      do ns = 1, nstreams
        if (f_fbri(1:1) /= 'x') &
         call define_hist_field(n_fbri,"fbri","1",tstr3Dc, tcstr,        &
             "ice vol frac. with dynamic sal, cat",                     &
             "none", c1, c0,       &
             ns, f_fbri)
      enddo ! ns

      end subroutine init_hist_bgc_3Dc

!=======================================================================
!
!BOP
!
! !IROUTINE: init_hist - initialize history files
!
! !INTERFACE:
!
      subroutine init_hist_bgc_3Db
!
! !DESCRIPTION:
!
! Initialize history files
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_calendar, only: nstreams
      use ice_constants, only: c0, secday
      use ice_history_shared, only: tstr3Db, tcstr, define_hist_field
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: ns
      
     ! !! 3D (vertical) ice biology variables !

      do ns = 1, nstreams
      
       if (f_upNO(1:1) /= 'x') &
         call define_hist_field(n_upNO,"upNO","mmol/m^3/d",tstr3Db, tcstr, &
             "Algal NO uptake rate",                           &
             "Positive flux is NO to N pool", secday, c0, &
             ns, f_upNO)

       if (f_upNH(1:1) /= 'x') &
         call define_hist_field(n_upNH,"upNH","mmol/m^3/d",tstr3Db, tcstr, &
             "Algal NH uptake rate",                           &
             "Positive flux is NH to N pool", secday, c0,&
             ns, f_upNH)
      enddo   !ns

      end subroutine init_hist_bgc_3Db

!=======================================================================
!
!BOP
!
! !IROUTINE: init_hist - initialize history files
!
! !INTERFACE:
!
      subroutine init_hist_bgc_4Db
!
! !DESCRIPTION:
!
! Initialize history files
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_calendar, only: nstreams
      use ice_constants, only: c0, c1, c100, secday
      use ice_history_shared, only: tstr4Db, tcstr, define_hist_field
!
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: ns
      
       ! biology vertical grid

      do ns = 1, nstreams
      
       if (f_zTin(1:1) /= 'x') &
            call define_hist_field(n_zTin,"zTinb","C",tstr4Db, tcstr, &
                "ice internal temperatures on bio grid", "interpolated to bio grid", c1, c0,  &
                ns, f_zTin)
      
       if (f_zphi(1:1) /= 'x') &
            call define_hist_field(n_zphi,"zphin","%",tstr4Db, tcstr, &
                "porosity", "brine volume fraction", c100, c0, &
                ns, f_zphi)
         
       if (f_iDi(1:1) /= 'x') &
             call define_hist_field(n_iDi,"iDin","m^2/d",tstr4Db, tcstr, &
                "interface diffusivity", "on bio interface grid", secday, c0, &
                ns, f_iDi)
      
       if (f_iki(1:1) /= 'x') &
            call define_hist_field(n_iki,"ikin","mm^2",tstr4Db, tcstr, &
                "permeability", "on bio interface grid", 1.0e6_dbl_kind, c0,&
                ns, f_iki)

       if (f_bgc_NO(1:1) /= 'x') &
            call define_hist_field(n_bgc_NO,"bgc_NO","mmol/m^3",tstr4Db, tcstr, &
                "bulk nitrate ", "on bio grid", c1, c0,     &
                ns, f_bgc_NO)
 
       if (f_bgc_NH(1:1) /= 'x') &
            call define_hist_field(n_bgc_NH,"bgc_NH","mmol/m^3",tstr4Db, tcstr, &
                "bulk ammonia/um ", "on bio grid", c1, c0,  &
                ns, f_bgc_NH)

       if (f_bgc_N(1:1) /= 'x') &
            call define_hist_field(n_bgc_N,"bgc_N","mmol/m^3",tstr4Db, tcstr, &
                "bulk algal N conc. ", "on bio grid", c1, c0, &
                ns, f_bgc_N)
     
       if (f_bgc_C(1:1) /= 'x') &
            call define_hist_field(n_bgc_C,"bgc_C","mmol/m^3",tstr4Db, tcstr, &
                "bulk algal carbon ", "on bio grid", c1, c0, &
                ns, f_bgc_C)

       if (f_bgc_chl(1:1) /= 'x') &
            call define_hist_field(n_bgc_chl,"bgc_chl","mg/m^3",tstr4Db, tcstr, &
                "bulk algal chlorophyll ", "on bio grid", c1, c0,&
                ns, f_bgc_chl)
      
       if (f_bgc_Sil(1:1) /= 'x') &
            call define_hist_field(n_bgc_Sil,"bgc_Sil","mmol/m^3",tstr4Db, tcstr, &
                "bulk silicate ", "on bio grid", c1, c0, &
                ns, f_bgc_Sil)
      
       if (f_bgc_DMSPp(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMSPp,"bgc_DMSPp","mmol/m^3",tstr4Db, tcstr, &
                "bulk algal DMSP ", "on bio grid", c1, c0,&
                ns, f_bgc_DMSPp)
      
       if (f_bgc_DMSPd(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMSPd,"bgc_DMSPd","mmol/m^3",tstr4Db, tcstr, &
                "bulk dissolved DMSP ", "on bio grid", c1, c0, &
                ns, f_bgc_DMSPd)
  
       if (f_bgc_DMS(1:1) /= 'x') &
            call define_hist_field(n_bgc_DMS,"bgc_DMS","mmol/m^3",tstr4Db, tcstr, &
                "bulk DMS gas ", "on bio grid", c1, c0, &
                ns, f_bgc_DMS)
     
       if (f_bgc_PON(1:1) /= 'x') &
            call define_hist_field(n_bgc_PON,"bgc_PON","mmol/m^3",tstr4Db, tcstr, &
                "other bulk nitrogen pool ", "on bio grid valid for (2:nblyr+1)", c1, c0, &
                ns, f_bgc_PON)
 
       if (f_bgc_S(1:1) /= 'x') &
            call define_hist_field(n_bgc_S,"bgc_S","ppt",tstr4Db, tcstr, &
                "bulk salinity", "on bio grid valid for (2:nblyr+1)", c1, c0, &
                ns, f_bgc_S)
     
       if (f_growN(1:1) /= 'x') &
            call define_hist_field(n_growN,"growN","d^-1",tstr4Db, tcstr, &
                "Specific algal growth rate", "on bio grid valid for (2:nblyr+1)", secday , c0,  &
                ns, f_growN)
      
       if (f_zfswin(1:1) /= 'x') &
            call define_hist_field(n_zfswin,"zfswin","W/m^2",tstr4Db, tcstr, &
                "internal ice PAR", "on bio grid", c1, c0, &
                ns, f_zfswin)
    

      enddo  !ns

      end subroutine init_hist_bgc_4Db

!=======================================================================
!
!BOP
!
! !IROUTINE: accum_hist - accumulate average ice quantities or snapshots
!
! !INTERFACE:
!
      subroutine accum_hist_bgc (iblk)
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
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_constants, only: c0, puny
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nblyr_hist
      use ice_flux, only: faero_atm, faero_ocn, sss
      use ice_history_shared, only: n2D, a2D, a3Dc, n3Dccum, a3Db, &
          n4Dscum, a4Db, &
          ncat_hist, accum_hist_field, nzblyr
      use ice_state, only: trcrn, trcr, aicen, vice, vicen, nt_aero, nt_fbri, &
          nt_bgc_N, nt_bgc_C, nt_bgc_chl, nt_bgc_NO, nt_bgc_NH, nt_bgc_Sil, &
          nt_bgc_DMSPd, nt_bgc_DMSPp, nt_bgc_DMS, nt_bgc_PON, nt_bgc_S, &
          nt_bgc_N_sk, nt_bgc_C_sk, nt_bgc_chl_sk, nt_bgc_Nit_sk, &
          nt_bgc_Am_sk, nt_bgc_Sil_sk, nt_bgc_DMSPp_sk, nt_bgc_DMSPd_sk, &
          nt_bgc_DMS_sk, nt_bgc_Nit_ml, nt_bgc_Am_ml, nt_bgc_Sil_ml, &
          nt_bgc_DMSP_ml, nt_bgc_DMS_ml
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           iblk                 ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
         i,j,k,n, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist) :: &
         workz 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist,ncat) :: &
         workzn 

      type (block) :: &
         this_block           ! block information for current block

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

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

         ! skeletal layer bgc
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
                        ocean_bio(:,:,nlt_bgc_NO,iblk), a2D)  
         if (f_bgc_Am_ml(1:1)/= 'x') &
             call accum_hist_field(n_bgc_Am_ml,iblk, &
                        ocean_bio(:,:,nlt_bgc_NH,iblk), a2D)  
         if (f_bgc_Sil_ml(1:1)/= 'x') &
             call accum_hist_field(n_bgc_Sil_ml,iblk, &
                        ocean_bio(:,:,nlt_bgc_Sil,iblk), a2D)  
         if (f_bgc_DMSP_ml(1:1)/= 'x') &
             call accum_hist_field(n_bgc_DMSP_ml,iblk, &
                        ocean_bio(:,:,nlt_bgc_DMSPp,iblk), a2D)  
         if (f_bgc_DMS_ml(1:1)/= 'x') &
             call accum_hist_field(n_bgc_DMS_ml,iblk, &
                        ocean_bio(:,:,nt_bgc_DMS,iblk), a2D)  
 
         ! zbgc
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
         if (f_chlnet  (1:1) /= 'x') &
             call accum_hist_field(n_chlnet, iblk, chl_net(:,:,iblk), a2D)
         if (f_PPnet  (1:1) /= 'x') &
             call accum_hist_field(n_PPnet,   iblk, PP_net(:,:,iblk), a2D)
         if (f_NOnet  (1:1) /= 'x') &
             call accum_hist_field(n_NOnet,   iblk, NO_net(:,:,iblk), a2D)
         if (f_grownet  (1:1) /= 'x') &
             call accum_hist_field(n_grownet, iblk, grow_net(:,:,iblk), a2D)
         if (f_hbri  (1:1) /= 'x') &
             call accum_hist_field(n_hbri, iblk, hbri(:,:,iblk), a2D)

         ! 3D category fields

         if (f_fbri   (1:1) /= 'x') &
             call accum_hist_field(n_fbri-n2D, iblk, ncat_hist, &
                                   trcrn(:,:,nt_fbri,1:ncat_hist,iblk), a3Dc)

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

      end subroutine accum_hist_bgc

!=======================================================================

      end module ice_history_bgc

!=======================================================================
