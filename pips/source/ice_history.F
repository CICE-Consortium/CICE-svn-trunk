c $Id: ice_history.F,v 1.21 2004/02/25 17:35:13 eclare Exp $ 
!=======================================================================
!
!BOP
!
! !MODULE: ice_history - ice model history and restart files
!
! !DESCRIPTION:
!
! Output files: netCDF data, Fortran unformatted dumps
!
! !REVISION HISTORY:
!
! authors Tony Craig, NCAR
!         Elizabeth C. Hunke, LANL
!         C. M. Bitz, UW
!         Bruce P. Briegleb, NCAR
!         William H. Lipscomb, LANL
!
! !INTERFACE:
!
      module ice_history
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain
      use ice_read_write
      use ice_fileunits
      use ice_work, only: work_g1, work_gr, worka
!
!EOP
!
      implicit none

      logical (kind=log_kind) ::
     &   restart ! if true, initialize using restart file instead of defaults
     &,  hist_avg  ! if true, write averaged data instead of snapshots

      character (len=char_len) :: 
     &   dump_dir          ! output directory for restart dump
     &,  dump_file         ! output file for restart dump
     &,  history_file      ! output file for history 
     &,  incond_file       ! output file for snapshot initial conditions
     &,  runtype           ! initial, continue, or branch

      character (len=char_len_long) :: 
     &   restart_dir         ! directory name for restart dump
     &,  history_dir         ! directory name for history file
     &,  incond_dir          ! directory for snapshot initial conditions

      character (len=char_len_long) :: 
     &   runid             ! identifier for the run
     &,  pointer_file      ! input pointer file for restarts

      !---------------------------------------------------------------
      ! Instructions for adding a field:
      !     Here:
      ! (1) Increase avgsiz
      ! (2) Add to logical flags
      ! (3) Add to namelist (here and also in ice_in)
      ! (4) Assign an index number
      !     In init_hist:
      ! (5) Add to list of fields, vname
      ! (6) Add to field descriptions, vdesc
      ! (7) Add to field units, vunit
      ! (8) Load iout array with the logical flag
      ! (9) Specify unit conversion factor if necessary
      ! (10)Increment the field in ice_write_hist
      !
      ! Note: Ice area and volume fields defined for ncat <= 10
      !---------------------------------------------------------------

      !---------------------------------------------------------------
      ! primary info for the history file
      !---------------------------------------------------------------

      integer (kind=int_kind), parameter ::
     &   avgsiz = 91       ! number of fields that can be written

      real (kind=dbl_kind) ::
     &   avgct    ! average sample counter
     &,  aa(ilo:ihi,jlo:jhi,avgsiz) ! field accumulations and averages
     &,  cona(avgsiz)               ! multiplicative conversion factor
     &,  conb(avgsiz)               ! additive conversion factor

      logical (kind=log_kind) ::
     &   iout(avgsiz)      ! true if field is written to output file

      character (len=16) ::
     &   vname(avgsiz)     ! variable names
     &,  vunit(avgsiz)     ! variable units

      character (len=40) :: 
     &   vdesc(avgsiz)     ! variable descriptions

      !---------------------------------------------------------------
      ! logical flags: write to output file if true
      !---------------------------------------------------------------

      logical (kind=log_kind) ::
     &     f_hi,        f_hs
     &,    f_Tsfc,      f_aice
     &,    f_uvel,      f_vvel
     &,    f_fswdn,     f_flwdn
     &,    f_snow,      f_snow_ai
     &,    f_rain,      f_rain_ai
     &,    f_sst,       f_sss
     &,    f_uocn,      f_vocn
     &,    f_frzmlt
     &,    f_fswabs,    f_fswabs_ai    
     &,    f_albsni
     &,    f_flat,      f_flat_ai      
     &,    f_fsens,     f_fsens_ai
     &,    f_flwup,     f_flwup_ai
     &,    f_evap,      f_evap_ai
     &,    f_Tair,      f_Tref
     &,    f_congel,    f_frazil
     &,    f_snoice,    f_meltt
     &,    f_meltb,     f_meltl
     &,    f_fresh,     f_fresh_ai  
     &,    f_fsalt,     f_fsalt_ai    
     &,    f_fhnet,     f_fhnet_ai
     &,    f_fswthru,   f_fswthru_ai
     &,    f_strairx,   f_strairy
     &,    f_strtltx,   f_strtlty
     &,    f_strcorx,   f_strcory
     &,    f_strocnx,   f_strocny
     &,    f_strintx,   f_strinty
     &,    f_strength  
     &,    f_divu,      f_shear
     &,    f_sig1,      f_sig2
     &,    f_dvidtt,    f_dvidtd
     &,    f_daidtt,    f_daidtd
     &,    f_mlt_onset, f_frz_onset
     &,    f_dardg1dt,  f_dardg2dt
     &,    f_dvirdgdt,  f_opening
     &,    f_aice1,     f_aice2
     &,    f_aice3,     f_aice4
     &,    f_aice5,     f_aice6
     &,    f_aice7,     f_aice8
     &,    f_aice9,     f_aice10
     &,    f_vice1,     f_vice2
     &,    f_vice3,     f_vice4
     &,    f_vice5,     f_vice6
     &,    f_vice7,     f_vice8
     &,    f_vice9,     f_vice10
     &,    f_hisnap,    f_aisnap

      !---------------------------------------------------------------
      ! namelist variables (same as logical flags)
      !---------------------------------------------------------------

      namelist / icefields_nml /
     &     f_hi,        f_hs
     &,    f_Tsfc,      f_aice
     &,    f_uvel,      f_vvel
     &,    f_fswdn,     f_flwdn
     &,    f_snow,      f_snow_ai      
     &,    f_rain,      f_rain_ai
     &,    f_sst,       f_sss
     &,    f_uocn,      f_vocn
     &,    f_frzmlt   
     &,    f_fswabs,    f_fswabs_ai    
     &,    f_albsni
     &,    f_flat,      f_flat_ai
     &,    f_fsens,     f_fsens_ai
     &,    f_flwup,     f_flwup_ai
     &,    f_evap,      f_evap_ai
     &,    f_Tair,      f_Tref     
     &,    f_congel,    f_frazil
     &,    f_snoice,    f_meltt
     &,    f_meltb,     f_meltl
     &,    f_fresh,     f_fresh_ai  
     &,    f_fsalt,     f_fsalt_ai
     &,    f_fhnet,     f_fhnet_ai
     &,    f_fswthru,   f_fswthru_ai
     &,    f_strairx,   f_strairy
     &,    f_strtltx,   f_strtlty
     &,    f_strcorx,   f_strcory
     &,    f_strocnx,   f_strocny
     &,    f_strintx,   f_strinty
     &,    f_strength
     &,    f_divu,      f_shear
     &,    f_sig1,      f_sig2
     &,    f_dvidtt,    f_dvidtd
     &,    f_daidtt,    f_daidtd
     &,    f_mlt_onset, f_frz_onset
     &,    f_dardg1dt,  f_dardg2dt
     &,    f_dvirdgdt,  f_opening
     &,    f_aice1,     f_aice2
     &,    f_aice3,     f_aice4
     &,    f_aice5,     f_aice6
     &,    f_aice7,     f_aice8
     &,    f_aice9,     f_aice10
     &,    f_vice1,     f_vice2
     &,    f_vice3,     f_vice4
     &,    f_vice5,     f_vice6
     &,    f_vice7,     f_vice8
     &,    f_vice9,     f_vice10
     &,    f_hisnap,    f_aisnap

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), parameter ::
     &     n_hi         = 1
     &,    n_hs         = 2
     &,    n_Tsfc       = 3
     &,    n_aice       = 4
     &,    n_uvel       = 5
     &,    n_vvel       = 6
     &,    n_fswdn      = 7
     &,    n_flwdn      = 8
     &,    n_snow       = 9
     &,    n_snow_ai    = 10
     &,    n_rain       = 11
     &,    n_rain_ai    = 12
     &,    n_sst        = 13
     &,    n_sss        = 14
     &,    n_uocn       = 15
     &,    n_vocn       = 16
     &,    n_frzmlt     = 17
     &,    n_fswabs     = 18
     &,    n_fswabs_ai  = 19
     &,    n_albsni     = 20
     &,    n_flat       = 21
     &,    n_flat_ai    = 22
     &,    n_fsens      = 23
     &,    n_fsens_ai   = 24
     &,    n_flwup      = 25
     &,    n_flwup_ai   = 26
     &,    n_evap       = 27
     &,    n_evap_ai    = 28
     &,    n_Tair       = 29
     &,    n_Tref       = 30
     &,    n_congel     = 31
     &,    n_frazil     = 32
     &,    n_snoice     = 33
     &,    n_meltt      = 34
     &,    n_meltb      = 35
     &,    n_meltl      = 36
     &,    n_fresh      = 37
     &,    n_fresh_ai   = 38
     &,    n_fsalt      = 39
     &,    n_fsalt_ai   = 40
     &,    n_fhnet      = 41
     &,    n_fhnet_ai   = 42
     &,    n_fswthru    = 43
     &,    n_fswthru_ai = 44
     &,    n_strairx    = 45
     &,    n_strairy    = 46
     &,    n_strtltx    = 47
     &,    n_strtlty    = 48
     &,    n_strcorx    = 49
     &,    n_strcory    = 50
     &,    n_strocnx    = 51
     &,    n_strocny    = 52
     &,    n_strintx    = 53
     &,    n_strinty    = 54
     &,    n_strength   = 55
     &,    n_divu       = 56
     &,    n_shear      = 57
     &,    n_sig1       = 58
     &,    n_sig2       = 59
     &,    n_dvidtt     = 60
     &,    n_dvidtd     = 61
     &,    n_daidtt     = 62
     &,    n_daidtd     = 63
     &,    n_mlt_onset  = 64
     &,    n_dardg1dt   = 65
     &,    n_dardg2dt   = 66
     &,    n_dvirdgdt   = 67
     &,    n_opening    = 68
     &,    n_aice1      = 69
     &,    n_aice2      = 70
     &,    n_aice3      = 71
     &,    n_aice4      = 72
     &,    n_aice5      = 73
     &,    n_aice6      = 74
     &,    n_aice7      = 75
     &,    n_aice8      = 76
     &,    n_aice9      = 77
     &,    n_aice10     = 78
     &,    n_vice1      = 79
     &,    n_vice2      = 80
     &,    n_vice3      = 81
     &,    n_vice4      = 82
     &,    n_vice5      = 83
     &,    n_vice6      = 84
     &,    n_vice7      = 85
     &,    n_vice8      = 86
     &,    n_vice9      = 87
     &,    n_vice10     = 88
     &,    n_hisnap     = 89
     &,    n_aisnap     = 90
     &,    n_frz_onset  = 91

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
      subroutine init_hist
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
!
! !USES:
!
      use ice_constants
      use ice_calendar
      use ice_flux, only: mlt_onset, frz_onset
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: n, k
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag

      !---------------------------------------------------------------
      ! field names
      !---------------------------------------------------------------

      vname(n_hi        ) = 'hi'
      vname(n_hs        ) = 'hs'
      vname(n_Tsfc      ) = 'Tsfc'
      vname(n_aice      ) = 'aice'
      vname(n_uvel      ) = 'uvel'
      vname(n_vvel      ) = 'vvel'
      vname(n_fswdn     ) = 'fswdn'
      vname(n_flwdn     ) = 'flwdn'
      vname(n_snow      ) = 'snow'
      vname(n_snow_ai   ) = 'snow_ai'
      vname(n_rain      ) = 'rain'
      vname(n_rain_ai   ) = 'rain_ai'
      vname(n_sst       ) = 'sst'
      vname(n_sss       ) = 'sss'
      vname(n_uocn      ) = 'uocn'
      vname(n_vocn      ) = 'vocn'
      vname(n_frzmlt    ) = 'frzmlt'
      vname(n_fswabs    ) = 'fswabs'
      vname(n_fswabs_ai ) = 'fswabs_ai'
      vname(n_albsni    ) = 'albsni'  
      vname(n_flat      ) = 'flat'  
      vname(n_flat_ai   ) = 'flat_ai'  
      vname(n_fsens     ) = 'fsens'  
      vname(n_fsens_ai  ) = 'fsens_ai'  
      vname(n_flwup     ) = 'flwup'  
      vname(n_flwup_ai  ) = 'flwup_ai'  
      vname(n_evap      ) = 'evap'  
      vname(n_evap_ai   ) = 'evap_ai'  
      vname(n_Tair      ) = 'Tair'
      vname(n_Tref      ) = 'Tref'
      vname(n_congel    ) = 'congel'  
      vname(n_frazil    ) = 'frazil'  
      vname(n_snoice    ) = 'snoice'  
      vname(n_meltt     ) = 'meltt'  
      vname(n_meltb     ) = 'meltb'  
      vname(n_meltl     ) = 'meltl'  
      vname(n_fresh     ) = 'fresh' 
      vname(n_fresh_ai  ) = 'fresh_ai' 
      vname(n_fsalt     ) = 'fsalt'  
      vname(n_fsalt_ai  ) = 'fsalt_ai'  
      vname(n_fhnet     ) = 'fhnet'  
      vname(n_fhnet_ai  ) = 'fhnet_ai'
      vname(n_fswthru   ) = 'fswthru ' 
      vname(n_fswthru_ai) = 'fswthru_ai'
      vname(n_strairx   ) = 'strairx'
      vname(n_strairy   ) = 'strairy'
      vname(n_strtltx   ) = 'strtltx'
      vname(n_strtlty   ) = 'strtlty'
      vname(n_strcorx   ) = 'strcorx'
      vname(n_strcory   ) = 'strcory'
      vname(n_strocnx   ) = 'strocnx'
      vname(n_strocny   ) = 'strocny '
      vname(n_strintx   ) = 'strintx'
      vname(n_strinty   ) = 'strinty'
      vname(n_strength  ) = 'strength'
      vname(n_divu      ) = 'divu'
      vname(n_shear     ) = 'shear'
      vname(n_sig1      ) = 'sig1'
      vname(n_sig2      ) = 'sig2'
      vname(n_dvidtt    ) = 'dvidtt'
      vname(n_dvidtd    ) = 'dvidtd' 
      vname(n_daidtt    ) = 'daidtt'
      vname(n_daidtd    ) = 'daidtd' 
      vname(n_mlt_onset ) = 'mlt_onset'
      vname(n_frz_onset ) = 'frz_onset'
      vname(n_dardg1dt  ) = 'dardg1dt' 
      vname(n_dardg2dt  ) = 'dardg2dt'
      vname(n_dvirdgdt  ) = 'dvirdgdt'
      vname(n_opening   ) = 'opening'
      vname(n_aice1     ) = 'aice1'
      vname(n_aice2     ) = 'aice2'
      vname(n_aice3     ) = 'aice3'
      vname(n_aice4     ) = 'aice4'
      vname(n_aice5     ) = 'aice5'
      vname(n_aice6     ) = 'aice6'
      vname(n_aice7     ) = 'aice7'
      vname(n_aice8     ) = 'aice8'
      vname(n_aice9     ) = 'aice9'
      vname(n_aice10    ) = 'aice10'
      vname(n_vice1     ) = 'vice1'
      vname(n_vice2     ) = 'vice2'
      vname(n_vice3     ) = 'vice3'
      vname(n_vice4     ) = 'vice4'
      vname(n_vice5     ) = 'vice5'
      vname(n_vice6     ) = 'vice6'
      vname(n_vice7     ) = 'vice7'
      vname(n_vice8     ) = 'vice8'
      vname(n_vice9     ) = 'vice9'
      vname(n_vice10    ) = 'vice10'
      vname(n_hisnap    ) = 'hisnap'
      vname(n_aisnap    ) = 'aisnap'

      !---------------------------------------------------------------
      ! field descriptions
      !---------------------------------------------------------------

      vdesc(n_hi        ) = 'grid box mean ice thickness'  
      vdesc(n_hs        ) = 'grid box mean snow thickness'
      vdesc(n_Tsfc      ) = 'snow/ice surface temperature'      
      vdesc(n_aice      ) = 'ice area  (aggregate)'
      vdesc(n_uvel      ) = 'ice velocity (x)'            
      vdesc(n_vvel      ) = 'ice velocity (y)'        
      vdesc(n_fswdn     ) = 'down solar flux'             
      vdesc(n_flwdn     ) = 'down longwave flux'
      vdesc(n_snow      ) = 'snowfall rate (cpl)'         
      vdesc(n_snow_ai   ) = 'snowfall rate'
      vdesc(n_rain      ) = 'rainfall rate (cpl)'         
      vdesc(n_rain_ai   ) = 'rainfall rate'
      vdesc(n_sst       ) = 'sea surface temperature'     
      vdesc(n_sss       ) = 'sea surface salinity'
      vdesc(n_uocn      ) = 'ocean current (x)'           
      vdesc(n_vocn      ) = 'ocean current (y)'         
      vdesc(n_frzmlt    ) = 'freeze/melt potential'    
      vdesc(n_fswabs    ) = 'absorbed solar flux (cpl)'   
      vdesc(n_fswabs_ai ) = 'absorbed solar flux'      
      vdesc(n_albsni    ) = 'snw/ice broad band albedo'
      vdesc(n_flat      ) = 'latent heat flux (cpl)'      
      vdesc(n_flat_ai   ) = 'latent heat flux'
      vdesc(n_fsens     ) = 'sensible heat flux (cpl)'    
      vdesc(n_fsens_ai  ) = 'sensible heat flux'
      vdesc(n_flwup     ) = 'upward longwave flx (cpl)'   
      vdesc(n_flwup_ai  ) = 'upward longwave flux'
      vdesc(n_evap      ) = 'evaporative water flux (cpl)'
      vdesc(n_evap_ai   ) = 'evaporative water flux'   
      vdesc(n_Tair      ) = 'surface air temperature'     
      vdesc(n_Tref      ) = 'sfc reference temperature'
      vdesc(n_congel    ) = 'congelation ice growth'      
      vdesc(n_frazil    ) = 'frazil ice growth'
      vdesc(n_snoice    ) = 'snow-ice formation'          
      vdesc(n_meltt     ) = 'top ice melt'
      vdesc(n_meltb     ) = 'basal ice melt'              
      vdesc(n_meltl     ) = 'lateral ice melt'            
      vdesc(n_fresh     ) = 'frshwtr flx ice to ocn (cpl)'
      vdesc(n_fresh_ai  ) = 'fresh wtr flx ice to ocn'
      vdesc(n_fsalt     ) = 'salt flux ice to ocn (cpl)'  
      vdesc(n_fsalt_ai  ) = 'salt flux ice to ocean'   
      vdesc(n_fhnet     ) = 'heat flux ice to ocn (cpl)'  
      vdesc(n_fhnet_ai  ) = 'heat flux ice to ocean'        
      vdesc(n_fswthru   ) = 'SW thru ice to ocean (cpl)'  
      vdesc(n_fswthru_ai) = 'SW flux thru ice to ocean'
      vdesc(n_strairx   ) = 'atm/ice stress (x)'          
      vdesc(n_strairy   ) = 'atm/ice stress (y)'
      vdesc(n_strtltx   ) = 'sea sfc tilt stress (x)'     
      vdesc(n_strtlty   ) = 'sea sfc tilt stress (y)'
      vdesc(n_strcorx   ) = 'coriolis stress (x)'         
      vdesc(n_strcory   ) = 'coriolis stress (y)'
      vdesc(n_strocnx   ) = 'ocean/ice stress (x)'        
      vdesc(n_strocny   ) = 'ocean/ice stress (y)' 
      vdesc(n_strintx   ) = 'internal ice stress (x)'     
      vdesc(n_strinty   ) = 'internal ice stress (y)'  
      vdesc(n_strength  ) = 'compressive ice strength' 
      vdesc(n_divu      ) = 'strain rate (divergence)'    
      vdesc(n_shear     ) = 'strain rate (shear)'     
      vdesc(n_sig1      ) = 'norm. principal stress 1'    
      vdesc(n_sig2      ) = 'norm. principal stress 2' 
      vdesc(n_dvidtt    ) = 'volume tendency thermo'      
      vdesc(n_dvidtd    ) = 'volume tendency dynamics' 
      vdesc(n_daidtt    ) = 'area tendency thermo'        
      vdesc(n_daidtd    ) = 'area tendency dynamics' 
      vdesc(n_mlt_onset ) = 'melt onset date'            
      vdesc(n_dardg1dt  ) = 'ice area ridging rate'       
      vdesc(n_dardg2dt  ) = 'ridge area formation rate'
      vdesc(n_dvirdgdt  ) = 'ice volume ridging rate'     
      vdesc(n_opening   ) = 'lead area opening rate'
      vdesc(n_aice1     ) = 'ice area (category 1)'
      vdesc(n_aice2     ) = 'ice area (category 2)'
      vdesc(n_aice3     ) = 'ice area (category 3)'
      vdesc(n_aice4     ) = 'ice area (category 4)'
      vdesc(n_aice5     ) = 'ice area (category 5)'
      vdesc(n_aice6     ) = 'ice area (category 6)'
      vdesc(n_aice7     ) = 'ice area (category 7)'
      vdesc(n_aice8     ) = 'ice area (category 8)'
      vdesc(n_aice9     ) = 'ice area (category 9)'
      vdesc(n_aice10    ) = 'ice area (category 10)'
      vdesc(n_vice1     ) = 'ice volume (category 1)'
      vdesc(n_vice2     ) = 'ice volume (category 2)'
      vdesc(n_vice3     ) = 'ice volume (category 3)'
      vdesc(n_vice4     ) = 'ice volume (category 4)'
      vdesc(n_vice5     ) = 'ice volume (category 5)'
      vdesc(n_vice6     ) = 'ice volume (category 6)'
      vdesc(n_vice7     ) = 'ice volume (category 7)'
      vdesc(n_vice8     ) = 'ice volume (category 8)'
      vdesc(n_vice9     ) = 'ice volume (category 9)'
      vdesc(n_vice10    ) = 'ice volume (category 10)'
      vdesc(n_hisnap    ) = 'ice volume snapshot'         
      vdesc(n_aisnap    ) = 'ice area snapshot' 
      vdesc(n_frz_onset ) = 'freeze onset date'

      !---------------------------------------------------------------
      ! field units
      !---------------------------------------------------------------

      vunit(n_hi        ) = 'm'
      vunit(n_hs        ) = 'm'
      vunit(n_Tsfc      ) = 'C'
      vunit(n_aice      ) = ' '
      vunit(n_uvel      ) = 'm/s'
      vunit(n_vvel      ) = 'm/s'
ccsm      vunit(n_aice      ) = '%'
ccsm      vunit(n_uvel      ) = 'cm/s'
ccsm      vunit(n_vvel      ) = 'cm/s'
      vunit(n_fswdn     ) = 'W/m^2'
      vunit(n_flwdn     ) = 'W/m^2'
      vunit(n_snow      ) = 'cm/day'
      vunit(n_snow_ai   ) = 'cm/day'
      vunit(n_rain      ) = 'cm/day'
      vunit(n_rain_ai   ) = 'cm/day'
      vunit(n_sst       ) = 'C'
      vunit(n_sss       ) = 'psu'
      vunit(n_uocn      ) = 'm/s'
      vunit(n_vocn      ) = 'm/s'
ccsm      vunit(n_uocn      ) = 'cm/s'
ccsm      vunit(n_vocn      ) = 'cm/s'
      vunit(n_frzmlt    ) = 'W/m^2'
      vunit(n_fswabs    ) = 'W/m^2'
      vunit(n_fswabs_ai ) = 'W/m^2'
      vunit(n_albsni    ) = '%'
      vunit(n_flat      ) = 'W/m^2'
      vunit(n_flat_ai   ) = 'W/m^2'
      vunit(n_fsens     ) = 'W/m^2'
      vunit(n_fsens_ai  ) = 'W/m^2'
      vunit(n_flwup     ) = 'W/m^2'
      vunit(n_flwup_ai  ) = 'W/m^2'
      vunit(n_evap      ) = 'cm/day'
      vunit(n_evap_ai   ) = 'cm/day'
      vunit(n_Tair      ) = 'C'
      vunit(n_Tref      ) = 'C'
      vunit(n_congel    ) = 'cm/day'
      vunit(n_frazil    ) = 'cm/day'
      vunit(n_snoice    ) = 'cm/day'
      vunit(n_meltt     ) = 'cm/day'
      vunit(n_meltb     ) = 'cm/day'
      vunit(n_meltl     ) = 'cm/day'
      vunit(n_fresh     ) = 'cm/day'
      vunit(n_fresh_ai  ) = 'cm/day'
      vunit(n_fsalt     ) = 'kg/m^2 s'
      vunit(n_fsalt_ai  ) = 'kg/m^2 s'
      vunit(n_fhnet     ) = 'W/m^2'
      vunit(n_fhnet_ai  ) = 'W/m^2'
      vunit(n_fswthru   ) = 'W/m^2'
      vunit(n_fswthru_ai) = 'W/m^2'
      vunit(n_strairx   ) = 'N/m^2'
      vunit(n_strairy   ) = 'N/m^2'
      vunit(n_strtltx   ) = 'N/m^2'
      vunit(n_strtlty   ) = 'N/m^2'
      vunit(n_strcorx   ) = 'N/m^2'
      vunit(n_strcory   ) = 'N/m^2'
      vunit(n_strocnx   ) = 'N/m^2'
      vunit(n_strocny   ) = 'N/m^2'
      vunit(n_strintx   ) = 'N/m^2'
      vunit(n_strinty   ) = 'N/m^2'
      vunit(n_strength  ) = 'N/m'
      vunit(n_divu      ) = '%/day'
      vunit(n_shear     ) = '%/day'
      vunit(n_sig1      ) = ' '
      vunit(n_sig2      ) = ' '
      vunit(n_dvidtt    ) = 'cm/day'
      vunit(n_dvidtd    ) = 'cm/day'
      vunit(n_daidtt    ) = '%/day'
      vunit(n_daidtd    ) = '%/day'
      vunit(n_mlt_onset ) = 'day of year'
      vunit(n_frz_onset ) = 'day of year'
      vunit(n_dardg1dt  ) = '%/day'
      vunit(n_dardg2dt  ) = '%/day'
      vunit(n_dvirdgdt  ) = 'cm/day'
      vunit(n_opening   ) = '%/day'
      vunit(n_aice1     ) = ' '
      vunit(n_aice2     ) = ' '
      vunit(n_aice3     ) = ' '
      vunit(n_aice4     ) = ' '
      vunit(n_aice5     ) = ' '
      vunit(n_aice6     ) = ' '
      vunit(n_aice7     ) = ' '
      vunit(n_aice8     ) = ' '
      vunit(n_aice9     ) = ' '
      vunit(n_aice10    ) = ' '
ccsm      vunit(n_aice1     ) = '%'
ccsm      vunit(n_aice2     ) = '%'
ccsm      vunit(n_aice3     ) = '%'
ccsm      vunit(n_aice4     ) = '%'
ccsm      vunit(n_aice5     ) = '%'
ccsm      vunit(n_aice6     ) = '%'
ccsm      vunit(n_aice7     ) = '%'
ccsm      vunit(n_aice8     ) = '%'
ccsm      vunit(n_aice9     ) = '%'
ccsm      vunit(n_aice10    ) = '%'
      vunit(n_vice1     ) = 'm'
      vunit(n_vice2     ) = 'm'
      vunit(n_vice3     ) = 'm'
      vunit(n_vice4     ) = 'm'
      vunit(n_vice5     ) = 'm'
      vunit(n_vice6     ) = 'm'
      vunit(n_vice7     ) = 'm'
      vunit(n_vice8     ) = 'm'
      vunit(n_vice9     ) = 'm'
      vunit(n_vice10    ) = 'm'
      vunit(n_hisnap    ) = 'm'
      vunit(n_aisnap    ) = ' ' 
ccsm      vunit(n_aisnap    ) = '%' 

      !-----------------------------------------------------------------
      ! read namelist
      !-----------------------------------------------------------------
      nml_error = -1
      open (nu_nml, file='ice_in', status='old')
   10 continue       !*** keep reading until right namelist is found
      read(nu_nml, nml=icefields_nml,err=10,end=20)
      close(nu_nml)
      nml_error = 0
   20 continue

      !-----------------------------------------------------------------
      ! fill iout array with namelist values
      !-----------------------------------------------------------------

      iout=.true.  ! all fields are written by default

      iout(n_hi        ) = f_hi
      iout(n_hs        ) = f_hs
      iout(n_Tsfc      ) = f_Tsfc  
      iout(n_aice      ) = f_aice  
      iout(n_uvel      ) = f_uvel  
      iout(n_vvel      ) = f_vvel  
      iout(n_fswdn     ) = f_fswdn  
      iout(n_flwdn     ) = f_flwdn  
      iout(n_snow      ) = f_snow
      iout(n_snow_ai   ) = f_snow_ai
      iout(n_rain      ) = f_rain  
      iout(n_rain_ai   ) = f_rain_ai  
      iout(n_sst       ) = f_sst  
      iout(n_sss       ) = f_sss  
      iout(n_uocn      ) = f_uocn  
      iout(n_vocn      ) = f_vocn  
      iout(n_frzmlt    ) = f_frzmlt  
      iout(n_fswabs    ) = f_fswabs  
      iout(n_fswabs_ai ) = f_fswabs_ai  
      iout(n_albsni    ) = f_albsni  
      iout(n_flat      ) = f_flat  
      iout(n_flat_ai   ) = f_flat_ai  
      iout(n_fsens     ) = f_fsens  
      iout(n_fsens_ai  ) = f_fsens_ai  
      iout(n_flwup     ) = f_flwup  
      iout(n_flwup_ai  ) = f_flwup_ai  
      iout(n_evap      ) = f_evap  
      iout(n_evap_ai   ) = f_evap_ai  
      iout(n_Tair      ) = f_Tair
      iout(n_Tref      ) = f_Tref
      iout(n_congel    ) = f_congel  
      iout(n_frazil    ) = f_frazil  
      iout(n_snoice    ) = f_snoice  
      iout(n_meltt     ) = f_meltt  
      iout(n_meltb     ) = f_meltb  
      iout(n_meltl     ) = f_meltl  
      iout(n_fresh     ) = f_fresh 
      iout(n_fresh_ai  ) = f_fresh_ai 
      iout(n_fsalt     ) = f_fsalt  
      iout(n_fsalt_ai  ) = f_fsalt_ai  
      iout(n_fhnet     ) = f_fhnet  
      iout(n_fhnet_ai  ) = f_fhnet_ai  
      iout(n_fswthru   ) = f_fswthru  
      iout(n_fswthru_ai) = f_fswthru_ai  
      iout(n_strairx   ) = f_strairx  
      iout(n_strairy   ) = f_strairy  
      iout(n_strtltx   ) = f_strtltx  
      iout(n_strtlty   ) = f_strtlty  
      iout(n_strcorx   ) = f_strcorx  
      iout(n_strcory   ) = f_strcory  
      iout(n_strocnx   ) = f_strocnx
      iout(n_strocny   ) = f_strocny 
      iout(n_strintx   ) = f_strintx  
      iout(n_strinty   ) = f_strinty  
      iout(n_strength  ) = f_strength
      iout(n_divu      ) = f_divu  
      iout(n_shear     ) = f_shear  
      iout(n_sig1      ) = f_sig1  
      iout(n_sig2      ) = f_sig2  
      iout(n_dvidtt    ) = f_dvidtt
      iout(n_dvidtd    ) = f_dvidtd  
      iout(n_daidtt    ) = f_daidtt
      iout(n_daidtd    ) = f_daidtd  
      iout(n_mlt_onset ) = f_mlt_onset  
      iout(n_dardg1dt  ) = f_dardg1dt  
      iout(n_dardg2dt  ) = f_dardg2dt  
      iout(n_dvirdgdt  ) = f_dvirdgdt 
      iout(n_opening   ) = f_opening
      iout(n_aice1     ) = f_aice1  
      iout(n_aice2     ) = f_aice2  
      iout(n_aice3     ) = f_aice3  
      iout(n_aice4     ) = f_aice4  
      iout(n_aice5     ) = f_aice5  
      iout(n_aice6     ) = f_aice6  
      iout(n_aice7     ) = f_aice7  
      iout(n_aice8     ) = f_aice8  
      iout(n_aice9     ) = f_aice9  
      iout(n_aice10    ) = f_aice10  
      iout(n_vice1     ) = f_vice1  
      iout(n_vice2     ) = f_vice2  
      iout(n_vice3     ) = f_vice3  
      iout(n_vice4     ) = f_vice4  
      iout(n_vice5     ) = f_vice5  
      iout(n_vice6     ) = f_vice6  
      iout(n_vice7     ) = f_vice7  
      iout(n_vice8     ) = f_vice8  
      iout(n_vice9     ) = f_vice9  
      iout(n_vice10    ) = f_vice10  
      iout(n_hisnap    ) = f_hisnap
      iout(n_aisnap    ) = f_aisnap
      iout(n_frz_onset ) = f_frz_onset  

      if (my_task.eq.master_task) then
        write(nu_diag,*) ''
        write(nu_diag,*) 'The following variables will be ',
     &                   'written to the history tape: '
        write(nu_diag,*) ' description                           units',
     &       '     netcdf variable'
         do n=1,avgsiz
            if (iout(n)) write(nu_diag,100) vdesc(n), vunit(n), vname(n)
         enddo
         write(nu_diag,*) ' '
      endif
  100 format (1x,a40,2x,a10,2x,a10)

      !-----------------------------------------------------------------
      ! initialize the history arrays
      !-----------------------------------------------------------------
      aa = c0
      avgct = c0

      do k=1,avgsiz
         cona(k) = c1   ! multiply by 1.
         conb(k) = c0   ! add 0.
      enddo

      if (restart .and. yday .ge. c2) then
         mlt_onset = 999. ! restarting midyear gives erroneous onset dates
         frz_onset = 999. ! restarting midyear gives erroneous onset dates
      else
         mlt_onset = c0
         frz_onset = c0
      endif

      !---------------------------------------------------------------
      ! set conversion factors   
      !---------------------------------------------------------------

!ccsm      cona(n_aice  ) = c100                   ! aice  fraction to %
!ccsm      cona(n_aice1 ) = c100                   ! aice1 fraction to %
!ccsm      cona(n_aice2 ) = c100                   ! aice2 fraction to %
!ccsm      cona(n_aice3 ) = c100                   ! aice3 fraction to %
!ccsm      cona(n_aice4 ) = c100                   ! aice4 fraction to %
!ccsm      cona(n_aice5 ) = c100                   ! aice5 fraction to %
!ccsm      cona(n_aice6 ) = c100                   ! aice6 fraction to %
!ccsm      cona(n_aice7 ) = c100                   ! aice7 fraction to %
!ccsm      cona(n_aice8 ) = c100                   ! aice8 fraction to %
!ccsm      cona(n_aice9 ) = c100                   ! aice9 fraction to %
!ccsm      cona(n_aice10) = c100                   ! aice10 fraction to %
!ccsm      cona(n_aisnap) = c100                   ! aisnap fraction to %
!ccsm      cona(n_uvel  ) = m_to_cm                ! u m/s to cm/s
!ccsm      cona(n_vvel  ) = m_to_cm                ! v m/s to cm/s

      cona(n_snow   ) = mps_to_cmpdy/rhofresh  ! snow kg/m2/s to cm/day
      cona(n_snow_ai) = mps_to_cmpdy/rhofresh  ! snow kg/m2/s to cm/day
      cona(n_rain   ) = mps_to_cmpdy/rhofresh  ! rain kg/m2/s to cm/day
      cona(n_rain_ai) = mps_to_cmpdy/rhofresh  ! rain kg/m2/s to cm/day

!ccsm      cona(n_uocn    ) = m_to_cm                ! uocn m/s to cm/s
!ccsm      cona(n_vocn    ) = m_to_cm                ! vocn m/s to cm/s

      cona(n_albsni ) = c100              ! avg of spectral albedos to %
      cona(n_evap   ) = mps_to_cmpdy/rhofresh   ! evap kg/m2/s to cm/day 
      cona(n_evap_ai) = mps_to_cmpdy/rhofresh   ! evap kg/m2/s to cm/day 
      conb(n_Tair   ) = -tffresh                ! Tair K to C
      conb(n_Tref   ) = -tffresh                ! Tref K to C

      cona(n_congel ) = mps_to_cmpdy/dt  ! congel m per step to cm/day
      cona(n_frazil ) = mps_to_cmpdy/dt  ! frazil m per step to cm/day
      cona(n_snoice ) = mps_to_cmpdy/dt  ! snoice m per step to cm/day
      cona(n_meltt  ) = mps_to_cmpdy/dt  ! meltt  m per step to cm/day
      cona(n_meltb  ) = mps_to_cmpdy/dt  ! meltb  m per step to cm/day
      cona(n_meltl  ) = mps_to_cmpdy/dt  ! meltl  m per step to cm/day
      cona(n_fresh  ) = mps_to_cmpdy/rhofresh ! frshwtr flx kg/m2/s to cm/day 
      cona(n_fresh_ai)= mps_to_cmpdy/rhofresh ! frshwtr flx kg/m2/s to cm/day 

      cona(n_divu  ) = secday*c100      ! divu from 1/s to %/day
      cona(n_shear ) = secday*c100      ! shear from 1/s to %/day

      cona(n_dvidtt) = mps_to_cmpdy     ! dvidtt m/s to cm/day
      cona(n_dvidtd) = mps_to_cmpdy     ! dvidtd m/s to cm/day
      cona(n_daidtt) = secday*c100      ! daidtt frac/s to %/day
      cona(n_daidtd) = secday*c100      ! daidtd frac/s to %/day

      cona(n_dardg1dt) = secday*c100    ! dardg1dt frac/s to %/day
      cona(n_dardg2dt) = secday*c100    ! dardg2dt frac/s to %/day
      cona(n_dvirdgdt) = mps_to_cmpdy   ! dvirdgdt m/s to cm/day
      cona(n_opening)  = secday*c100    ! opening  frac/s to %/day

      end subroutine init_hist

!=======================================================================
!
!BOP
!
! !IROUTINE: ice_write_hist - write average ice quantities or snapshots
!
! !INTERFACE:
!
      subroutine ice_write_hist
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
      use ice_flux
      use ice_albedo
      use ice_mechred_cice
      use ice_grid
      use ice_calendar
      use ice_state
      use ice_dyn_evp
      use ice_constants
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i,j,k,n

      real (kind=dbl_kind) :: 
     &     ravgct               ! 1/avgct
     &,    ai                   ! aice_init

      real (kind=dbl_kind), parameter ::
     &     big = 1.0e30                ! big number to mask out land 

      ! ice vol. tendency for history, due to dynamics
      do j=jlo,jhi
        do i=ilo,ihi
           dvidtd(i,j) = (vice(i,j)-dvidtd(i,j))/dt
           daidtd(i,j) = (aice(i,j)-daidtd(i,j))/dt
        enddo
      enddo

      !---------------------------------------------------------------
      ! increment step counter
      !---------------------------------------------------------------

      if (.not. hist_avg) then  ! write snapshots
        aa = c0
        avgct = c1
      else                      ! write averages over time histfreq
        avgct = avgct + c1
      endif

      !---------------------------------------------------------------
      ! increment fields
      !---------------------------------------------------------------

!lipscomb  Split into multiple loops for improved performance.

      do j=jlo,jhi
       do i=ilo,ihi
        if (tmask(i,j)) then  
         ai = aice_init(i,j)

         aa(i,j,n_hi      ) = aa(i,j,n_hi      ) + vice(i,j)
         aa(i,j,n_hs      ) = aa(i,j,n_hs      ) + vsno(i,j)
         aa(i,j,n_Tsfc    ) = aa(i,j,n_Tsfc    ) + Tsfc(i,j)
         aa(i,j,n_aice    ) = aa(i,j,n_aice    ) + aice(i,j)

         aa(i,j,n_uvel    ) = aa(i,j,n_uvel    ) + uvel(i,j)
         aa(i,j,n_vvel    ) = aa(i,j,n_vvel    ) + vvel(i,j)
         
         aa(i,j,n_fswdn   ) = aa(i,j,n_fswdn   ) + fsw(i,j)
         aa(i,j,n_flwdn   ) = aa(i,j,n_flwdn   ) + flw(i,j)
         aa(i,j,n_snow    ) = aa(i,j,n_snow    ) + fsnow(i,j)
         aa(i,j,n_snow_ai ) = aa(i,j,n_snow_ai ) + fsnow(i,j)*ai
         aa(i,j,n_rain    ) = aa(i,j,n_rain    ) + frain(i,j)
         aa(i,j,n_rain_ai ) = aa(i,j,n_rain_ai ) + frain(i,j)*ai
         aa(i,j,n_sst     ) = aa(i,j,n_sst     ) + sst(i,j)
         aa(i,j,n_sss     ) = aa(i,j,n_sss     ) + sss(i,j)
         aa(i,j,n_uocn    ) = aa(i,j,n_uocn    ) + uocn(i,j)
         aa(i,j,n_vocn    ) = aa(i,j,n_vocn    ) + vocn(i,j)
         aa(i,j,n_frzmlt  ) = aa(i,j,n_frzmlt  ) + frzmlt(i,j)
         
         aa(i,j,n_fswabs  ) = aa(i,j,n_fswabs  ) + fswabs(i,j)
         aa(i,j,n_fswabs_ai)= aa(i,j,n_fswabs_ai)+ fswabs(i,j)*ai
         aa(i,j,n_albsni  ) = aa(i,j,n_albsni  ) + alvdr(i,j)*awtvdr
     &                                           + alidr(i,j)*awtidr
     &                                           + alvdf(i,j)*awtvdf
     &                                           + alidf(i,j)*awtidf
         aa(i,j,n_flat    ) = aa(i,j,n_flat    ) + flat(i,j)
         aa(i,j,n_flat_ai ) = aa(i,j,n_flat_ai ) + flat(i,j)*ai
         aa(i,j,n_fsens   ) = aa(i,j,n_fsens   ) + fsens(i,j)
         aa(i,j,n_fsens_ai) = aa(i,j,n_fsens_ai) + fsens(i,j)*ai
         aa(i,j,n_flwup   ) = aa(i,j,n_flwup   ) + flwout(i,j)
         aa(i,j,n_flwup_ai) = aa(i,j,n_flwup_ai) + flwout(i,j)*ai
         aa(i,j,n_evap    ) = aa(i,j,n_evap    ) + evap(i,j)
         aa(i,j,n_evap_ai ) = aa(i,j,n_evap_ai ) + evap(i,j)*ai
         aa(i,j,n_Tair    ) = aa(i,j,n_Tair    ) + Tair(i,j)
         aa(i,j,n_Tref    ) = aa(i,j,n_Tref    ) + Tref(i,j) 

        endif   ! tmask(i,j)
       enddo
      enddo

      do j=jlo,jhi
       do i=ilo,ihi
        if (tmask(i,j)) then  
         ai = aice_init(i,j)

         aa(i,j,n_congel  ) = aa(i,j,n_congel  ) + congel(i,j)
         aa(i,j,n_frazil  ) = aa(i,j,n_frazil  ) + frazil(i,j)
         aa(i,j,n_snoice  ) = aa(i,j,n_snoice  ) + snoice(i,j)
         aa(i,j,n_meltt   ) = aa(i,j,n_meltt   ) + meltt(i,j)
         aa(i,j,n_meltb   ) = aa(i,j,n_meltb   ) + meltb(i,j)
         aa(i,j,n_meltl   ) = aa(i,j,n_meltl   ) + meltl(i,j)
         aa(i,j,n_fresh   ) = aa(i,j,n_fresh   ) + fresh_hist(i,j)
         aa(i,j,n_fresh_ai) = aa(i,j,n_fresh_ai) + fresh_hist(i,j)*ai
         aa(i,j,n_fsalt   ) = aa(i,j,n_fsalt   ) + fsalt_hist(i,j)
         aa(i,j,n_fsalt_ai) = aa(i,j,n_fsalt_ai) + fsalt_hist(i,j)*ai
         aa(i,j,n_fhnet   ) = aa(i,j,n_fhnet   ) + fhnet_hist(i,j)
         aa(i,j,n_fhnet_ai) = aa(i,j,n_fhnet_ai) + fhnet_hist(i,j)*ai
         aa(i,j,n_fswthru ) = aa(i,j,n_fswthru ) + fswthru_hist(i,j)
         aa(i,j,n_fswthru_ai)=aa(i,j,n_fswthru_ai)+fswthru_hist(i,j)*ai

         aa(i,j,n_strairx ) = aa(i,j,n_strairx ) + strairx(i,j)
         aa(i,j,n_strairy ) = aa(i,j,n_strairy ) + strairy(i,j)
         aa(i,j,n_strtltx ) = aa(i,j,n_strtltx ) + strtltx(i,j)
         aa(i,j,n_strtlty ) = aa(i,j,n_strtlty ) + strtlty(i,j)
         aa(i,j,n_strcorx ) = aa(i,j,n_strcorx ) + fm(i,j)*vvel(i,j)
         aa(i,j,n_strcory ) = aa(i,j,n_strcory ) - fm(i,j)*uvel(i,j)
         aa(i,j,n_strocnx ) = aa(i,j,n_strocnx ) + strocnx(i,j)
         aa(i,j,n_strocny ) = aa(i,j,n_strocny ) + strocny(i,j)
         aa(i,j,n_strintx ) = aa(i,j,n_strintx ) + strintx(i,j)
         aa(i,j,n_strinty ) = aa(i,j,n_strinty ) + strinty(i,j)
         aa(i,j,n_strength) = aa(i,j,n_strength) + strength(i,j)

! The following fields (divu, shear, sig1, and sig2) will be smeared 
!  if averaged over more than a few days.
! Snapshots may be more useful (see below).

c         aa(i,j,n_divu    ) = aa(i,j,n_divu    ) + divu(i,j)
c         aa(i,j,n_shear   ) = aa(i,j,n_shear   ) + shear(i,j)
c         aa(i,j,n_sig1    ) = aa(i,j,n_sig1    ) + sig1(i,j)
c         aa(i,j,n_sig2    ) = aa(i,j,n_sig2    ) + sig2(i,j)

         aa(i,j,n_dvidtt  ) = aa(i,j,n_dvidtt  ) + dvidtt(i,j)
         aa(i,j,n_dvidtd  ) = aa(i,j,n_dvidtd  ) + dvidtd(i,j)
         aa(i,j,n_daidtt  ) = aa(i,j,n_daidtt  ) + daidtt(i,j)
         aa(i,j,n_daidtd  ) = aa(i,j,n_daidtd  ) + daidtd(i,j)
         aa(i,j,n_dardg1dt) = aa(i,j,n_dardg1dt) + dardg1dt(i,j)
         aa(i,j,n_dardg2dt) = aa(i,j,n_dardg2dt) + dardg2dt(i,j)
         aa(i,j,n_dvirdgdt) = aa(i,j,n_dvirdgdt) + dvirdgdt(i,j)
         aa(i,j,n_opening ) = aa(i,j,n_opening ) + opening(i,j)

        endif   ! tmask(i,j)
       enddo
      enddo

      do n=1,ncat
       do j=jlo,jhi
        do i=ilo,ihi
         if (tmask(i,j)) then
            ! assume consecutive indices
            aa(i,j,n_aice1+n-1) = aa(i,j,n_aice1+n-1) + aicen(i,j,n)
            aa(i,j,n_vice1+n-1) = aa(i,j,n_vice1+n-1) + vicen(i,j,n)
         endif
        enddo
       enddo
      enddo

      !---------------------------------------------------------------
      ! Write output files at prescribed intervals
      !---------------------------------------------------------------

      if (write_history .or. write_ic) then

        ravgct = c1/avgct
        do k=1,avgsiz
         do j=jlo,jhi
          do i=ilo,ihi
             if (.not. tmask(i,j)) then ! mask out land points
                aa(i,j,k) = big
             else                       ! convert units
                aa(i,j,k) = cona(k)*aa(i,j,k)*ravgct + conb(k) 
             endif
          enddo
         enddo  
        enddo

      ! snapshots
        call principal_stress   ! to compute sig1 and sig2
        do j=jlo,jhi
         do i=ilo,ihi
            if (.not. tmask(i,j)) then ! mask out land points
               aa(i,j,n_divu )     = big
               aa(i,j,n_shear)     = big
               aa(i,j,n_sig1 )     = big
               aa(i,j,n_sig2 )     = big
               aa(i,j,n_mlt_onset) = big
               aa(i,j,n_frz_onset) = big
               aa(i,j,n_hisnap   ) = big
               aa(i,j,n_aisnap   ) = big
            else
               aa(i,j,n_divu)      = divu(i,j)*cona(n_divu)
               aa(i,j,n_shear)     = shear(i,j)*cona(n_shear)
               aa(i,j,n_sig1)      = sig1(i,j)*cona(n_sig1)  
               aa(i,j,n_sig2)      = sig2(i,j)*cona(n_sig2) 
               aa(i,j,n_mlt_onset) = mlt_onset(i,j)
               aa(i,j,n_frz_onset) = frz_onset(i,j)
               aa(i,j,n_hisnap   ) = vice(i,j)
               aa(i,j,n_aisnap   ) = aice(i,j)
            endif
         enddo
        enddo

        call icecdf
        aa = c0
        avgct = c0

      endif  ! write_history or write_ic

      if (new_year) then
        do j=jlo,jhi
          do i=ilo,ihi
            mlt_onset(i,j) = mlt_onset(i,j) * mask_s(i,j) ! reset NH Jan 1
            frz_onset(i,j) = frz_onset(i,j) * mask_n(i,j) ! reset SH Jan 1
          enddo
        enddo
      endif
      if ((yday.gt.181._dbl_kind) .and.
     &    (yday.lt.181._dbl_kind+dt/secday)) then
        do j=jlo,jhi
          do i=ilo,ihi
            mlt_onset(i,j) = mlt_onset(i,j) * mask_n(i,j) ! reset SH Jul 1
            frz_onset(i,j) = frz_onset(i,j) * mask_s(i,j) ! reset NH Jan 1
          enddo
        enddo
      endif

      end subroutine ice_write_hist

!=======================================================================
!
!BOP
!
! !IROUTINE: icecdf - write netCDF history file
!
! !INTERFACE:
!
      subroutine icecdf
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
      use ice_model_size
      use ice_constants
      use ice_mpi_internal
      use ice_grid
      use ice_calendar
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      include "netcdf.inc"

      integer (kind=int_kind) :: i,j,n
     &,  ncid,status,imtid,jmtid,timid,varid
     &,  length
      integer (kind=int_kind), dimension(3) :: dimid,start,count
      real (kind=real_kind) :: ltime
      character (char_len) :: ncfile,title

      integer (kind=int_kind) :: iyear, imonth, iday

      if (my_task.eq.master_task) then
        allocate (work_gr(imt_global,jmt_global))

        ltime=time/int(secday)

        ! construct filename
        if (write_ic) then
           write(ncfile,'(a,a,i8.8,a)') 
     &      history_file(1:lenstr(history_file)),'.',idate,'.incond.nc'
        else

         iyear = nyr
         imonth = month
         iday = mday
         if (new_year) then
           iyear = nyr - 1
           imonth = 12
           iday = 31
         elseif (new_month) then
           imonth = month - 1
           iday = daymo(imonth)
         elseif (new_day) then
           iday = iday - 1
         endif

         if (histfreq.eq.'1') then ! instantaneous, write every dt
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)') 
     &      history_file(1:lenstr(history_file)),'.',iyear,'-',
     &       imonth,'-',iday,'-',sec,'.nc'

         elseif (hist_avg) then    ! write averaged data

          if (histfreq.eq.'d'.or.histfreq.eq.'D') then     ! daily
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a)') 
     &      history_file(1:lenstr(history_file)),'_davg.',iyear,'-',
     &       imonth,'-',iday,'.nc'
          elseif (histfreq.eq.'m'.or.histfreq.eq.'M') then ! monthly
           write(ncfile,'(a,a,i4.4,a,i2.2,a)') 
     &      history_file(1:lenstr(history_file)),'_mavg.',iyear,'-',
     &       imonth,'.nc'
          elseif (histfreq.eq.'w'.or.histfreq.eq.'W') then ! weekly
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a)') 
     &      history_file(1:lenstr(history_file)),'_wavg.',iyear,'-',
     &       imonth,'-',iday,'.nc'
          elseif (histfreq.eq.'y'.or.histfreq.eq.'Y') then ! yearly
           write(ncfile,'(a,a,i4.4,a)')
     &     history_file(1:lenstr(history_file)),'_yavg',iyear,'.nc'
          endif

         else                      ! instantaneous, histfreq > dt
           write(ncfile,'(a,a,i8.8,a)') 
     &      history_file(1:lenstr(history_file)),'.',idate,'.nc'
         endif
        endif

        ! create file
        status = nf_create(ncfile, nf_clobber, ncid)

      !-----------------------------------------------------------------
      ! define dimensions
      !-----------------------------------------------------------------

        status = nf_def_dim(ncid,'lon',imt_global,imtid)
        status = nf_def_dim(ncid,'lat',jmt_global,jmtid)
        status = nf_def_dim(ncid,'time',nf_unlimited,timid)

      !-----------------------------------------------------------------
      ! define coordinate variables
      !-----------------------------------------------------------------

        dimid(1) = imtid
        dimid(2) = jmtid

        status = nf_def_var(ncid,'tlon',nf_float,2,dimid,varid)
        status = nf_put_att_text(ncid,varid,'long_name',
     $                           21,'grid center longitude')
        status = nf_put_att_text(ncid,varid,'units',
     $                           12,'degrees_east')

        status = nf_def_var(ncid,'tlat',nf_float,2,dimid,varid)
        status = nf_put_att_text(ncid,varid,'long_name',
     $                           20,'grid center latitude')
        status = nf_put_att_text(ncid,varid,'units',
     $                           13,'degrees_north')

        dimid(1) = timid
        status = nf_def_var(ncid,'time',nf_float,1,dimid,varid)
        status = nf_put_att_text(ncid,varid,'long_name',
     $                           10,'model time')
        status = nf_put_att_text(ncid,varid,'units',
     $                  30,'days since 0000-00-00 00:00:00')  ! for now
        status = nf_put_att_text(ncid,varid,'calendar',
     $                           6,'noleap')
        if (histfreq .eq. '1' .or. .not.hist_avg) then
        status = nf_put_att_text(ncid,varid,'time_rep',
     $                          13,'instantaneous')
        else
        status = nf_put_att_text(ncid,varid,'time_rep',
     $                           8,'averaged')
        endif

      !-----------------------------------------------------------------
      ! define variables and variable attributes
      !-----------------------------------------------------------------

        dimid(1) = imtid
        dimid(2) = jmtid
        dimid(3) = timid

        status  = nf_def_var(ncid, 'hmask', nf_float, 
     $                       2, dimid, varid)
        status = nf_put_att_text(ncid,varid,
     $             'units',4,'flag')
        status = nf_put_att_text(ncid,varid,
     $             'long_name',15,'ocean grid mask')
        status = nf_put_att_text(ncid,varid,'coordinates',
     $                             5,'tlon tlat')

        status  = nf_def_var(ncid, 'tarea', nf_float, 
     $                       2, dimid, varid)
        status = nf_put_att_text(ncid,varid,
     $             'units',3,'m^2')
        status = nf_put_att_text(ncid,varid,
     $             'long_name',11,'t grid area')
        status = nf_put_att_text(ncid,varid,'coordinates',
     $                             5,'tlon tlat')

        status  = nf_def_var(ncid, 'rotang', nf_float, 
     $                       2, dimid, varid)
        status = nf_put_att_text(ncid,varid,
     $             'units',7,'radians')
        status = nf_put_att_text(ncid,varid,
     $             'long_name',19,'grid rotation angle')
        status = nf_put_att_text(ncid,varid,'coordinates',
     $                             5,'tlon tlat')

        do n=1,avgsiz
          if (iout(n)) then
            status  = nf_def_var(ncid, vname(n), nf_float, 
     $                         3, dimid, varid)
            length = lenstr(vunit(n))
            status = nf_put_att_text(ncid,varid,
     $             'units',length,vunit(n))
            length = lenstr(vdesc(n))
            status = nf_put_att_text(ncid,varid,
     $             'long_name',length,vdesc(n))
            status = nf_put_att_text(ncid,varid,'coordinates',
     $                             8,'i j time')
            status = nf_put_att_real(ncid,varid,'missing_value',
     $                             nf_float,1,spval)

            status = nf_put_att_real(ncid,varid,'_fillvalue',
     $                             nf_float,1,spval)
            if (histfreq .eq. '1' .or. .not.hist_avg
     &          .or. n.eq.n_divu .or. n.eq.n_shear     ! snapshots
     &          .or. n.eq.n_sig1 .or. n.eq.n_sig2 
     &          .or. n.eq.n_mlt_onset .or. n.eq.n_frz_onset
     &          .or. n.eq.n_hisnap .or. n.eq.n_aisnap) then
            status = nf_put_att_text(ncid,varid,'time_rep',
     $                              13,'instantaneous')
            else
            status = nf_put_att_text(ncid,varid,'time_rep',
     $                               8,'averaged')
            endif

          endif
        enddo

      !-----------------------------------------------------------------
      ! global attributes
      !-----------------------------------------------------------------
      ! ... the user should change these to something useful ...
      !-----------------------------------------------------------------
        title  = 'sea ice model output for cice'
        length = lenstr(title)
        status = nf_put_att_text(ncid,nf_global,'title',length,title)
        title = 'sea ice model: LANL sea ice model (cice)'
        length = lenstr(title)
        status = nf_put_att_text(ncid,nf_global,'source',length,title)
        title = 'log file ice.log.'                      ! for now
        length = lenstr(title)
        status = nf_put_att_text(ncid,nf_global,'history',length,title)
        title = 'of LANL cice'
        length = lenstr(title)
        status = 
     &       nf_put_att_text(ncid,nf_global,'conventions',length,title)

      !-----------------------------------------------------------------
      ! end define mode
      !-----------------------------------------------------------------

        status = nf_enddef(ncid)

      !-----------------------------------------------------------------
      ! write coordinate variables
      !-----------------------------------------------------------------

        work_gr= TLON_G*rad_to_deg  ! single precision
        status = nf_inq_varid(ncid,'tlon',varid)
        status = nf_put_var_real(ncid,varid,work_gr)

        work_gr= TLAT_G*rad_to_deg  ! single precision
        status = nf_inq_varid(ncid,'tlat',varid)
        status = nf_put_var_real(ncid,varid,work_gr)

        status = nf_inq_varid(ncid,'time',varid)
        status = nf_put_vara_real(ncid,varid,1,1,ltime)

      endif          ! master_task

      !-----------------------------------------------------------------
      ! write grid mask, area and rotation angle
      !-----------------------------------------------------------------

      worka(:,:) = uvm(ilo:ihi,jlo:jhi)
      call global_gather(work_g1,worka)
      if (my_task.eq.master_task) then
        work_gr=work_g1
        status = nf_inq_varid(ncid,'hmask',varid)
        status = nf_put_var_real(ncid,varid,work_gr)
      endif

      worka(:,:) = tarea(ilo:ihi,jlo:jhi)
      call global_gather(work_g1,worka)
      if (my_task.eq.master_task) then
        work_gr=work_g1
        status = nf_inq_varid(ncid,'tarea',varid)
        status = nf_put_var_real(ncid,varid,work_gr)
      endif

      call global_gather(work_g1,angle)
      if (my_task.eq.master_task) then
        work_gr=work_g1
        status = nf_inq_varid(ncid,'rotang',varid)
        status = nf_put_var_real(ncid,varid,work_gr)
      endif

      !-----------------------------------------------------------------
      ! write variable data 
      !-----------------------------------------------------------------

      start=1
      count(1)=imt_global
      count(2)=jmt_global
      count(3)=1
      do n=1,avgsiz
        if (iout(n)) then
          worka=aa(:,:,n)
          call global_gather(work_g1,worka)
          if (my_task.eq.master_task) then
            work_gr=work_g1
            status = nf_inq_varid(ncid,vname(n),varid)
            status = nf_put_vara_real(ncid,varid,start,count,work_gr)
 
          endif
        endif
      enddo

      if (my_task.eq.master_task) deallocate (work_gr)

      !-----------------------------------------------------------------
      ! close output dataset
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) then
        status = nf_close(ncid)
        write(nu_diag,*) ' '
        write(nu_diag,*) 'Finished writing ',ncfile(1:lenstr(ncfile))
        write(nu_diag,*) ' '
      endif

      end subroutine icecdf

!=======================================================================
!
!BOP
!
! !IROUTINE: integer function lenstr(label) - compute length string
!
! !INTERFACE:
!
      integer function lenstr(label)
!
! !DESCRIPTION:
!
! Compute length of string by finding first non-blank
! character from the right.
!
! !REVISION HISTORY:
!
! author:   ?
!
! !INPUT/OUTPUT PARAMETERS:
!
      character*(*) label
!
!EOP
!
      integer (kind=int_kind) :: 
     &   length  ! length of character string 
     &,  n       ! loop index

      length = len(label)
      do n=length,1,-1
        if( label(n:n) .ne. ' ' ) exit 
      enddo
      lenstr = n

      end function lenstr

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: dumpfile - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine dumpfile
!
! !DESCRIPTION:
!
! Dumps all values needed for a restart 
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_model_size
      use ice_flux
      use ice_grid
      use ice_calendar
      use ice_state
      use ice_dyn_evp
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, k, n
      character(len=80) :: filename
      logical (kind=log_kind) :: gather

      ! write pointer (just filename)
      write(filename,'(a,a,i8.8)') 
     &  dump_file(1:lenstr(dump_file)),'.',idate
      if (my_task.eq.master_task) then
        open(nu_rst_pointer,file=pointer_file)
        write(nu_rst_pointer,'(a80)') filename
        close(nu_rst_pointer)
      endif

      ! write dump (dir/filename)
      write(filename,'(a,a,a,i8.8)') 
     &        dump_dir(1:lenstr(dump_dir)),
     &        dump_file(1:lenstr(dump_file)),'.',idate
      call ice_open(nu_dump,filename,0)
      if (my_task.eq.master_task) then
        write (nu_dump) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
        write(nu_diag,*) 'Restart written ',istep1,time,time_forc
      endif

      gather = .true.

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------
      do n=1,ncat
        worka(ilo:ihi,jlo:jhi)=aicen(ilo:ihi,jlo:jhi,n)
        call ice_write(nu_dump,0,worka,'ruf8',gather)
        worka(ilo:ihi,jlo:jhi)=vicen(ilo:ihi,jlo:jhi,n)
        call ice_write(nu_dump,0,worka,'ruf8',gather)
        worka(ilo:ihi,jlo:jhi)=vsnon(ilo:ihi,jlo:jhi,n)
        call ice_write(nu_dump,0,worka,'ruf8',gather)
        worka(ilo:ihi,jlo:jhi)=Tsfcn(ilo:ihi,jlo:jhi,n)
        call ice_write(nu_dump,0,worka,'ruf8',gather)
        worka(ilo:ihi,jlo:jhi)=esnon(ilo:ihi,jlo:jhi,n)
        call ice_write(nu_dump,0,worka,'ruf8',gather)
      enddo

      do k=1,ntilay
        worka(ilo:ihi,jlo:jhi)=eicen(ilo:ihi,jlo:jhi,k)
        call ice_write(nu_dump,0,worka,'ruf8',gather)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      worka(ilo:ihi,jlo:jhi)=uvel(ilo:ihi,jlo:jhi)
      call ice_write(nu_dump,0,worka,'ruf8',gather)
      worka(ilo:ihi,jlo:jhi)=vvel(ilo:ihi,jlo:jhi)
      call ice_write(nu_dump,0,worka,'ruf8',gather)

      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux
      !-----------------------------------------------------------------
      worka(ilo:ihi,jlo:jhi)=fresh(ilo:ihi,jlo:jhi)      
      call ice_write(nu_dump,0,worka,'ruf8',gather)           
      worka(ilo:ihi,jlo:jhi)=fsalt(ilo:ihi,jlo:jhi)      
      call ice_write(nu_dump,0,worka,'ruf8',gather)          
      worka(ilo:ihi,jlo:jhi)=fhnet(ilo:ihi,jlo:jhi)      
      call ice_write(nu_dump,0,worka,'ruf8',gather)          

      !-----------------------------------------------------------------
      ! ice strength 
      !-----------------------------------------------------------------
      worka(ilo:ihi,jlo:jhi)=strength(ilo:ihi,jlo:jhi)   
      call ice_write(nu_dump,0,worka,'ruf8',gather)           

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      worka(ilo:ihi,jlo:jhi)=strocnxT(ilo:ihi,jlo:jhi)
      call ice_write(nu_dump,0,worka,'ruf8',gather)
      worka(ilo:ihi,jlo:jhi)=strocnyT(ilo:ihi,jlo:jhi)
      call ice_write(nu_dump,0,worka,'ruf8',gather)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
         worka(ilo:ihi,jlo:jhi)=stressp_1(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)
         worka(ilo:ihi,jlo:jhi)=stressm_1(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)
         worka(ilo:ihi,jlo:jhi)=stress12_1(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)

         worka(ilo:ihi,jlo:jhi)=stressp_2(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)
         worka(ilo:ihi,jlo:jhi)=stressm_2(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)
         worka(ilo:ihi,jlo:jhi)=stress12_2(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)

         worka(ilo:ihi,jlo:jhi)=stressp_3(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)
         worka(ilo:ihi,jlo:jhi)=stressm_3(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)
         worka(ilo:ihi,jlo:jhi)=stress12_3(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)

         worka(ilo:ihi,jlo:jhi)=stressp_4(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)
         worka(ilo:ihi,jlo:jhi)=stressm_4(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)
         worka(ilo:ihi,jlo:jhi)=stress12_4(ilo:ihi,jlo:jhi)
         call ice_write(nu_dump,0,worka,'ruf8',gather)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      do j=jlo,jhi
       do i=ilo,ihi
        worka(i,j) = 0.
        if (iceumask(i,j)) worka(i,j) = 1.
       enddo
      enddo
      call ice_write(nu_dump,0,worka,'ruf8',gather)

      ! for mixed layer model
      worka(ilo:ihi,jlo:jhi)=sst(ilo:ihi,jlo:jhi)
      call ice_write(nu_dump,0,worka,'ruf8',gather)
      worka(ilo:ihi,jlo:jhi)=frzmlt(ilo:ihi,jlo:jhi)
      call ice_write(nu_dump,0,worka,'ruf8',gather)

      if (my_task.eq.master_task) close(nu_dump)

      end subroutine dumpfile

!=======================================================================
!BOP
!
! !IROUTINE: restartfile  - restarts from a dumpfile
!
! !INTERFACE:
!
      subroutine restartfile
!
! !DESCRIPTION:
!
! Restarts from a dump
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_model_size
      use ice_flux
      use ice_mpi_internal
      use ice_grid
      use ice_calendar
      use ice_state
      use ice_dyn_evp
      use ice_itd
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, k, n
      character(len=80) :: filename, filename0
      logical (kind=log_kind) :: scatter, diag

      if (my_task.eq.master_task) then
        open(nu_rst_pointer,file=pointer_file)
        read(nu_rst_pointer,'(a80)') filename0
        close(nu_rst_pointer)
        write(filename,'(a,a)') 
     &        dump_dir(1:lenstr(dump_dir)),
     &        filename0(1:lenstr(filename0))
      endif
      call ice_open(nu_restart,filename,0)
      if (my_task.eq.master_task) then
        read (nu_restart) istep0,time,time_forc
        write(nu_diag,*) 'read ',pointer_file(1:lenstr(pointer_file))
        write(nu_diag,*) 'dump read at istep=',istep0,time,time_forc
      endif
      call ice_bcast_iscalar(istep0)
      istep1=istep0
      call ice_bcast_rscalar(time)
      call ice_bcast_rscalar(time_forc)

      scatter = .true.  ! distribute data to MPI subdomains
      diag = .true.     ! write min/max diagnostics for field

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------
      do n=1,ncat
        if (my_task.eq.master_task) 
     &   write(nu_diag,*) 'cat ',n,
     &               ' min/max area, vol ice, vol snow, Tsfc, esnon'
        call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
        aicen(ilo:ihi,jlo:jhi,n)=worka(ilo:ihi,jlo:jhi)
        call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
        vicen(ilo:ihi,jlo:jhi,n)=worka(ilo:ihi,jlo:jhi)
        call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
        vsnon(ilo:ihi,jlo:jhi,n)=worka(ilo:ihi,jlo:jhi)
        call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
        Tsfcn(ilo:ihi,jlo:jhi,n)=worka(ilo:ihi,jlo:jhi)
        call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
        esnon(ilo:ihi,jlo:jhi,n)=worka(ilo:ihi,jlo:jhi)
      enddo

      if (my_task.eq.master_task)
     & write(nu_diag,*) 'min/max eicen for each layer'
      do k=1,ntilay
        call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
        eicen(ilo:ihi,jlo:jhi,k)=worka(ilo:ihi,jlo:jhi)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write(nu_diag,*) 'min/max velocity components'
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      uvel(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      vvel(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux       
      !-----------------------------------------------------------------
      if (my_task.eq.master_task)
     & write(nu_diag,*) 'min/max fresh water and heat flux components'
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      fresh(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      fsalt(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      fhnet(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! ice strength 
      !-----------------------------------------------------------------
      if (my_task.eq.master_task)
     & write(nu_diag,*) 'min/max ice strength'
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      strength(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write(nu_diag,*) 'min/max ocean stress components'
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      strocnxT(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      strocnyT(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) write(nu_diag,*) 
     & 'min/max internal stress components (3 per triangle)'

         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stressp_1(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stressm_1(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stress12_1(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)

         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stressp_2(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stressm_2(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stress12_2(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)

         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stressp_3(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stressm_3(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stress12_3(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)

         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stressp_4(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stressm_4(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
         call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
         stress12_4(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write(nu_diag,*) 'ice mask for dynamics'
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      do j=jlo,jhi
       do i=ilo,ihi
        iceumask(i,j) = .false.
        if (worka(i,j).gt.p5) iceumask(i,j) = .true.
       enddo
      enddo

      ! for mixed layer model
      if (my_task.eq.master_task) 
     & write(nu_diag,*) 'min/max sst, frzmlt'
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      sst(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)
      call ice_read(nu_restart,0,worka,'ruf8',scatter,diag)
      frzmlt(ilo:ihi,jlo:jhi)=worka(ilo:ihi,jlo:jhi)

      if (my_task.eq.master_task) close(nu_restart)

      !-----------------------------------------------------------------
      ! update boundary conditions
      !-----------------------------------------------------------------

      call bound_state          ! aicen, vicen, etc.
      call bound(uvel)      
      call bound(vvel)      

      !-----------------------------------------------------------------
      ! ensure ice is binned in correct categories 
      ! should not be necessary unless restarting from a run with
      ! different category boundaries
      !-----------------------------------------------------------------
      ! These calls do not allow for exact restart capability, so they
      ! have been commented out.  This may be a problem for a run that
      ! uses different category boundaries.
      !-----------------------------------------------------------------
c      call zap_small_areas
c      if (ncat > 1) call rebin

c!!       call cleanup_restart

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------
      call aggregate

      call bound_aggregate

      end subroutine restartfile

!=======================================================================
!BOP
!
! !IROUTINE: cleanup_restart  - fix up restart data
!
! !INTERFACE:
!
      subroutine cleanup_restart
!
! !DESCRIPTION:
!
! Fixes up restart data from incompatible run
! (in this case, ice energy did not have 4 layers in all categories,
! and snow energy did not exist)
! 
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_model_size
      use ice_flux
      use ice_mpi_internal
      use ice_grid
      use ice_calendar
      use ice_state
      use ice_itd
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, k, n
      character(len=80) :: filename
      logical (kind=log_kind) :: scatter

      do j=jlo,jhi
      do i=ilo,ihi
        ! for now, until we get esnon in the restart file
        do n=1,ncat 
          esnon(i,j,n) = -rhos*(Lfresh-cp_ice*min(c0,Tsfcn(i,j,n)))
     &                        *vsnon(i,j,n)
        enddo
        ! categories 1 and 2 have only 2 layers in the restart file
        do k=20,9,-1
          eicen(i,j,k) = eicen(i,j,k-4)
        enddo
        eicen(i,j,8) = eicen(i,j,4)*0.5
        eicen(i,j,7) = eicen(i,j,4)*0.5
        eicen(i,j,6) = eicen(i,j,3)*0.5
        eicen(i,j,5) = eicen(i,j,3)*0.5
        eicen(i,j,4) = eicen(i,j,2)*0.5
        eicen(i,j,3) = eicen(i,j,2)*0.5
        eicen(i,j,2) = eicen(i,j,1)*0.5
        eicen(i,j,1) = eicen(i,j,1)*0.5
      enddo
      enddo
      
      call zap_small_areas
      if (ncat > 1) call rebin

      end subroutine cleanup_restart

!=======================================================================

      end module ice_history

!=======================================================================
