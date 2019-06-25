c $Id: $ 
c=======================================================================
!---! output files: netCDF data, Fortran unformatted dumps
!---!
!---! authors Tony Craig, Bruce P. Briegleb, NCAR
!---!         Elizabeth C. Hunke, LANL
!---!         Cecilia M. Bitz, UW
c=======================================================================

      module ice_history

      use ice_kinds_mod
      use ice_domain
      use ice_read_write

      implicit none

      logical (kind=log_kind) ::
     &   restart ! if true, initialize using restart file instead of defaults
     &,  hist_avg  ! if true, write averaged data instead of snapshots

      character (len=char_len) :: 
     &   dump_file         ! output file for restart dump
     &,  restrt_file       !  input file for restarting
     &,  history_file      ! output file for history 

      ! extra arrays used only for holding history information
      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   congel    ! basal ice growth         (m)
     &,  frazil    ! frazil ice growth        (m)
     &,  meltt     ! top ice melt             (m)
     &,  meltb     ! basal ice melt           (m)
     &,  sabs      ! ice solar absorbed flx (W/m**2)
     &,  daidtd    ! ice area tendency advect. (s^-1)
     &,  dvidtt    ! ice vol. tendency thermo. (m/s)
     &,  dvidtd    ! ice vol. tendency advect. (m/s)

      !---------------------------------------------------------------
      ! primary info for the history file
      !---------------------------------------------------------------

      integer (kind=int_kind), parameter ::
     &   avgsiz = 52  ! number of fields to be written

      real (kind=dbl_kind) ::
     &   avgct    ! average sample counter
     &,  aa(ilo:ihi,jlo:jhi,avgsiz) ! field accumulations and averages
     &,  cona(avgsiz)
     &,  conb(avgsiz)

      logical (kind=log_kind) ::
     &   iout(avgsiz)     ! logical if field should be written

      character (len=16) ::
     &   vname(avgsiz)     ! variable names
     &,  vunit(avgsiz)     ! variable units
      character (len=40) :: 
     &   vdesc(avgsiz)     ! variable descriptions

      data vname / 
     &     'hi',      'hs',      'Tsfc',    'aice',
     &     'aice1',   'aice2',   'aice3',   'aice4',   'aice5',
     &     'uvel',    'vvel',    'Fswdn',   'Flwdn',   
     &     'snow',    'rain',    
     &     'sst',     'sss',     'uocn',    'vocn',    'Focnht',
     &     'Fswabs',  'albsni',  'Flat',    'Fsens',   
     &     'Flwup',   'evap',    'Tair',    
     &     'congel',  'frazil',  'snoice',  'meltb',
     &     'meltt',              'Ffrshw',  'Fioht',  
     &     'straix',  'straiy',  'strtltx', 'strtlty',     
     &     'strcorx', 'strcory', 'stroix',  'stroiy',    
     &     'strintx', 'strinty', 'strength','divu',    
     &     'shear',              'sig1',    'sig2',    
     &     'dvidtt',  'dvidtd',             'daidtd' /

      data vdesc / 
     &  'grid box mean ice thickness'  
     & ,'grid box mean snow thickness'
     & ,'snow/ice surface temperature'      
     & ,'ice area  (aggregate)'    ,'ice area (category 1)' 
     & ,'ice area (category 2)'    ,'ice area (category 3)' 
     & ,'ice area (category 4)'    ,'ice area (category 5)' 
     & ,'ice velocity (x)'         ,'ice velocity (y)'        
     & ,'down solar flux'          ,'down longwave flux'
     & ,'snow fall rate'           ,'rain fall rate'
     & ,'sea surface temperature'  ,'sea surface salinity'
     & ,'ocean current (x)'        ,'ocean current (y)'         
     & ,'freeze/melt potential'    
     & ,'absorbed solar flx'       ,'snw/ice broad band albedo'
     & ,'latent heat flux'         ,'sensible heat flux'
     & ,'outgoing long wave flux'  ,'evaporative water flux'   
     & ,'air       temperature'
     & ,'congelation ice growth'   ,'frazil ice growth'
     & ,'snow-ice formation'       ,'basal ice melt'
     & ,'top ice melt'            
     & ,'fresh wat flx ice to ocn' ,'heat flux ice to ocean'       
     & ,'atm/ice stress (x)'       ,'atm/ice stress (y)'
     & ,'tilt (x)'                 ,'tilt (y)'
     & ,'coriolis stress (x)'      ,'coriolis stress (y)'
     & ,'ocean/ice stress (x)'     ,'ocean/ice stress (y)' 
     & ,'internal ice stress (x)'  ,'internal ice stress (y)'  
     & ,'compressive ice strength' ,'strain rate (divergence)' 
     & ,'strain rate (shear)'     
     & ,'norm. principal stress 1' ,'norm. principal stress 2' 
     & ,'volume tendency thermo'   ,'volume tendency dynamics' 
     &                             ,'area tendency dynamics' /

      data vunit / 
c     &     'm',       'm',       'C',       '%',
c     &     '%',       '%',       '%',       '%',     '%',
c     &     'cm/s',    'cm/s',    'W/m^2',   'W/m^2',   
     &     'm',       'm',       'C',       ' ',
     &     ' ',       ' ',       ' ',       ' ',     ' ',
     &     ' m/s',    ' m/s',    'W/m^2',   'W/m^2',   
     &     'cm/day',  'cm/day',  
     &     'C',       'psu',     'cm/s',    'cm/s',  'W/m^2',   
     &     'W/m^2',   '%',       'W/m^2',   'W/m^2',   
     &     'W/m^2',   'cm/day',  'C',    
     &     'cm/day',  'cm/day',  'cm/day',  'cm/day',
     &     'cm/day',             'cm/day',  'W/m^2',  
     &     'N/m^2',   'N/m^2',   'N/m^2',   'N/m^2',     
     &     'N/m^2',   'N/m^2',   'N/m^2',   'N/m^2',     
     &     'N/m^2',   'N/m^2',   'N/m',     '%/day',    
     &     '%/day',              ' ',       ' ',    
     &     'cm/day',  'cm/day',             'cm/day' /

      logical
     &     f_hi,      f_hs,      f_Tsfc,    f_aice,
     &     f_aice1,   f_aice2,   f_aice3,   f_aice4,   f_aice5,
     &     f_uvel,    f_vvel,    f_Fswdn,   f_Flwdn,   
     &     f_snow,    f_rain,    
     &     f_sst,     f_sss,     f_uocn,    f_vocn,    f_Focnht,  
     &     f_Fswabs,  f_albsni,  f_Flat,    f_Fsens,   
     &     f_Flwup,   f_evap,    f_Tair,    
     &     f_congel,  f_frazil,  f_snoice,  f_meltb,
     &     f_meltt,              f_Ffrshw,  f_Fioht,  
     &     f_straix,  f_straiy,  f_strtltx, f_strtlty,     
     &     f_strcorx, f_strcory, f_stroix,  f_stroiy,    
     &     f_strintx, f_strinty, f_strength,f_divu,    
     &     f_shear,              f_sig1,    f_sig2,    
     &     f_dvidtt,  f_dvidtd,             f_daidtd

      namelist / icefields_nml /
     &     f_hi,      f_hs,      f_Tsfc,    f_aice,
     &     f_aice1,   f_aice2,   f_aice3,   f_aice4,   f_aice5,
     &     f_uvel,    f_vvel,    f_Fswdn,   f_Flwdn,   
     &     f_snow,    f_rain,    
     &     f_sst,     f_sss,     f_uocn,    f_vocn,    f_Focnht,  
     &     f_Fswabs,  f_albsni,  f_Flat,    f_Fsens,   
     &     f_Flwup,   f_evap,    f_Tair,    
     &     f_congel,  f_frazil,  f_snoice,  f_meltb,
     &     f_meltt,              f_Ffrshw,  f_Fioht,  
     &     f_straix,  f_straiy,  f_strtltx, f_strtlty,     
     &     f_strcorx, f_strcory, f_stroix,  f_stroiy,    
     &     f_strintx, f_strinty, f_strength,f_divu,    
     &     f_shear,              f_sig1,    f_sig2,    
     &     f_dvidtt,  f_dvidtd,             f_daidtd

      equivalence (iout(1),f_hi)
      equivalence (iout(2),f_hs)
      equivalence (iout(3),f_Tsfc)
      equivalence (iout(4),f_aice)
      equivalence (iout(5),f_aice1)
      equivalence (iout(6),f_aice2)
      equivalence (iout(7),f_aice3)
      equivalence (iout(8),f_aice4)
      equivalence (iout(9),f_aice5)
      equivalence (iout(10),f_uvel)
      equivalence (iout(11),f_vvel)
      equivalence (iout(12),f_Fswdn)
      equivalence (iout(13),f_Flwdn)
      equivalence (iout(14),f_snow)
      equivalence (iout(15),f_rain)
      equivalence (iout(16),f_sst)
      equivalence (iout(17),f_sss)
      equivalence (iout(18),f_uocn)
      equivalence (iout(19),f_vocn)
      equivalence (iout(20),f_focnht)
      equivalence (iout(21),f_fswabs)
      equivalence (iout(22),f_albsni)
      equivalence (iout(23),f_flat)
      equivalence (iout(24),f_fsens)
      equivalence (iout(25),f_flwup)
      equivalence (iout(26),f_evap)
      equivalence (iout(27),f_tair)
      equivalence (iout(28),f_congel)
      equivalence (iout(29),f_frazil)
      equivalence (iout(30),f_snoice)
      equivalence (iout(31),f_meltb)
      equivalence (iout(32),f_meltt)
      equivalence (iout(33),f_ffrshw)
      equivalence (iout(34),f_fioht)
      equivalence (iout(35),f_straix)
      equivalence (iout(36),f_straiy)
      equivalence (iout(37),f_strtltx)
      equivalence (iout(38),f_strtlty)
      equivalence (iout(39),f_strcorx)
      equivalence (iout(40),f_strcory)
      equivalence (iout(41),f_stroix)
      equivalence (iout(42),f_stroiy)
      equivalence (iout(43),f_strintx)
      equivalence (iout(44),f_strinty)
      equivalence (iout(45),f_strength)
      equivalence (iout(46),f_divu)
      equivalence (iout(47),f_shear)
      equivalence (iout(48),f_sig1)
      equivalence (iout(49),f_sig2)
      equivalence (iout(50),f_dvidtt)
      equivalence (iout(51),f_dvidtd)
      equivalence (iout(52),f_daidtd)

c=======================================================================

      contains

c=======================================================================

      subroutine init_hist

      use ice_constants
      use ice_calendar

      integer (kind=int_kind) :: n, k
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag

      iout=.true.

      !-----------------------------------------------------------------
      ! read namelist
      !-----------------------------------------------------------------
        nml_error = -1
        open (21, file='ice_in', status='old')
   10   continue  !*** keep reading until right namelist is found
        read(21, nml=icefields_nml,err=10,end=20)
        close(21)
        nml_error = 0
   20   continue

      if (my_task.eq.master_task) then
        write (6,*) 'the following variables will be',
     &    ' written to the history tape: '
        write (6,*) ' description                             units',
     &    '     netcdf variable'
        do n=1,avgsiz
           if (iout(n))   write (6,1000) vdesc(n), vunit(n), vname(n)
        enddo
        write (6,*) ' '
        call flush(6)
      endif
1000  format (1x,a40,2x,a10,2x,a10)

      !-----------------------------------------------------------------
      ! initialize the history arrays
      !-----------------------------------------------------------------
      aa = c0
      avgct = c0

      do k=1,avgsiz
        cona(k) = c1
        conb(k) = c0
      enddo

      ! set up modified factors   
c      cona( 4) = c100                    ! ai  fraction to %
c      cona( 5) = c100                    ! ai1 fraction to %
c      cona( 6) = c100                    ! ai2 fraction to %
c      cona( 7) = c100                    ! ai3 fraction to %
c      cona( 8) = c100                    ! ai4 fraction to %
c      cona( 9) = c100                    ! ai5 fraction to %
c      cona(10) = m_to_cm                 ! u m/s to cm/s
c      cona(11) = m_to_cm                 ! v m/s to cm/s

      cona(14) = mps_to_cmpdy/rhofresh   ! snow kg/m2/s to cm/day
      cona(15) = mps_to_cmpdy/rhofresh   ! rain kg/m2/s to cm/day

c      cona(18) = m_to_cm                 ! uocn m/s to cm/s
c      cona(19) = m_to_cm                 ! vocn m/s to cm/s

      cona(22) = c100                    ! avg of spectral 4 albedos to %
      cona(26) = mps_to_cmpdy/rhofresh   ! evap kg/m2/s to cm/day 
      conb(27) = -tffresh                ! tair k to c

      cona(28) = mps_to_cmpdy/dt         ! congel m per step to cm/day
      cona(29) = mps_to_cmpdy/dt         ! frazil m per step to cm/day
      cona(30) = mps_to_cmpdy/dt         ! snoice m per step to cm/day
      cona(31) = mps_to_cmpdy/dt         ! meltb  m per step to cm/day
      cona(32) = mps_to_cmpdy/dt         ! meltt  m per step to cm/day
      cona(33) = mps_to_cmpdy/rhofresh   ! fresh water flux kg/m2/s to cm/day 

      cona(46) = secday*c100             ! divu from 1/s to %/day
      cona(47) = secday*c100             ! shear from 1/s to %/day

      cona(50) = mps_to_cmpdy            ! dvidtt m/s to cm/day
      cona(51) = mps_to_cmpdy            ! dvidtd m/s to cm/day
      cona(52) = secday*c100             ! daidtd frac/s to %/day

      end subroutine init_hist

c=======================================================================

      subroutine ice_write_hist

!---!-------------------------------------------------------------------
!---! write average ice quantities or snapshots
!---!-------------------------------------------------------------------

      use ice_flux
      use ice_albedo
      use ice_mechred_cice
      use ice_grid
      use ice_calendar
      use ice_state
      use ice_dyn_evp
      use ice_itd

      integer (kind=int_kind) :: i,j,k

      ! ice vol. tendency for history, due to dynamics
      do j=jlo,jhi
        do i=ilo,ihi
           dvidtd(i,j) = (vice(i,j)-dvidtd(i,j))/dt
           daidtd(i,j) = (aice(i,j)-daidtd(i,j))/dt
        enddo
      enddo

      if (.not. hist_avg) then  ! write snapshots
        aa = c0
        avgct = c1
      else                      ! write averages over time histfreq
        avgct = avgct + c1
      endif

      do j=jlo,jhi
       do i=ilo,ihi
        if (iout( 1)) aa(i,j, 1) = aa(i,j, 1) + vice(i,j)
        if (iout( 2)) aa(i,j, 2) = aa(i,j, 2) + vsno(i,j)
        if (iout( 3)) aa(i,j, 3) = aa(i,j, 3) + Tsfc(i,j)
        if (iout( 4)) aa(i,j, 4) = aa(i,j, 4) + aice(i,j)

        if (ncat .ge. 5 ) then
          if (iout( 5)) aa(i,j, 5) = aa(i,j, 5) + aicen(i,j,1)
          if (iout( 6)) aa(i,j, 6) = aa(i,j, 6) + aicen(i,j,2)
          if (iout( 7)) aa(i,j, 7) = aa(i,j, 7) + aicen(i,j,3)
          if (iout( 8)) aa(i,j, 8) = aa(i,j, 8) + aicen(i,j,4)
          if (iout( 9)) aa(i,j, 9) = aa(i,j, 9) + aicen(i,j,5)
        else
          if (iout( 5)) aa(i,j, 5) = 0.0
          if (iout( 6)) aa(i,j, 6) = 0.0
          if (iout( 7)) aa(i,j, 7) = 0.0
          if (iout( 8)) aa(i,j, 8) = 0.0
          if (iout( 9)) aa(i,j, 9) = 0.0
        endif 

        if (iout(10)) aa(i,j,10) = aa(i,j,10) + u(i,j)
        if (iout(11)) aa(i,j,11) = aa(i,j,11) + v(i,j)

        if (iout(12)) aa(i,j,12) = aa(i,j,12) + Fsw(i,j)
        if (iout(13)) aa(i,j,13) = aa(i,j,13) + Flw(i,j)
        if (iout(14)) aa(i,j,14) = aa(i,j,14) + Fsnow(i,j)
        if (iout(15)) aa(i,j,15) = aa(i,j,15) + Frain(i,j)

        if (iout(16)) aa(i,j,16) = aa(i,j,16) + sst(i,j)
        if (iout(17)) aa(i,j,17) = aa(i,j,17) + sss(i,j)
        if (iout(18)) aa(i,j,18) = aa(i,j,18) + uocn(i,j)
        if (iout(19)) aa(i,j,19) = aa(i,j,19) + vocn(i,j)
        if (iout(20)) aa(i,j,20) = aa(i,j,20) + frzmlt(i,j)

        if (iout(21)) aa(i,j,21) = aa(i,j,21) + sabs(i,j)
        if (iout(22)) aa(i,j,22) = aa(i,j,22) + 
     &                             alvdr(i,j)*awtvdr +
     &                             alidr(i,j)*awtidr +
     &                             alvdf(i,j)*awtvdf +
     &                             alidf(i,j)*awtidf
        if (iout(23)) aa(i,j,23) = aa(i,j,23) + Flatent(i,j)
        if (iout(24)) aa(i,j,24) = aa(i,j,24) + Fsensible(i,j)
        if (iout(25)) aa(i,j,25) = aa(i,j,25) + Flwout(i,j)
        if (iout(26)) aa(i,j,26) = aa(i,j,26) + evap(i,j)
        if (iout(27)) aa(i,j,27) = aa(i,j,27) + Tair(i,j)

        if (iout(28)) aa(i,j,28) = aa(i,j,28) + congel(i,j)
        if (iout(29)) aa(i,j,29) = aa(i,j,29) + frazil(i,j)
        if (iout(30)) aa(i,j,30) = aa(i,j,30) + snoice(i,j)
        if (iout(31)) aa(i,j,31) = aa(i,j,31) + meltb(i,j)
        if (iout(32)) aa(i,j,32) = aa(i,j,32) + meltt(i,j)
        if (iout(33)) aa(i,j,33) = aa(i,j,33) + Frsh_hist(i,j)
        if (iout(34)) aa(i,j,34) = aa(i,j,34) + Foht_hist(i,j)
        if (iout(35)) aa(i,j,35) = aa(i,j,35) + strairx(i,j)
        if (iout(36)) aa(i,j,36) = aa(i,j,36) + strairy(i,j)
        if (iout(37)) aa(i,j,37) = aa(i,j,37) + strtltx(i,j)
        if (iout(38)) aa(i,j,38) = aa(i,j,38) + strtlty(i,j)
        if (iout(39)) aa(i,j,39) = aa(i,j,39) + fm(i,j)*v(i,j)
        if (iout(40)) aa(i,j,40) = aa(i,j,40) - fm(i,j)*u(i,j)
        if (iout(41)) aa(i,j,41) = aa(i,j,41) + strocnx(i,j)
        if (iout(42)) aa(i,j,42) = aa(i,j,42) + strocny(i,j)
        if (iout(43)) aa(i,j,43) = aa(i,j,43) + strintx(i,j)
        if (iout(44)) aa(i,j,44) = aa(i,j,44) + strinty(i,j)
        if (iout(45)) aa(i,j,45) = aa(i,j,45) + strength(i,j)
        if (iout(46)) aa(i,j,46) = aa(i,j,46) + divu(i,j)
        if (iout(47)) aa(i,j,47) = aa(i,j,47) + shear(i,j)
        if (iout(50)) aa(i,j,50) = aa(i,j,50) + dvidtt(i,j)
        if (iout(51)) aa(i,j,51) = aa(i,j,51) + dvidtd(i,j)
        if (iout(52)) aa(i,j,52) = aa(i,j,52) + daidtd(i,j)
       enddo
      enddo

      if (write_history .or. write_ic) then

        do k=1,avgsiz
         do j=jlo,jhi
          do i=ilo,ihi
           if (.not. tmask(i,j)) then  ! mask out land points
            aa(i,j,k) = 1.0e30
           elseif (iout(k)) then       ! convert units
            aa(i,j,k) = (cona(k)*aa(i,j,k)/avgct + conb(k))
           endif
          enddo
         enddo  
        enddo

      ! snapshots 
      call principal_stress
      do j=jlo,jhi
       do i=ilo,ihi
        if (iout(48)) aa(i,j,48) = sig1(i,j) 
        if (iout(49)) aa(i,j,49) = sig2(i,j) 
       enddo
      enddo

        call icecdf
        aa = c0
        avgct = c0

      endif  ! write_history or write_ic

      end subroutine ice_write_hist

c=======================================================================

      subroutine icecdf

      use ice_model_size
      use ice_constants
      use ice_mpi_internal
      use ice_grid
      use ice_calendar

      include "netcdf.inc"

      integer (kind=int_kind) :: i,j,n
     &,  ncid,status,imtid,jmtid,timid,varid
     &,  length
      integer (kind=int_kind), dimension(3) :: dimid,start,count
      real (kind=real_kind) :: ltime
      character (char_len) :: ncfile,title
      real (kind=real_kind), parameter :: spval=1.0e30

      real (kind=dbl_kind),dimension(ilo:ihi,jlo:jhi) :: work
      real (kind=dbl_kind),dimension(imt_global,jmt_global) :: arg
      real (kind=real_kind),dimension(imt_global,jmt_global) :: ar

      integer (kind=int_kind) :: iyear, imonth, iday

      if (my_task.eq.master_task) then

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
     &          .or. n.eq.48 .or. n.eq.49) then ! sig1, sig2 snapshots
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
        title  = 'sea ice model output for ccsm2-si'
        length = lenstr(title)
        status = nf_put_att_text(ncid,nf_global,'title',length,title)
        title = 'sea ice model: ncar ccsm2-sea ice'
        length = lenstr(title)
        status = nf_put_att_text(ncid,nf_global,'source',length,title)
        title = 'log file ice.log.'                      ! for now
        length = lenstr(title)
        status = nf_put_att_text(ncid,nf_global,'history',length,title)
        title = 'of ccsm pcwg'
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

        ar=tlon_g*rad_to_deg  ! single precision
        status = nf_inq_varid(ncid,'tlon',varid)
        status = nf_put_var_real(ncid,varid,ar)

        ar=tlat_g*rad_to_deg  ! single precision
        status = nf_inq_varid(ncid,'tlat',varid)
        status = nf_put_var_real(ncid,varid,ar)

        status = nf_inq_varid(ncid,'time',varid)
        status = nf_put_vara_real(ncid,varid,1,1,ltime)

      endif          ! master_task

      !-----------------------------------------------------------------
      ! write grid mask, area and rotation angle
      !-----------------------------------------------------------------

      work(:,:) = uvm(ilo:ihi,jlo:jhi)
      call global_gather(arg,work)
      if (my_task.eq.master_task) then
        ar=arg
        status = nf_inq_varid(ncid,'hmask',varid)
        status = nf_put_var_real(ncid,varid,ar)
      endif

      work(:,:) = tarea(ilo:ihi,jlo:jhi)
      call global_gather(arg,work)
      if (my_task.eq.master_task) then
        ar=arg
        status = nf_inq_varid(ncid,'tarea',varid)
        status = nf_put_var_real(ncid,varid,ar)
      endif

      call global_gather(arg,angle)
      if (my_task.eq.master_task) then
        ar=arg
        status = nf_inq_varid(ncid,'rotang',varid)
        status = nf_put_var_real(ncid,varid,ar)
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
          work=aa(:,:,n)
          call global_gather(arg,work)
          if (my_task.eq.master_task) then
            ar=arg
            status = nf_inq_varid(ncid,vname(n),varid)
            status = nf_put_vara_real(ncid,varid,start,count,ar)
          endif
        endif
      enddo

      !-----------------------------------------------------------------
      ! close output dataset
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) then
        status = nf_close(ncid)
        write(6,*) ' '
        write(6,*) 'Finished writing ',ncfile(1:lenstr(ncfile))
        write(6,*) ' '
      endif

      end subroutine icecdf

c=======================================================================

      integer function lenstr(label)

!---!-------------------------------------------------------------------
!---! compute length of string by finding first non-blank
!---! character from the right
!---!-------------------------------------------------------------------

      character*(*) label
      integer (kind=int_kind) :: 
     &   length  ! length of character string 
     &,  n       ! loop index

      length = len(label)
      do n=length,1,-1
        if( label(n:n) .ne. ' ' ) goto 10 
      enddo
 10   lenstr = n

      end function lenstr

c=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
c=======================================================================

      subroutine dumpfile

!---!-------------------------------------------------------------------
!---! Dumps all values needed for a restart 
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_flux
      use ice_grid
      use ice_calendar
      use ice_state
      use ice_dyn_evp

      integer (kind=int_kind) :: i, j, k, n
      character(len=80) :: filename
      logical (kind=log_kind) :: gather
      real (kind=dbl_kind) :: work(ilo:ihi,jlo:jhi)

      write(filename,'(a,a,i8.8)') 
     &  dump_file(1:lenstr(dump_file)),'.',idate
      if (my_task.eq.master_task) then
        open(50,file=restrt_file)
        write(50,*) filename
        close(50)
      endif
      call ice_open(50,filename,0)
      if (my_task.eq.master_task) then
        write (50) istep1,time,time_forc
        write(6,*) 'Writing ',filename(1:lenstr(filename))
        write(6,*) 'Restart written ',istep1,time,time_forc
      endif

      gather = .true.

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------
      do n=1,ncat
        work(ilo:ihi,jlo:jhi)=aicen(ilo:ihi,jlo:jhi,n)
        call ice_write(50,0,work,'ruf8',gather)
        work(ilo:ihi,jlo:jhi)=vicen(ilo:ihi,jlo:jhi,n)
        call ice_write(50,0,work,'ruf8',gather)
        work(ilo:ihi,jlo:jhi)=vsnon(ilo:ihi,jlo:jhi,n)
        call ice_write(50,0,work,'ruf8',gather)
        work(ilo:ihi,jlo:jhi)=Tsfcn(ilo:ihi,jlo:jhi,n)
        call ice_write(50,0,work,'ruf8',gather)
        work(ilo:ihi,jlo:jhi)=esnon(ilo:ihi,jlo:jhi,n)
        call ice_write(50,0,work,'ruf8',gather)
      enddo

      do k=1,ntilay
        work(ilo:ihi,jlo:jhi)=eicen(ilo:ihi,jlo:jhi,k)
        call ice_write(50,0,work,'ruf8',gather)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      work(ilo:ihi,jlo:jhi)=u(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8',gather)
      work(ilo:ihi,jlo:jhi)=v(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8',gather)

      !-----------------------------------------------------------------
      ! fresh water and heat flux
      !-----------------------------------------------------------------
      work(ilo:ihi,jlo:jhi)=Fresh(ilo:ihi,jlo:jhi)      
      call ice_write(50,0,work,'ruf8',gather)           
      work(ilo:ihi,jlo:jhi)=Fhnet(ilo:ihi,jlo:jhi)      
      call ice_write(50,0,work,'ruf8',gather)          

      !-----------------------------------------------------------------
      ! ice strength 
      !-----------------------------------------------------------------
      work(ilo:ihi,jlo:jhi)=strength(ilo:ihi,jlo:jhi)   
      call ice_write(50,0,work,'ruf8',gather)           

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      work(ilo:ihi,jlo:jhi)=strocnxT(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8',gather)
      work(ilo:ihi,jlo:jhi)=strocnyT(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8',gather)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      do k=1,4
      work(ilo:ihi,jlo:jhi)=stressp(k,ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8',gather)
      work(ilo:ihi,jlo:jhi)=stressm(k,ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8',gather)
      work(ilo:ihi,jlo:jhi)=stress12(k,ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8',gather)
      enddo

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      do j=jlo,jhi
       do i=ilo,ihi
        work(i,j) = 0.
        if (iceumask(i,j)) work(i,j) = 1.
       enddo
      enddo
      call ice_write(50,0,work,'ruf8',gather)

!#ifndef coupled
      ! for mixed layer model
      work(ilo:ihi,jlo:jhi)=sst(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8',gather)
      work(ilo:ihi,jlo:jhi)=Fhocn(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8',gather)
      work(ilo:ihi,jlo:jhi)=frzmlt(ilo:ihi,jlo:jhi)
      call ice_write(50,0,work,'ruf8',gather)
!#endif

      if (my_task.eq.master_task) close(50)

      end subroutine dumpfile

c=======================================================================

      subroutine restartfile

!---!-------------------------------------------------------------------
!---! Restarts from a dump
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_flux
      use ice_mpi_internal
      use ice_grid
      use ice_calendar
      use ice_state
      use ice_dyn_evp
      use ice_itd

      integer (kind=int_kind) :: i, j, k, n
      character(len=80) :: filename
      logical (kind=log_kind) :: scatter, diag
      real (kind=dbl_kind) :: work(ilo:ihi,jlo:jhi)

      if (my_task.eq.master_task) then
        open(50,file=restrt_file)
        read(50,*) filename
        close(50)
      endif
      call ice_open(50,filename,0)
      if (my_task.eq.master_task) then
        read (50) istep0,time,time_forc
        write(6,*) 'read ',restrt_file(1:lenstr(restrt_file))
        write(6,*) 'dump read at istep=',istep0,time,time_forc
      endif
      call ice_bcast_iscalar(istep0)
      istep1=istep0
      call ice_bcast_rscalar(time)
      call ice_bcast_rscalar(time_forc)

      scatter = .true.
      diag = .true.

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------
      do n=1,ncat
        if (my_task.eq.master_task) 
     &   write (6,*) 'cat ',n,
     &               ' min/max area, vol ice, vol snow, Tsfc, esnon'
        call ice_read(50,0,work,'ruf8',scatter,diag)
        aicen(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
        call ice_read(50,0,work,'ruf8',scatter,diag)
        vicen(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
        call ice_read(50,0,work,'ruf8',scatter,diag)
        vsnon(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
        call ice_read(50,0,work,'ruf8',scatter,diag)
        Tsfcn(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
        call ice_read(50,0,work,'ruf8',scatter,diag)
        esnon(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
      enddo

      if (my_task.eq.master_task)
     & write (6,*) 'min/max eicen for each layer'
      do k=1,ntilay
        call ice_read(50,0,work,'ruf8',scatter,diag)
        eicen(ilo:ihi,jlo:jhi,k)=work(ilo:ihi,jlo:jhi)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write (6,*) 'min/max velocity components'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      u(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      v(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! fresh water and heat flux       
      !-----------------------------------------------------------------
      if (my_task.eq.master_task)
     & write (6,*) 'min/max fresh water and heat flux components'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      Fresh(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      Fhnet(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! ice strength 
      !-----------------------------------------------------------------
      if (my_task.eq.master_task)
     & write (6,*) 'min/max ice strength'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      strength(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write (6,*) 'min/max ocean stress components'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      strocnxT(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      strocnyT(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write (6,*) 'min/max internal stress components (3 per triangle)'
      do k=1,4
      call ice_read(50,0,work,'ruf8',scatter,diag)
      stressp(k,ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      stressm(k,ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      stress12(k,ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      enddo

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write (6,*) 'ice mask for dynamics'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      do j=jlo,jhi
       do i=ilo,ihi
        iceumask(i,j) = .false.
        if (work(i,j).gt.p5) iceumask(i,j) = .true.
       enddo
      enddo

!#ifndef coupled
      ! for mixed layer model
      if (my_task.eq.master_task) 
     & write (6,*) 'min/max sst, Fhocn, frzmlt'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      sst(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      Fhocn(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      frzmlt(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
!#endif

      if (my_task.eq.master_task) close(50)

      !-----------------------------------------------------------------
      ! update boundary conditions
      !-----------------------------------------------------------------
      do n = 1,ncat
        call bound(aicen(:,:,n))
        call bound(vicen(:,:,n))
        call bound(vsnon(:,:,n))
        call bound(Tsfcn(:,:,n))
      end do
      do k = 1,ntilay
        call bound(eicen(:,:,k))
      end do
      call bound(u)      
      call bound(v)      

      !-----------------------------------------------------------------
      ! ensure ice is binned in correct categories 
      ! should not be necessary unless restarting from a run with
      ! different category boundaries
      !-----------------------------------------------------------------
      ! These calls do not allow for exact restart capability, so they
      ! have been commented out.  This may be a problem for a run that
      ! uses different category boundaries.
c      do j=jlo,jhi
c        do i=ilo,ihi
c          call to_column(i,j) 
c          call rebin(i,j)
c          call from_column(i,j) 
c        enddo
c      enddo
c      call distr_check

c!!       call cleanup_restart

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------
      call aggregate 

      end subroutine restartfile

c=======================================================================

      subroutine restartfile_CCSM

!---!-------------------------------------------------------------------
!---! Restarts from a dump
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_flux
      use ice_mpi_internal
      use ice_grid
      use ice_calendar
      use ice_state
      use ice_dyn_evp
      use ice_itd

      integer (kind=int_kind) :: i, j, k, n
      character(len=80) :: filename
      logical (kind=log_kind) :: scatter, diag
      real (kind=dbl_kind) :: work(ilo:ihi,jlo:jhi)

      if (my_task.eq.master_task) then
        open(50,file=restrt_file)
        read(50,*) filename
        close(50)
      endif
      call ice_open(50,filename,0)
      if (my_task.eq.master_task) then
        read (50) istep0,time,time_forc
        write(6,*) 'read ',restrt_file(1:lenstr(restrt_file))
        write(6,*) 'dump read at istep=',istep0,time,time_forc
      endif
      call ice_bcast_iscalar(istep0)
      istep1=istep0
      call ice_bcast_rscalar(time)
      call ice_bcast_rscalar(time_forc)

      scatter = .true.
      diag = .true.

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------
      do n=1,ncat
        if (my_task.eq.master_task) 
     &   write (6,*) 'cat ',n,
     &               ' min/max area, vol ice, vol snow, Tsfc, esnon'
        call ice_read(50,0,work,'ruf8',scatter,diag)
        aicen(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
        call ice_read(50,0,work,'ruf8',scatter,diag)
        vicen(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
        call ice_read(50,0,work,'ruf8',scatter,diag)
        vsnon(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
        call ice_read(50,0,work,'ruf8',scatter,diag)
        Tsfcn(ilo:ihi,jlo:jhi,n)=work(ilo:ihi,jlo:jhi)
      enddo

      if (my_task.eq.master_task)
     & write (6,*) 'min/max eicen for each layer'
c      do k=1,ntilay
      do k=1,16
        call ice_read(50,0,work,'ruf8',scatter,diag)
        eicen(ilo:ihi,jlo:jhi,k)=work(ilo:ihi,jlo:jhi)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write (6,*) 'min/max velocity components'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      u(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      v(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! fresh water and heat flux       
      !-----------------------------------------------------------------
      if (my_task.eq.master_task)
     & write (6,*) 'min/max fresh water and heat flux components'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      Fresh(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      Fhnet(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! ice strength 
      !-----------------------------------------------------------------
      if (my_task.eq.master_task)
     & write (6,*) 'min/max ice strength'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      strength(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write (6,*) 'min/max ocean stress components'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      strocnxT(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      strocnyT(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write (6,*) 'min/max internal stress components (3 per triangle)'
      do k=1,4
      call ice_read(50,0,work,'ruf8',scatter,diag)
      stressp(k,ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      stressm(k,ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      stress12(k,ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      enddo

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) 
     & write (6,*) 'ice mask for dynamics'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      do j=jlo,jhi
       do i=ilo,ihi
        iceumask(i,j) = .false.
        if (work(i,j).gt.p5) iceumask(i,j) = .true.
       enddo
      enddo

       goto 999  !!! if starting stand-alone

#ifndef coupled
      ! for mixed layer model
      if (my_task.eq.master_task) 
     & write (6,*) 'min/max sst, Fhocn, frzmlt'
      call ice_read(50,0,work,'ruf8',scatter,diag)
      sst(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      Fhocn(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(50,0,work,'ruf8',scatter,diag)
      frzmlt(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
#endif

999    continue          

      if (my_task.eq.master_task) close(50)

      !-----------------------------------------------------------------
      ! update boundary conditions
      !-----------------------------------------------------------------
      do n = 1,ncat
        call bound(aicen(:,:,n))
        call bound(vicen(:,:,n))
        call bound(vsnon(:,:,n))
        call bound(Tsfcn(:,:,n))
      end do
      do k = 1,ntilay
        call bound(eicen(:,:,k))
      end do
      call bound(u)      
      call bound(v)      

      !-----------------------------------------------------------------
      ! ensure ice is binned in correct categories 
      ! should not be necessary unless restarting from a run with
      ! different category boundaries
      !-----------------------------------------------------------------
      ! These calls do not allow for exact restart capability, so they
      ! have been commented out.  This may be a problem for a run that
      ! uses different category boundaries.
c      do j=jlo,jhi
c        do i=ilo,ihi
c          call to_column(i,j) 
c          call rebin(i,j)
c          call from_column(i,j) 
c        enddo
c      enddo
c      call distr_check

       call cleanup_restart

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------
      call aggregate 

      end subroutine restartfile_CCSM

c=======================================================================

      subroutine cleanup_restart

!---!-------------------------------------------------------------------
!---! Restarts from a dump
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_flux
      use ice_mpi_internal
      use ice_grid
      use ice_calendar
      use ice_state
      use ice_dyn_evp
      use ice_itd

      integer (kind=int_kind) :: i, j, k, n
      character(len=80) :: filename
      logical (kind=log_kind) :: scatter
      real (kind=dbl_kind) :: work(ilo:ihi,jlo:jhi)

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

      do j=jlo,jhi
        do i=ilo,ihi
          call to_column(i,j) 
          call rebin(i,j)
          call from_column(i,j) 
        enddo
      enddo
      call distr_check

      end subroutine cleanup_restart

c=======================================================================

      end module ice_history

c=======================================================================
