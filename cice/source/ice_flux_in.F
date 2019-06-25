c $Id: ice_flux_in.F,v 1.14 2004/03/01 16:56:23 eclare Exp $
c=======================================================================
!
!BOP
!
! !MODULE: ice_flux_in - reads and interpolates input forcing data
!
! !DESCRIPTION:
!
! Reads and interpolates forcing data for atmospheric and oceanic quantities.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! !INTERFACE:
!
      module ice_flux_in
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain
      use ice_constants
      use ice_flux
      use ice_calendar
      use ice_read_write
      use ice_fileunits
!
!EOP
!
      implicit none

      integer (kind=int_kind) ::
     &   ycycle              ! number of years in forcing cycle
     &,  fyear_init          ! first year of data in forcing cycle
     &,  fyear_final         ! last year in cycle

      integer (kind=int_kind) :: 
     &   fyear                  ! current year in forcing cycle

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::
     &    cldf                ! cloud fraction

      character (char_len_long) ::        ! input data file names
     &   height_file
     &,   uwind_file
     &,   vwind_file
     &,    potT_file
     &,    tair_file
     &,   humid_file
     &,    rhoa_file
     &,     fsw_file
     &,     flw_file
     &,    rain_file
     &,     sst_file
     &,     sss_file

      real (kind=dbl_kind) :: 
     &     c1intp, c2intp     ! cubic interpolation coefficients  
     &,    ftime              ! forcing time (for restart)

      integer (kind=int_kind) :: 
     &     oldrecnum = 0      ! old record number (save between steps)

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,2)  ::
     &      fsw_data      ! field values at 2 temporal data points
     &,    cldf_data
     &,   fsnow_data
     &,    Tair_data
     &,    uatm_data
     &,    vatm_data
     &,      Qa_data
     &,    rhoa_data
     &,    potT_data
     &,    zlvl_data
     &,     flw_data
     &,     sst_data
     &,     sss_data

      character(char_len_long) :: 
     &   atm_data_dir     ! top directory for atmospheric data
     &,  ocn_data_dir     ! top directory for ocean data

c=======================================================================

      contains

c=======================================================================
!
!BOP
!
! !IROUTINE: init_getflux - initialize input forcing files
!
! !INTERFACE:
!
      subroutine init_getflux
!
! !DESCRIPTION:
!
! Determines the current and final years of the forcing cycle based on
! namelist input, initializes the forcing data filenames, and initializes
! surface temperature and salinity from data. \\
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_history, only: restart
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      fyear_final = fyear_init + ycycle - 1 ! last year in forcing cycle
      fyear = fyear_init + mod(nyr-year_init,ycycle)  ! current year

      if (my_task.eq.master_task) then
      write (nu_diag,*) ' Initial forcing data year = ',fyear_init
      write (nu_diag,*) ' Current forcing data year = ',fyear
      write (nu_diag,*) ' Final   forcing data year = ',fyear_final
      endif

!ech-later      call Tair_limit

      call NCAR_files(fyear)       ! data for individual years

c      call sss_clim                      ! climatology (12-month avg)
c      if ( .not. (restart)) call sst_ic  ! not interpolated but depends on sss

      end subroutine init_getflux

c=======================================================================
!
!BOP
!
! !IROUTINE: getflux - Get forcing data and interpolate as necessary
!
! !INTERFACE:
!
      subroutine getflux
!
! !DESCRIPTION:
!
! Get forcing data and interpolate as necessary! 
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      fyear = fyear_init + mod(nyr-year_init,ycycle)  ! current year

      ftime = time         ! forcing time
      time_forc = ftime    ! for restarting

    !-------------------------------------------------------------------
    ! Read and interpolate annual climatologies of SSS and SST.
    ! Restore model SST to data.
    !-------------------------------------------------------------------

c      call sss_sst_restore

    !-------------------------------------------------------------------
    ! Read and interpolate atmospheric data
    !-------------------------------------------------------------------
      call NCAR_bulk_dat
      call prepare_forcing

      end subroutine getflux

c=======================================================================
!
!BOP
!
! !IROUTINE: read_data - Read data needed for interpolation
!
! !INTERFACE:
!
      subroutine read_data (flag, recd, yr, imx, ixx, ipx,
     &                      maxrec, data_file, field_data)
!
! !DESCRIPTION:
!
! If data is at the beginning of a one-year record, get data from
!  the previous year.
! If data is at the end of a one-year record, get data from the 
!  following year.
! If no earlier data exists (beginning of fyear\_init), then \\
!  (1) For monthly data, get data from the end of fyear\_final. \\
!  (2) For more frequent data, let the imx value equal the 
!      first value of the year. \\
! If no later data exists (end of fyear\_final), then \\
!  (1) For monthly data, get data from the beginning of fyear\_init. \\
!  (2) For more frequent data, let the ipp and, if necessary, ipx 
!      values equal the last value of the year. \\
! In other words, we assume persistence when daily or 6-hourly
!   data is missing, and we assume periodicity when monthly data
!   is missing. \\
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! !USES:
!
      use ice_diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
!
      logical (kind=log_kind), intent(in) :: flag

      integer (kind=int_kind), intent(in) :: 
     &  recd                    ! baseline record number
     &, yr                      ! year of forcing data
     &, imx, ixx, ipx           ! record numbers of 3 data values
                                ! relative to recd
     &, maxrec                  ! maximum record value

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,2), intent(out) ::
     &  field_data              ! 2 values needed for interpolation
!
!EOP
!
      character (char_len_long) ::
     &  data_file               ! data file to be read 

      integer (kind=int_kind) :: 
     &  nbits                ! = 32 for single precision, 64 for double
     &, nrec                 ! record number to read
     &, n2, n4               ! like imx and ipx, but 
                             ! adjusted at beginning and end of data
     &, arg                  ! value of time argument in field_data

      logical (kind=log_kind) :: scatter, diag

      nbits = 64                ! double precision data

      scatter = .true.     ! scatter data to processors
      diag = .false.       ! write diagnostic info

      if (istep1 > check_step) diag = .true.  !! debugging

c      if (my_task==master_task .and. (diag)) then
c         write(nu_diag,*) '  ',data_file
c      endif

      if (flag) then
      !-----------------------------------------------------------------
      ! Initialize record counters
      ! (n2, n4 will change only at the very beginning or end of
      !  a forcing cycle.)
      !-----------------------------------------------------------------
      n2 = imx
      n4 = ipx
      arg = 0
      
      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

      if (imx .ne. 99) then
      ! currently in first half of data interval
        if (ixx<=1) then
          if (yr > fyear_init) then  ! get data from previous year
            call file_year (data_file, yr-1)
          else                   ! yr = fyear_init, no prior data exists
            if (maxrec > 12) then  ! extrapolate from first record
              if (ixx==1) n2 = ixx
            else                 ! go to end of fyear_final
              call file_year (data_file, fyear_final)
            endif
          endif                  ! yr > fyear_init
        endif                    ! ixx <= 1
        call ice_open (nu_forcing, data_file, nbits)

        arg = 1
        nrec = recd + n2
        call ice_read (nu_forcing, nrec, field_data(:,:,arg),
     &                'rda8', scatter, diag)
 
        if (ixx==1 .and. my_task==master_task) close(nu_forcing)
      endif  ! imx ne 99

      ! always read ixx data from data file for current year
      call file_year (data_file, yr)
      call ice_open (nu_forcing, data_file, nbits)

      arg = arg + 1
      nrec = recd + ixx
      call ice_read (nu_forcing, nrec, field_data(:,:,arg),
     &                'rda8', scatter, diag)

      if (ipx .ne. 99) then 
      ! currently in latter half of data interval
        if (ixx==maxrec) then
          if (yr < fyear_final) then ! get data from following year
            if (my_task.eq.master_task) close(nu_forcing)
            call file_year (data_file, yr+1)
            call ice_open (nu_forcing, data_file, nbits)
          else                     ! yr = fyear_final, no more data exists
            if (maxrec > 12) then  ! extrapolate from ixx
              n4 = ixx
            else                   ! go to beginning of fyear_init
              if (my_task.eq.master_task) close(nu_forcing)
              call file_year (data_file, fyear_init)
              call ice_open (nu_forcing, data_file, nbits)
            endif
          endif                    ! yr < fyear_final
        endif                      ! ixx = maxrec

        arg = arg + 1
        nrec = recd + n4
        call ice_read (nu_forcing, nrec, field_data(:,:,arg),
     &                'rda8', scatter, diag)
      endif  ! ipx ne 99

      if (my_task.eq.master_task) close(nu_forcing)
      endif  ! flag

      end subroutine read_data

c=======================================================================
!
!BOP
!
! !IROUTINE: read_clim_data - read annual climatological data
!
! !INTERFACE:
!
      subroutine read_clim_data (readflag, recd, imx, ixx, ipx,
     &                           data_file, field_data)
!
! !DESCRIPTION:
!
! Read data needed for interpolation, as in read\_data.
! Assume a one-year cycle of climatological data, so that there is
!  no need to get data from other years or to extrapolate data beyond
!  the forcing time period.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! !USES:
!
      use ice_diagnostics  !! debugging
!
! !INPUT/OUTPUT PARAMETERS:
!
      logical (kind=log_kind),intent(in) :: readflag

      integer (kind=int_kind), intent(in) :: 
     &  recd                ! baseline record number
     &, imx,ixx,ipx         ! record numbers of 3 data values
                            ! relative to recd

      character (char_len_long), intent(in) ::  data_file

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,2), intent(out) ::
     &  field_data         ! 2 values needed for interpolation
!
!EOP
!
      integer (kind=int_kind) :: 
     &  nbits              ! = 32 for single precision, 64 for double
     &, nrec               ! record number to read
     &, arg                ! value of time argument in field_data

      logical (kind=log_kind) :: scatter, diag

      nbits = 64                ! double precision data

      scatter = .true.     ! scatter data to processors
      diag = .false.       ! write diagnostic info

      if (istep1 > check_step) diag = .true.  !! debugging

      if (my_task==master_task .and. (diag)) 
     &  write(nu_diag,*) '  ',data_file

      if (readflag) then
      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------
      
      call ice_open (nu_forcing, data_file, nbits)

      arg = 0
      if (imx .ne. 99) then
        arg = 1
        nrec = recd + imx
        call ice_read (nu_forcing, nrec, field_data(:,:,arg),
     &                 'rda8', scatter, diag)
      endif

      arg = arg + 1
      nrec = recd + ixx
      call ice_read (nu_forcing, nrec, field_data(:,:,arg),
     &                 'rda8', scatter, diag)

      if (ipx .ne. 99) then
        arg = arg + 1
        nrec = recd + ipx
        call ice_read (nu_forcing, nrec, field_data(:,:,arg),
     &                 'rda8', scatter, diag)
      endif

      if (my_task.eq.master_task) close (nu_forcing)
      endif  ! readflag

      end subroutine read_clim_data

c=======================================================================
!
!BOP
!
! !IROUTINE: interp_coeff_monthly - Compute monthly data interpolation coefficients
!
! !INTERFACE:
!
      subroutine interp_coeff_monthly (recslot)
!
! !DESCRIPTION:
!
! Compute coefficients for interpolating monthly data to current time step.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &    recslot         ! slot (1 or 2) for current record
!
!EOP
!
      real (kind=dbl_kind) ::
     &    tt               ! seconds elapsed in current year
     &,   t1, t2           ! seconds elapsed at month midpoint

      real (kind=int_kind) :: 
     &    daymid(0:13)       ! month mid-points

      daymid(1:13) = 14.   ! time frame ends 0 sec into day 15
      daymid(0)    = -17.   ! Dec 15, 0 sec

      ! make time cyclic
      tt = mod(ftime/secday,c365)
 
      ! Find neighboring times
      
      if (recslot==2) then      ! first half of month
        t2 = daycal(month) + daymid(month)   ! midpoint, current month
        if (month == 1) then
          t1 = daymid(0)                 ! Dec 15 (0 sec)
        else
          t1 = daycal(month-1) + daymid(month-1) ! midpoint, previous month
        endif
      else                      ! second half of month
        t1 = daycal(month) + daymid(month)    ! midpoint, current month
        t2 = daycal(month+1) + daymid(month+1)! day 15 of next month (0 sec)
      endif
 
      ! Compute coefficients
      c1intp = (t2 - tt) / (t2 - t1)
      c2intp =  c1 - c1intp

      end subroutine interp_coeff_monthly

c=======================================================================
!
!BOP
!
! !IROUTINE: interp_coeff
!
! !INTERFACE:
!
      subroutine interp_coeff (recnum, recslot, secint)
!
! !DESCRIPTION:
!
! Compute coefficients for interpolating data to current time step.
! Works for any data interval that divides evenly into a 365-day 
!  year (daily, 6-hourly, etc.)
! Use interp\_coef\_monthly for monthly data.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! !USES:
!
      integer (kind=int_kind), intent(in) :: 
     &    recnum          ! record number for current data value
     &,   recslot         ! spline slot for current record 

      real (kind=dbl_kind), intent(in) ::
     &    secint                    ! seconds in data interval
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind), parameter ::
     &    secyr = c365 * secday     ! seconds in a 365-day year 

      real (kind=dbl_kind) ::
     &    tt               ! seconds elapsed in current year
     &,   t1, t2           ! seconds elapsed at data points

      tt = mod(ftime,secyr)
 
      ! Find neighboring times
      if (recslot==2) then          ! current record goes in slot 2 (NCEP)
         t2 = real(recnum)*secint
         t1 = t2 - secint           !  - 1 interval
      else                          ! recslot = 1
         t1 = real(recnum)*secint
         t2 = t1 + secint           !  + 1 interval
      endif
 
      ! Compute coefficients
      c1intp =  abs((t2 - tt) / (t2 - t1))
      c2intp =  c1 - c1intp

      end subroutine interp_coeff

c=======================================================================
!
!BOP
!
! !IROUTINE: file_year - 
!
! !INTERFACE:
!
      subroutine file_year (data_file, yr)
!
! !DESCRIPTION:
!
! Construct the correct name of the atmospheric data file
! to be read, given the year and assuming the naming convention
! that filenames end with 'yyyy.dat'.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (char_len_long), intent(inout) ::  data_file
!
!EOP
!
      integer (kind=int_kind), intent(in) :: yr

      character (char_len_long) :: tmpname

      integer (kind=int_kind) :: i

      i = index(data_file,'.dat') - 5
      tmpname = data_file
      write(data_file,'(a,i4.4,a)') tmpname(1:i), yr, '.dat'

      end subroutine file_year

c=======================================================================
!
!BOP
!
! !IROUTINE: NCAR_files - construct filenames for NCAR_bulk atmospheric data
!
! !INTERFACE:
!
      subroutine NCAR_files (yr)
!
! !DESCRIPTION:
!
! This subroutine is based on the LANL file naming conventions.
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file\_year will insert the correct year.
! 
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: 
     &     yr                   ! current forcing year
!
!EOP
!
      fsw_file =
     &     trim(atm_data_dir)//'ISCCPM/MONTHLY/RADFLX/swdn.1996.dat'
      call file_year(fsw_file,yr)

      flw_file =
     &     trim(atm_data_dir)//'ISCCPM/MONTHLY/RADFLX/cldf.1996.dat'
      call file_year(flw_file,yr)

      rain_file =
     &     trim(atm_data_dir)//'MXA/MONTHLY/PRECIP/prec.1996.dat'
      call file_year(rain_file,yr)

      uwind_file =
     &     trim(atm_data_dir)//'NCEP/4XDAILY/STATES/u_10.1996.dat'
      call file_year(uwind_file,yr)

      vwind_file =
     &     trim(atm_data_dir)//'NCEP/4XDAILY/STATES/v_10.1996.dat'
      call file_year(vwind_file,yr)

      tair_file =
     &     trim(atm_data_dir)//'NCEP/4XDAILY/STATES/t_10.1996.dat'
      call file_year(tair_file,yr)

      humid_file =
     &     trim(atm_data_dir)//'NCEP/4XDAILY/STATES/q_10.1996.dat'
      call file_year(humid_file,yr)

      rhoa_file =
     &     trim(atm_data_dir)//'NCEP/4XDAILY/STATES/dn10.1996.dat'
      call file_year(rhoa_file,yr)

      if (my_task.eq.master_task) then
         write (nu_diag,*) ''
         write (nu_diag,*) 'Forcing data year = ', fyear         
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) fsw_file
         write (nu_diag,*) flw_file
         write (nu_diag,*) rain_file
         write (nu_diag,*) uwind_file
         write (nu_diag,*) vwind_file
         write (nu_diag,*) tair_file
         write (nu_diag,*) humid_file
         write (nu_diag,*) rhoa_file
      endif                     ! master_task

      end subroutine NCAR_files

c=======================================================================
!
!BOP
!
! !IROUTINE: NCAR_bulk_dat - read NCAR_bulk atmospheric data
!
! !INTERFACE:
!
      subroutine NCAR_bulk_dat
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: 
     &    i, j
     &,   imx,ixx,ipx         ! record numbers for neighboring months
     &,   recnum              ! record number
     &,   maxrec              ! maximum record number
     &,   recslot             ! spline slot for current record
     &,   midmonth            ! middle day of month

      real (kind=dbl_kind) ::
     &    sec6hr              ! number of seconds in 6 hours

      logical (kind=log_kind) :: readm, read6

    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------
      
      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      imx  = mod(month+maxrec-2,maxrec) + 1
      ipx  = mod(month,         maxrec) + 1
      if (mday >= midmonth) imx = 99  ! other two points will be used
      if (mday <  midmonth) ipx = 99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values 
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

      call read_data (readm, 0, fyear, imx, month, ipx, 
     &      maxrec, fsw_file, fsw_data)
      call read_data (readm, 0, fyear, imx, month, ipx, 
     &      maxrec, flw_file, cldf_data)
      call read_data (readm, 0, fyear, imx, month, ipx, 
     &      maxrec, rain_file, fsnow_data)

      call interpolate_data (fsw_data, fsw)
      call interpolate_data (cldf_data, cldf)
      call interpolate_data (fsnow_data, fsnow)

    !-------------------------------------------------------------------
    ! 6-hourly data
    ! 
    ! Assume that the 6-hourly value is located at the end of the
    !  6-hour period.  This is the convention for NCEP reanalysis data.
    !  E.g. record 1 gives conditions at 6 am GMT on 1 January.
    !-------------------------------------------------------------------

      sec6hr = secday/c4        ! seconds in 6 hours
      maxrec = 1460             ! 365*4

      ! current record number
      recnum = 4*int(yday) - 3 + int(real(sec)/sec6hr)

      ! Compute record numbers for surrounding data (2 on each side)

      imx = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
!      ipx = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      recslot = 2
      ipx = 99
      call interp_coeff (recnum, recslot, sec6hr)

      ! Read
      read6 = .false.
      if (istep==1 .or. oldrecnum .ne. recnum) read6 = .true.

      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec,
     &                tair_file, Tair_data)
      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec,
     &                uwind_file, uatm_data)
      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec,
     &                vwind_file, vatm_data)
      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec,
     &                rhoa_file, rhoa_data)
      call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec,
     &                humid_file, Qa_data)

      ! Interpolate
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (uatm_data, uatm)
      call interpolate_data (vatm_data, vatm)
      call interpolate_data (rhoa_data, rhoa)
      call interpolate_data (Qa_data, Qa)

      ! Save record number
      oldrecnum = recnum

      end subroutine NCAR_bulk_dat

c=======================================================================
!
!BOP
!
! !IROUTINE: prepare_forcing - finish manipulating forcing
!
! !INTERFACE:
!
      subroutine prepare_forcing
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_grid, only: ANGLET, t2ugrid
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: workx, worky

      do j=jlo,jhi
       do i=ilo,ihi

      !-----------------------------------------------------------------
      ! Fix interpolation errors
      !-----------------------------------------------------------------
        fsw (i,j)  = max(fsw(i,j),c0)
        cldf(i,j)  = max(min(cldf(i,j),c1),c0)
        fsnow(i,j) = max(fsnow(i,j),c0)   
        rhoa(i,j)  = max(rhoa(i,j),c0)
        Qa  (i,j)  = max(Qa(i,j),c0)

      !-----------------------------------------------------------------
      ! Correct known biases
      !
      ! NOTE: These correction factors may not be appropriate for other
      !       data sets.
      !-----------------------------------------------------------------    
        Qa(i,j)  = Qa(i,j)  * 0.94_dbl_kind ! as in the dummy atm (latm)
        fsw(i,j) = fsw(i,j) * 0.92_dbl_kind ! as in the dummy atm (latm)

      !-----------------------------------------------------------------
      ! Compute other fields needed by model
      !-----------------------------------------------------------------

        zlvl(i,j) = c10
        wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s
        potT(i,j) = Tair(i,j)

        ! divide shortwave into spectral bands
        swvdr(i,j) = fsw(i,j)*(.28)         ! visible direct
        swvdf(i,j) = fsw(i,j)*(.24)         ! visible diffuse
        swidr(i,j) = fsw(i,j)*(.31)         ! near IR direct
        swidf(i,j) = fsw(i,j)*(.17)         ! near IR diffuse
                                            ! as in the dummy atm (latm)

        ! longwave as in Parkinson and Washington (1979)
        flw(i,j) = stefan_boltzmann*Tair(i,j)**4   !!! downward longwave !!!
     &   * (c1-
     &   0.261_dbl_kind*exp(-7.77e-4_dbl_kind*(Tffresh - Tair(i,j))**2))
     &   * (c1 + 0.275_dbl_kind*cldf(i,j))

        ! determine whether precip is rain or snow
        fsnow(i,j) = fsnow(i,j)/2.592e+06_dbl_kind  ! mm/month -> kg/m^2 s
        frain(i,j) = c0                     
!        if (Tair(i,j) >= Tffresh) then
        if (Tair(i,j) .ge. Tffresh) then
            frain(i,j) = fsnow(i,j)
            fsnow(i,j) = c0
        endif

      !-----------------------------------------------------------------
      ! rotate zonal/meridional vectors to local coordinates
      ! Vector fields come in on T grid, but are oriented geographically
      ! need to rotate to pop-grid FIRST using ANGLET
      ! then interpolate to the U-cell centers  (otherwise we
      ! interpolate across the pole)
      ! use ANGLET which is on the T grid !
      ! atmo variables are needed in T cell centers in subroutine stability,
      ! and are interpolated to the U grid later as necessary
      !-----------------------------------------------------------------
        workx      = uatm(i,j)                ! wind velocity, m/s
        worky      = vatm(i,j) 
        uatm (i,j) = workx*cos(ANGLET(i,j))   ! convert to POP grid
     1             + worky*sin(ANGLET(i,j))   ! note uatm, vatm, wind
        vatm (i,j) = worky*cos(ANGLET(i,j))   !  are on the T-grid here
     1             - workx*sin(ANGLET(i,j))

       enddo
      enddo

      end subroutine prepare_forcing

c=======================================================================
!
!BOP
!
! !IROUTINE: interpolate_data 
!
! !INTERFACE:
!
      subroutine interpolate_data (field_data, field)
!
! !DESCRIPTION:
!
! Linear interpolation
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,2), intent(in) ::
     &  field_data    ! 2 values used for interpolation

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(out) ::
     &  field         ! interpolated field
!
!EOP
!
      integer (kind=int_kind) :: i,j

      do j=jlo,jhi
       do i=ilo,ihi
          field(i,j) = c1intp * field_data(i,j,1) 
     &               + c2intp * field_data(i,j,2) 
       enddo
      enddo

      end subroutine interpolate_data

c=======================================================================
!
!BOP
!
! !IROUTINE: sss_sst_restore - interpolates sss, sst; restores sst
!
! !INTERFACE:
!
      subroutine sss_sst_restore
!
! !DESCRIPTION:
!
! Interpolate monthly sss, sst data to timestep.
! Restore sst computed by ice model to data.
! 
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_ocean
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: 
     &    i, j
     &,   imx,ipx         ! record numbers for neighboring months
     &,   maxrec              ! maximum record number
     &,   recslot             ! spline slot for current record
     &,   midmonth            ! middle day of month

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: 
     &    sstdat              ! data value toward which SST is restored

      real (kind=dbl_kind) :: Trestore

      logical (kind=log_kind) :: readm

      sss_file = trim(ocn_data_dir)//'sss_Lev.mm'
      sst_file = trim(ocn_data_dir)//'sst_Lev.mm'

c      Trestore = dt       ! STRONG restoring 
c      Trestore = 1.296e6  ! 15 days
c      Trestore = 2.592e6  ! 30 days
      Trestore = 7.776e6  ! 90 days

      if (my_task.eq.master_task .and. istep.eq.1) then
         write (nu_diag,*) ''
         write (nu_diag,*) 'SSS and SST data interpolated to timestep:'
         write (nu_diag,*) sss_file
         write (nu_diag,*) sst_file
         write (nu_diag,*) 
     &        'SST restoring timescale = ',Trestore/secday,' days'
      endif

    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------
      
      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      imx  = mod(month+maxrec-2,maxrec) + 1
      ipx  = mod(month,         maxrec) + 1
      if (mday >= midmonth) imx = 99  ! other two points will be used
      if (mday <  midmonth) ipx = 99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values 
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

      call read_clim_data (readm, 0, imx, month, ipx,
     &                     sss_file, sss_data)
      call read_clim_data (readm, 0, imx, month, ipx,
     &                     sst_file, sst_data)

    ! Interpolate
      call interpolate_data (sss_data, sss)
      call interpolate_data (sst_data, sstdat)

    !-------------------------------------------------------------------
    ! Restore sst to data
    !-------------------------------------------------------------------
      do j = jlo,jhi
       do i = ilo,ihi
        if (sss(i,j).lt.c0) sss(i,j) = c0
        sst(i,j) = sst(i,j) + (sstdat(i,j)-sst(i,j))*dt/Trestore 
       enddo
      enddo

      call complete_getflux_ocn

      end subroutine sss_sst_restore

c=======================================================================
!
!BOP
!
! !IROUTINE: complete_getflux_ocn - compute remaining ocean forcing fields
!
! !INTERFACE:
!
      subroutine complete_getflux_ocn
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j

    !-------------------------------------------------------------------
    ! Compute freezing temperature based on salinity
    ! Make sure SST is not below Tf.
    !-------------------------------------------------------------------
      do j=jlo,jhi
       do i=ilo,ihi
         Tf(i,j) = -depressT*sss(i,j)   ! freezing temp, deg C
         sst(i,j) = max (sst(i,j), Tf(i,j))
       enddo
      enddo

    !-------------------------------------------------------------------
    ! interpolate ocean dynamics variables from T-cell centers to 
    ! U-cell centers   
    !-------------------------------------------------------------------
c      call t2ugrid(uocn)
c      call t2ugrid(vocn)
c      call t2ugrid(ss_tltx)
c      call t2ugrid(ss_tlty)

      end subroutine complete_getflux_ocn

c=======================================================================
!
!BOP
!
! !IROUTINE: sss_clim - annual mean climatology for Levitus sss
!
! !INTERFACE:
!
      subroutine sss_clim
!
! !DESCRIPTION:
!
! Creates annual mean climatology for Levitus sss from 12-month climatology.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_work, only:  worka
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, nbits, k

      logical (kind=log_kind) :: scatter, diag

      scatter = .true.
      diag = .false.

      nbits = 64                ! double precision data

      sss_file = trim(ocn_data_dir)//'sss_Lev.mm'

      if (my_task.eq.master_task) then
         write (nu_diag,*) ''
         write (nu_diag,*) 'SSS climatology computed from:'
         write (nu_diag,*) sss_file
      endif

      if (my_task.eq.master_task) 
     &     call ice_open (nu_forcing, sss_file, nbits)

    !-------------------------------------------------------------------
    ! create surface salinity climatology from monthly data
    !-------------------------------------------------------------------

      do j = jlo,jhi
       do i = ilo,ihi
        worka(i,j) = c0
       enddo
      enddo

      do k = 1,12   ! loop over 12 months

         call ice_read (nu_forcing, k, worka, 'rda8', scatter, diag)
         do j = jlo,jhi
          do i = ilo,ihi
             sss(i,j) = worka(i,j) + sss(i,j)
          enddo
         enddo

      enddo  ! k

      do j = jlo,jhi
       do i = ilo,ihi
        sss(i,j) = sss(i,j) / c12       ! annual average salinity
        if (sss(i,j).lt.c0) sss(i,j) = c0
       enddo
      enddo

      ! close file
      if (my_task.eq.master_task) close(nu_forcing)

      call complete_getflux_ocn

      end subroutine sss_clim

c=======================================================================
!
!BOP
!
! !IROUTINE: sst_ic - sst initial condition
!
! !INTERFACE:
!
      subroutine sst_ic
!
! !DESCRIPTION:
!
! Reads sst data for current month, and adjusts sst based on freezing 
! temperature.  Does not interpolate.
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, nbits, k

      logical (kind=log_kind) :: scatter, diag

      scatter = .true.
      diag = .false.

      nbits = 64                ! double precision data

      !-----------------------------------------------------------------
      ! Shea, Trenberth and Reynolds SST
      !-----------------------------------------------------------------
!      sst_file = '/n/filer3/climate/eclare/sst.sss/sst.mm.100x116.da'
      sst_file = trim(ocn_data_dir)//'sst_Lev.mm'

      if (my_task.eq.master_task) then
         write (nu_diag,*) ''
         write (nu_diag,*) 'SST initial condition:'
         write (nu_diag,*) sst_file
      endif

      if (my_task.eq.master_task) 
     &     call ice_open (nu_forcing, sst_file, nbits)

      call ice_read (nu_forcing, month, sst, 'rda8', scatter, diag)

      if (my_task.eq.master_task) close(nu_forcing)

      do j=jlo,jhi
       do i=ilo,ihi
         ! Make sure sst is not less than Tf
         sst(i,j) = max(sst(i,j),Tf(i,j))
       enddo
      enddo

      end subroutine sst_ic

c=======================================================================

      end module ice_flux_in

c=======================================================================
