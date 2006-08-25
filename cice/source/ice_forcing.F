!=======================================================================
!
!BOP
!
! !MODULE: ice_forcing - reads and interpolates input forcing data
!
! !DESCRIPTION:
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Block structure added by WHL
! 2005: ECMWF option added by WHL
! 2006: LY option added by ECH; module name changed from ice_flux_in
!
! !INTERFACE:
!
      module ice_forcing
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size
      use ice_communicate, only: my_task, master_task
      use ice_constants
      use ice_calendar, only: istep, istep1, time, time_forc, year_init,
     &                        sec, mday, month, nyr, yday, daycal
      use ice_fileunits
      use ice_exit
!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) ::
     &   ycycle              ! number of years in forcing cycle
     &,  fyear_init          ! first year of data in forcing cycle
     &,  fyear               ! current year in forcing cycle
     &,  fyear_final         ! last year in cycle

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
     &     c1intp, c2intp     ! interpolation coefficients
     &,    ftime              ! forcing time (for restart)

      integer (kind=int_kind) ::
     &     oldrecnum = 0      ! old record number (save between steps)
     &,    oldrecslot = 1     ! old record slot (save between steps)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) ::
     &    cldf                ! cloud fraction

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks) ::
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

      character(char_len) :: 
     &   atm_data_type    ! 'default', 'ncar', 'ecmwf', or 'LYq'
     &,  sss_data_type    ! 'default', 'clim', or 'ncar'
     &,  sst_data_type    ! 'default', 'clim', or 'ncar'
     &,  precip_units     ! 'mm_per_month', 'mm_per_sec', 'mks'
 
      character(char_len_long) :: 
     &   atm_data_dir     ! top directory for atmospheric data
     &,  ocn_data_dir     ! top directory for ocean data
     &,  oceanmixed_file  ! netCDF file name for ocean forcing data

      integer (kind=int_kind), parameter :: 
     &   nfld = 8    ! number of fields to search for in netCDF file

      real (kind=dbl_kind), dimension (nx_block,ny_block,nfld,12) :: 
     &   ocn_frc_m   ! ocn data for 12 months

      logical (kind=log_kind) ::
     &   restore_sst                 ! restore sst if true

      integer (kind=int_kind) ::
     &   trestore                    ! restoring time scale (days)

      real (kind=dbl_kind) :: 
     &   trest                       ! restoring time scale (sec)

      logical (kind=log_kind) ::
     &   dbug             ! prints debugging output if true

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: init_forcing_atmo - initialize atmospheric forcing
!
! !INTERFACE:
!
      subroutine init_forcing_atmo
!
! !DESCRIPTION:
!
! Determine the current and final year of the forcing cycle based on
! namelist input; initialize the forcing data filenames.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      fyear       = fyear_init + mod(nyr-1,ycycle) ! current year
      fyear_final = fyear_init + ycycle - 1 ! last year in forcing cycle

      if (trim(atm_data_type) /= 'default' .and.
     &                    my_task == master_task) then
         write (nu_diag,*) ' Initial forcing data year = ',fyear_init
         write (nu_diag,*) ' Final   forcing data year = ',fyear_final
      endif

    !-------------------------------------------------------------------
    ! Get filenames for input forcing data     
    !-------------------------------------------------------------------

      ! default forcing values from init_flux_atm
      if (trim(atm_data_type) == 'ncar') then
         call NCAR_files(fyear)
      elseif (trim(atm_data_type) == 'ecmwf') then
         call ecmwf_files(fyear)    
      elseif (trim(atm_data_type) == 'LYq') then
         call LY_files(fyear)
      endif

      end subroutine init_forcing_atmo

!=======================================================================
!BOP
!
! !IROUTINE: init_forcing_ocn - initialize sss and sst
!
! !INTERFACE:
!
      subroutine init_forcing_ocn(dt)
!
! !DESCRIPTION:
!
! Set sea surface salinity to annual mean value using a 12-month
!  climatology.
! Read sst data for current month, and adjust sst based on freezing
! temperature.  No interpolation in time.

! Note: SST is subsequently prognosed if CICE is run with a mixed layer
! ocean (oceanmixed\_ice = T).  SSS is not prognosed by CICE. 
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_domain, only: nblocks
      use ice_flux, only: sss, sst, Tf
      use ice_work, only:  work1
      use ice_read_write
      use ice_therm_vertical, only: ustar_scale
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::
     &     dt                   ! time step
!
!EOP
!
      integer (kind=int_kind) ::
     &     i, j, iblk           ! horizontal indices
     &,    k                    ! month index
     &,    nbits

      nbits = 64                ! double precision data

      if (restore_sst) then
         if (trestore == 0) then
            trest = dt          ! use data instantaneously
         else
            trest = real(trestore,kind=dbl_kind) * secday ! seconds
         endif
      endif

      if (trim(sss_data_type) == 'clim') then

         sss_file = trim(ocn_data_dir)//'sss_Lev.mm'
!!!         sss_file = trim(ocn_data_dir)//'sss_12.r'

         if (my_task == master_task) then
            write (nu_diag,*) ''
            write (nu_diag,*) 'SSS climatology computed from:'
            write (nu_diag,*) trim(sss_file)
         endif

         if (my_task == master_task)
     &        call ice_open (nu_forcing, sss_file, nbits)

    !-------------------------------------------------------------------
    ! create surface salinity climatology from monthly data
    !-------------------------------------------------------------------

         sss(:,:,:) = c0

         do k = 1,12            ! loop over 12 months
            call ice_read (nu_forcing, k, work1, 'rda8', dbug)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  sss(i,j,iblk) = sss(i,j,iblk) + work1(i,j,iblk)
               enddo
               enddo
            enddo
         enddo                  ! k

         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sss(i,j,iblk) = sss(i,j,iblk) / c12   ! annual average
               sss(i,j,iblk) = max(sss(i,j,iblk),c0)
            enddo
            enddo
         enddo

         ! close file
         if (my_task == master_task) close(nu_forcing)

      endif                     ! sss_data_type

      if (trim(sst_data_type) == 'clim') then

         if (nx_global == 320) then ! gx1
            sst_file = trim(ocn_data_dir)//'sst_clim_hurrell.dat'
         else                   ! gx3
            sst_file = trim(ocn_data_dir)//'sst_Lev.mm'
!!!            sst_file = trim(ocn_data_dir)//'sst_12.r'
         endif

         if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'Initial SST file:', trim(sst_file)
         endif

         if (my_task == master_task)
     &        call ice_open (nu_forcing, sst_file, nbits)

         call ice_read (nu_forcing, month, sst, 'rda8', dbug)

         if (my_task == master_task) close(nu_forcing)

         ! Make sure sst is not less than freezing temperature Tf
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               Tf (i,j,iblk) = -depressT*sss(i,j,iblk) ! deg C
               sst(i,j,iblk) = max(sst(i,j,iblk),Tf(i,j,iblk))
            enddo
            enddo
         enddo

      endif                     ! init_sst_data

      if (trim(sst_data_type) == 'ncar' .or. 
     &    trim(sss_data_type) == 'ncar') then
!lipscomb - not implemented yet
!!!         call getforcing_ocn_ncar_init
         call abort_ice ('ncar ocean forcing not implemented yet')
      endif

      ! set ustar_scale for case of zero currents
      ! default value of c1 (for nonzero currents) set in init_thermo_vertical

      if (trim(sst_data_type) /= 'ncar' .or. 
     &    trim(sss_data_type) /= 'ncar') then
         ustar_scale = c10            ! for zero currents
      endif

      end subroutine init_forcing_ocn

!=======================================================================
!BOP
!
! !IROUTINE: get_forcing_atmo - Get atmospheric forcing data and interpolate
!
! !INTERFACE:
!
      subroutine get_forcing_atmo
!
! !DESCRIPTION:
!
! Get atmospheric forcing data and interpolate as necessary
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_domain
      use ice_blocks
      use ice_flux
      use ice_state
      use ice_grid, only: ANGLET, hm
!
!EOP
!
      integer (kind=int_kind) ::
     &   iblk                 ! block index

      fyear = fyear_init + mod(nyr-1,ycycle)  ! current year
      if (trim(atm_data_type) /= 'default' .and. istep <= 1
     &             .and. my_task == master_task) then
         write (nu_diag,*) ' Current forcing data year = ',fyear
      endif

      ftime = time         ! forcing time
      time_forc = ftime    ! for restarting

    !-------------------------------------------------------------------
    ! Read and interpolate atmospheric data
    !-------------------------------------------------------------------

      if (trim(atm_data_type) == 'ncar') then
         call ncar_data
      elseif (trim(atm_data_type) == 'ecmwf') then
         call ecmwf_data
      elseif (trim(atm_data_type) == 'LYq') then
         call LY_data
!     else    ! default values set in init_flux
      endif

      do iblk = 1, nblocks

    !-------------------------------------------------------------------
    ! Convert forcing data to fields needed by ice model
    !-------------------------------------------------------------------

         call prepare_forcing (nx_block,  ny_block,
     &                         nghost,
     &                         hm    (:,:,iblk),
     &                         Tair  (:,:,iblk),
     &                         fsw   (:,:,iblk),     
     &                         cldf  (:,:,iblk),
     &                         flw   (:,:,iblk),
     &                         frain (:,:,iblk),
     &                         fsnow (:,:,iblk),
     &                         Qa    (:,:,iblk),
     &                         rhoa  (:,:,iblk),
     &                         uatm  (:,:,iblk),
     &                         vatm  (:,:,iblk),
     &                         zlvl  (:,:,iblk),
     &                         wind  (:,:,iblk),
     &                         swvdr (:,:,iblk),
     &                         swvdf (:,:,iblk),
     &                         swidr (:,:,iblk),
     &                         swidf (:,:,iblk),
     &                         potT  (:,:,iblk),
     &                         ANGLET(:,:,iblk), 
     &                         trcr  (:,:,1,iblk),
     &                         sst   (:,:,iblk),
     &                         aice  (:,:,iblk) )

      enddo                     ! iblk

      end subroutine get_forcing_atmo

!=======================================================================
!BOP
!
! !IROUTINE: get_forcing_ocn - interpolate sss, sst; restore sst
!
! !INTERFACE:
!
      subroutine get_forcing_ocn (dt)
!
! !DESCRIPTION:
!
! Interpolate monthly sss, sst data to timestep.
! Restore prognostic sst to data.
! Interpolate fields from U grid to T grid if necessary.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! !USES:
!
      use ice_domain, only: nblocks
      use ice_ocean
      use ice_flux, only: Tf, sss, sst, uocn, vocn, ss_tltx, ss_tlty
      use ice_grid, only: t2ugrid_vector
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::
     &   dt      ! time step
!
!EOP
!
      integer (kind=int_kind) ::
     &    i, j, iblk      ! horizontal indices
     &,   ixm,ixp         ! record numbers for neighboring months
     &,   maxrec          ! maximum record number
     &,   recslot         ! spline slot for current record
     &,   midmonth        ! middle day of month

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) ::
     &    sstdat              ! data value toward which SST is restored

      logical (kind=log_kind) :: readm

      if (my_task == master_task .and. istep == 1) then
         if (trim(sss_data_type)=='clim') then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'SSS interpolated to data:'
            write (nu_diag,*) trim(sss_file)
         endif
         if (trim(sst_data_type)=='clim') then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'SST restored to data:'
            write (nu_diag,*) trim(sst_file)
            write (nu_diag,*)
     &        'SST restoring timescale (days) =', trestore
         endif
!lipscomb - not implemented yet
         if (trim(sss_data_type)=='ncar' .or.
     &       trim(sst_data_type)=='ncar') then
            call abort_ice('ncar ocean forcing not implemented yet')
         endif
      endif                     ! my_task, istep

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

      if (trim(sss_data_type)=='clim' .or. 
     &    trim(sst_data_type)=='clim') then

         midmonth = 15          ! data is given on 15th of every month
!!!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

         ! Compute record numbers for surrounding months
         maxrec = 12
         ixm  = mod(month+maxrec-2,maxrec) + 1
         ixp  = mod(month,         maxrec) + 1
         if (mday >= midmonth) ixm = 99 ! other two points will be used
         if (mday <  midmonth) ixp = 99

         ! Determine whether interpolation will use values 1:2 or 2:3
         ! recslot = 2 means we use values 1:2, with the current value (2)
         !  in the second slot
         ! recslot = 1 means we use values 2:3, with the current value (2)
         !  in the first slot
         recslot = 1            ! latter half of month
         if (mday < midmonth) recslot = 2 ! first half of month

         ! Find interpolation coefficients
         call interp_coeff_monthly (recslot)

         readm = .false.
         if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

      endif   ! sss/sst_data_type

    !-------------------------------------------------------------------
    ! Read two monthly SSS values and interpolate.
    ! Note: SSS is restored instantaneously to data.
    !-------------------------------------------------------------------

      if (trim(sss_data_type)=='clim') then
         call read_clim_data (readm, 0, ixm, month, ixp,
     &                        sss_file, sss_data)
         call interpolate_data (sss_data, sss)
         sss(:,:,:) = max(sss(:,:,:), c0)
      endif

    !-------------------------------------------------------------------
    ! Read two monthly SST values and interpolate.
    ! Restore toward interpolated value.
    ! Make sure SST is not below Tf.
    !-------------------------------------------------------------------

      if (trim(sst_data_type)=='clim') then
         call read_clim_data (readm, 0, ixm, month, ixp,
     &                        sst_file, sst_data)
         call interpolate_data (sst_data, sstdat)

         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = sst(i,j,iblk) 
     &                   + (sstdat(i,j,iblk)-sst(i,j,iblk))*dt/trest
               Tf (i,j,iblk) = -depressT*sss(i,j,iblk) ! deg C
               sst(i,j,iblk) =  max (sst(i,j,iblk), Tf(i,j,iblk))
            enddo
            enddo
         enddo
      endif

    !-------------------------------------------------------------------
    ! Interpolate ocean dynamics variables from T-cell centers to
    !  U-cell centers.
    ! Typically these fields are set to zero for standalone ice runs.
    !-------------------------------------------------------------------

!      call t2ugrid_vector(uocn)
!      call t2ugrid_vector(vocn)
!      call t2ugrid_vector(ss_tltx)
!      call t2ugrid_vector(ss_tlty)

      end subroutine get_forcing_ocn

!=======================================================================
!
!BOP
!
! !IROUTINE: read_data - Read data needed for interpolation
!
! !INTERFACE:
!
      subroutine read_data (flag, recd, yr, ixm, ixx, ixp,
     &                      maxrec, data_file, field_data)
!
! !DESCRIPTION:
!
! If data is at the beginning of a one-year record, get data from
!  the previous year.
! If data is at the end of a one-year record, get data from the
!  following year.
! If no earlier data exists (beginning of fyear_init), then
!  (1) For monthly data, get data from the end of fyear_final.
!  (2) For more frequent data, let the ixm value equal the
!      first value of the year.
! If no later data exists (end of fyear_final), then
!  (1) For monthly data, get data from the beginning of fyear_init.
!  (2) For more frequent data, let the ipx value
!      equal the last value of the year.
! In other words, we assume persistence when daily or 6-hourly
!   data is missing, and we assume periodicity when monthly data
!   is missing.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_read_write
      use ice_diagnostics, only: check_step
!
! !INPUT/OUTPUT PARAMETERS:
!
      logical (kind=log_kind), intent(in) :: flag

      integer (kind=int_kind), intent(in) ::
     &   recd                    ! baseline record number
     &,  yr                      ! year of forcing data
     &,  ixm, ixx, ixp           ! record numbers of 3 data values
                                 ! relative to recd
     &,  maxrec                  ! maximum record value

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks),
     &   intent(out) ::
     &   field_data              ! 2 values needed for interpolation
!
!EOP
!
      character (char_len_long) ::
     &   data_file               ! data file to be read

      integer (kind=int_kind) ::
     &   nbits                ! = 32 for single precision, 64 for double
     &,  nrec                 ! record number to read
     &,  n2, n4               ! like ixm and ixp, but
                              ! adjusted at beginning and end of data
     &,  arg                  ! value of time argument in field_data

      nbits = 64              ! double precision data

      if (istep1 > check_step) dbug = .true.  !! debugging

      if (my_task==master_task .and. (dbug)) then
         write(nu_diag,*) '  ', trim(data_file)
      endif

      if (flag) then

      !-----------------------------------------------------------------
      ! Initialize record counters
      ! (n2, n4 will change only at the very beginning or end of
      !  a forcing cycle.)
      !-----------------------------------------------------------------
         n2 = ixm
         n4 = ixp
         arg = 0

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         if (ixm /= 99) then
      ! currently in first half of data interval
            if (ixx <= 1) then
               if (yr > fyear_init) then ! get data from previous year
                  call file_year (data_file, yr-1)
               else             ! yr = fyear_init, no prior data exists
                  if (maxrec > 12) then ! extrapolate from first record
                     if (ixx == 1) n2 = ixx
                  else          ! go to end of fyear_final
                     call file_year (data_file, fyear_final)
                  endif
               endif            ! yr > fyear_init
            endif               ! ixx <= 1
            call ice_open (nu_forcing, data_file, nbits)

            arg = 1
            nrec = recd + n2
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:),
     &                     'rda8', dbug)
            if (ixx==1 .and. my_task == master_task) close(nu_forcing)
         endif                  ! ixm ne 99

      ! always read ixx data from data file for current year
         call file_year (data_file, yr)
         call ice_open (nu_forcing, data_file, nbits)

         arg = arg + 1
         nrec = recd + ixx
         call ice_read (nu_forcing, nrec, field_data(:,:,arg,:),
     &                  'rda8', dbug)

         if (ixp /= 99) then
         ! currently in latter half of data interval
            if (ixx==maxrec) then
               if (yr < fyear_final) then ! get data from following year
                  if (my_task == master_task) close(nu_forcing)
                  call file_year (data_file, yr+1)
                  call ice_open (nu_forcing, data_file, nbits)
               else             ! yr = fyear_final, no more data exists
                  if (maxrec > 12) then ! extrapolate from ixx
                     n4 = ixx
                  else          ! go to beginning of fyear_init
                     if (my_task == master_task) close(nu_forcing)
                     call file_year (data_file, fyear_init)
                     call ice_open (nu_forcing, data_file, nbits)
                  endif
               endif            ! yr < fyear_final
            endif               ! ixx = maxrec

            arg = arg + 1
            nrec = recd + n4
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:),
     &                     'rda8', dbug)
         endif                  ! ixp /= 99

         if (my_task == master_task) close(nu_forcing)

      endif                     ! flag

      end subroutine read_data

!=======================================================================
!
!BOP
!
! !IROUTINE: read_clim_data - read annual climatological data
!
! !INTERFACE:
!
      subroutine read_clim_data (readflag, recd, ixm, ixx, ixp,
     &                           data_file, field_data)
!
! !DESCRIPTION:
!
! Read data needed for interpolation, as in read_data.
! Assume a one-year cycle of climatological data, so that there is
!  no need to get data from other years or to extrapolate data beyond
!  the forcing time period.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_read_write
      use ice_diagnostics, only: check_step
!
! !INPUT/OUTPUT PARAMETERS:
!
      logical (kind=log_kind),intent(in) :: readflag

      integer (kind=int_kind), intent(in) ::
     &  recd                ! baseline record number
     &, ixm,ixx,ixp         ! record numbers of 3 data values
                            ! relative to recd

      character (char_len_long), intent(in) ::  data_file

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks),
     &  intent(out) ::
     &  field_data         ! 2 values needed for interpolation
!
!EOP
!
      integer (kind=int_kind) ::
     &  nbits              ! = 32 for single precision, 64 for double
     &, nrec               ! record number to read
     &, arg                ! value of time argument in field_data

      nbits = 64                ! double precision data

      if (istep1 > check_step) dbug = .true.  !! debugging

      if (my_task==master_task .and. (dbug))
     &  write(nu_diag,*) '  ', trim(data_file)

      if (readflag) then

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         call ice_open (nu_forcing, data_file, nbits)

         arg = 0
         if (ixm /= 99) then
            arg = 1
            nrec = recd + ixm
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:),
     &                     'rda8', dbug)
         endif

         arg = arg + 1
         nrec = recd + ixx
         call ice_read (nu_forcing, nrec, field_data(:,:,arg,:),
     &                  'rda8', dbug)

         if (ixp /= 99) then
            arg = arg + 1
            nrec = recd + ixp
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:),
     &                     'rda8', dbug)
         endif

         if (my_task == master_task) close (nu_forcing)
      endif                     ! readflag

      end subroutine read_clim_data

!=======================================================================
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
! authors: same as module
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

      real (kind=dbl_kind) ::
     &    daymid(0:13)     ! month mid-points

      daymid(1:13) = 14._dbl_kind   ! time frame ends 0 sec into day 15
      daymid(0)    = -17._dbl_kind  ! Dec 15, 0 sec

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

!=======================================================================
!
!BOP
!
! !IROUTINE: interp_coeff
!
! !INTERFACE:
!
      subroutine interp_coeff (recnum, recslot, secint, dataloc)
!
! !DESCRIPTION:
!
! Compute coefficients for interpolating data to current time step.
! Works for any data interval that divides evenly into a 365-day
!  year (daily, 6-hourly, etc.)
! Use interp_coef_monthly for monthly data.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      integer (kind=int_kind), intent(in) ::
     &    recnum          ! record number for current data value
     &,   recslot         ! spline slot for current record
     &,   dataloc         ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval

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
     &,   rcnum            ! recnum => dbl_kind

      tt = mod(ftime,secyr)

      ! Find neighboring times
      rcnum = real(recnum,kind=dbl_kind)
      if (recslot==2) then           ! current record goes in slot 2
         if (dataloc==1) then        ! data located at middle of interval
            t2 = (rcnum-p5)*secint
         else                        !  data located at end of interval
            t2 = rcnum*secint
         endif
         t1 = t2 - secint            !  - 1 interval
      else                           ! recslot = 1
         if (dataloc==1) then        ! data located at middle of interval
            t1 = (rcnum-p5)*secint
         else                        
            t1 = rcnum*secint       ! data located at end of interval
         endif
         t2 = t1 + secint            !  + 1 interval
      endif

      ! Compute coefficients
      c1intp =  abs((t2 - tt) / (t2 - t1))
      c2intp =  c1 - c1intp

      end subroutine interp_coeff
!=======================================================================
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
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain, only: nblocks
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks),
     &  intent(in) ::
     &  field_data    ! 2 values used for interpolation

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks),
     &  intent(out) ::
     &  field         ! interpolated field
!
!EOP
!
      integer (kind=int_kind) :: i,j, iblk

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            field(i,j,iblk) = c1intp * field_data(i,j,1,iblk)
     &                      + c2intp * field_data(i,j,2,iblk)
         enddo
         enddo
      enddo

      end subroutine interpolate_data

!=======================================================================
!
!BOP
!
! !IROUTINE: file_year - construct name of atmospheric data file
!
! !INTERFACE:
!
      subroutine file_year (data_file, yr)
!
! !DESCRIPTION:
!
! Construct the correct name of the atmospheric data file
! to be read, given the year and assuming the naming convention
! that filenames end with 'yyyy.dat' or 'yyyy.r'.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
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

      if (trim(atm_data_type) == 'ecmwf') then ! NPS/ECMWF naming convention
         i = index(data_file,'.r') - 5
         tmpname = data_file
         write(data_file,'(a,i4.4,a)') tmpname(1:i), yr, '.r'
      else                                     ! LANL/NCAR naming convention
         i = index(data_file,'.dat') - 5
         tmpname = data_file
         write(data_file,'(a,i4.4,a)') tmpname(1:i), yr, '.dat'
      endif

      end subroutine file_year

!=======================================================================
!
!BOP
!
! !IROUTINE: ncar_files - construct filenames for NCAR bulk atmospheric data
!
! !INTERFACE:
!
      subroutine ncar_files (yr)
!
! !DESCRIPTION:
!
! Construct filenames based on the LANL naming conventions for NCAR data.
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file\_year will insert the correct year.
!
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
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

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Forcing data year =', fyear
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(fsw_file)
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
         write (nu_diag,*) trim(rhoa_file)
      endif                     ! master_task

      end subroutine ncar_files

!=======================================================================
!BOP
!
! !IROUTINE: ecmwf_files - construct filenames for ECMWF atmospheric data
!
! !INTERFACE:
!
      subroutine ecmwf_files (yr)
!
! !DESCRIPTION:
!
! Construct filenames based on naming conventions used by Wieslaw Maslowski
!  for reading (mostly) ECMWF atmospheric data. 
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file\_year will insert the correct year.
!
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL (based on ncar_files)
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
     &     trim(atm_data_dir)//'sol_2002.r'
      call file_year(fsw_file,yr)

      flw_file =
     &     trim(atm_data_dir)//'flo_2002.r'
      call file_year(flw_file,yr)

      rain_file =
     &     trim(atm_data_dir)//'prec_lanl_12.r'
      ! Comment out the file_year call if rain file is from climatology
!!!      call file_year(rain_file,yr)

      uwind_file =
     &     trim(atm_data_dir)//'ucmp_2002.r'
      call file_year(uwind_file,yr)

      vwind_file =
     &     trim(atm_data_dir)//'vcmp_2002.r'
      call file_year(vwind_file,yr)

      tair_file =
     &     trim(atm_data_dir)//'tair_2002.r'
      call file_year(tair_file,yr)

      humid_file =
     &     trim(atm_data_dir)//'qa_2002.r'
      call file_year(humid_file,yr)

      rhoa_file =
     &     trim(atm_data_dir)//'rhoa_ncar85-88_12.r'
      ! Comment out the file_year call if rhoa file is from climatology
!!!      call file_year(rhoa_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Forcing data year = ', fyear
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(fsw_file)
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
         write (nu_diag,*) trim(rhoa_file)
      endif                     ! master_task

      end subroutine ecmwf_files

!=======================================================================
!
!BOP
!
! !IROUTINE: LY_files - construct filenames for Large and Yeager data
!     note:  includes AOMIP (OMIP) cldf climatology 
!
! !INTERFACE:
!
      subroutine LY_files (yr)
!
! !DESCRIPTION:
!
! Construct filenames based on the LANL naming conventions for NCAR data.
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file_year will insert the correct year.
!
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
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
      flw_file =
     &     trim(atm_data_dir)//'MONTHLY/cldf.omip.dat'

      rain_file =
     &     trim(atm_data_dir)//'MONTHLY/prec.nmyr.dat'

      uwind_file =
     &     trim(atm_data_dir)//'4XDAILY/u_10.1996.dat'
      call file_year(uwind_file,yr)

      vwind_file =
     &     trim(atm_data_dir)//'4XDAILY/v_10.1996.dat'
      call file_year(vwind_file,yr)

      tair_file =
     &     trim(atm_data_dir)//'4XDAILY/t_10.1996.dat'
      call file_year(tair_file,yr)

      humid_file =
     &     trim(atm_data_dir)//'4XDAILY/q_10.1996.dat'
      call file_year(humid_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Forcing data year = ', fyear         
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
      endif                     ! master_task

      end subroutine LY_files

c=======================================================================
!
!BOP
!
! !IROUTINE: ncar_data - read NCAR bulk atmospheric data
!
! !INTERFACE:
!
      subroutine ncar_data
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) ::
     &    i, j
     &,   ixm,ixx,ixp     ! record numbers for neighboring months
     &,   recnum          ! record number
     &,   maxrec          ! maximum record number
     &,   recslot         ! spline slot for current record
     &,   dataloc         ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
     &,   midmonth        ! middle day of month

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
      ixm  = mod(month+maxrec-2,maxrec) + 1
      ixp  = mod(month,         maxrec) + 1
      if (mday >= midmonth) ixm = 99  ! other two points will be used
      if (mday <  midmonth) ixp = 99

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

      call read_data (readm, 0, fyear, ixm, month, ixp,
     &                maxrec, fsw_file, fsw_data)
      call read_data (readm, 0, fyear, ixm, month, ixp,
     &                maxrec, flw_file, cldf_data)
      call read_data (readm, 0, fyear, ixm, month, ixp,
     &                maxrec, rain_file, fsnow_data)

      ! Interpolate to current time step
      call interpolate_data (fsw_data,   fsw)
      call interpolate_data (cldf_data,  cldf)
      call interpolate_data (fsnow_data, fsnow)

    !-------------------------------------------------------------------
    ! 6-hourly data
    !
    ! Assume that the 6-hourly value is located at the end of the
    !  6-hour period.  This is the convention for NCEP reanalysis data.
    !  E.g. record 1 gives conditions at 6 am GMT on 1 January.
    !-------------------------------------------------------------------

      dataloc = 2               ! data located at end of interval
      sec6hr = secday/c4        ! seconds in 6 hours
      maxrec = 1460             ! 365*4

      ! current record number
      recnum = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec6hr)

      ! Compute record numbers for surrounding data

      ixm = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
!      ixp = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record always goes in slot 2.

      recslot = 2
      ixp = 99
      call interp_coeff (recnum, recslot, sec6hr, dataloc)

      ! Read
      read6 = .false.
      if (istep==1 .or. oldrecnum /= recnum) read6 = .true.

      call read_data (read6, 0, fyear, ixm, ixx, ixp,
     &                maxrec, tair_file, Tair_data)
      call read_data (read6, 0, fyear, ixm, ixx, ixp,
     &                maxrec, uwind_file, uatm_data)
      call read_data (read6, 0, fyear, ixm, ixx, ixp,
     &                maxrec, vwind_file, vatm_data)
      call read_data (read6, 0, fyear, ixm, ixx, ixp,
     &                maxrec, rhoa_file, rhoa_data)
      call read_data (read6, 0, fyear, ixm, ixx, ixp,
     &                maxrec, humid_file, Qa_data)

      ! Interpolate
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (uatm_data, uatm)
      call interpolate_data (vatm_data, vatm)
      call interpolate_data (rhoa_data, rhoa)
      call interpolate_data (Qa_data,   Qa)

      ! Save record number for next time step
      oldrecnum = recnum

      end subroutine ncar_data

!=======================================================================
!BOP
!
! !IROUTINE: ecmwf_data - read ECMWF atmospheric data
!
! !INTERFACE:
!
      subroutine ECMWF_data
!
! !DESCRIPTION:
!
! Read ECMWF atmospheric data.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Wieslaw Maslowski, NPS
!          Based on ncar_data
!
! !USES:
!
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) ::
     &    i, j
     &,   ixm,ixx,ixp     ! record numbers for neighboring months
     &,   recnum          ! record number
     &,   maxrec          ! maximum record number
     &,   recslot         ! spline slot for current record
     &,   dataloc         ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
     &,   midmonth        ! middle day of month

      logical (kind=log_kind) :: readm, readd

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
      ixm  = mod(month+maxrec-2,maxrec) + 1
      ixp  = mod(month,         maxrec) + 1
      if (mday >= midmonth) ixm = 99  ! other two points will be used
      if (mday <  midmonth) ixp = 99

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

      call read_clim_data (readm, 0,  ixm, month, ixp,
     &     rhoa_file, rhoa_data)
      call read_clim_data (readm, 0,  ixm, month, ixp,
     &     rain_file, fsnow_data)

      ! Interpolate to current time step
      call interpolate_data (fsnow_data, fsnow)
      call interpolate_data (rhoa_data, rhoa)  

    !-------------------------------------------------------------------
    ! Daily data
    !
    ! Assume that the daily value is located in the middle of the
    !  24-hour period.
    !-------------------------------------------------------------------

      dataloc = 1          ! data located in middle of interval
      maxrec = 365         ! days in a year (no leap years)

      ! current record number
      recnum = int(yday)   ! current record number

      ! Compute record numbers for surrounding data

      ixm = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
      ixp = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      if (real(sec,kind=dbl_kind) < p5*secday-puny) then  ! first half of day
         recslot = 2
         ixp = 99
      else                             ! second half of day
         recslot = 1
         ixm = 99
      endif

      call interp_coeff (recnum, recslot, secday, dataloc)

      ! Read new data at midpoint of day

      readd = .false.
      if (istep==1 .or. (recslot==1 .and. oldrecslot==2))
     &     readd = .true.

      call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec,
     &                tair_file, Tair_data)
      call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec,
     &                uwind_file, uatm_data)
      call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec,
     &                vwind_file, vatm_data)
      call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec, 
     &                fsw_file, fsw_data)
      call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec,
     &                flw_file, flw_data)
      call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec,
     &                humid_file, Qa_data)

      ! Interpolate
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (uatm_data, uatm)
      call interpolate_data (vatm_data, vatm)
      call interpolate_data ( fsw_data, fsw)
      call interpolate_data ( flw_data, flw)
      call interpolate_data (  Qa_data, Qa)

      ! Save recslot for next time step
      oldrecslot = recslot

      end subroutine ecmwf_data

!=======================================================================
!
!BOP
!
! !IROUTINE: LY_data - read Large and Yeager atmospheric data
!        note:  also uses AOMIP protocol, in part
!
! !INTERFACE:
!
      subroutine LY_data
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_global_reductions
      use ice_domain, only: nblocks, distrb_info
      use ice_flux 
      use ice_grid, only: hm, tlon, tlat, tmask, umask
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: 
     &    i, j
     &,   imx,ixx,ipx     ! record numbers for neighboring months
     &,   recnum          ! record number
     &,   maxrec          ! maximum record number
     &,   recslot         ! spline slot for current record
     &,   midmonth        ! middle day of month
     &,   dataloc         ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
     &,   iblk            ! block index

      real (kind=dbl_kind) ::
     &    sec6hr              ! number of seconds in 6 hours
     &,   vmin, vmax

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

      call read_clim_data (readm, 0, imx, month, ipx, 
     &       flw_file, cldf_data)
      call read_clim_data (readm, 0, imx, month, ipx, 
     &       rain_file, fsnow_data)

      call interpolate_data (cldf_data, cldf)
      call interpolate_data (fsnow_data, fsnow)  ! units mm/s = kg/m^2/s

    !-------------------------------------------------------------------
    ! 6-hourly data
    ! 
    ! Assume that the 6-hourly value is located at the end of the
    !  6-hour period.  This is the convention for NCEP reanalysis data.
    !  E.g. record 1 gives conditions at 6 am GMT on 1 January.
    !-------------------------------------------------------------------

      dataloc = 2               ! data located at end of interval
      sec6hr = secday/c4        ! seconds in 6 hours
      maxrec = 1460             ! 365*4

      ! current record number
      recnum = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec6hr)

      ! Compute record numbers for surrounding data (2 on each side)

      imx = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
!      ipx = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      recslot = 2
      ipx = 99
      call interp_coeff (recnum, recslot, sec6hr, dataloc)

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
     &                humid_file, Qa_data)

      ! Interpolate
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (uatm_data, uatm)
      call interpolate_data (vatm_data, vatm)
      call interpolate_data (Qa_data, Qa)

      call Qa_fixLY(nx_block,  ny_block,
     &                         Tair (:,:,iblk),
     &                         Qa   (:,:,iblk))

      do iblk = 1, nblocks
        do j = 1, ny_block
          do i = 1, nx_block
            Qa  (i,j,iblk) = Qa  (i,j,iblk) * hm(i,j,iblk)
            Tair(i,j,iblk) = Tair(i,j,iblk) * hm(i,j,iblk)
            uatm(i,j,iblk) = uatm(i,j,iblk) * hm(i,j,iblk)
            vatm(i,j,iblk) = vatm(i,j,iblk) * hm(i,j,iblk)
          enddo
        enddo

      ! AOMIP
      call compute_shortwave(nx_block,  ny_block, nghost,
     &                       TLON (:,:,iblk),
     &                       TLAT (:,:,iblk),
     &                       hm   (:,:,iblk),
     &                       Qa   (:,:,iblk),
     &                       cldf (:,:,iblk),
     &                       fsw  (:,:,iblk))

      enddo  ! iblk

      ! Save record number
      oldrecnum = recnum

         if (dbug) then
           if (my_task == master_task) write (nu_diag,*) 'LY_bulk_data'
           vmin = global_minval(fsw
     &                         ,distrb_info,field_loc_center,tmask)
           vmax = global_maxval(fsw
     &                         ,distrb_info,field_loc_center,tmask)
           if (my_task.eq.master_task) 
     &         write (nu_diag,*) 'fsw',vmin,vmax
           vmin = global_minval(cldf
     &                         ,distrb_info,field_loc_center,tmask)
           vmax = global_maxval(cldf
     &                         ,distrb_info,field_loc_center,tmask)
           if (my_task.eq.master_task) 
     &         write (nu_diag,*) 'cldf',vmin,vmax
           vmin =global_minval(fsnow
     &                         ,distrb_info,field_loc_center,tmask)
           vmax =global_maxval(fsnow
     &                         ,distrb_info,field_loc_center,tmask)
           if (my_task.eq.master_task) 
     &         write (nu_diag,*) 'fsnow',vmin,vmax
           vmin = global_minval(Tair
     &                         ,distrb_info,field_loc_center,tmask)
           vmax = global_maxval(Tair
     &                         ,distrb_info,field_loc_center,tmask)
           if (my_task.eq.master_task) 
     &         write (nu_diag,*) 'Tair',vmin,vmax
           vmin = global_minval(uatm
     &                         ,distrb_info,field_loc_NEcorner,umask)
           vmax = global_maxval(uatm
     &                         ,distrb_info,field_loc_NEcorner,umask)
           if (my_task.eq.master_task) 
     &         write (nu_diag,*) 'uatm',vmin,vmax
           vmin = global_minval(vatm
     &                         ,distrb_info,field_loc_NEcorner,umask)
           vmax = global_maxval(vatm
     &                         ,distrb_info,field_loc_NEcorner,umask)
           if (my_task.eq.master_task) 
     &         write (nu_diag,*) 'vatm',vmin,vmax
           vmin = global_minval(Qa
     &                         ,distrb_info,field_loc_center,tmask)
           vmax = global_maxval(Qa
     &                         ,distrb_info,field_loc_center,tmask)
           if (my_task.eq.master_task) 
     &         write (nu_diag,*) 'Qa',vmin,vmax

        endif                   ! dbug

      end subroutine LY_data

c=======================================================================

      subroutine compute_shortwave(nx_block,  ny_block, nghost,
     &                             TLON, TLAT, hm, Qa, cldf, fsw)

!---!-------------------------------------------------------------------
!---! AOMIP shortwave forcing
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block ! block dimensions
     &,  nghost             ! number of ghost cells

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::
     &   TLON, TLAT         ! longitude, latitude
     &,  Qa                 ! specific humidity
     &,  cldf               ! cloud fraction
     &,  hm                 ! land mask

      real (kind=dbl_kind), dimension(nx_block,ny_block), 
     &   intent(inout) ::
     &   fsw                ! shortwave

      real (kind=dbl_kind) ::
     &   hour_angle
     &,  solar_time
     &,  declin
     &,  cosZ
     &,  e, d
     &,  sw0
     &,  deg2rad

      integer (kind=int_kind) ::
     &   i, j
     &,  ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do j=jlo,jhi
       do i=ilo,ihi
        deg2rad = pi/180._dbl_kind
        solar_time = mod(real(sec,kind=dbl_kind),secday)/3600._dbl_kind
     &             + c12*sin(p5*TLON(i,j))
        hour_angle = (c12 - solar_time)*pi/c12
        declin = 23.44_dbl_kind*cos((172._dbl_kind-yday)
     &           * c2*pi/c365)*deg2rad
        cosZ = sin(TLAT(i,j))*sin(declin)            
     &       + cos(TLAT(i,j))*cos(declin)*cos(hour_angle)
        cosZ = max(cosZ,c0)
        e = 1.e5*Qa(i,j)/(0.622_dbl_kind + 0.378_dbl_kind*Qa(i,j))
        d = (cosZ+2.7_dbl_kind)*e*1.e-5_dbl_kind+1.085_dbl_kind*cosZ+p1
        sw0 = 1353._dbl_kind*cosZ**2/d
        sw0 = max(sw0,c0)

        ! total downward shortwave for cice
        Fsw(i,j) = sw0*(c1-p6*cldf(i,j)**3) 
        Fsw(i,j) = Fsw(i,j)*hm(i,j)
       enddo
      enddo

      end subroutine compute_shortwave

c=======================================================================

      subroutine Qa_fixLY(nx_block, ny_block, Tair, Qa)

      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block ! block dimensions

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::
     &   Tair               ! air temperature

      real (kind=dbl_kind), dimension(nx_block,ny_block), 
     &   intent(inout) ::
     &   Qa                 ! specific humidity

      real (kind=dbl_kind), dimension(nx_block, ny_block) ::
     &    work

      work = Tair - Tffresh
      work = c2 + (0.7859_dbl_kind + 0.03477_dbl_kind*work)
     &               /(c1 + 0.00412_dbl_kind*work) ! 2+ converts ea mb -> Pa
     &          + 0.00422_dbl_kind*work            ! for ice
      ! vapor pressure
      work = (c10**work)      ! saturated 
      work = max(work,puny)   ! puny over land to prevent division by zero
      ! specific humidity
      work = 0.622_dbl_kind*work/(1.e5 - 0.378_dbl_kind*work)

      Qa = min(Qa, work)

      end subroutine Qa_fixLY

c=======================================================================
!
!BOP
!
! !IROUTINE: prepare_forcing - finish manipulating forcing
!
! !INTERFACE:
!
      subroutine prepare_forcing (nx_block, ny_block,
     &                            nghost,   hm,
     &                            Tair,     fsw,     
     &                            cldf,     flw,
     &                            frain,    fsnow,
     &                            Qa,       rhoa,
     &                            uatm,     vatm,
     &                            zlvl,     wind,
     &                            swvdr,    swvdf,
     &                            swidr,    swidf,
     &                            potT,     ANGLET,
     &                            Tsfc,     sst,
     &                            aice)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block ! block dimensions
     &,  nghost             ! number of ghost cells

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::
     &   Tair        ! air temperature  (K)
     &,  ANGLET      ! ANGLE converted to T-cells
     &,  Tsfc        ! ice skin temperature
     &,  sst         ! sea surface temperature
     &,  aice        ! ice area fraction
     &,  hm          ! land mask
     
      real (kind=dbl_kind), dimension(nx_block,ny_block),
     &   intent(inout) ::
     &   fsw         ! incoming shortwave radiation (W/m^2)
     &,  cldf        ! cloud fraction
     &,  frain       ! rainfall rate (kg/m^2 s)
     &,  fsnow       ! snowfall rate (kg/m^2 s)
     &,  Qa          ! specific humidity (kg/kg)
     &,  rhoa        ! air density (kg/m^3)
     &,  uatm        ! wind speed (m/s)
     &,  vatm
     &,  zlvl        ! atm level height (m)
     &,  wind        ! wind speed (m/s)
     &,  flw         ! incoming longwave radiation (W/m^2)
     &,  swvdr       ! sw down, visible, direct  (W/m^2)
     &,  swvdf       ! sw down, visible, diffuse (W/m^2)
     &,  swidr       ! sw down, near IR, direct  (W/m^2)
     &,  swidf       ! sw down, near IR, diffuse (W/m^2)
     &,  potT        ! air potential temperature  (K)

      ! as in the dummy atm (latm)
      real (kind=dbl_kind), parameter ::
     &   frcvdr = 0.28_dbl_kind ! frac of incoming sw in vis direct band
     &,  frcvdf = 0.24_dbl_kind ! frac of incoming sw in vis diffuse band
     &,  frcidr = 0.31_dbl_kind ! frac of incoming sw in near IR direct band
     &,  frcidf = 0.17_dbl_kind ! frac of incoming sw in near IR diffuse band
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j
     &,  ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      real (kind=dbl_kind) :: workx, worky
     &,  fcc, sstk, rtea, ptem, qlwm

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do j = jlo, jhi
      do i = ilo, ihi

      !-----------------------------------------------------------------
      ! make sure interpolated values are physically realistic
      !-----------------------------------------------------------------
         cldf (i,j) = max(min(cldf(i,j),c1),c0)
         fsw  (i,j) = max(fsw(i,j),c0)
         fsnow(i,j) = max(fsnow(i,j),c0)
         rhoa (i,j) = max(rhoa(i,j),c0)
         Qa   (i,j) = max(Qa(i,j),c0)

      enddo                     ! i
      enddo                     ! j

      !-----------------------------------------------------------------
      ! calculations specific to datasets
      !-----------------------------------------------------------------

      if (trim(atm_data_type) == 'ncar') then
         do j = jlo, jhi
         do i = ilo, ihi

      !-----------------------------------------------------------------
      ! correct known biases in NCAR data (as in CCSM latm)
      !-----------------------------------------------------------------

            Qa (i,j) = Qa (i,j) * 0.94_dbl_kind
            fsw(i,j) = fsw(i,j) * 0.92_dbl_kind

      !-----------------------------------------------------------------
      ! compute downward longwave as in Parkinson and Washington (1979)
      ! (for now)
      !-----------------------------------------------------------------

            flw(i,j) = stefan_boltzmann*Tair(i,j)**4
     &               * (c1 - 0.261_dbl_kind
     &                *exp(-7.77e-4_dbl_kind*(Tffresh - Tair(i,j))**2))
     &               * (c1 + 0.275_dbl_kind*cldf(i,j))

      ! precip is in mm/month; converted to mks below

         enddo
         enddo

      elseif (trim(atm_data_type) == 'ecmwf') then

      !-----------------------------------------------------------------
      ! The following assumes that the input Qa is really dew point temp
      ! (deg K) and need to be converted to specific humidity (kg/kg).
      ! Cf. ice_atmo module.
      !-----------------------------------------------------------------

         Qa (i,j) = (qqqocn/rhoa(i,j)) * exp(-TTTocn/Qa(i,j))
         Qa (i,j) = max(Qa(i,j),c0)

      !-----------------------------------------------------------------
      ! compute downward longwave as in Parkinson and Washington (1979)
      ! (for now)
      !-----------------------------------------------------------------

         flw(i,j) = stefan_boltzmann*Tair(i,j)**4
     &            * (c1 - 0.261_dbl_kind
     &             *exp(-7.77e-4_dbl_kind*(Tffresh - Tair(i,j))**2))
     &            * (c1 + 0.275_dbl_kind*cldf(i,j))

      ! precip is in mm/month; converted to mks below

      elseif (trim(atm_data_type) == 'LYq') then

      !-----------------------------------------------------------------
      ! longwave, Rosati and Miyakoda, JPO 18, p. 1607 (1988) - sort of 
      !-----------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi
            fcc = c1 - 0.8_dbl_kind * cldf(i,j)
            sstk = (Tsfc(i,j) * aice(i,j)
     &           + sst(i,j) * (c1 - aice(i,j))) + Tffresh
            rtea = sqrt(c1000*Qa(i,j) / 
     &            (0.622_dbl_kind+0.378_dbl_kind*Qa(i,j)))
            ptem = Tair(i,j)    ! get this from stability?
            qlwm = ptem * ptem * ptem 
     &           * ( ptem*(0.39_dbl_kind-0.05_dbl_kind*rtea)*fcc 
     &                                + c4*(sstk-ptem) )
            flw(i,j) = emissivity*stefan_boltzmann * ( sstk**4 - qlwm )
            flw(i,j) = flw(i,j) * hm(i,j) ! land mask
         enddo
         enddo

      endif                     ! atm_data_type

      !-----------------------------------------------------------------
      ! Compute other fields needed by model
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi

         zlvl(i,j) = c10
         wind(i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2)
         potT(i,j) = Tair(i,j)

        ! divide shortwave into spectral bands
         swvdr(i,j) = fsw(i,j)*frcvdr        ! visible direct
         swvdf(i,j) = fsw(i,j)*frcvdf        ! visible diffuse
         swidr(i,j) = fsw(i,j)*frcidr        ! near IR direct
         swidf(i,j) = fsw(i,j)*frcidf        ! near IR diffuse
                 
        ! determine whether precip is rain or snow

        ! determine whether precip is rain or snow
        if (trim(precip_units) == 'mm_per_month') then
          fsnow(i,j) = fsnow(i,j)/2.592e+06_dbl_kind  ! mm/month -> kg/m^2 s
!       elseif (trim(precip_units) == 'mm_per_sec' .or.
!               trim(precip_units) == 'mks') then
!         no change:  mm/sec = kg/m^2 s
        endif
        frain(i,j) = c0                     
        if (Tair(i,j) >= Tffresh) then
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
         workx      = uatm(i,j) ! wind velocity, m/s
         worky      = vatm(i,j)
         uatm (i,j) = workx*cos(ANGLET(i,j)) ! convert to POP grid
     &              + worky*sin(ANGLET(i,j)) ! note uatm, vatm, wind
         vatm (i,j) = worky*cos(ANGLET(i,j)) !  are on the T-grid here
     &              - workx*sin(ANGLET(i,j))
      enddo                     ! i
      enddo                     ! j

      end subroutine prepare_forcing

!=======================================================================

      end module ice_forcing

!=======================================================================
