module ice_comp_mct

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: ice_comp_mct
!
! !DESCRIPTION:
! This interface 
!
! !USES:

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use shr_inputInfo_mod
  use shr_sys_mod

  use mct_mod
  use seq_flds_mod
  use seq_flds_indices
  use seq_cdata_mod

  use eshr_timemgr_mod

  use perf_mod

  use ice_flux
  use ice_state
  use ice_domain_size
  use ice_domain
  use ice_blocks
  use ice_grid
  use ice_constants
  use ice_calendar
  use ice_timers
  use ice_kinds_mod
  use ice_init
  use ice_boundary
  use ice_prescribed_mod
  use ice_scam
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ice_init_mct
  public :: ice_run_mct
  public :: ice_final_mct
  SAVE
  private                              ! By default make data private
!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Jacob Sewall
!
!EOP
! !PRIVATE MEMBER FUNCTIONS:
  private :: ice_export_mct
  private :: ice_import_mct
  private :: ice_SetGSMap_mct
  private :: ice_domain_mct

!
! !PRIVATE VARIABLES
  integer, dimension(:), allocatable ::  perm  ! permutation array to reorder points

!=======================================================================

contains

!=======================================================================
!BOP
!
! !IROUTINE: ice_init_mct
!
! !INTERFACE:
  subroutine ice_init_mct( cdata_i, x2i_i, i2x_i, SyncClock, NLFilename )
!
! !DESCRIPTION:
! Initialize thermodynamic ice model and obtain relevant atmospheric model
! arrays back from driver 
!
! !USES:

    use CICE_InitMod
    use ice_restart, only: runid, runtype, restart_dir
    use ice_history, only: history_dir, history_file
!
! !ARGUMENTS:
    type(seq_cdata)          , intent(inout) :: cdata_i
    type(mct_aVect)          , intent(inout) :: x2i_i, i2x_i
    type(eshr_timemgr_clockType), intent(in) :: SyncClock  ! Synchronization clock
    character(len=*), optional  , intent(in) :: NLFilename ! Namelist filename
!
! !LOCAL VARIABLES:
!
    integer                               :: ICEID	
    integer                               :: mpicom_ice
    type(mct_gsMap)             , pointer :: gsMap_ice
    type(mct_gGrid)             , pointer :: dom_i
    type(shr_inputInfo_initType), pointer :: CCSMInit   ! Input init object
    integer                               :: lsize

    character(len=256) :: drvarchdir         ! driver archive directory
    integer            :: start_ymd          ! Start date (YYYYMMDD)
    integer            :: start_tod          ! start time of day (s)
    integer            :: ref_ymd            ! Reference date (YYYYMMDD)
    integer            :: ref_tod            ! reference time of day (s)
    integer            :: iyear              ! yyyy

! !REVISION HISTORY:
! Author: Jacob Sewall
!EOP
!-----------------------------------------------------------------------

    !=============================================================
    ! Set cdata pointers
    !=============================================================

    call seq_cdata_setptrs(cdata_i, ID=ICEID, mpicom=mpicom_ice, &
         gsMap=gsMap_ice, dom=dom_i, CCSMInit=CCSMInit)

    !=============================================================
    ! use CCSMInit to determine type of run
    !=============================================================

    ! Preset single column values

    single_column = .false.
    scmlat = -999.
    scmlon = -999.

    call shr_inputInfo_initGetData( CCSMInit, case_name=runid   ,  &  
                                    single_column=single_column ,  &
             	                    scmlat=scmlat,scmlon=scmlon)

    if (      shr_inputInfo_initIsStartup(  CCSMInit ) )then
       runtype = "initial"
    else if ( shr_inputInfo_initIsContinue( CCSMInit ) )then
       runtype = "continue"
    else if ( shr_inputInfo_initIsBranch(   CCSMInit ) )then
       runtype = "branch"
    end if

    !=============================================================
    ! Initialize cice because grid information is needed for
    ! creation of GSMap_ice.  cice_init also sets time manager info
    !=============================================================

    call t_startf ('cice_init')
    call cice_init( mpicom_ice )
    call t_stopf ('cice_init')

    ! use SyncClock to reset calendar information on initial start
    ! - the following logic duplicates the logic for the concurrent system - 
    ! cice_init is called then init_cpl is called where the start date is received
    ! from the flux coupler
    ! - in the following calculation for the variable time, iyear is used
    ! rather than (iyear-1) (as in init_cpl) since the sequential system permits you 
    ! to start with year "0" not year "1"
    ! - on restart run 
    !   - istep0, time and time_forc are read from restart file
    !   - istep1 is set to istep0
    !   - idate is determined from time via the call to calendar (see below)
    ! - on initial run 
    !   - iyear, month and mday obtained from sync clock
    !   - time determined from iyear, month and mday
    !   - istep0 and istep1 are set to 0 

    if (runtype == 'initial') then
       call eshr_timemgr_clockGet(                                     &
            SyncClock, start_ymd=start_ymd, start_tod=start_tod,       &
            ref_ymd=ref_ymd, ref_tod=ref_tod)

       if (ref_ymd /= start_ymd .or. ref_tod /= start_tod) then
          write(nu_diag,*) 'ice_comp_mct: ref_ymd ',ref_ymd, &
                           ' must equal start_ymd ',start_ymd
          write(nu_diag,*) 'ice_comp_mct: ref_ymd ',ref_tod, &
                           ' must equal start_ymd ',start_tod
          call shr_sys_abort()
       end if

       write(nu_diag,*) '(ice_init_mct) idate from sync clock = ', &
                        start_ymd
       write(nu_diag,*) '(ice_init_mct)   tod from sync clock = ', &
                        start_tod
       write(nu_diag,*) &
             '(ice_init_mct) resetting idate to match sync clock'

       idate = start_ymd
       iyear = (idate/10000)                     ! integer year of basedate
       month = (idate-iyear*10000)/100           ! integer month of basedate
       mday  =  idate-iyear*10000-month*100-1    ! day of month of basedate
       time  = (((iyear)*daycal(13)+daycal(month)+mday)*secday) &
             + start_tod
       call shr_sys_flush(nu_diag)
    end if
    call calendar(time)     ! update calendar info

    !=============================================================
    ! Initialize MCT attribute vectors and indices
    !=============================================================

    call t_startf ('cice_mct_init')

    ! Initialize ice gsMap

    call ice_SetGSMap_mct( mpicom_ice, ICEID, GSMap_ice ) 	
    lsize = mct_gsMap_lsize(gsMap_ice, mpicom_ice)

    ! Initialize mct ice domain (needs ice initialization info)

    call ice_domain_mct( lsize, gsMap_ice, dom_i )

    ! Inialize mct attribute vectors

    call mct_aVect_init(x2i_i, rList=seq_flds_x2i_fields, lsize=lsize)
    call mct_aVect_zero(x2i_i)

    call mct_aVect_init(i2x_i, rList=seq_flds_i2x_fields, lsize=lsize) 
    call mct_aVect_zero(i2x_i)

    !============================================================= 
    ! send intial state to driver
    !=============================================================

    call ice_export_mct (i2x_i)  !Send initial state to driver
    call t_stopf ('cice_mct_init')

  end subroutine ice_init_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ice_run_mct
!
! !INTERFACE:
  subroutine ice_run_mct( cdata_i, x2i_i, i2x_i, SyncClock )
!
! !DESCRIPTION:
! Run thermodynamic CICE
!
! !USES:
    use CICE_RunMod
    use ice_history
    use ice_restart
    use ice_diagnostics

    use eshr_timemgr_mod, only: eshr_timemgr_clockIsOnLastStep,  &
                                eshr_timemgr_clockAlarmIsOnRes,  &
                                eshr_timemgr_clockDateInSync,    &
                                eshr_timeMgr_clockGet
!
! !ARGUMENTS:
    type(seq_cdata), intent(inout) :: cdata_i
    type(mct_aVect), intent(inout) :: x2i_i
    type(mct_aVect), intent(inout) :: i2x_i
    type(eshr_timemgr_clockType), intent(IN) :: SyncClock  ! Synchronization clock

! !LOCAL VARIABLES:
    integer :: k             ! index
    logical :: rstwr         ! .true. ==> write a restart file
    integer :: ymd           ! Current date (YYYYMMDD)
    integer :: tod           ! Current time of day (sec)
    integer :: yr_sync       ! Sync current year
    integer :: mon_sync      ! Sync current month
    integer :: day_sync      ! Sync current day
    integer :: tod_sync      ! Sync current time of day (sec)
    integer :: ymd_sync      ! Current year of sync clock
    character(len=char_len_long) :: fname
    character(len=*), parameter  :: SubName = "ice_run_mct"
!
! !REVISION HISTORY:
! Author: Jacob Sewall, Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! run thermodynamic sea ice
    !-------------------------------------------------------------------
    
    ! Get import state
    
    call t_startf ('cice_import')
    call ice_import_mct( x2i_i )
    call t_stopf ('cice_import')
 
    ! Update time

    istep  = istep  + 1    ! update time step counters
    istep1 = istep1 + 1
    time = time + dt       ! determine the time and date
    call calendar(time)    ! at the end of the timestep
    
    ! If prescribed ice

    if(prescribed_ice) then  ! read prescribed ice
       call ice_prescribed_run(idate, sec)
    endif
    
     ! Step through the first part of the thermodynamics

    call t_startf ('cice_therm1')
    call step_therm1(dt)
    call t_stopf ('cice_therm1')
    
    ! Send export state to driver (this matches cpl6 logic)
    
    call t_startf ('cice_export')
    call ice_export_mct ( i2x_i )
    call t_stopf ('cice_export')
    
    ! Step through the second part of the thermodynamics and dynamics
    
    call t_startf ('cice_therm2_dyn')
    if (.not.prescribed_ice) then
       call step_therm2 (dt)  ! post-coupler thermodynamics
       do k = 1, ndyn_dt
          call step_dynamics (dyn_dt) ! dynamics, transport, ridging
       enddo
    endif ! not prescribed_ice
    call t_stopf ('cice_therm2_dyn')
    
    !-----------------------------------------------------------------
    ! write data
    !-----------------------------------------------------------------
    
    call ice_timer_start(timer_readwrite)  ! reading/writing
    
    call t_startf ('cice_diag')
    if (mod(istep,diagfreq) == 0) call runtime_diags(dt) ! log file
    call t_stopf ('cice_diag')
    
    call t_startf ('cice_hist')
    call ice_write_hist (dt)    ! history file
    call t_stopf ('cice_hist')

    rstwr = eshr_timemgr_clockAlarmIsOnRes( SyncClock )
    if (rstwr) then
       call eshr_timemgr_clockGet(SyncClock, year=yr_sync, &
          month=mon_sync, day=day_sync, CurrentTOD=tod_sync)
       fname = restart_filename(yr_sync, mon_sync, day_sync, tod_sync)
       write(nu_diag,*) &
          'ice_comp_mct: callinng dumpfile for restart filename= ',fname
       call dumpfile(fname)
    end if

    call ice_timer_stop(timer_readwrite)  ! reading/writing
    
    !--------------------------------------------------------------------
    ! Check that internal clock is in sync with master clock
    !--------------------------------------------------------------------

    tod = sec
    ymd = idate
    if ( .not. eshr_timemgr_clockDateInSync( SyncClock, ymd, tod ) )then
       call eshr_timemgr_clockGet( Syncclock, CurrentYMD=ymd_sync, &
          CurrentTOD=tod_sync )
       write(nu_diag,*)' cice ymd=',ymd     ,'  cice tod= ',tod
       write(nu_diag,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
       call shr_sys_abort( SubName// &
          ":: Internal sea-ice clock not in sync with Sync Clock")
    end if
   
  end subroutine ice_run_mct

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ice_final_mct
!
! !INTERFACE:
  subroutine ice_final_mct( )
!
! !DESCRIPTION:
! Finalize CICE
!
!
! !USES:
    use ice_exit
    use ice_fileunits
!
!------------------------------------------------------------------------------
!BOP
!
! !ARGUMENTS:
!
! !REVISION HISTORY:
! Author: Jacob Sewall
!
!EOP
!---------------------------------------------------------------------------

   !-------------------------------------------------------------------
   ! stop timers and print timer info
   !-------------------------------------------------------------------

    call ice_timer_stop(timer_total)        ! stop timing entire run
    call ice_timer_print_all(stats=.false.) ! print timing information
    
    if (nu_diag /= 6) close (nu_diag) ! diagnostic output
    
    ! do *NOT* call end_run from this subroutine.  For a serial
    ! run it is a moot point as end_run does nothing.  But for an
    ! MPI run end_run will kill MPI before the sequential driver is
    ! ready for that to happen.
      
  end subroutine ice_final_mct

!=================================================================================

  subroutine ice_SetGSMap_mct( mpicom_ice, ICEID, gsMap_ice )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    implicit none
    integer        , intent(in)    :: mpicom_ice
    integer        , intent(in)    :: ICEID
    type(mct_gsMap), intent(inout) :: gsMap_ice
    !
    ! Local variables
    !
    integer,allocatable :: gindex(:)
    integer     :: lat
    integer     :: lon
    integer     :: i, j, iblk, n, gi
    integer     :: lsize,gsize
    integer     :: ier
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    !-------------------------------------------------------------------

    ! Build the CICE grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used
    ! in SCRIP

    ! number the local grid

    lsize = block_size_x*block_size_y*nblocks
    gsize = nx_global*ny_global

    allocate(gindex(lsize),stat=ier)
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             lon = this_block%i_glob(i)
             lat = this_block%j_glob(j)
             gi = (lat-1)*nx_global + lon
             gindex(n) = gi
          enddo !i
       enddo    !j
    enddo        !iblk
    
    allocate(perm(lsize),stat=ier)

    ! reorder gindex to be in ascending order.

    ! initialize a permutation array

    call mct_indexset(perm)

    ! derive a permutation that puts gindex in ascending order
    ! the default for IndexSort is Ascending.

    call mct_indexsort(lsize,perm,gindex)

    ! Sort gindex in-place

    call mct_permute(gindex,perm,lsize)

    call mct_gsMap_init( gsMap_ice, gindex, mpicom_ice, ICEID, lsize, gsize )

    deallocate(gindex)

  end subroutine ice_SetGSMap_mct


!====================================================================================

  subroutine ice_export_mct( i2x_i )   

    !-----------------------------------------------------
    type(mct_aVect)   , intent(inout) :: i2x_i

    integer :: i, j, iblk, n, ij 
    integer :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    integer (kind=int_kind)                                :: icells ! number of ocean/ice cells
    integer (kind=int_kind), dimension (nx_block*ny_block) :: indxi  ! compressed indices in i
    integer (kind=int_kind), dimension (nx_block*ny_block) :: indxj  ! compressed indices in i

    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
        Tsrf  &      ! surface temperature
     ,  tauxa &      ! atmo/ice stress
     ,  tauya &
     ,  tauxo &      ! ice/ocean stress
     ,  tauyo &
     ,  sicthk &     ! needed for cam/som only 
     ,  ailohi       ! fractional ice area

    real (kind=dbl_kind) :: &
       gsum, workx, worky           ! tmps for converting grid

    type(block)        :: this_block                           ! block information for current block
    logical :: flag
    !-----------------------------------------------------

    flag=.false.

    !calculate ice thickness from aice and vice. Also
    !create Tsrf from the first tracer (trcr) in ice_state.F

    do iblk = 1, nblocks
       do j = 1, ny_block
       do i = 1, nx_block
             
          ! sea-ice thickness needed for cam-som only
          sicthk(i,j,iblk) = vice(i,j,iblk)/(aice(i,j,iblk)+puny)

          ! ice fraction
          ailohi(i,j,iblk) = aice(i,j,iblk)

          if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then   !??? ask Dave about this extra logic - trcr failed
          ! surface temperature
          Tsrf(i,j,iblk)  = Tffresh + trcr(i,j,1,iblk)             !Kelvin (original ???)
          
          ! wind stress  (on POP T-grid:  convert to lat-lon)
          workx = strairxT(i,j,iblk)                             ! N/m^2
          worky = strairyT(i,j,iblk)                             ! N/m^2
          tauxa(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                          - worky*sin(ANGLET(i,j,iblk))
          tauya(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                          + workx*sin(ANGLET(i,j,iblk))
          
          ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
          workx = -strocnxT(i,j,iblk)                            ! N/m^2
          worky = -strocnyT(i,j,iblk)                            ! N/m^2
          tauxo(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                          - worky*sin(ANGLET(i,j,iblk))
          tauyo(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                          + workx*sin(ANGLET(i,j,iblk))
          endif

       enddo
       enddo
    enddo

    do iblk = 1, nblocks
       do j = 1, ny_block
       do i = 1, nx_block
          if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
             flag = .true.
          endif
       end do
       end do
    end do
    if (flag) then
       do iblk = 1, nblocks
          do j = 1, ny_block
          do i = 1, nx_block
             if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                write(nu_diag,*) &
                     ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j,iblk)
                call shr_sys_flush(nu_diag)
             endif
          end do
          end do
       end do
    endif

    ! Fill export state i2x_i

     n=0
     do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            n = n+1

            !-------states-------------------- 
            i2x_i%rAttr(index_i2x_Si_sicthk,n)    = sicthk(i,j,iblk) ! (needed by CAM/SOM only)
            i2x_i%rAttr(index_i2x_Si_ifrac ,n)    = ailohi(i,j,iblk)   

            if ( tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then
               !-------states-------------------- 
               i2x_i%rAttr(index_i2x_Si_t     ,n)    = Tsrf(i,j,iblk)
               i2x_i%rAttr(index_i2x_Si_avsdr ,n)    = alvdr(i,j,iblk)
               i2x_i%rAttr(index_i2x_Si_anidr ,n)    = alidr(i,j,iblk)
               i2x_i%rAttr(index_i2x_Si_avsdf ,n)    = alvdf(i,j,iblk)
               i2x_i%rAttr(index_i2x_Si_anidf ,n)    = alidf(i,j,iblk)
               i2x_i%rAttr(index_i2x_Si_tref  ,n)    = Tref(i,j,iblk)
               i2x_i%rAttr(index_i2x_Si_qref  ,n)    = Qref(i,j,iblk)
            
               !--- a/i fluxes computed by ice
               i2x_i%rAttr(index_i2x_Faii_taux ,n)   = tauxa(i,j,iblk)    
               i2x_i%rAttr(index_i2x_Faii_tauy ,n)   = tauya(i,j,iblk)    
               i2x_i%rAttr(index_i2x_Faii_lat  ,n)   = flat(i,j,iblk)     
               i2x_i%rAttr(index_i2x_Faii_sen  ,n)   = fsens(i,j,iblk)    
               i2x_i%rAttr(index_i2x_Faii_lwup ,n)   = flwout(i,j,iblk)   
               i2x_i%rAttr(index_i2x_Faii_evap ,n)   = evap(i,j,iblk)     
               i2x_i%rAttr(index_i2x_Faii_swnet,n)   = fswabs(i,j,iblk)
            
               !--- i/o fluxes computed by ice
               i2x_i%rAttr(index_i2x_Fioi_melth,n)   = fhocn(i,j,iblk)
               i2x_i%rAttr(index_i2x_Fioi_swpen,n)   = fswthru(i,j,iblk) ! hf from melting          
               i2x_i%rAttr(index_i2x_Fioi_melth,n)   = fhocn(i,j,iblk)   
               i2x_i%rAttr(index_i2x_Fioi_meltw,n)   = fresh(i,j,iblk)   ! h2o flux from melting    ???
               i2x_i%rAttr(index_i2x_Fioi_salt ,n)   = fsalt(i,j,iblk)   ! salt flux from melting   ???
               i2x_i%rAttr(index_i2x_Fioi_taux ,n)   = tauxo(i,j,iblk)   ! stress : i/o zonal       ???
               i2x_i%rAttr(index_i2x_Fioi_tauy ,n)   = tauyo(i,j,iblk)   ! stress : i/o meridional  ???
            end if
         enddo    !i
         enddo    !j
     enddo        !iblk


    ! permute before using the Rearrange call.

    call mct_aVect_permute(i2x_i,perm)

  end subroutine ice_export_mct

!====================================================================================

  subroutine ice_import_mct( x2i_i )

    !-----------------------------------------------------

    implicit none
    type(mct_aVect)   , intent(inout) :: x2i_i

    integer :: i, j, iblk, n
    integer :: ilo, ihi, jlo, jhi !beginning and end of physical domain
    type(block) :: this_block      ! block information for current block

    real (kind=dbl_kind) :: &
         gsum, workx, worky
    !-----------------------------------------------------

    ! unpermute

    call mct_aVect_unpermute(x2i_i, perm)

    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2 which is what CICE requires.
    ! Note also that the read in below includes only values needed
    ! by the thermodynamic component of CICE.  Variables uocn, vocn,
    ! ss_tltx, and ss_tlty are excluded. Also, because the SOM and
    ! DOM don't  compute SSS.   SSS is not read in and is left at
    ! the initilized value (see ice_flux.F init_coupler_flux) of
    ! 34 ppt

    !-----------------------------------------------------------------
    ! Zero stuff while waiting, only filling in active cells.
    !-----------------------------------------------------------------
    
    zlvl   (:,:,:) = c0
    uatm   (:,:,:) = c0
    vatm   (:,:,:) = c0
    potT   (:,:,:) = c0
    Tair   (:,:,:) = c0
    Qa     (:,:,:) = c0
    rhoa   (:,:,:) = c0
    swvdr  (:,:,:) = c0
    swvdf  (:,:,:) = c0
    swidr  (:,:,:) = c0
    swidf  (:,:,:) = c0
    flw    (:,:,:) = c0
    frain  (:,:,:) = c0
    fsnow  (:,:,:) = c0
    uocn   (:,:,:) = c0
    vocn   (:,:,:) = c0
    sst    (:,:,:) = c0
    ss_tltx(:,:,:) = c0
    ss_tlty(:,:,:) = c0
    frzmlt (:,:,:) = c0

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo, jhi
       do i = ilo, ihi

          n = n+1
          !--- ocn states--
          sst  (i,j,iblk)   = x2i_i%rAttr(index_x2i_So_t   ,n) 
          sss  (i,j,iblk)   = x2i_i%rAttr(index_x2i_So_s   ,n)
          uocn (i,j,iblk)   = x2i_i%rAttr(index_x2i_So_u   ,n)
          vocn (i,j,iblk)   = x2i_i%rAttr(index_x2i_So_v   ,n)
          
          !--- atm states-
          zlvl (i,j,iblk)   = x2i_i%rAttr(index_x2i_Sa_z   ,n)
          uatm (i,j,iblk)   = x2i_i%rAttr(index_x2i_Sa_u   ,n)
          vatm (i,j,iblk)   = x2i_i%rAttr(index_x2i_Sa_v   ,n)
          potT (i,j,iblk)   = x2i_i%rAttr(index_x2i_Sa_ptem,n)
          Tair (i,j,iblk)   = x2i_i%rAttr(index_x2i_Sa_tbot,n)
          Qa   (i,j,iblk)   = x2i_i%rAttr(index_x2i_Sa_shum,n)
          rhoa (i,j,iblk)   = x2i_i%rAttr(index_x2i_Sa_dens,n)
          
          !--- ocn states--
          ss_tltx(i,j,iblk) = x2i_i%rAttr(index_x2i_So_dhdx ,n)
          ss_tlty(i,j,iblk) = x2i_i%rAttr(index_x2i_So_dhdy ,n)
          frzmlt (i,j,iblk) = x2i_i%rAttr(index_x2i_Fioo_q  ,n)
          
          !--- atm fluxes--
          swvdr(i,j,iblk)   = x2i_i%rAttr(index_x2i_Faxa_swvdr,n)
          swidr(i,j,iblk)   = x2i_i%rAttr(index_x2i_Faxa_swndr,n)
          swvdf(i,j,iblk)   = x2i_i%rAttr(index_x2i_Faxa_swvdf,n)
          swidf(i,j,iblk)   = x2i_i%rAttr(index_x2i_Faxa_swndf,n)
          flw  (i,j,iblk)   = x2i_i%rAttr(index_x2i_Faxa_lwdn ,n)
          frain(i,j,iblk)   = x2i_i%rAttr(index_x2i_Faxa_rain ,n)
          fsnow(i,j,iblk)   = x2i_i%rAttr(index_x2i_Faxa_snow ,n)

         enddo    !i
         enddo    !j

     enddo        !iblk

     !-------------------------------------------------------
     ! Update ghost cells for imported quantities
     !-------------------------------------------------------

     call update_ghost_cells(sst    , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(sss    , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(uocn   , bndy_info, field_loc_center, &
                                                 field_type_vector)
     call update_ghost_cells(vocn   , bndy_info, field_loc_center, &
                                                 field_type_vector)
     call update_ghost_cells(zlvl   , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(uatm   , bndy_info, field_loc_center, &
                                                 field_type_vector)
     call update_ghost_cells(vatm   , bndy_info, field_loc_center, &
                                                 field_type_vector)
     call update_ghost_cells(potT   , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(Tair   , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(Qa     , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(rhoa   , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(ss_tltx, bndy_info, field_loc_center, &
                                                 field_type_vector)
     call update_ghost_cells(ss_tlty, bndy_info, field_loc_center, &
                                                 field_type_vector)
     call update_ghost_cells(frzmlt , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(swvdr  , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(swidr  , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(swvdf  , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(swidf  , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(flw    , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(frain  , bndy_info, field_loc_center, &
                                                 field_type_scalar)
     call update_ghost_cells(fsnow  , bndy_info, field_loc_center, &
                                                 field_type_scalar)

      !-----------------------------------------------------------------
      ! rotate zonal/meridional vectors to local coordinates
      ! compute data derived quantities
      !-----------------------------------------------------------------

      ! Vector fields come in on T grid, but are oriented geographically
      ! need to rotate to pop-grid FIRST using ANGLET
      ! then interpolate to the U-cell centers  (otherwise we
      ! interpolate across the pole)
      ! use ANGLET which is on the T grid !

      do iblk = 1, nblocks

       do j = 1,ny_block
       do i = 1,nx_block

          ! ocean
          workx      = uocn  (i,j,iblk) ! currents, m/s 
          worky      = vocn  (i,j,iblk)
          uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                         + worky*sin(ANGLET(i,j,iblk))
          vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                         - workx*sin(ANGLET(i,j,iblk))

          workx      = ss_tltx  (i,j,iblk)           ! sea sfc tilt, m/m
          worky      = ss_tlty  (i,j,iblk)
          ss_tltx(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                            + worky*sin(ANGLET(i,j,iblk))
          ss_tlty(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                            - workx*sin(ANGLET(i,j,iblk))

          sst(i,j,iblk) = sst(i,j,iblk) - Tffresh       ! sea sfc temp (C)
          Tf (i,j,iblk) = -1.8_dbl_kind                 ! hardwired for NCOM
!         Tf (i,j,iblk) = -depressT*sss(i,j,iblk)       ! freezing temp (C)
!         Tf (i,j,iblk) = -depressT*max(sss(i,j,iblk),ice_ref_salinity)

       enddo
       enddo
      enddo

      ! Interpolate ocean dynamics variables from T-cell centers to 
      ! U-cell centers.

      call t2ugrid_vector(uocn)
      call t2ugrid_vector(vocn)
      call t2ugrid_vector(ss_tltx)
      call t2ugrid_vector(ss_tlty)

      ! Atmosphere variables are needed in T cell centers in
      ! subroutine stability and are interpolated to the U grid
      ! later as necessary.

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
  
         ! atmosphere
         workx      = uatm(i,j,iblk) ! wind velocity, m/s
         worky      = vatm(i,j,iblk) 
         uatm (i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                         + worky*sin(ANGLET(i,j,iblk))   ! note uatm, vatm, wind
         vatm (i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & ! are on the T-grid here
                         - workx*sin(ANGLET(i,j,iblk))

         wind (i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
         fsw  (i,j,iblk) = swvdr(i,j,iblk) + swvdf(i,j,iblk) &
                         + swidr(i,j,iblk) + swidf(i,j,iblk)
         enddo
         enddo
      enddo

    call mct_aVect_permute(x2i_i, perm)

   end subroutine ice_import_mct

!=======================================================================

  subroutine ice_domain_mct( lsize, gsMap_i, dom_i )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_i
    type(mct_ggrid), intent(inout) :: dom_i     
    !
    ! Local Variables
    !
    integer :: i, j, iblk, n, gi           ! indices
    integer :: ilo, ihi, jlo, jhi          ! beginning and end of physical domain
    real(dbl_kind), pointer :: work_dom(:) ! temporary
    real(dbl_kind), pointer :: data(:)     ! temporary
    integer       , pointer :: idata(:)    ! temporary
    type(block)             :: this_block  ! block information for current block
    !-------------------------------------------------------------------
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (ocean), 0 (non-ocean)
    !
    call mct_gGrid_init( GGrid=dom_i, CoordChars="lat:lon", OtherChars="area:aream:mask:frac", lsize=lsize )
    !  
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_i, my_task, idata)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_i,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_i,"mask",data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"frac",data,lsize) 
    !
    ! Fill in correct values for domain components
    !
    allocate(work_dom(lsize)) 
    work_dom(:) = 0.0_dbl_kind

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          data(n) = TLON(i,j,iblk)*rad_to_deg 
       enddo    !i
       enddo    !j
    enddo       !iblk
    call mct_gGrid_importRattr(dom_i,"lon",data,lsize) 

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          data(n) = TLAT(i,j,iblk)*rad_to_deg 
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"lat",data,lsize) 

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          data(n) = tarea(i,j,iblk)/(radius*radius)
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"area",data,lsize) 

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          data(n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"mask",data,lsize) 

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          data(n) = 1._dbl_kind
       enddo   !i
       enddo   !j
    enddo      !iblk
    call mct_gGrid_importRattr(dom_i,"frac",data,lsize) 

    ! Permute dom_i to have ascending order

    call mct_gGrid_permute(dom_i, perm)

    deallocate(data)
    deallocate(idata)
    deallocate(work_dom)

  end subroutine ice_domain_mct

!=======================================================================
! BOP
!
! !ROUTINE: restart_filename
!
! !INTERFACE:
  character(len=char_len_long) function restart_filename( yr_spec, mon_spec, day_spec, sec_spec )
!
! !DESCRIPTION: Create a filename from a filename specifyer. Interpret filename specifyer 
! string with: 
! %c for case, 
! %t for optional number argument sent into function
! %y for year
! %m for month
! %d for day
! %s for second
! %% for the "%" character
! If the filename specifyer has spaces " ", they will be trimmed out
! of the resulting filename.
!
! !USES:
    use ice_restart, only: runid
!
! !INPUT/OUTPUT PARAMETERS:
  integer         , intent(in)  :: yr_spec         ! Simulation year
  integer         , intent(in)  :: mon_spec        ! Simulation month
  integer         , intent(in)  :: day_spec        ! Simulation day
  integer         , intent(in)  :: sec_spec        ! Seconds into current simulation day
!
! EOP
!
  integer             :: i, n      ! Loop variables
  integer             :: year      ! Simulation year
  integer             :: month     ! Simulation month
  integer             :: day       ! Simulation day
  integer             :: ncsec     ! Seconds into current simulation day
  character(len=char_len_long) :: string    ! Temporary character string 
  character(len=char_len_long) :: format    ! Format character string 
  character(len=char_len_long) :: filename_spec = '%c.cice.r.%y-%m-%d-%s' ! ice restarts

  !-----------------------------------------------------------------
  ! Determine year, month, day and sec to put in filename
  !-----------------------------------------------------------------

  if ( len_trim(filename_spec) == 0 )then
     call shr_sys_abort ('restart_filename: filename specifier is empty')
  end if
  if ( index(trim(filename_spec)," ") /= 0 )then
     call shr_sys_abort ('restart_filename: filename specifier can not contain a space:'//trim(filename_spec))
  end if

  year  = yr_spec
  month = mon_spec
  day   = day_spec
  ncsec = sec_spec

  ! Go through each character in the filename specifyer and interpret if special string

  i = 1
  restart_filename = ''
  do while ( i <= len_trim(filename_spec) )
     if ( filename_spec(i:i) == "%" )then
        i = i + 1
        select case( filename_spec(i:i) )
        case( 'c' )   ! runid
           string = trim(runid)
        case( 'y' )   ! year
           if ( year > 99999   ) then
              format = '(i6.6)'
           else if ( year > 9999    ) then
              format = '(i5.5)'
           else
              format = '(i4.4)'
           end if
           write(string,format) year
        case( 'm' )   ! month
           write(string,'(i2.2)') month
        case( 'd' )   ! day
           write(string,'(i2.2)') day
        case( 's' )   ! second
           write(string,'(i5.5)') ncsec
        case( '%' )   ! percent character
           string = "%"
        case default
           call shr_sys_abort ('restart_filename: Invalid expansion character: '//filename_spec(i:i))
        end select
     else
        n = index( filename_spec(i:), "%" )
        if ( n == 0 ) n = len_trim( filename_spec(i:) ) + 1
        if ( n == 0 ) exit 
        string = filename_spec(i:n+i-2)
        i = n + i - 2
     end if
     if ( len_trim(restart_filename) == 0 )then
        restart_filename = trim(string)
     else
        if ( (len_trim(restart_filename)+len_trim(string)) >= char_len_long )then
           call shr_sys_abort ('restart_filename Resultant filename too long')
        end if
        restart_filename = trim(restart_filename) // trim(string)
     end if
     i = i + 1
  end do
  if ( len_trim(restart_filename) == 0 )then
     call shr_sys_abort ('restart_filename: Resulting filename is empty')
  end if

end function restart_filename

end module ice_comp_mct

