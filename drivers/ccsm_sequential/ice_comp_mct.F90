
module ice_comp_mct

#if (defined COUP_CAM)

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: ice_comp_mct
!
! !DESCRIPTION:
! This interface 
!
! !USES:
  use seq_mct_mod
  use seq_init_mct
  use seq_flds_mod
  use seq_flds_indices
!!  use shr_kind_mod,   only : r8 => shr_kind_r8
!CICE Uses
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

  use eshr_timemgr_mod
  use shr_inputInfo_mod
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
  integer, dimension(:), allocatable ::   &
     	 perm     ! permutation array to reorder points
#if !(defined SPMD)
  integer, parameter :: mpicom = 1
#endif

!=======================================================================

contains

!=======================================================================
!BOP
!
! !IROUTINE: ice_init_mct
!
! !INTERFACE:
  subroutine ice_init_mct( ICEID, mpicom_ice, gsMap_ice, dom_i, x2i_i, i2x_i, &
                           CCSMInit, SyncClock, NLFilename )
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
    integer        , intent(in)    :: ICEID     
    integer        , intent(in)    :: mpicom_ice
    type(mct_gsMap), intent(inout) :: GSMap_ice
    type(mct_gGrid), intent(inout) :: dom_i
    type(mct_aVect), intent(inout) :: x2i_i, i2x_i

    type(shr_inputInfo_initType), intent(IN) :: CCSMInit   ! Input init object
    type(eshr_timemgr_clockType), intent(IN) :: SyncClock  ! Synchronization clock
    character(len=*), optional  , intent(in) :: NLFilename ! Namelist filename

!
! !LOCAL VARIABLES:
!
    character(len=256) :: drvarchdir! driver archive directory
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
    ! use CCSMInit to set namelist parameters needed when linking
    ! with CAM.
    !=============================================================

    call shr_inputInfo_initGetData( CCSMInit, case_name=runid )    

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
    call cice_init()
    call t_stopf ('cice_init')

    !=============================================================
    ! Initialize MCT attribute vectors and indices
    !=============================================================

    call t_startf ('cice_mct_init')
    call ice_SetGSMap_mct( mpicom_ice, ICEID, GSMap_ice ) 	

    call ice_domain_mct( mpicom_ice, GSMap_ice, dom_i )

    call mct_aVect_init(x2i_i, rList=seq_flds_x2i_fields,    &
         lsize=MCT_GSMap_lsize(GSMap_ice, mpicom_ice))
    call mct_aVect_zero(x2i_i)

    call mct_aVect_init(i2x_i, rList=seq_flds_i2x_fields,    &
         lsize=MCT_GSMap_lsize(GSMap_ice, mpicom_ice))
    call mct_aVect_zero(i2x_i)

    !=============================================================
    ! use SyncClock to reset calendar information when linking
    ! with CAM on initial start.
    !=============================================================
    if (runtype == "initial") then
!jsewall Should this use ref_ymd and ref tod to reset calendar instead?
       call eshr_timemgr_clockGet(                                          &
            SyncClock, start_ymd=start_ymd, start_tod=start_tod,            &
            ref_ymd=ref_ymd, ref_tod=ref_tod)

       write(nu_diag,*) '(ice_init_mct) idate from sync clock = ',start_ymd
       write(nu_diag,*) '(ice_init_mct) resetting idate to match sync clock'
       idate = start_ymd
       iyear = (idate/10000)                     ! integer year of basedate
       month = (idate-iyear*10000)/100           ! integer month of basedate
       mday  =  idate-iyear*10000-month*100-1    ! day of month of basedate
       !jsewall CAM starts with year "0" not year "1"
       time  = (((iyear)*daycal(13)+daycal(month)+mday)*secday) + start_tod

!jsewall
    else
     istep  = istep  + 1    ! update time step counters on restart
     istep1 = istep1 + 1
     time = time + dt       ! determine the time and date for first step
                            ! of a continuation run.
    endif

    call calendar(time)                                       ! update calendar info
    time_forc = time
!jsewall
    write_restart = 0  !set write_restart to 0 as the update of the timestep
                       !from the value read in on the restart file causes CICE
                       !to write a restart file on the first timestep of a restart
                       !run.

    call broadcast_scalar(idate,        master_task)
    call broadcast_scalar(time,         master_task)
    call broadcast_scalar(time_forc,    master_task)
    call broadcast_scalar(write_restart,master_task)
    call broadcast_scalar(year_init,    master_task) !for restart file names
                                                     !calendar is already set

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
  subroutine ice_run_mct( x2i_i, i2x_i, SyncClock )
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
                                eshr_timemgr_clockDateInSync

!
! !ARGUMENTS:
    type(mct_aVect), intent(inout) :: x2i_i
    type(mct_aVect), intent(inout) :: i2x_i
    type(eshr_timemgr_clockType), intent(IN) :: SyncClock  ! Synchronization clock

! !LOCAL VARIABLES:
    character(len=*), parameter :: SubName = "ice_run_mct"
    logical :: rstwr ! .true. ==> write a restart file
    integer :: ymd   ! Current date (YYYYMMDD)
    integer :: tod   ! Current time of day (sec)

!
! !REVISION HISTORY:
! Author: Jacob Sewall
!
!EOP
!---------------------------------------------------------------------------

    !--------------------------------------------------------------------
    ! Check that internal clock is in sync with master clock
    !--------------------------------------------------------------------
    tod = sec
    ymd = idate
    if ( .not. eshr_timemgr_clockDateInSync( SyncClock, ymd, tod ) )then
       call shr_sys_abort( SubName//":: Internal sea-ice clock not in sync with "// &
                           "Master Synchronization clock" )
    end if
   
    !-------------------------------------------------------------------
    ! run thermodynamic sea ice
    !-------------------------------------------------------------------
    
    ! Get import state
    
    call t_startf ('cice_import')
    call ice_import_mct( x2i_i )
    call t_stopf ('cice_import')
 
    !set restart flag
    rstwr = eshr_timemgr_clockAlarmIsOnRes( SyncClock )

    ! Step through the first part of the thermodynamics
    call t_startf ('cice_therm1')
    call step_therm1(dt)
    call t_stopf ('cice_therm1')
    
    ! Send export state to driver
    
    call t_startf ('cice_export')
    call ice_export_mct ( i2x_i )
    call t_stopf ('cice_export')
    
    ! Step through the second part of the thermodynamics
    
    call t_startf ('cice_therm2')
    call step_therm2(dt)
    call t_stopf ('cice_therm2')
    
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
!jsewall        
    !----------------------------------------------------------------
    ! If rstwr is true write a restart.  Restart files are written for
    ! the last timestep that CICE completed.  istep is then incremented
    ! on restart.  This matches the other components of sequential CAM.
    !----------------------------------------------------------------
    
    if (rstwr) write_restart=1  !rstwr is the restart flag from the sequential driver

    if (write_restart == 1) call dumpfile ! dumps for restarting

    call ice_timer_stop(timer_readwrite)  ! reading/writing
    
    !----------------------------------------------------------------
    ! Update time and step information at the *end* of the step
    ! which is where the sequential driver does it.  Otherwise CICE
    ! is one step ahead of the other models.
    !----------------------------------------------------------------
    
    istep  = istep  + 1    ! update time step counters
    istep1 = istep1 + 1
    time = time + dt       ! determine the time and date
    call calendar(time)    ! at the end of the timestep

    
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

   ! not sure if we need to write a final history/output file here???
   !-------------------------------------------------------------------
   ! stop timers and print timer info
   !-------------------------------------------------------------------

      call ice_timer_stop(timer_total)        ! stop timing entire run
      call ice_timer_print_all(stats=.false.) ! print timing information

      if (nu_diag /= 6) close (nu_diag) ! diagnostic output

!jsewall do *NOT* call end_run from this subroutine.  For a serial
!jsewall run it is a moot point as end_run does nothing.  But for an
!jsewall MPI run end_run will kill MPI before the sequential driver is
!jsewall ready for that to happen.

  end subroutine ice_final_mct

!=================================================================================

  subroutine ice_SetGSMap_mct( mpicom_ice, ICEID, gsMap_ice )

    !-------------------------------------------------------------------

    implicit none
    integer        , intent(in)  :: mpicom_ice
    integer        , intent(in)  :: ICEID
    type(mct_gsMap), intent(inout) :: gsMap_ice

    integer,allocatable :: gindex(:)
    integer :: lat
    integer :: lon
    integer :: i, j, iblk, n, gi
    integer :: lsize,gsize
    integer :: ier
    integer :: ilo, ihi, jlo, jhi !beginning and end of physical domain

    type(block) :: this_block      ! block information for current block

    !-------------------------------------------------------------------

    ! Build the CICE grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used
    ! in SCRIP
 
    ! number the local grid

!CICE
    lsize = block_size_x*block_size_y*nblocks
    gsize = nx_global*ny_global

    allocate(gindex(lsize),stat=ier)

!    allocate(lats(lsize))
!    allocate(lons(lsize))

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

!CICE
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
!jsewall
    integer (kind=int_kind) :: icells ! number of ocean/ice cells
    integer (kind=int_kind), dimension (nx_block*ny_block) :: indxi ! compressed indices in i
    integer (kind=int_kind), dimension (nx_block*ny_block) :: indxj ! compressed indices in i

    real(dbl_kind):: sicthk(nx_block,ny_block,max_blocks) !aggregate ice thickness(m)
    real(dbl_kind):: ifrac(nx_block,ny_block,max_blocks)  !ice fraction wrt gridcell
    real(dbl_kind):: Tsrf(nx_block,ny_block,max_blocks)   !surface temp K
    type(block) :: this_block      ! block information for current block
    !-----------------------------------------------------

    !calculate ice thickness from aice and vice. Also
    !create Tsrf from the first tracer (trcr) in ice_state.F

     do iblk = 1, nblocks
!jsewall
         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block

         !----ice thickness-----
         sicthk(i,j,iblk) = vice(i,j,iblk)/(aice(i,j,iblk)+puny)

         if(aice(i,j,iblk).gt.c1) then
         write(6,*) 'JSEWALL aice = ', aice(i,j,iblk), 'i=', i, 'j=', j
         endif
!jsewall
         !---------Initialize Tsrf to Tffresh to fill array----
         Tsrf(i,j,iblk) = Tffresh 

         !-----set icells to update Tsrf over ice cells--------
         if (tmask(i,j,iblk)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif                  ! tmask

         enddo    !i
         enddo    !j

         !--- loop over only ice cells to match ice_itd.F changes---
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            !-----surface temperature----
!jsewall
             Tsrf(i,j,iblk) = Tffresh + trcr(i,j,1,iblk) !Kelvin
          enddo   !icells
     enddo        !iblk

    !fill export state i2x_i
    !jsewall - leave all fluxes as positive down.  CAM
    !jsewall still expects varied signs but the driver
    !jsewall uses positive down convention so CAM makes
    !jsewall the change in atm_comp_mct.F90.  The comment
    !jsewall "minus" below indicates those that will need
    !jsewall to be changed if CAM/the driver changes again.  

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
         i2x_i%rAttr(index_i2x_Si_t     ,n)    = Tsrf(i,j,iblk)
         i2x_i%rAttr(index_i2x_Si_tref  ,n)    = Tref(i,j,iblk)
         i2x_i%rAttr(index_i2x_Si_qref  ,n)    = Qref(i,j,iblk)
         i2x_i%rAttr(index_i2x_Si_avsdr ,n)    = alvdr(i,j,iblk)
         i2x_i%rAttr(index_i2x_Si_anidr ,n)    = alidr(i,j,iblk)
         i2x_i%rAttr(index_i2x_Si_avsdf ,n)    = alvdf(i,j,iblk)
         i2x_i%rAttr(index_i2x_Si_anidf ,n)    = alidf(i,j,iblk)
         i2x_i%rAttr(index_i2x_Si_ifrac ,n)    = aice(i,j,iblk)
         i2x_i%rAttr(index_i2x_Si_sicthk,n)    = sicthk(i,j,iblk)

         !-------fluxes-------------------
         i2x_i%rAttr(index_i2x_Faii_taux ,n)    = strairxT(i,j,iblk) !minus
         i2x_i%rAttr(index_i2x_Faii_tauy ,n)    = strairyT(i,j,iblk) !minus
         i2x_i%rAttr(index_i2x_Faii_lat  ,n)    = flat(i,j,iblk)     !minus
         i2x_i%rAttr(index_i2x_Faii_sen  ,n)    = fsens(i,j,iblk)    !minus
         i2x_i%rAttr(index_i2x_Faii_lwup ,n)    = flwout(i,j,iblk)   !minus
         i2x_i%rAttr(index_i2x_Faii_evap ,n)    = evap(i,j,iblk)     !minus
         i2x_i%rAttr(index_i2x_Faii_swnet,n)    = fswabs(i,j,iblk)
         i2x_i%rAttr(index_i2x_Fioi_melth,n)    = fhocn(i,j,iblk)

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

      sst    (:,:,:) = c0
      frzmlt (:,:,:) = c0


!jsewall Check sign conventions, some fluxes may need to be changed
!jsewall Should all be positive down coming in, but keep checking as
!jsewall  CAM keeps changing.   

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
         sst  (i,j,iblk) =  x2i_i%rAttr(index_x2i_So_t      ,n) 

         !--- atm states-
         zlvl (i,j,iblk) =  x2i_i%rAttr(index_x2i_Sa_z      ,n)
         uatm (i,j,iblk) =  x2i_i%rAttr(index_x2i_Sa_u      ,n)
         vatm (i,j,iblk) =  x2i_i%rAttr(index_x2i_Sa_v      ,n)
         potT (i,j,iblk) =  x2i_i%rAttr(index_x2i_Sa_ptem   ,n)
         Tair (i,j,iblk) =  x2i_i%rAttr(index_x2i_Sa_tbot   ,n)
         Qa   (i,j,iblk) =  x2i_i%rAttr(index_x2i_Sa_shum   ,n)
         rhoa (i,j,iblk) =  x2i_i%rAttr(index_x2i_Sa_dens   ,n)

         !--- ocn states--
         frzmlt (i,j,iblk) =  x2i_i%rAttr(index_x2i_Fioo_q  ,n)

         !--- atm fluxes--
         swvdr(i,j,iblk) =  x2i_i%rAttr(index_x2i_Faxa_swvdr,n)
         swidr(i,j,iblk) =  x2i_i%rAttr(index_x2i_Faxa_swndr,n)
         swvdf(i,j,iblk) =  x2i_i%rAttr(index_x2i_Faxa_swvdf,n)
         swidf(i,j,iblk) =  x2i_i%rAttr(index_x2i_Faxa_swndf,n)
         flw  (i,j,iblk) =  x2i_i%rAttr(index_x2i_Faxa_lwdn ,n)
         frain(i,j,iblk) =  x2i_i%rAttr(index_x2i_Faxa_rain ,n)
         fsnow(i,j,iblk) =  x2i_i%rAttr(index_x2i_Faxa_snow ,n)

         enddo    !i
         enddo    !j

     enddo        !iblk

      !-------------------------------------------------------
      ! Update ghost cells for imported quantities
      !-------------------------------------------------------

     call update_ghost_cells(sst    , bndy_info, field_loc_center, &
                                                 field_type_scalar)
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


      ! Determine derived quantities for required fields

      do iblk = 1, nblocks
       do j = 1, ny_block
       do i = 1, nx_block
         !ocean
         sst(i,j,iblk) = sst(i,j,iblk) - Tffresh    ! sea sfc temp (C)
                                                    ! presumes import is Kelvin
         Tf (i,j,iblk) = -1.8_dbl_kind                ! hardwired for seq CAM
         !atmosphere
         wind (i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
         fsw  (i,j,iblk) = swvdr(i,j,iblk) + swvdf(i,j,iblk)   &
                         + swidr(i,j,iblk) + swidf(i,j,iblk)

       enddo                    ! i
       enddo                    ! j
      enddo                     ! iblk


   end subroutine ice_import_mct

!=======================================================================

  subroutine ice_domain_mct( mpicom_ice, gsMap_ice, dom_i )


    implicit none
    integer        , intent(in)    :: mpicom_ice
    type(mct_gsMap), intent(inout) :: gsMap_ice
    type(mct_ggrid), intent(out)   :: dom_i 

    integer :: i, j, iblk, n, gi
    integer :: lsize
    integer :: ilo, ihi, jlo, jhi !beginning and end of physical domain

    real(dbl_kind), pointer:: work_dom(:)  !temporary
    real(dbl_kind), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:)     ! temporary

    type(block) :: this_block      ! block information for current block

    !-------------------------------------------------------------------
    !
    ! Initialize mct domain type
    !
    call mct_gGrid_init( GGrid=dom_i, CoordChars="lat:lon",  &
                         OtherChars="area:mask:maxfrac",     &
                         lsize=mct_gsMap_lsize(gsMap_ice, mpicom_ice) )

    !
    ! Allocate memory
    !
    lsize = mct_gGrid_lsize(dom_i)
  
    allocate(data(lsize))
    allocate(idata(lsize))
    !
    ! Determine cam domain 
    ! Numbering scheme is: West to East and South to North to South pole
    !
    ! Initialize attribute vector with special value
    !
    idata(:) = -999
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)

    data(:) = -9999.0_dbl_kind  ! generic special value       
    call mct_gGrid_importRAttr(dom_i,"lat" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"lon" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"area",data,lsize) 

    data(:) = 0.0_dbl_kind  ! generic special value   
    call mct_gGrid_importRAttr(dom_i,"mask"   ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"maxfrac",data,lsize) 

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
               work_dom(n) = TLON(i,j,iblk) ! in radians
               data(n) = work_dom(n)*rad_to_deg !convert to degrees
!jsewall
!               write(6,*) 'n= ',n,' iblk= ',iblk,' i= ',i,' j= ',j,' TLON= ',work_dom(n)
            enddo !i
         enddo    !j
     enddo        !iblk

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
               work_dom(n) = TLAT(i,j,iblk) ! in radians
               data(n) = work_dom(n)*rad_to_deg !convert to degrees
!jsewall               write(6,*) 'n= ',n,' iblk= ',iblk,' i= ',i,' j= ',j,' TLAT= ',work_dom(n)
            enddo !i
         enddo    !j
     enddo        !iblk

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
               work_dom(n) = tarea(i,j,iblk) ! in m2
               data(n) = work_dom(n)*m2_to_km2 !convert to km2
            enddo !i
         enddo    !j
     enddo        !iblk

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
               work_dom(n) = hm(i,j,iblk) !1=ocn, 0 = land
               !convert to 0 = ocn
               if (work_dom(n) == c1) then 
                  data(n) = c0
               elseif (work_dom(n) == c0) then
                  data(n) = c1
               else
                  write(6,*) 'ERROR: ice land mask not 0 or 1'
               endif
            enddo !i
         enddo    !j
     enddo        !iblk

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

               work_dom(n) = CAMFRAC(i,j,iblk) !From CAM so 1 = land
                                               !0 = ocean. To match
                                               !maxfrac in ocn model
                                               !must set data(n) to
                                               !1.0 - CAMFRAC

               data(n) = 1._dbl_kind - work_dom(n)
            enddo !i
         enddo    !j
     enddo        !iblk

    call mct_gGrid_importRattr(dom_i,"maxfrac",data,lsize) 

    ! Permute dom_i to have ascending order

    call mct_gGrid_permute(dom_i, perm)

    deallocate(data)
    deallocate(idata)
    deallocate(work_dom)

  end subroutine ice_domain_mct

!---------------------------------------------------------------

#endif

end module ice_comp_mct

