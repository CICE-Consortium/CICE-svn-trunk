!=======================================================================
!
!BOP
!
! !MODULE: CICE_InitMod - performs CICE initialization
!
! !DESCRIPTION:
!
!  This module contains the CICE initialization routine that sets model
!  parameters and initializes the grid and CICE state variables.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
!  authors Elizabeth C. Hunke, LANL
!          William H. Lipscomb, LANL
!          Philip W. Jones, LANL
!
! 2006: Converted to free form source (F90) by Elizabeth Hunke
! 2008: E. Hunke moved ESMF code to its own driver
!
! !INTERFACE:
!
      module CICE_InitMod
!
! !USES:
!
      use ice_aerosol
      use ice_age
      use ice_calendar
      use ice_communicate
      use ice_diagnostics
      use ice_domain
      use ice_dyn_eap
      use ice_dyn_evp
      use ice_exit
      use ice_fileunits
      use ice_firstyear
      use ice_flux
      use ice_forcing
      use ice_grid
      use ice_history
      use ice_restart
      use ice_init
      use ice_itd
      use ice_kinds_mod
      use ice_mechred
      use ice_meltpond_cesm
      use ice_meltpond_lvl
      use ice_meltpond_topo
      use ice_algae
      use ice_ocean
      use ice_orbital
      use ice_lvl
      use ice_restoring
      use ice_shortwave
      use ice_therm_itd
      use ice_therm_vertical
      use ice_therm_oned
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use ice_work
      use ice_zbgc
      use ice_zsalinity
      use ice_domain_size, only: nltrcr
#ifdef popcice
      use drv_forcing, only: sst_sss
#endif

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: CICE_Initialize, cice_init

!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: CICE_Initialize - initialize CICE model
!
! !DESCRIPTION:
!
!  Initialize the basic state, grid and all necessary parameters for
!  running the CICE model.  Return the initial state in routine
!  export state.
!  Note: This initialization driver is designed for standalone and
!        CCSM-coupled applications.  For other
!        applications (e.g., standalone CAM), this driver would be
!        replaced by a different driver that calls subroutine cice_init,
!        where most of the work is done.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine CICE_Initialize
!
!EOP
!BOC
!
   !--------------------------------------------------------------------
   ! model initialization
   !--------------------------------------------------------------------

      call cice_init
!
!EOC
!
      end subroutine CICE_Initialize

!=======================================================================
!BOP
!
! !ROUTINE: cice_init - initialize CICE model
!
! !DESCRIPTION:
!
!  Initialize CICE model.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine cice_init
!
!EOP
!

      call init_communicate     ! initial setup for message passing
      call init_fileunits       ! unit numbers
      call input_data           ! namelist variables
      if (trim(runtype) == 'bering') &
         call check_finished_file  ! quit if finished file already exists
      call init_zbgc            ! vertical biogeochemistry namelist
      call init_work            ! work arrays

      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution
      call init_ice_timers      ! initialize all timers
      call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! grid variables

      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file

      if (kdyn == 2) then
         call init_eap (dyn_dt) ! define eap dynmaics parameters, variables
      else                      ! for both kdyn = 0 or 1
         call init_evp (dyn_dt) ! define evp dynamics parameters, variables
      endif

      call init_coupler_flux    ! initialize fluxes exchanged with coupler
#ifdef popcice
      call sst_sss              ! POP data for CICE initialization
#endif 
      call init_thermo_vertical ! initialize vertical thermodynamics
      call init_itd             ! initialize ice thickness distribution
      call calendar(time)       ! determine the initial date

      call init_forcing_ocn(dt) ! initialize sss and sst from data
      call init_state           ! initialize the ice state
      call init_transport       ! initialize horizontal transport
      call ice_HaloRestore_init ! restored boundary conditions

      if (restart_ext) then           ! read extended grid
         if (trim(runtype) == 'continue' .or. trim(runtype) == 'bering') then 
            ! start from core restart file
            call restartfile_ext()       ! given by pointer in ice_in
            call calendar(time)          ! update time parameters
         else if (restart) then          ! ice_ic = core restart file
            call restartfile_ext(ice_ic) !  or 'default' or 'none'
         endif         
      else                            ! read physical grid
         if (trim(runtype) == 'continue' .or. trim(runtype) == 'bering') then 
            ! start from core restart file
            call restartfile()           ! given by pointer in ice_in
            call calendar(time)          ! update time parameters
         else if (restart) then          ! ice_ic = core restart file
            call restartfile (ice_ic)    !  or 'default' or 'none'
         endif         
      endif         

      ! tracers
      if (tr_iage)      call init_age            ! ice age tracer
      if (tr_FY)        call init_FY             ! first-year area tracer
      if (tr_lvl)       call init_lvl            ! level ice tracer
      if (tr_pond_cesm) call init_meltponds_cesm ! CESM melt ponds
      if (tr_pond_lvl)  call init_meltponds_lvl  ! level-ice melt ponds
      if (tr_pond_topo) call init_meltponds_topo ! topographic melt ponds
      if (tr_aero)      call init_aerosol        ! ice aerosol
      if (tr_bgc_S)     call init_zsalinity(sss) ! salinity on bio-grid
      if (tr_bgc_NO .or. solve_bgc) call init_bgc! layer biogeochemistry

      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables

      ! Initialize shortwave components using swdn from previous timestep 
      ! if restarting. These components will be scaled to current forcing 
      ! in prep_radiation.
      if (trim(runtype) == 'continue' .or. trim(runtype) == 'bering' .or. restart) &
         call init_shortwave    ! initialize radiative transfer

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date
         call calendar(time)    ! at the end of the first timestep

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------

      call init_forcing_atmo    ! initialize atmospheric forcing (standalone)
#ifdef oned
      call init_therm_oned
#endif

#ifndef coupled
      call get_forcing_atmo     ! atmospheric forcing from data
      call get_forcing_ocn(dt)  ! ocean forcing from data
!      if (tr_aero) call faero_data          ! aerosols
      if (tr_aero) call faero_default ! aerosols
      if (tr_bgc_NO .or. tr_bgc_Sil) call get_forcing_bgc
#endif

      if (runtype == 'initial' .and. .not. restart) &
         call init_shortwave    ! initialize radiative transfer using current swdn

      call init_flux_atm        ! initialize atmosphere fluxes sent to coupler
      call init_flux_ocn        ! initialize ocean fluxes sent to coupler

      if (write_ic) call ice_write_hist(dt) ! write initial conditions 

      end subroutine cice_init

!=======================================================================

      subroutine check_finished_file()

        character(len=char_len_long) :: filename

        logical :: lexist = .false.

        if (my_task == master_task) then
           
           filename = trim(restart_dir)//"finished"
           
           inquire(file=filename, exist=lexist)
           
           if (lexist) then
              
              call abort_ice("Found already finished file - quitting")
              
           end if

        endif

      end subroutine check_finished_file

!=======================================================================

      end module CICE_InitMod

!=======================================================================
