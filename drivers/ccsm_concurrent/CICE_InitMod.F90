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
!  SVN:$Id: CICE_InitMod.F90 52 2007-01-30 18:04:24Z eclare $
!
!  authors Elizabeth C. Hunke, LANL
!          William H. Lipscomb, LANL
!          Philip W. Jones, LANL
!
! 2006: Converted to free form source (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module CICE_InitMod
!
! !USES:
!
      use ice_age
      use ice_calendar
      use ice_communicate
      use ice_coupling
      use ice_diagnostics
      use ice_domain
      use ice_dyn_evp
      use ice_exit
      use ice_fileunits
      use ice_flux
      use ice_grid
      use ice_history
      use ice_restart
      use ice_init
      use ice_itd
      use ice_kinds_mod
      use ice_mechred
      use ice_meltpond
      use ice_orbital
      use ice_shortwave
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use ice_work
      use shr_msg_mod           ! for CCSM coupled runs only
      use ice_prescribed_mod

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: CICE_Init

!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: CICE_Init - initialize CICE model
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
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
!     local temporary variables

   !--------------------------------------------------------------------
   ! CCSM-specific stuff to redirect stdin,stdout
   !--------------------------------------------------------------------

      call shr_msg_chdir('ice')! change cwd

   !--------------------------------------------------------------------
   ! model initialization
   !--------------------------------------------------------------------

      call init_communicate     ! initial setup for message passing
      if (my_task == master_task) call shr_msg_dirio('ice')    ! redirect stdin/stdout
      call init_fileunits       ! unit numbers
      call input_data           ! namelist variables
      call init_work            ! work arrays

      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution
      call init_ice_timers      ! initialize all timers
      call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! grid variables

      call init_transport       ! initialize horizontal transport
      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file
      call init_evp (dt)        ! define evp dynamics parameters, variables
      call init_coupler_flux    ! initialize fluxes exchanged with coupler
      call init_thermo_vertical ! initialize vertical thermodynamics
      if (trim(shortwave) == 'dEdd') then
         call init_orbit        ! initialize orbital parameters
      endif
      call init_itd             ! initialize ice thickness distribution
      call calendar(time)       ! determine the initial date
      call init_state           ! initialize the ice state
      call ice_prescribed_init

      if (runtype /= 'continue') then
         ! for non-continuation run, determine if should read restart file
         if (trim(inic_file) /= 'default' .and. trim(inic_file) /= 'none') then
            call restartfile(inic_file)      
         end if
      else	
         ! for continuation run, always start for restart pointer file
         call restartfile()
         call calendar(time)       ! For restart runs, use initial time
      end if

      ! tracers
      if (tr_iage) call init_age        ! ice age tracer
      if (tr_pond) call init_meltponds  ! melt ponds

      call init_shortwave       ! initialize radiative transfer
      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables

      write_ic = .true.        ! write initial conditions
      if(.not.prescribed_ice) call ice_write_hist(dt)
      write_ic = .false.

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------

      call init_cpl             ! initialize message passing (CCSM only)

      end subroutine CICE_Init

!=======================================================================

      end module CICE_InitMod

!=======================================================================
