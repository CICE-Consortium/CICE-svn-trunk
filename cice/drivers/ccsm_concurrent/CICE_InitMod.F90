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
      use ice_calendar
      use ice_communicate
      use ice_coupling
      use ice_diagnostics
      use ice_domain
      use ice_dyn_evp
      use ice_exit
      use ice_fileunits
      use ice_flux
      use ice_forcing
      use ice_grid
      use ice_history
      use ice_restart
      use ice_init
      use ice_itd
      use ice_kinds_mod
      use ice_mechred
      use ice_meltpond
      use ice_ocean
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
!        CCSM-coupled applications, with or without ESMF.  For other
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
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!BOC
!
   !--------------------------------------------------------------------
   ! CCSM-specific stuff to redirect stdin,stdout
   !--------------------------------------------------------------------

      call shr_msg_chdir('ice')! change cwd

   !--------------------------------------------------------------------
   ! model initialization
   !--------------------------------------------------------------------

      call cice_init

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------

      call init_cpl             ! initialize message passing (CCSM only)

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
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
!     local temporary variables

      call init_communicate     ! initial setup for message passing
      if (my_task == master_task) call shr_msg_dirio('ice')    ! redirect stdin/stdout
      call input_data           ! namelist variables
      call init_fileunits       ! fileunits
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
      if (shortwave == 'dEdd') then
         call init_orbit        ! initialize orbital parameters
         call init_dEdd         ! initialize delta-Eddington scheme
      endif
      call init_itd             ! initialize ice thickness distribution
      call calendar(time)       ! determine the initial date
      call init_state           ! initialize the ice state
      call ice_prescribed_init
      if (restart) call restartfile      ! start from restart file

      if (kpond == 1) call init_meltponds
      call init_shortwave
      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables

      write_ic = .true.        ! write initial conditions
      if(.not.prescribed_ice) call ice_write_hist(dt)
      write_ic = .false.

      end subroutine cice_init

!=======================================================================

      end module CICE_InitMod

!=======================================================================
