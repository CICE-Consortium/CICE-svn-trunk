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
! !INTERFACE:
!
      module CICE_InitMod
!
! !USES:
!
#ifdef USE_ESMF
      use esmf_mod
#endif
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
      use ice_ocean
      use ice_orbital
      use ice_shortwave
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_upwind
      use ice_transport_remap
      use ice_work

#ifdef CCSM
      use shr_msg_mod           ! for CCSM coupled runs only
      use ice_prescribed_mod
#endif

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: CICE_Initialize

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

      subroutine CICE_Initialize(CICE_Comp,  importState, exportState,
     &                           synchClock, errorCode)
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!

   !--------------------------------------------------------------------
   ! Argument types depend on whether the model is using ESMF.
   !--------------------------------------------------------------------

#ifdef USE_ESMF

      type (ESMF_GridComp), intent(inout) ::
     &     CICE_Comp            ! defined ESMF component for CICE

      type (ESMF_State), intent(inout) ::
     &     importState          ! CICE import state
     &,    exportState          ! CICE export state

      type (ESMF_Clock), intent(inout) ::
     &     synchClock           ! ESMF clock to check init time

      integer (int_kind), intent(inout) ::
     &     errorCode            ! returns an error code if any init fails

#else
! declare as integer dummy arguments

      integer (int_kind) , intent(inout) ::
     &     CICE_Comp            ! dummy argument
     &,    importState          ! dummy argument
     &,    exportState          ! dummy argument
     &,    synchClock           ! dummy argument
     &,    errorCode            ! dummy argument

#endif
!
!EOP
!BOC
!
   !--------------------------------------------------------------------
   !  local variables
   !--------------------------------------------------------------------


   !--------------------------------------------------------------------
   !  initialize return flag
   !--------------------------------------------------------------------

      errorCode = ESMF_SUCCESS

#ifdef CCSM
   !--------------------------------------------------------------------
   ! CCSM-specific stuff to redirect stdin,stdout
   !--------------------------------------------------------------------

      call shr_msg_chdir('ice')! change cwd

#endif


   !--------------------------------------------------------------------
   ! model initialization
   !--------------------------------------------------------------------

      call cice_init

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------

#ifdef CCSM
      call init_cpl             ! initialize message passing (CCSM only)
#else
      call init_forcing_atmo    ! initialize atmospheric forcing (standalone)
#endif

#ifdef USE_ESMF
   !--------------------------------------------------------------------
   !  initialize and fill the export state with initial fields
   !--------------------------------------------------------------------

      call CICE_CoupledInit(importState, exportState, errorCode)

      if (errorCode /= ESMF_Success) then
         write(nu_diag,*)
     &        '(ice) CICE_Initialize: error filling export state'
         return
      endif
#endif

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
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      call init_communicate     ! initial setup for message passing
      call input_data           ! namelist variables
      call init_work            ! work arrays

      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution
      call init_ice_timers      ! initialize all timers
      call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! grid variables

      if (advection == 'remap')
     &     call init_remap      ! grid variables for remapping transport
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

      call init_forcing_ocn(dt) ! initialize sss and sst from data
      call init_state           ! initialize the ice state
#ifdef CCSM
      if(prescribed_ice) then  ! initialize prescribed ice
         call ice_prescribed_init
      endif
#endif
      if (restart) call restartfile      ! start from restart file
      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables

      end subroutine cice_init

!=======================================================================

      end module CICE_InitMod

!=======================================================================
