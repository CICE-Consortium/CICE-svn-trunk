!=======================================================================
!
!BOP
!
! !MODULE: CICE_RunMod - contains main run method for CICE
!
! !DESCRIPTION:
!
!  Contains main driver routine for time stepping of CICE.
!
! !REVISION HISTORY:
!  SVN:$Id: CICE_RunMod.F90 52 2007-01-30 18:04:24Z eclare $
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2006 ECH: moved exit timeLoop to prevent execution of unnecessary timestep
! 2006 ECH: Streamlined for efficiency 
! 2006 ECH: Converted to free source form (F90)
! 2007 BPB: Modified Delta-Eddington shortwave interface
!
! !INTERFACE:
!

      module CICE_RunMod
!
! !USES:
!
      use ice_step_mod
      use ice_age
      use ice_calendar
      use ice_diagnostics
      use ice_flux
      use ice_timers
      use ice_prescribed_mod
      use ice_coupling
      use ice_exit
      use ice_history	
      use ice_restart 	
      use ice_meltpond
      use ice_shortwave

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: CICE_Run
!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: CICE_Run - advances CICE model forward in time
!
! !DESCRIPTION:
!
!  This is the main driver routine for advancing CICE forward in time.
!  Accepts forcing fields at beginning of time step and returns the
!  model state at the end of the time step.
!
!  The module is divided into three parts that are called independently:
!  step_therm1, step_therm2, and step_dynamics.  The thermodynamics is
!  split up such that the fields needed for coupling are computed in
!  step_therm1, and the rest of the work is done in step_therm2.
!
! !REVISION HISTORY:
!
!  author Elizabeth C. Hunke, LANL
!         Philip W. Jones, LANL
!         William H. Lipscomb, LANL
!
! !INTERFACE:
!

      subroutine CICE_Run

!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!BOC
!
   !--------------------------------------------------------------------
   !  local variables
   !--------------------------------------------------------------------

      real (dbl_kind) :: &
           coupledInterval      ! time (seconds) for each coupling interval

      integer (kind=int_kind) :: k

   !--------------------------------------------------------------------
   !  initialize error code and step timer
   !--------------------------------------------------------------------

      call ice_timer_start(timer_step)   ! start timing entire run

   !--------------------------------------------------------------------
   ! timestep loop
   !--------------------------------------------------------------------

      timeLoop: do

         call ice_timer_start(timer_sndrcv)   ! timing between send-recv

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date
         call calendar(time)    ! at the end of the timestep

         call ice_timer_stop(timer_sndrcv)   ! timing between send-recv

         call from_coupler      ! get updated info from CCSM coupler

         if (stop_now >= 1) exit timeLoop

         call ice_timer_start(timer_rcvsnd)   ! timing between send-recv

         if ((istep == 1) .and. (trim(runtype) == 'startup') .and. &
            (trim(shortwave) == 'dEdd')) call init_dEdd

         call init_mass_diags   ! diagnostics per timestep

         if(prescribed_ice) then  ! read prescribed ice
            call ice_prescribed_run(idate, sec)
         endif

         call init_flux_atm

      !-----------------------------------------------------------------
      ! Pre-thermo radiation
      !-----------------------------------------------------------------

         call step_rad1 (dt)

      !-----------------------------------------------------------------
      ! thermodynamics
      !-----------------------------------------------------------------

      ! Used to be pre-coupler thermo, mostly just fluxes
         call step_therm1 (dt) 

      ! Used to be post-coupled thermo, mostly the itd.
         if (.not.prescribed_ice) then
            call step_therm2 (dt)  
         endif

      !-----------------------------------------------------------------
      ! dynamics, transport, riding
      !-----------------------------------------------------------------

         if (.not.prescribed_ice) then
            do k = 1, ndyn_dt
               call step_dynamics (dyn_dt)
            enddo
         endif ! not prescribed_ice

      !-----------------------------------------------------------------
      ! Pre-coupler radiation. For CCSM shortwave, this is just the
      ! albedos, but for Delta-Eddington, this is the full radiative
      ! calculation.
      !-----------------------------------------------------------------

         call step_rad2 (dt)

         call ice_timer_stop(timer_rcvsnd)   ! timing between recv-send

         call to_coupler        ! collect/send data to CCSM coupler

         call ice_timer_start(timer_sndrcv)   ! timing between send-recv

      !-----------------------------------------------------------------
      ! write data
      !-----------------------------------------------------------------

         call ice_timer_start(timer_readwrite)  ! reading/writing

         if (mod(istep,diagfreq) == 0) call runtime_diags(dt) ! log file

         call ice_write_hist (dt)    ! history file

         if (write_restart == 1) then
            call dumpfile ! dumps for restarting
            if (tr_iage) call write_restart_age
            if (tr_pond) call write_restart_pond
            if (trim(shortwave) == 'dEdd') call write_restart_dEdd
         endif

         call ice_timer_stop(timer_readwrite)  ! reading/writing

         call ice_timer_stop(timer_sndrcv)   ! timing between send-recv

      enddo timeLoop

   !--------------------------------------------------------------------
   ! end of timestep loop
   !--------------------------------------------------------------------

      call ice_timer_stop(timer_step)   ! end timestepping loop timer     
!
!EOC
!
      end subroutine CICE_Run

!=======================================================================

      end module CICE_RunMod

!=======================================================================
