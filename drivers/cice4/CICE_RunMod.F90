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
!  SVN:$Id$
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
#ifdef USE_ESMF
      use esmf_mod
#endif
      use ice_age
      use ice_atmo
      use ice_calendar
      use ice_communicate
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
      use ice_state
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use ice_work

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: CICE_Run, step_therm1, step_therm2, step_dynamics, ice_step
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

      subroutine CICE_Run(CICE_Comp,  importState, exportState, &
                          synchClock, errorCode)

!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:

#ifdef USE_ESMF

      type (ESMF_GridComp), intent(inout) :: &
           CICE_Comp            ! defined ESMF component for CICE

      type (ESMF_State), intent(in) :: &
           importState          ! CICE import state (forcing data)

      type (ESMF_State), intent(out) :: &
           exportState          ! CICE export state (surface fields)

      type (ESMF_Clock), intent(inout) :: &
           synchClock           ! ESMF clock to check init time

      integer (int_kind), intent(inout) :: &
           errorCode            ! returns an error code if any init fails

#else
! declare as integer dummy arguments

      integer (int_kind), intent(inout), optional :: &
           CICE_Comp        , & ! dummy argument
           importState      , & ! dummy argument
           exportState      , & ! dummy argument
           synchClock       , & ! dummy argument
           errorCode            ! dummy argument

#endif

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

#ifdef ESMF  
      errorCode = ESMF_SUCCESS

   !--------------------------------------------------------------------
   !  check clock to make sure models agree on time info
   !--------------------------------------------------------------------
   !
   !   call CICE_CoupledCheckTime(synchClock, 'start', errorCode)
   !
   !   if (errorCode /= ESMF_SUCCESS) then
   !      write (nu_diag,*) &
   !         '(ice) CICE_Run: error in clock synchronization'
   !      return
   !   endif
   !
   !UNCOMMENT THESE TWO SECTIONS FOR COUPLED MODEL SIMULATIONS WITH ESMF
   !--------------------------------------------------------------------
   !  reset timer to stop when requested by input clock
   !--------------------------------------------------------------------
   !
   ! Use appropriate ESMF Clock query functions to determine the
   ! coupling interval in seconds.  Because this depends on how the
   ! coupled model driver is using an ESMF clock, this is left to
   ! the user of the coupled model.
   !
   !coupledInterval = ?
   !
   ! reset number of time steps based on coupling interval
   !
   !npt = nint(coupledInterval/dt)
   !
   !--------------------------------------------------------------------
   !  extract data from import state
   !--------------------------------------------------------------------
   
      call CICE_CoupledExtractImport(importState, errorCode)
   
      if (errorCode /= ESMF_SUCCESS) then
         write(nu_diag,*) &
              '(ice) CICE_Run: error extracting data from import state'
         return
      endif
#endif

   !--------------------------------------------------------------------
   ! timestep loop
   !--------------------------------------------------------------------

      timeLoop: do

         call ice_step

#ifdef USE_ESMF
         call CICE_CoupledAccumulateExport(errorCode)
#endif

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date
         call calendar(time)    ! at the end of the timestep

         if (stop_now >= 1) exit timeLoop

#ifndef coupled
         call ice_timer_start(timer_couple)  ! atm/ocn coupling
         call get_forcing_atmo     ! atmospheric forcing from data
         call get_forcing_ocn(dt)  ! ocean forcing from data
         call ice_timer_stop(timer_couple)   ! atm/ocn coupling
#endif

         call init_flux_atm     ! initialize atmosphere fluxes sent to coupler
         call init_flux_ocn     ! initialize ocean fluxes sent to coupler

      enddo timeLoop

   !--------------------------------------------------------------------
   ! end of timestep loop
   !--------------------------------------------------------------------

      call ice_timer_stop(timer_step)   ! end timestepping loop timer     


#ifdef USE_ESMF
   !--------------------------------------------------------------------
   !  check clock to make sure models agree on time info
   !--------------------------------------------------------------------
   !
   !   call CICE_CoupledCheckTime(synchClock, 'stop', errorCode)
   !
   !   if (errorCode /= ESMF_SUCCESS) then
   !      write(nu_diag,*)
   !          '(ice) CICE_Run: error in stop time synchronization'
   !      return
   !   endif
   !
   !--------------------------------------------------------------------
   !  fill the export state
   !--------------------------------------------------------------------

      call CICE_CoupledFillExport(exportState, errorCode)

      if (errorCode /= ESMF_SUCCESS) then
         write(nu_diag,*) &
              '(ice) CICE_Run: error filling export state'
         return
      endif
#endif
!
!EOC
!

      end subroutine CICE_Run

!=======================================================================

      subroutine ice_step

      use ice_restoring, only: restore_ice, ice_HaloRestore

      integer (kind=int_kind) :: k

      !-----------------------------------------------------------------
      ! restoring on grid boundaries
      !-----------------------------------------------------------------

         if (restore_ice) call ice_HaloRestore

      !-----------------------------------------------------------------
      ! initialize diagnostics
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics/history
         call init_mass_diags   ! diagnostics per timestep
         call ice_timer_stop(timer_diags)   ! diagnostics/history

      !-----------------------------------------------------------------
      ! thermodynamics
      !-----------------------------------------------------------------

         call step_therm1 (dt)  ! pre-coupler thermodynamics

         call step_therm2 (dt)  ! post-coupler thermodynamics

      !-----------------------------------------------------------------
      ! dynamics, transport, ridging
      !-----------------------------------------------------------------

         do k = 1, ndyn_dt
            call step_dynamics (dyn_dt)
         enddo

      !-----------------------------------------------------------------
      ! get ready for coupling
      !-----------------------------------------------------------------

         call coupling_prep

      !-----------------------------------------------------------------
      ! write data
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics
         if (mod(istep,diagfreq) == 0) call runtime_diags(dt) ! log file
         call ice_timer_stop(timer_diags)   ! diagnostics

         call ice_timer_start(timer_hist)   ! history
         call ice_write_hist (dt)           ! history file
         call ice_timer_stop(timer_hist)    ! history

         call ice_timer_start(timer_readwrite)  ! reading/writing
         if (write_restart == 1) then
            call dumpfile ! core variables for restarting
            if (tr_iage) call write_restart_age
!            if (tr_pond) call write_restart_ponds
         endif
         call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine ice_step
    
!=======================================================================
!BOP
!
! !ROUTINE: step_therm1 - step pre-coupler thermodynamics
!
! !DESCRIPTION:
!
! Driver for updating ice and snow internal temperatures and
! computing thermodynamic growth rates and coupler fluxes.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!
! !INTERFACE:

      subroutine step_therm1 (dt)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         iblk        , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n           , & ! thickness category index
         il1, il2    , & ! ice layer indices for eice
         sl1, sl2        ! snow layer indices for esno

      integer (kind=int_kind), save :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), save :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

! 2D coupler variables (computed for each category, then aggregated)

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsensn      , & ! surface downward sensible heat     (W/m^2)
         fswabsn     , & ! shortwave absorbed by ice          (W/m^2)
         flwoutn     , & ! upward LW at surface               (W/m^2)
         evapn       , & ! flux of vapor, atmos to ice   (kg m-2 s-1)
         freshn      , & ! flux of water, ice to ocean     (kg/m^2/s)
         fsaltn      , & ! flux of salt, ice to ocean      (kg/m^2/s)
         fhocnn      , & ! fbot corrected for leftover energy (W/m^2)
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         strairxn    , & ! air/ice zonal  stress,             (N/m^2)
         strairyn    , & ! air/ice meridional stress,         (N/m^2)
         Trefn       , & ! air tmp reference level                (K)
         Qrefn           ! air sp hum reference level         (kg/kg)

      ! other local variables
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         Tbot        , & ! ice bottom surface temperature (deg C)
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         fswintn         ! SW absorbed in ice interior, below surface (W m-2)

      ! Local variables to keep track of melt for ponds
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         melts_old, &
         meltt_old, &
         melts_tmp, &
         meltt_tmp

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr) :: &
         Iswabsn         ! SW radiation absorbed in ice layers (W m-2)

      ! snow variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsn             ! snow horizontal fraction
      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) :: &
         rhosnwn     , & ! snow density (kg/m3)
         rsnwn       , & ! snow grain radius (micro-meters)
         Sswabsn         ! SW radiation absorbed in snow layers (W m-2)

      ! pond variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fpn         , & ! pond fraction
         hpn             ! pond depth (m)

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      call ice_timer_start(timer_column)  ! column physics
      call ice_timer_start(timer_thermo)  ! thermodynamics

      call init_history_therm    ! initialize thermo history variables

      l_stop = .false.

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Save the ice area passed to the coupler (so that history fields
      !  can be made consistent with coupler fields).
      ! Save the initial ice area and volume in each category.
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            aice_init (i,j,  iblk) = aice (i,j,  iblk)
         enddo
         enddo

         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen_init(i,j,n,iblk) = aicen(i,j,n,iblk)
            vicen_init(i,j,n,iblk) = vicen(i,j,n,iblk)
         enddo
         enddo
         enddo

      !-----------------------------------------------------------------
      ! Adjust frzmlt to account for ice-ocean heat fluxes since last
      !  call to coupler.
      ! Compute lateral and bottom heat fluxes.
      !-----------------------------------------------------------------

         call frzmlt_bottom_lateral                                      &
                                (nx_block,           ny_block,           &
                                 ilo, ihi,           jlo, jhi,           &
                                 dt,                                     &
                                 aice  (:,:,  iblk), frzmlt(:,:,  iblk), &
                                 eicen (:,:,:,iblk), esnon (:,:,:,iblk), &
                                 sst   (:,:,  iblk), Tf    (:,:,  iblk), &
                                 strocnxT(:,:,iblk), strocnyT(:,:,iblk), &
                                 Tbot,               fbot,               &
                                 rside (:,:,  iblk) )


      !-----------------------------------------------------------------
      ! Compute cosine of solar zenith angle.
      ! This is used by the delta-Eddington shortwave module.
      ! Albedos are aggregated in merge_fluxes only for cells w/ coszen > 0.
      ! For basic shortwave, simply set coszen to a constant between 0 and 1.
      !-----------------------------------------------------------------

         call ice_timer_start(timer_sw)

         if (trim(shortwave) == 'dEdd') then ! delta Eddington

            ! identify ice-ocean cells
            icells = 0
            do j = 1, ny_block
            do i = 1, nx_block
               if (tmask(i,j,iblk)) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif
            enddo               ! i
            enddo               ! j

            call compute_coszen (nx_block,         ny_block,       &
                                 icells,                           &
                                 indxi,            indxj,          &
                                 tlat  (:,:,iblk), tlon(:,:,iblk), &
                                 coszen(:,:,iblk), dt)

         else                     ! basic (ccsm3) shortwave
            coszen(:,:,iblk) = p5 ! sun above the horizon
         endif

         call ice_timer_stop(timer_sw)

         do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify cells with nonzero ice area
      !-----------------------------------------------------------------
           
            icells = 0
            do j = jlo, jhi
            do i = ilo, ihi
               if (aicen(i,j,n,iblk) > puny) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif
            enddo               ! i
            enddo               ! j

      !-----------------------------------------------------------------
      ! Solar radiation: albedo and absorbed shortwave
      !-----------------------------------------------------------------

            call ice_timer_start(timer_sw)

            if (calc_Tsfc) then

               if (trim(shortwave) == 'dEdd') then   ! delta Eddington

      ! note that rhoswn, rsnw, fp, hp and Sswabs ARE NOT dimensioned with ncat
      ! BPB 19 Dec 2006

                  ! set snow properties
                  call shortwave_dEdd_set_snow(nx_block, ny_block,        &
                                 icells,                                  &
                                 indxi,               indxj,              &
                                 aicen(:,:,n,iblk),   vsnon(:,:,n,iblk),  &
                                 trcrn(:,:,nt_Tsfc,n,iblk), fsn,          &
                                 rhosnwn,             rsnwn)

                  if (.not. tr_pond) then

                  ! set pond properties
                  call shortwave_dEdd_set_pond(nx_block, ny_block,         &
                                 icells,                                   &
                                 indxi,               indxj,               &
                                 aicen(:,:,n,iblk),                        &
                                 trcrn(:,:,nt_Tsfc,n,iblk),                &
                                 fsn,                 fpn,                 &
                                 hpn)

                  else
                  
                     fpn(:,:) = apondn(:,:,n,iblk)
                     hpn(:,:) = hpondn(:,:,n,iblk)

                  endif

                  call shortwave_dEdd(nx_block,     ny_block,            &
                                 icells,                                 &
                                 indxi,             indxj,               &
                                 coszen(:,:, iblk),                      &
                                 aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                                 vsnon(:,:,n,iblk), fsn,                 &
                                 rhosnwn,           rsnwn,               &
                                 fpn,               hpn,                 &
                                 swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                                 swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                                 alvdrn(:,:,n,iblk),alvdfn(:,:,n,iblk),  &
                                 alidrn(:,:,n,iblk),alidfn(:,:,n,iblk),  &
                                 fswsfcn,           fswintn,             &
                                 fswthrun,          Sswabsn,             &
                                 Iswabsn)

               else
                  Sswabsn(:,:,:) = c0

                 call shortwave_ccsm3(nx_block,     ny_block,            &
                                 icells,                                 &
                                 indxi,             indxj,               &
                                 aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                                 vsnon(:,:,n,iblk),                      &
                                 trcrn(:,:,nt_Tsfc,n,iblk),              &
                                 swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                                 swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                                 alvdrn(:,:,n,iblk),alidrn(:,:,n,iblk),  &
                                 alvdfn(:,:,n,iblk),alidfn(:,:,n,iblk),  &
                                 fswsfcn,           fswintn,             &
                                 fswthrun,          Iswabsn,             &
                                 apondn(:,:,n,iblk),hpondn(:,:,n,iblk))

               endif

            else    ! .not. calc_Tsfc

               ! Initialize for safety
               alvdrn(:,:,n,iblk)  = c0
               alidrn(:,:,n,iblk)  = c0
               alvdfn(:,:,n,iblk)  = c0
               alidfn(:,:,n,iblk)  = c0
               fswsfcn(:,:) = c0
               fswintn(:,:) = c0
               fswthrun(:,:)= c0
               Iswabsn(:,:,:) = c0
               Sswabsn(:,:,:) = c0

            endif    ! calc_Tsfc

            call ice_timer_stop(timer_sw)

            if (calc_Tsfc .or. calc_strair) then 

      !-----------------------------------------------------------------
      ! Atmosphere boundary layer calculation; compute coefficients
      ! for sensible and latent heat fluxes.
      !
      ! NOTE: The wind stress is computed here for later use if 
      !       calc_strair = .true.   Otherwise, the wind stress
      !       components are set to the data values.
      !-----------------------------------------------------------------

               if (trim(atmbndy) == 'constant') then
                   call atmo_boundary_const &
                                   (nx_block,      ny_block,        &
                                    'ice',          icells,         &
                                    indxi,          indxj,          &
                                    uatm(:,:,iblk), vatm(:,:,iblk), &
                                    wind(:,:,iblk), rhoa(:,:,iblk), &
                                    strairxn,       strairyn,       &
                                    lhcoef,         shcoef)
               else ! default
                   call atmo_boundary_layer & 
                                  (nx_block,       ny_block,       &
                                   'ice',          icells,         &
                                   indxi,          indxj,          &
                                   trcrn(:,:,nt_Tsfc,n,iblk),      &
                                   potT(:,:,iblk),                 &
                                   uatm(:,:,iblk), vatm(:,:,iblk), &
                                   wind(:,:,iblk), zlvl(:,:,iblk), &
                                   Qa  (:,:,iblk), rhoa(:,:,iblk), &
                                   strairxn,       strairyn,       &
                                   Trefn,          Qrefn,          &
                                   worka,          workb,          &
                                   lhcoef,         shcoef)
               endif ! atmbndy

            else

               ! Initialize for safety
               Trefn (:,:)  = c0
               Qrefn (:,:)  = c0
               lhcoef(:,:)  = c0
               shcoef(:,:)  = c0

            endif   ! calc_Tsfc or calc_strair

            if (.not.(calc_strair)) then
               ! Set to data values (on T points)
               strairxn(:,:) = strax(:,:,iblk)
               strairyn(:,:) = stray(:,:,iblk)
            endif

      !-----------------------------------------------------------------
      ! If not calculating surface temperature and fluxes, set surface 
      ! fluxes (flatn, fsurfn, and fcondtopn) to be used in
      ! thickness_changes
      !-----------------------------------------------------------------
 
            if (.not.(calc_Tsfc)) then

               call set_sfcflux(nx_block,  ny_block,  &
                                n,         iblk,      &
                                icells,               & 
                                indxi,     indxj,     &
                                aicen    (:,:,n,iblk),&
                                flatn    (:,:,n,iblk),&
                                fsurfn   (:,:,n,iblk),&
                                fcondtopn(:,:,n,iblk) )

            endif

      !-----------------------------------------------------------------
      ! Update ice age
      ! This is further adjusted for freezing in the thermodynamics.
      ! Melting does not alter the ice age.
      !-----------------------------------------------------------------

            if (tr_iage) then
               call increment_age (nx_block, ny_block,      &
                                   dt, icells,              &
                                   indxi, indxj,            &
                                   trcrn(:,:,nt_iage,n,iblk))
            endif

      !-----------------------------------------------------------------
      ! Vertical thermodynamics: Heat conduction, growth and melting.
      !----------------------------------------------------------------- 

            il1 = ilyr1(n)
            il2 = ilyrn(n)
            sl1 = slyr1(n)
            sl2 = slyrn(n)

            melts_old = melts(:,:,iblk)
            meltt_old = meltt(:,:,iblk)

!lipscomb - Temporary subroutine to compute fsurfn and fcondtopn.
!           Used for testing the option calc_Tsfc = F

            if (.not. calc_Tsfc) then

               call explicit_calc_Tsfc (nx_block,          ny_block,     &
                                        my_task,           icells,       &      
                                        indxi,             indxj,        &
                                        trcrn(:,:,nt_Tsfc,n,iblk),       &
                                        aicen(:,:,n,iblk),               &
                                        vicen(:,:,n,iblk),               &
                                        vsnon(:,:,n,iblk),               &
                                        eicen  (:,:,il1:il2,iblk),       &
                                        esnon  (:,:,sl1:sl2,iblk),       &
                                        rhoa(:,:,iblk),    flw(:,:,iblk),&
                                        potT(:,:,iblk),    Qa (:,:,iblk),&
                                        shcoef,            lhcoef,       &
                                        fswsfcn,           flwoutn,      &
                                        fsensn,                          & 
                                        flatn(:,:,n,iblk),               &
                                        fsurfn(:,:,n,iblk),              &
                                        fcondtopn(:,:,n,iblk))

            endif

            call thermo_vertical                                       &
                            (nx_block,            ny_block,            &
                             dt,                  icells,              &
                             indxi,               indxj,               &
                             aicen(:,:,n,iblk),                        &
                             trcrn(:,:,:,n,iblk),                      &
                             vicen(:,:,n,iblk),   vsnon(:,:,n,iblk),   &
                             eicen  (:,:,il1:il2,iblk),                &
                             esnon  (:,:,sl1:sl2,iblk),                &
                             flw    (:,:,iblk),   potT (:,:,iblk),     &
                             Qa     (:,:,iblk),   rhoa (:,:,iblk),     &
                             fsnow  (:,:,iblk),                        &
                             fbot,                Tbot,                &
                             lhcoef,              shcoef,              &
                             fswsfcn,             fswintn,             &
                             fswthrun,                                 &
                             Sswabsn,             Iswabsn,             &
                             fsurfn(:,:,n,iblk),  fcondtopn(:,:,n,iblk),&
                             fsensn,              flatn(:,:,n,iblk),   &
                             fswabsn,             flwoutn,             &
                             evapn,               freshn,              &
                             fsaltn,              fhocnn,              &
                             meltt   (:,:,iblk),  melts   (:,:,iblk),  &
                             meltb   (:,:,iblk),                       &
                             congel  (:,:,iblk),  snoice  (:,:,iblk),  &
                             mlt_onset(:,:,iblk), frz_onset(:,:,iblk), &
                             yday,                l_stop,              &
                             istop,               jstop)

         if (l_stop) then
            write (nu_diag,*) 'istep1, my_task, iblk =', &
                               istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0) &
                 write(nu_diag,*) 'Global i and j:', &
                                  this_block%i_glob(istop), &
                                  this_block%j_glob(jstop) 
                 write(nu_diag,*) 'Lat, Lon:', &
                                 TLAT(istop,jstop,iblk)*rad_to_deg, &
                                 TLON(istop,jstop,iblk)*rad_to_deg
            call abort_ice ('ice: Vertical thermo error')
         endif

      !-----------------------------------------------------------------
      ! Melt ponds
      !-----------------------------------------------------------------

         if (tr_pond) then

            melts_tmp = melts(:,:,iblk) - melts_old
            meltt_tmp = meltt(:,:,iblk) - meltt_old

            call compute_ponds(nx_block, ny_block,                      &
                               ilo, ihi, jlo, jhi,                      &
                               meltt_tmp, melts_tmp, frain(:,:,iblk),   &
                               aicen (:,:,n,iblk), vicen (:,:,n,iblk),  &
                               vsnon (:,:,n,iblk), trcrn (:,:,:,n,iblk),&
                               apondn(:,:,n,iblk), hpondn(:,:,n,iblk))

         endif

      !-----------------------------------------------------------------
      ! Increment area-weighted fluxes.
      !-----------------------------------------------------------------

         call merge_fluxes (nx_block,           ny_block,             &
                            icells,                                   &
                            indxi,              indxj,                &
                            aicen_init(:,:,n,iblk),                   &
                            flw(:,:,iblk),      coszen(:,:,iblk),     &
                            alvdrn(:,:,n,iblk), alidrn(:,:,n,iblk),   &
                            alvdfn(:,:,n,iblk), alidfn(:,:,n,iblk),   &
                            strairxn,           strairyn,             &
                            fsurfn(:,:,n,iblk), fcondtopn(:,:,n,iblk),&
                            fsensn,             flatn(:,:,n,iblk),    &
                            fswabsn,            flwoutn,              &
                            evapn,                                    &
                            Trefn,              Qrefn,                &
                            freshn,             fsaltn,               &
                            fhocnn,             fswthrun,             &
                            alvdr   (:,:,iblk), alidr     (:,:,iblk), &
                            alvdf   (:,:,iblk), alidf     (:,:,iblk), &
                            strairxT(:,:,iblk), strairyT  (:,:,iblk), &
                            fsurf   (:,:,iblk), fcondtop  (:,:,iblk), &
                            fsens   (:,:,iblk), flat      (:,:,iblk), &
                            fswabs  (:,:,iblk), flwout    (:,:,iblk), &
                            evap    (:,:,iblk),                       &
                            Tref    (:,:,iblk), Qref      (:,:,iblk), &
                            fresh   (:,:,iblk), fsalt     (:,:,iblk), &
                            fhocn   (:,:,iblk), fswthru   (:,:,iblk))

         enddo                  ! ncat

      enddo                      ! iblk

      call ice_timer_stop(timer_thermo) ! thermodynamics
      call ice_timer_stop(timer_column) ! column physics

      end subroutine step_therm1

!=======================================================================
!BOP
!
! !ROUTINE: step_therm2 - step post-coupler thermodynamics
!
! !DESCRIPTION:
!
!-----------------------------------------------------------------------
! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! NOTE: Ocean fluxes are initialized here.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:

      subroutine step_therm2 (dt)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
!lipscomb - delete hicen later?
!      real (kind=dbl_kind), &
!         dimension (nx_block,ny_block,ncat,max_blocks) :: &

      integer (kind=int_kind) :: &
         iblk        , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, n

      integer (kind=int_kind), save :: &
         icells          ! number of ice/ocean cells 

      integer (kind=int_kind), dimension(nx_block*ny_block), save :: &
         indxi, indxj    ! indirect indices for ice/ocean cells

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts

      call ice_timer_start(timer_column)  ! column physics
      call ice_timer_start(timer_thermo)  ! thermodynamics

      l_stop = .false.
      
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Let rain drain through to the ocean.
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            fresh     (i,j,iblk) = fresh(i,j,iblk)       &
                                 + frain(i,j,iblk)*aice(i,j,iblk)
         enddo
         enddo

      !-----------------------------------------------------------------
      ! Given thermodynamic growth rates, transport ice between
      ! thickness categories.
      !-----------------------------------------------------------------

         call ice_timer_start(timer_catconv)    ! category conversions

         if (kitd == 1) then
      !-----------------------------------------------------------------
      ! Compute fractional ice area in each grid cell.
      !-----------------------------------------------------------------
            call aggregate_area (nx_block,          ny_block, &
                                 aicen(:,:,:,iblk),           &
                                 aice (:,:,  iblk), aice0(:,:,iblk))

      !-----------------------------------------------------------------
      ! Identify grid cells with ice.
      !-----------------------------------------------------------------

            icells = 0
            do j = jlo,jhi
            do i = ilo,ihi
               if (aice(i,j,iblk) > puny) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif
            enddo
            enddo

            if (icells > 0) then

            call linear_itd (nx_block, ny_block,       &
                             icells, indxi, indxj,     &
                             trcr_depend,    &
                             aicen_init(:,:,:,iblk),   &
                             vicen_init(:,:,:,iblk),   &
                             aicen     (:,:,:,iblk),   &
                             trcrn     (:,:,:,:,iblk), & 
                             vicen     (:,:,:,iblk),   &
                             vsnon     (:,:,:,iblk),   &
                             eicen     (:,:,:,iblk),   &
                             esnon     (:,:,:,iblk),   &
                             aice      (:,:,  iblk),   &
                             aice0     (:,:,  iblk),   &
                             l_stop,                   &
                             istop,    jstop)

            if (l_stop) then
               write (nu_diag,*) 'istep1, my_task, iblk =', &
                                  istep1, my_task, iblk
               write (nu_diag,*) 'Global block:', this_block%block_id
               if (istop > 0 .and. jstop > 0) &
                    write(nu_diag,*) 'Global i and j:', &
                                     this_block%i_glob(istop), &
                                     this_block%j_glob(jstop) 
               call abort_ice ('ice: Linear ITD error')
            endif

            endif

         endif  ! kitd = 1

         call ice_timer_stop(timer_catconv)    ! category conversions

      !-----------------------------------------------------------------
      ! Add frazil ice growing in leads.
      !-----------------------------------------------------------------

         ! identify ice-ocean cells
         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

         call add_new_ice (nx_block,              ny_block, &
                           icells,                          &
                           indxi,                 indxj,    &
                           tmask     (:,:,  iblk), dt,      &
                           aicen     (:,:,:,iblk),          &
                           trcrn     (:,:,:,:,iblk),        &
                           vicen     (:,:,:,iblk),          &
                           eicen     (:,:,:,iblk),          &
                           aice0     (:,:,  iblk),          &
                           aice      (:,:,  iblk),          &
                           frzmlt    (:,:,  iblk),          &
                           frazil    (:,:,  iblk),          &
                           frz_onset (:,:,  iblk), yday,    &
                           fresh     (:,:,  iblk),          &
                           fsalt     (:,:,  iblk),          &
                           Tf        (:,:,  iblk), l_stop,  &
                           istop                 , jstop)

         if (l_stop) then
            write (nu_diag,*) 'istep1, my_task, iblk =', &
                               istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0) &
                 write(nu_diag,*) 'Global i and j:', &
                                  this_block%i_glob(istop), &
                                  this_block%j_glob(jstop) 
            call abort_ice ('ice: add_new_ice error')
         endif

      !-----------------------------------------------------------------
      ! Melt ice laterally.
      !-----------------------------------------------------------------
         call lateral_melt (nx_block, ny_block,     &
                            ilo, ihi, jlo, jhi,     &
                            dt,                     &
                            fresh     (:,:,  iblk), &
                            fsalt     (:,:,  iblk), &    
                            fhocn     (:,:,  iblk), &
                            rside     (:,:,  iblk), &
                            meltl     (:,:,  iblk), &
                            aicen     (:,:,:,iblk), &
                            vicen     (:,:,:,iblk), &
                            vsnon     (:,:,:,iblk), &
                            eicen     (:,:,:,iblk), &
                            esnon     (:,:,:,iblk) )

      !-----------------------------------------------------------------
      ! For the special case of a single category, adjust the area and
      ! volume (assuming that half the volume change decreases the
      ! thickness, and the other half decreases the area).  
      !-----------------------------------------------------------------

!NOTE - this does not work - hicen_init is not defined - ECH

!         if (ncat==1) &
!              call reduce_area (nx_block, ny_block,     &
!                                ilo, ihi, jlo, jhi,     &
!                                tmask     (:,:,  iblk), &
!                                aicen     (:,:,:,iblk), &
!                                vicen     (:,:,:,iblk), &
!                                hicen_init(:,:,1,iblk), &
!                                hicen     (:,:,1,iblk)) 

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

         call cleanup_itd (nx_block,             ny_block,             &
                           ilo, ihi,             jlo, jhi,             &
                           dt,                                         &
                           aicen   (:,:,:,iblk), trcrn (:,:,:,:,iblk), &
                           vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
                           eicen   (:,:,:,iblk), esnon (:,:,  :,iblk), &
                           aice0   (:,:,  iblk), aice      (:,:,iblk), &
                           trcr_depend,                                &
                           fresh   (:,:,  iblk), fsalt     (:,:,iblk), &
                           fhocn   (:,:,  iblk),                       &
                           heat_capacity,        l_stop,               &
                           istop,                jstop)

         if (l_stop) then
            write (nu_diag,*) 'istep1, my_task, iblk =', &
                               istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0) &
                 write(nu_diag,*) 'Global i and j:', &
                                  this_block%i_glob(istop), &
                                  this_block%j_glob(jstop) 
            call abort_ice ('ice: ITD cleanup error')
         endif

      enddo                     ! iblk

      !-------------------------------------------------------------------
      ! Ghost cell updates for state variables.
      !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)
      call bound_state (aicen, trcrn, &
                        vicen, vsnon, &
                        eicen, esnon)
      call ice_timer_stop(timer_bound)

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Aggregate the updated state variables (includes ghost cells). 
      !----------------------------------------------------------------- 
 
         call aggregate (nx_block,          ny_block,             &
                         aicen(:,:,:,iblk), trcrn(:,:,:,:,iblk),  &
                         vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
                         eicen(:,:,:,iblk), esnon(:,:,  :,iblk),  &
                         aice (:,:,  iblk), trcr (:,:,:,  iblk),  &
                         vice (:,:,  iblk), vsno (:,:,    iblk),  &
                         eice (:,:,  iblk), esno (:,:,    iblk),  &
                         aice0(:,:,  iblk), tmask(:,:,    iblk),  &
                         trcr_depend) 

      !-----------------------------------------------------------------
      ! Compute thermodynamic area and volume tendencies.
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            daidtt(i,j,iblk) = (aice(i,j,iblk) - daidtt(i,j,iblk)) / dt
            dvidtt(i,j,iblk) = (vice(i,j,iblk) - dvidtt(i,j,iblk)) / dt
         enddo
         enddo

      enddo                     ! iblk

      call ice_timer_stop(timer_thermo)  ! thermodynamics
      call ice_timer_stop(timer_column)  ! column physics

      end subroutine step_therm2

!=======================================================================
!BOP
!
! !ROUTINE: step_dynamics - step ice dynamics, transport, and ridging
!
! !DESCRIPTION:
!
! Run one time step of dynamics, horizontal transport, and ridging.
! NOTE: The evp and transport modules include boundary updates, so
!       they cannot be done inside a single block loop.  Ridging
!       and cleanup, on the other hand, are single-column operations. 
!       They are called with argument lists inside block loops
!       to increase modularity.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!
! !INTERFACE:

      subroutine step_dynamics (dt)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: & 
         iblk        , & ! block index 
         i,j         , & ! horizontal indices
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      integer (kind=int_kind), save :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), save :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts

      call init_history_dyn     ! initialize dynamic history variables

      !-----------------------------------------------------------------
      ! Elastic-viscous-plastic ice dynamics
      !-----------------------------------------------------------------

      if (kdyn == 1) call evp (dt)

      !-----------------------------------------------------------------
      ! Horizontal ice transport
      !-----------------------------------------------------------------

      if (advection == 'upwind') then
         call transport_upwind (dt)    ! upwind
      else
         call transport_remap (dt)     ! incremental remapping
      endif

      !-----------------------------------------------------------------
      ! Ridging
      !-----------------------------------------------------------------

      call ice_timer_start(timer_column)
      call ice_timer_start(timer_ridge)

      l_stop = .false.

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk), iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Identify ice-ocean cells.
      ! Note:  We can not define icells here using aice>puny because
      !        aice has not yet been updated since the transport (and
      !        it may be out of whack, which the ridging helps fix).-ECH
      !-----------------------------------------------------------------
           
         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

         if (icells > 0) then

         call ridge_ice (nx_block,             ny_block,                 &
                         dt,                   icells,                   &
                         indxi,                indxj,                    &
!!                         Delt    (:,:,  iblk), divu      (:,:,  iblk), &
                         rdg_conv(:,:,  iblk), rdg_shear (:,:,  iblk),   &
                         aicen   (:,:,:,iblk), trcrn     (:,:,:,:,iblk), &
                         vicen   (:,:,:,iblk), vsnon     (:,:,:,iblk),   &
                         eicen   (:,:,:,iblk), esnon     (:,:,:,iblk),   &
                         aice0   (:,:,  iblk),                           &
                         trcr_depend,          l_stop,                   &
                         istop,                jstop,                    &   
                         dardg1dt(:,:,iblk),   dardg2dt  (:,:,iblk),     &
                         dvirdgdt(:,:,iblk),   opening   (:,:,iblk),     &
                         fresh   (:,:,iblk),   fhocn     (:,:,iblk))

         if (l_stop) then
            write (nu_diag,*) 'istep1, my_task, iblk =', &
                               istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0) &
                 write(nu_diag,*) 'Global i and j:', &
                                  this_block%i_glob(istop), &
                                  this_block%j_glob(jstop) 
            call abort_ice ('ice: Ridging error')
         endif

         endif

      enddo                     ! iblk

      call ice_timer_stop(timer_ridge)

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk), iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

         call cleanup_itd (nx_block,             ny_block,             &
                           ilo, ihi,             jlo, jhi,             &
                           dt,                                         &
                           aicen   (:,:,:,iblk), trcrn (:,:,:,:,iblk), &
                           vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
                           eicen   (:,:,:,iblk), esnon (:,:,  :,iblk), &
                           aice0   (:,:,  iblk), aice      (:,:,iblk), &
                           trcr_depend,                                &
                           fresh   (:,:,  iblk), fsalt     (:,:,iblk), &
                           fhocn   (:,:,  iblk),                       &
                           heat_capacity,        l_stop,               &
                           istop,                jstop)

         if (l_stop) then
            write (nu_diag,*) 'istep1, my_task, iblk =', &
                               istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0) &
                 write(nu_diag,*) 'Global i and j:', &
                                  this_block%i_glob(istop), &
                                  this_block%j_glob(jstop) 
            call abort_ice ('ice: ITD cleanup error')
         endif

      enddo                     ! iblk

      !-------------------------------------------------------------------
      ! Ghost cell updates for state variables.
      !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)
      call bound_state (aicen, trcrn, &
                        vicen, vsnon, &
                        eicen, esnon)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! Aggregate the updated state variables (includes ghost cells). 
      !----------------------------------------------------------------- 
 
      do iblk = 1, nblocks
         call aggregate (nx_block,          ny_block,             &
                         aicen(:,:,:,iblk), trcrn(:,:,:,:,iblk),  &
                         vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
                         eicen(:,:,:,iblk), esnon(:,:,  :,iblk),  &
                         aice (:,:,  iblk), trcr (:,:,:,  iblk),  &
                         vice (:,:,  iblk), vsno (:,:,    iblk),  &
                         eice (:,:,  iblk), esno (:,:,    iblk),  &
                         aice0(:,:,  iblk), tmask(:,:,    iblk),  &
                         trcr_depend) 
      enddo

      call ice_timer_stop(timer_column)

      end subroutine step_dynamics

!=======================================================================
!BOP
!
! !ROUTINE: coupling_prep
!
! !DESCRIPTION:
!
! Prepare for coupling
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!
! !INTERFACE:

      subroutine coupling_prep
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: & 
         iblk        , & ! block index 
         i,j         , & ! horizontal indices
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      call ice_timer_start(timer_column)

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk), iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

         if (oceanmixed_ice) &
            call ocean_mixed_layer (dt)   ! ocean surface fluxes and sst

      !----------------------------------------------------------------
      ! Store grid box mean fluxes before scaling by aice_init
      !----------------------------------------------------------------

         fresh_gbm  (:,:,iblk) = fresh  (:,:,iblk)
         fsalt_gbm  (:,:,iblk) = fsalt  (:,:,iblk)
         fhocn_gbm  (:,:,iblk) = fhocn  (:,:,iblk)
         fswthru_gbm(:,:,iblk) = fswthru(:,:,iblk)

      !-----------------------------------------------------------------
      ! Divide fluxes by ice area 
      !  - the CCSM coupler assumes fluxes are per unit ice area
      !  - also needed for global budget in diagnostics
      !-----------------------------------------------------------------

         call scale_fluxes (nx_block,            ny_block,           &
                            tmask    (:,:,iblk),                     &
                            aice_init(:,:,iblk), Tf      (:,:,iblk), &
                            Tair     (:,:,iblk), Qa      (:,:,iblk), &
                            strairxT (:,:,iblk), strairyT(:,:,iblk), &
                            fsens    (:,:,iblk), flat    (:,:,iblk), &
                            fswabs   (:,:,iblk), flwout  (:,:,iblk), &
                            evap     (:,:,iblk),                     &
                            Tref     (:,:,iblk), Qref    (:,:,iblk), &
                            fresh    (:,:,iblk), fsalt   (:,:,iblk), &
                            fhocn    (:,:,iblk), fswthru (:,:,iblk), &
                            alvdr    (:,:,iblk), alidr   (:,:,iblk), &
                            alvdf    (:,:,iblk), alidf   (:,:,iblk))

         if (.not. calc_Tsfc) then

       !---------------------------------------------------------------
       ! If surface fluxes were provided, conserve these fluxes at ice 
       ! free points by passing to ocean. 
       !---------------------------------------------------------------

            call sfcflux_to_ocn & 
                         (nx_block,              ny_block,             &
                          tmask   (:,:,iblk),    aice_init(:,:,iblk),  &
                          fsurfn_f (:,:,:,iblk), flatn_f(:,:,:,iblk),  &
                          fresh    (:,:,iblk),   fhocn    (:,:,iblk))
         endif                 

      enddo

      call ice_timer_stop(timer_column)

      end subroutine coupling_prep

!=======================================================================
!BOP
!
! !ROUTINE: set_sfcflux - set surface fluxes from forcing fields
!
! !DESCRIPTION:
!
! If model is not calculating surface temperature, set the surface
! flux values using values read in from forcing data or supplied via
! coupling (stored in ice_flux).
!
! If CICE is running in NEMO environment, convert fluxes from GBM values 
! to per unit ice area values. If model is not running in NEMO environment, 
! the forcing is supplied as per unit ice area values.
!
! !REVISION HISTORY:
!
! authors Alison McLaren, Met Office
!
! !INTERFACE:
!
      subroutine set_sfcflux (nx_block,  ny_block, &
                              n,         iblk,     &
                              icells,              & 
                              indxi,     indxj,    &
                              aicen,               &
                              flatn,               &
                              fsurfn,              &
                              fcondtopn)
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         n,                  & ! thickness category index
         iblk,               & ! block index
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      ! ice state variables
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen           ! concentration of ice

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         flatn       , & ! latent heat flux   (W/m^2) 
         fsurfn      , & ! net flux to top surface, not including fcondtopn
         fcondtopn       ! downward cond flux at top surface (W m-2)

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij              ! horizontal indices, combine i and j loops

#ifdef CICE_IN_NEMO
      real (kind=dbl_kind)  :: &
         raicen          ! 1/aicen
#endif

      logical (kind=log_kind) :: &
         extreme_flag    ! flag for extreme forcing values

      logical (kind=log_kind), parameter :: & 
         extreme_test=.true. ! test and write out extreme forcing data
!
!EOP
!
#ifdef CICE_IN_NEMO

!----------------------------------------------------------------------
! Convert fluxes from GBM values to per ice area values when 
! running in NEMO environment.  (When in standalone mode, fluxes
! are input as per ice area.)
!----------------------------------------------------------------------

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            raicen        = c1 / aicen(i,j)
            fsurfn(i,j)   = fsurfn_f(i,j,n,iblk)*raicen
            fcondtopn(i,j)= fcondtopn_f(i,j,n,iblk)*raicen
            flatn(i,j)    = flatn_f(i,j,n,iblk)*raicen

         enddo
#else
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            fsurfn(i,j)   = fsurfn_f(i,j,n,iblk)
            fcondtopn(i,j)= fcondtopn_f(i,j,n,iblk)
            flatn(i,j)    = flatn_f(i,j,n,iblk)

         enddo
#endif

!----------------------------------------------------------------
! Flag up any extreme fluxes
!---------------------------------------------------------------

         if (extreme_test) then
            extreme_flag = .false.

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)         

               if (fcondtopn(i,j) < -100.0_dbl_kind & 
                     .or. fcondtopn(i,j) > 20.0_dbl_kind) then
                  extreme_flag = .true.
               endif

               if (fsurfn(i,j) < -100.0_dbl_kind & 
                    .or. fsurfn(i,j) > 80.0_dbl_kind) then
                  extreme_flag = .true.
               endif

               if (flatn(i,j) < -20.0_dbl_kind & 
                     .or. flatn(i,j) > 20.0_dbl_kind) then
                  extreme_flag = .true.
               endif

            enddo  ! ij

            if (extreme_flag) then
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)         

                  if (fcondtopn(i,j) < -100.0_dbl_kind & 
                       .or. fcondtopn(i,j) > 20.0_dbl_kind) then
                     write(nu_diag,*) & 
                       'Extreme forcing: -100 > fcondtopn > 20'
                     write(nu_diag,*) & 
                       'i,j,n,iblk,aicen,fcondtopn = ', & 
                       i,j,n,iblk,aicen(i,j),fcondtopn(i,j)
                  endif

                  if (fsurfn(i,j) < -100.0_dbl_kind & 
                       .or. fsurfn(i,j) > 80.0_dbl_kind) then
                     write(nu_diag,*) & 
                       'Extreme forcing: -100 > fsurfn > 40'
                     write(nu_diag,*) & 
                       'i,j,n,iblk,aicen,fsurfn = ', & 
                        i,j,n,iblk,aicen(i,j),fsurfn(i,j)
                  endif

                  if (flatn(i,j) < -20.0_dbl_kind & 
                       .or. flatn(i,j) > 20.0_dbl_kind) then
                     write(nu_diag,*) & 
                       'Extreme forcing: -20 > flatn > 20'
                     write(nu_diag,*) & 
                       'i,j,n,iblk,aicen,flatn = ', & 
                        i,j,n,iblk,aicen(i,j),flatn(i,j)
                  endif

               enddo  ! ij
      
            endif  ! extreme_flag
         endif     ! extreme_test    

      end subroutine set_sfcflux 

!=======================================================================
!BOP
!
! !IROUTINE: sfcflux_to_ocn
!
! !DESCRIPTION:
!
! If surface heat fluxes are provided to CICE instead of CICE calculating
! them internally (i.e. .not. calc_Tsfc), then these heat fluxes can 
! be provided at points which do not have ice.  (This is could be due to
! the heat fluxes being calculated on a lower resolution grid or the
! heat fluxes not recalculated at every CICE timestep.)  At ice free points, 
! conserve energy and water by passing these fluxes to the ocean.
!
! !INTERFACE:
!
       subroutine sfcflux_to_ocn(nx_block,   ny_block,     &
                                 tmask,      aice,         &
                                 fsurfn_f,   flatn_f,      &
                                 fresh,      fhocn)
!
! !REVISION HISTORY:
!
! authors: A. McLaren, Met Office
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block  ! block dimensions

      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          tmask       ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(in):: &
          aice        ! initial ice concentration

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
          intent(in) :: &
          fsurfn_f, & ! net surface heat flux (provided as forcing)
          flatn_f     ! latent heat flux (provided as forcing)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &
          fresh        , & ! fresh water flux to ocean         (kg/m2/s)
          fhocn            ! actual ocn/ice heat flx           (W/m**2)
!
!EOP
!
#ifdef CICE_IN_NEMO
      integer (kind=int_kind) :: &
          i, j, n    ! horizontal indices
      
      real (kind=dbl_kind)    :: &
          rLsub            ! 1/Lsub

      rLsub = c1 / Lsub

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j) .and. aice(i,j) <= puny) then
               fhocn(i,j)      = fhocn(i,j)              &
                            + fsurfn_f(i,j,n) + flatn_f(i,j,n)
               fresh(i,j)      = fresh(i,j)              &
                                 + flatn_f(i,j,n) * rLsub
            endif
         enddo   ! i
         enddo   ! j
      enddo      ! n

#endif 
      end subroutine sfcflux_to_ocn

!=======================================================================
!BOP
!
! !ROUTINE: explicit_calc_Tsfc - temporary subroutine to compute sfc fluxes
!
! !DESCRIPTION:
!
! Compute fsurfn and fcondtopn, given temperature, thickness, and 
!  conductivity of top ice or snow layer.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine explicit_calc_Tsfc (nx_block,      ny_block, &
                                     my_task,       icells,   &
                                     indxi,         indxj,    &
                                     Tsfcn,         aicen,    &
                                     vicen,         vsnon,    &
                                     eicen,         esnon,    &
                                     rhoa,          flw,      &
                                     potT,          Qa,       &
                                     shcoef,        lhcoef,   &
                                     fswsfcn,       flwoutn,  &
                                     fsensn,        flatn,    &
                                     fsurfn,        fcondtopn)
!
! !USES:
!
      use ice_therm_vertical, only: hs_min, betak, kimin
! 
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         my_task           , & ! process ID
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension(icells), intent(in) :: &
         indxi, indxj    ! compressed indices for ice cells

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         Tsfcn     ! temperature of ice/snow top surface  (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr), &
         intent(in) :: &
         eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nslyr), &
         intent(in) :: &
         esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat
         
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn       ! conductive flux to top surface
!
!EOP
!
      integer (kind=int_kind) :: &
         isolve    ! number of cells with temps not converged (same as icells)

      integer (kind=int_kind), dimension(icells) :: &
         indxii, indxjj,&  ! compressed indices for cells not converged
         indxij            ! compressed 1D index for cells not converged

      real (kind=dbl_kind), dimension (icells) :: &
         Tsf         , & ! surface temperature
         khis        , & ! 2*k/h for top ice or snow layer
         Tis             ! temperature of top ice or snow layer

      real (kind=dbl_kind), dimension (icells) :: &
         dfsens_dT   , & ! deriv of fsens wrt Tsf (W m-2 deg-1)
         dflat_dT    , & ! deriv of flat wrt Tsf (W m-2 deg-1)
         dflwout_dT  , & ! deriv of flwout wrt Tsf (W m-2 deg-1)
         dfsurf_dT       ! derivative of fsurf wrt Tsf

      integer :: i, j, ij, k

      real (kind=dbl_kind) ::  &
         dTsf         , & ! change in Tsf
         aa1, bb1, cc1, & !
         hslyr, hilyr , & ! snow and ice layer thickness
         qsn, qin     , & ! snow and ice layer enthalpy
         kilyr        , & ! ice layer conductivity
         khmax        , & ! max allowed value of kh
         ci

      logical (kind=log_kind) ::   &
         l_snow           ! true if hsno > hs_min

!lipscomb - for testing - remove later
      if (my_task == mtest) then
         i = itest
         j = jtest
         print*, ''
         print*, 'Beginning Tsf calc, my_task, i, j =', mtest, i ,j
      endif

      ! Initialize fluxes

      fsurfn   (:,:) = c0
      fcondtopn(:,:) = c0
      flwoutn  (:,:) = c0
      fsensn   (:,:) = c0
      flatn    (:,:) = c0

      ! initialize surface temperature
      do ij = 1, icells	
         i = indxi(ij)
         j = indxj(ij)
         Tsf(ij) = Tsfcn(i,j)
      enddo

      !-----------------------------------------------------------------
      ! Initialize isolve and related indices to be identical to icells
      ! and related indices
      !-----------------------------------------------------------------

      isolve = icells
      do ij = 1, icells	
         indxii(ij) = indxi(ij)
         indxjj(ij) = indxj(ij)
         indxij(ij) = ij
      enddo

      ! Compute temperature of top layer and conductivity at upper interface.

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         k = 1   ! top layer of ice or snow

      ! Check if snow layer thickness hsno > hs_min

         hslyr = vsnon(i,j) / aicen(i,j)
         if (hslyr*nslyr > hs_min) then
            l_snow = .true.
         else
            l_snow = .false.
         endif

      ! Temperature and heat capacity of top layer

         if (l_snow) then

            ! Compute enthalpy of top ice layer
            qsn = esnon(i,j,k)*real(nslyr,kind=dbl_kind) / vsnon(i,j) 

            ! Compute snow temperature from enthalpy
            Tis(ij) = (Lfresh + qsn/rhos)/cp_ice
            Tis(ij) = min(Tis(ij), c0)

         else

            ! Compute enthalpy of top ice layer
            qin = eicen(i,j,k)*real(nilyr,kind=dbl_kind) / vicen(i,j)  

            ! Compute ice temperature from enthalpy using quadratic formula

            if (l_brine) then
               aa1 = cp_ice
               bb1 = (cp_ocn-cp_ice)*Tmlt(k) - qin/rhoi - Lfresh 
               cc1 = Lfresh * Tmlt(k)
               Tis(ij) =  (-bb1 - sqrt(bb1*bb1 - c4*aa1*cc1)) /  &
                          (c2*aa1)
               Tis(ij) = min(Tis(ij),Tmlt(k))
            else                ! fresh ice
               Tis(ij) = (Lfresh + qin/rhoi) / cp_ice
               Tis(ij) = min(Tis(ij),c0)
            endif

            ! Compute heat capacity of the ice layer
            if (l_brine) then
               ci = cp_ice - Lfresh*Tmlt(k) /  (Tis(ij)*Tis(ij))
            else
               ci = cp_ice
            endif
         endif

      ! Conductivity (k/h) at upper interface
      ! Limit to satisfy diffusive CFL condition

         if (l_snow) then
            khis(ij) = c2*ksno / hslyr
            khmax = rhos*cp_ice*hslyr / dt
         else
            k = 1

            kilyr = kice + betak*salin(k)/min(-puny,Tis(ij))
            kilyr = max (kilyr, kimin)

            hilyr = vicen(i,j) / aicen(i,j) 
            khis(ij) = c2*kilyr / hilyr
            khmax = rhoi*ci*hilyr / dt
         endif

         khis(ij) = min(khis(ij), khmax)

      enddo   ! ij

!lipscomb - Remove later

      if (my_task==mtest) then
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            if (i==itest .and. j==jtest) then
               print*, ''
               print*, 'kh =', khis(ij)
               print*, 'Tis =', Tis(ij)
               print*, ''
               print*, 'fswsfc =', fswsfcn(i,j)
               print*, 'flwdn =', emissivity*flw(i,j)
               print*, 'potT =', potT(i,j) - Tffresh
            endif
         enddo
      endif  ! my_task

      !-----------------------------------------------------------------
      ! Compute radiative and turbulent fluxes and their derivatives
      ! with respect to Tsf.
      !-----------------------------------------------------------------

      call surface_fluxes (nx_block,    ny_block,           &
                           isolve,      icells,             &
                           indxii,      indxjj,   indxij,   &
                           Tsf,         fswsfcn,            &
                           rhoa,        flw,                &
                           potT,        Qa,                 &
                           shcoef,      lhcoef,             &
                           flwoutn,     fsensn,             &
                           flatn,       fsurfn,             &
                           dflwout_dT,  dfsens_dT,          &
                           dflat_dT,    dfsurf_dT)

!lipscomb - for testing - remove later

      if (my_task == mtest) then
         i = itest
         j = jtest
         print*, ''
         print*, 'Got surface fluxes'

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            if (i==itest .and. j==jtest) then
               print*, ''
               print*, 'flwout =', flwoutn(i,j)
               print*, 'fsens =', fsensn(i,j)
               print*, 'flat =', flatn(i,j)
               print*, 'fsurf =', fsurfn(i,j)
               print*, 'dfsurf_dT =', dfsurf_dT(ij)
            endif
         enddo
      endif

      !-----------------------------------------------------------------
      ! Solve for the new surface temperature and fluxes
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         dTsf = (fsurfn(i,j) - khis(ij)*(Tsf(ij) - Tis(ij))) /   &
                (khis(ij) - dfsurf_dT(ij))

         Tsf(ij) = Tsf(ij) + dTsf

         if (Tsf(ij) > c0) then
            dTsf = dTsf - Tsf(ij) 
            Tsf(ij) = c0
         endif

         Tsfcn(i,j) = Tsf(ij)   ! for output

         fsensn (i,j) = fsensn (i,j) + dTsf*dfsens_dT(ij)
         flatn  (i,j) = flatn  (i,j) + dTsf*dflat_dT(ij)
         flwoutn(i,j) = flwoutn(i,j) + dTsf*dflwout_dT(ij)
         fsurfn (i,j) = fsurfn (i,j) + dTsf*dfsurf_dT(ij)

         fcondtopn(i,j) = khis(ij) * (Tsf(ij) - Tis(ij))

      enddo

!lipscomb - for testing - remove later

      if (my_task == mtest) then
         i = itest
         j = jtest
         print*, ''
         print*, 'New temp and fluxes:'

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            if (i==itest .and. j==jtest) then
               print*, ''
               print*, 'Tsf =', Tsf(ij)
               print*, 'fct =', fcondtopn(i,j)
               print*, ''
               print*, 'flwout =', flwoutn(i,j)
               print*, 'fsens =', fsensn(i,j)
               print*, 'flat =', flatn(i,j)
               print*, 'fsurf =', fsurfn(i,j)
               print*, 'dfsurf_dT =', dfsurf_dT(ij)
            endif
         enddo
      endif

      end subroutine explicit_calc_Tsfc

!=======================================================================

      end module CICE_RunMod

!=======================================================================
