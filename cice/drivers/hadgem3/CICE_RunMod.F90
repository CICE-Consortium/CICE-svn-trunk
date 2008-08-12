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
!  SVN:$Id: CICE_RunMod.F90 131 2008-05-30 16:53:40Z eclare $
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2006 ECH: moved exit timeLoop to prevent execution of unnecessary timestep
! 2006 ECH: Streamlined for efficiency 
! 2006 ECH: Converted to free source form (F90)
! 2007 BPB: Modified Delta-Eddington shortwave interface
! 2008 ECH: moved ESMF code to its own driver
!
! !INTERFACE:
!

      module CICE_RunMod
!
! !USES:
!
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
      use ice_step_mod
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

      public :: CICE_Run, ice_step, step_therm1 
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
!EOP
!BOC
!
   !--------------------------------------------------------------------
   !  local variables
   !--------------------------------------------------------------------

      integer (kind=int_kind) :: k

   !--------------------------------------------------------------------
   !  initialize error code and step timer
   !--------------------------------------------------------------------

      call ice_timer_start(timer_step)   ! start timing entire run

#ifndef CICE_IN_NEMO
   !--------------------------------------------------------------------
   ! timestep loop
   !--------------------------------------------------------------------

      timeLoop: do
#endif

         call ice_step

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date
         call calendar(time)    ! at the end of the timestep

#ifndef CICE_IN_NEMO
         if (stop_now >= 1) exit timeLoop
#endif

#ifndef coupled
         call ice_timer_start(timer_couple)  ! atm/ocn coupling
         call get_forcing_atmo     ! atmospheric forcing from data
         call get_forcing_ocn(dt)  ! ocean forcing from data
         call ice_timer_stop(timer_couple)   ! atm/ocn coupling
#endif

         call init_flux_atm     ! initialize atmosphere fluxes sent to coupler
         call init_flux_ocn     ! initialize ocean fluxes sent to coupler

#ifndef CICE_IN_NEMO
      enddo timeLoop
#endif

   !--------------------------------------------------------------------
   ! end of timestep loop
   !--------------------------------------------------------------------

      call ice_timer_stop(timer_step)   ! end timestepping loop timer     
!
!EOC
!
      end subroutine CICE_Run

!=======================================================================
!BOP
!
! !ROUTINE: ice_step
!
! !DESCRIPTION:
!
!  Calls drivers for physics components, some initialization, and output
!
! !REVISION HISTORY:
!
!  author Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ice_step
!
!EOP
!BOC
!
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
      ! Scale radiation fields
      !-----------------------------------------------------------------

         call prep_radiation (dt)

      !-----------------------------------------------------------------
      ! thermodynamics
      !-----------------------------------------------------------------

         call step_therm1 (dt)  ! vertical thermodynamics

         call step_therm2 (dt)  ! ice thickness distribution thermo

      !-----------------------------------------------------------------
      ! dynamics, transport, ridging
      !-----------------------------------------------------------------

         do k = 1, ndyn_dt
            call step_dynamics (dyn_dt)
         enddo

      !-----------------------------------------------------------------
      ! albedo, shortwave radiation
      !-----------------------------------------------------------------

         call step_radiation (dt)

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
         strairxn    , & ! air/ice zonal  stress,             (N/m^2)
         strairyn    , & ! air/ice meridional stress,         (N/m^2)
         Trefn       , & ! air tmp reference level                (K)
         Qrefn           ! air sp hum reference level         (kg/kg)

      ! other local variables
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         Tbot        , & ! ice bottom surface temperature (deg C)
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

      ! Local variables to keep track of melt for ponds
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         melts_old, &
         meltt_old, &
         melts_tmp, &
         meltt_tmp

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      real (kind=dbl_kind) :: &
         raice           ! 1/aice

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

#ifdef CICE_IN_NEMO
       !---------------------------------------------------------------
       ! Scale frain and fsnow by ice concentration as these fields
       ! are supplied by NEMO multiplied by ice concentration
       !---------------------------------------------------------------
 
         do j = 1, ny_block
         do i = 1, nx_block

            if (aice_init(i,j,iblk) > puny) then
               raice           = c1 / aice_init(i,j,iblk)
               frain(i,j,iblk) = frain(i,j,iblk)*raice
               fsnow(i,j,iblk) = fsnow(i,j,iblk)*raice
            else
               frain(i,j,iblk) = c0
               fsnow(i,j,iblk) = c0
            endif

         enddo
         enddo
#endif

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

#ifndef CICE_IN_NEMO
! Do not do the following for CICE_IN_NEMO as wind stress is supplied on
! u grid and multipied by ice concentration and set directly in evp
!
               ! Set to data values (on T points)
               strairxn(:,:) = strax(:,:,iblk)
               strairyn(:,:) = stray(:,:,iblk)
#else
               ! Zero for safety
               strairxn(:,:) = c0
               strairyn(:,:) = c0
#endif
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

            if (.not.(calc_Tsfc)) then

               ! If not calculating surface temperature and fluxes, set 
               ! surface fluxes (flatn, fsurfn, and fcondtopn) to be used 
               ! in thickness_changes
 
               ! hadgem routine sets fluxes to default values in ice-only mode
               call set_sfcflux(nx_block,  ny_block,  &
                                n,         iblk,      &
                                icells,               & 
                                indxi,     indxj,     &
                                aicen    (:,:,n,iblk),&
                                flatn    (:,:,n,iblk),&
                                fsurfn   (:,:,n,iblk),&
                                fcondtopn(:,:,n,iblk) )

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
                             fswsfcn(:,:,n,iblk), fswintn(:,:,n,iblk), &
                             fswthrun(:,:,n,iblk),                     &
                             Sswabsn(:,:,sl1:sl2,iblk),                &
                             Iswabsn(:,:,il1:il2,iblk),                &
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

         if (tr_pond .and. trim(shortwave) == 'dEdd') then

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
#ifdef CICE_IN_NEMO
      !-----------------------------------------------------------------
      ! Note that for CICE_IN_NEMO, strairxT/yT = 0 as wind stress 
      ! supplied on u grid and set directly in evp
      !-----------------------------------------------------------------
#endif
         call merge_fluxes (nx_block,           ny_block,             &
                            icells,                                   &
                            indxi,              indxj,                &
                            aicen_init(:,:,n,iblk),                   &
                            flw(:,:,iblk),      coszen(:,:,iblk),     &
                            strairxn,           strairyn,             &
                            fsurfn(:,:,n,iblk), fcondtopn(:,:,n,iblk),&
                            fsensn,             flatn(:,:,n,iblk),    &
                            fswabsn,            flwoutn,              &
                            evapn,                                    &
                            Trefn,              Qrefn,                &
                            freshn,             fsaltn,               &
                            fhocnn,             fswthrun(:,:,n,iblk), &
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
         n           , & ! thickness category index
         i,j         , & ! horizontal indices
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      call ice_timer_start(timer_column)

      !-----------------------------------------------------------------
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

      if (oceanmixed_ice) &
         call ocean_mixed_layer (dt)   ! ocean surface fluxes and sst

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Aggregate albedos
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = c0
            alidf(i,j,iblk) = c0
            alvdr(i,j,iblk) = c0
            alidr(i,j,iblk) = c0
         enddo
         enddo
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = alvdf(i,j,iblk) &
               + alvdfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidf(i,j,iblk) = alidf(i,j,iblk) &
               + alidfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alvdr(i,j,iblk) = alvdr(i,j,iblk) &
               + alvdrn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidr(i,j,iblk) = alidr(i,j,iblk) &
               + alidrn(i,j,n,iblk)*aicen(i,j,n,iblk)
         enddo
         enddo
         enddo

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf_gbm  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_gbm  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_gbm  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_gbm  (i,j,iblk) = alidr  (i,j,iblk)
            fresh_gbm  (i,j,iblk) = fresh  (i,j,iblk)
            fsalt_gbm  (i,j,iblk) = fsalt  (i,j,iblk)
            fhocn_gbm  (i,j,iblk) = fhocn  (i,j,iblk)
            fswthru_gbm(i,j,iblk) = fswthru(i,j,iblk)

      !-----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !-----------------------------------------------------------------
            scale_factor(i,j,iblk) = &
                       swvdr(i,j,iblk)*(c1 - alvdr_gbm(i,j,iblk)) &
                     + swvdf(i,j,iblk)*(c1 - alvdf_gbm(i,j,iblk)) &
                     + swidr(i,j,iblk)*(c1 - alidr_gbm(i,j,iblk)) &
                     + swidf(i,j,iblk)*(c1 - alidf_gbm(i,j,iblk))

         enddo
         enddo

      !-----------------------------------------------------------------
      ! Divide fluxes by ice area 
      !  - the CCSM coupler assumes fluxes are per unit ice area
      !  - also needed for global budget in diagnostics
      !-----------------------------------------------------------------

         call scale_fluxes (nx_block,            ny_block,           &
                            tmask    (:,:,iblk),                     &
                            aice     (:,:,iblk), Tf      (:,:,iblk), &
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

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (scale_factor,     halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

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

      real (kind=dbl_kind)  :: &
         raicen          ! 1/aicen

      logical (kind=log_kind) :: &
         extreme_flag    ! flag for extreme forcing values

      logical (kind=log_kind), parameter :: & 
         extreme_test=.true. ! test and write out extreme forcing data
!
!EOP
!
         raicen        = c1
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

#ifdef CICE_IN_NEMO
!----------------------------------------------------------------------
! Convert fluxes from GBM values to per ice area values when 
! running in NEMO environment.  (When in standalone mode, fluxes
! are input as per ice area.)
!----------------------------------------------------------------------
            raicen        = c1 / aicen(i,j)
#endif
            fsurfn(i,j)   = fsurfn_f(i,j,n,iblk)*raicen
            fcondtopn(i,j)= fcondtopn_f(i,j,n,iblk)*raicen
            flatn(i,j)    = flatn_f(i,j,n,iblk)*raicen

         enddo

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

      end module CICE_RunMod

!=======================================================================
