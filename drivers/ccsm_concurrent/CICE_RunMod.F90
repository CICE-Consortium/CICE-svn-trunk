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
!
! !INTERFACE:
!

      module CICE_RunMod
!
! !USES:
!
      use ice_atmo
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
      use ice_state
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use ice_work
      use ice_prescribed_mod

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

         call ice_step

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

      subroutine ice_step

      integer (kind=int_kind) :: k

         call ice_timer_start(timer_rcvsnd)   ! timing between send-recv

         call init_mass_diags   ! diagnostics per timestep

         if(prescribed_ice) then  ! read prescribed ice
            call ice_prescribed_run(idate, sec)
         endif

      !-----------------------------------------------------------------
      ! thermodynamics
      !-----------------------------------------------------------------

         call step_therm1 (dt)  ! pre-coupler thermodynamics

         call ice_timer_stop(timer_rcvsnd)   ! timing between recv-send

         call to_coupler        ! collect/send data to CCSM coupler

         call ice_timer_start(timer_sndrcv)   ! timing between send-recv

         if (.not.prescribed_ice) then

         call step_therm2 (dt)  ! post-coupler thermodynamics

      !-----------------------------------------------------------------
      ! dynamics, transport, ridging
      !-----------------------------------------------------------------

         do k = 1, ndyn_dt
            call step_dynamics (dyn_dt)
         enddo

         endif ! not prescribed_ice


      !-----------------------------------------------------------------
      ! write data
      !-----------------------------------------------------------------

         call ice_timer_start(timer_readwrite)  ! reading/writing

         if (mod(istep,diagfreq) == 0) call runtime_diags(dt) ! log file

         call ice_write_hist (dt)    ! history file

         if (write_restart == 1) call dumpfile ! dumps for restarting

         call ice_timer_stop(timer_readwrite)  ! reading/writing

         call ice_timer_stop(timer_sndrcv)   ! timing between send-recv

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
         i, j        , & ! horizontal indices
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
         flatn       , & ! surface downward latent heat       (W/m^2)
         fswabsn     , & ! shortwave absorbed by ice          (W/m^2)
         flwoutn     , & ! upward LW at surface               (W/m^2)
         evapn       , & ! flux of vapor, atmos to ice   (kg m-2 s-1)
         freshn      , & ! flux of water, ice to ocean     (kg/m^2/s)
         fsaltn      , & ! flux of salt, ice to ocean      (kg/m^2/s)
         fhocnn      , & ! fbot corrected for leftover energy (W/m^2)
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         strairxn    , & ! air/ice zonal  strss,              (N/m^2)
         strairyn    , & ! air/ice merdnl strss,              (N/m^2)
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

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr) :: &
         Iswabsn         ! SW radiation absorbed in ice layers (W m-2)

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      call ice_timer_start(timer_column)  ! column physics
      call ice_timer_start(timer_thermo)  ! thermodynamics

      call init_history_therm    ! initialize thermo history variables

      if (oceanmixed_ice) &
           call ocean_mixed_layer (dt)   ! ocean surface fluxes and sst

!      call ice_timer_start(timer_tmp)  ! temporary timer

      call init_flux_atm         ! initialize atmosphere fluxes sent to coupler

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
                                 nghost,             dt,                 &
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
                                 coszen(:,:,iblk))

         else                     ! basic (ccsm3) shortwave
            coszen(:,:,iblk) = p5 ! sun above the horizon
         endif

         do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify cells with nonzero ice area
      !-----------------------------------------------------------------
           
            icells = 0
            do j = jlo, jhi
            do i = ilo, ihi
!               if (aicen(i,j,n,iblk) > puny .and. .not.tmask(i,j,iblk)) then
!                 print*,my_task,i,j,n,iblk,aicen(i,j,n,iblk),tmask(i,j,iblk)
!               endif
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

            if (trim(shortwave) == 'dEdd') then   ! delta Eddington

               call shortwave_dEdd(nx_block,        ny_block,            &
                                 icells,                                 &
                                 indxi,             indxj,               &
                                 tlat (:,:,  iblk), coszen(:,:, iblk),   &
                                 aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                                 vsnon(:,:,n,iblk), trcrn(:,:,1,n,iblk), &
                                 swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                                 swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                                 alvdrn(:,:,n,iblk),alidrn(:,:,n,iblk),  &
                                 alvdfn(:,:,n,iblk),alidfn(:,:,n,iblk),  &
                                 fswsfcn,           fswintn,             &
                                 fswthrun,          Iswabsn)

            else
               call shortwave_ccsm3(nx_block,       ny_block,            &
                                 icells,                                 &
                                 indxi,             indxj,               &
                                 aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                                 vsnon(:,:,n,iblk), trcrn(:,:,1,n,iblk), &
                                 swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                                 swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                                 alvdrn(:,:,n,iblk),alidrn(:,:,n,iblk),  &
                                 alvdfn(:,:,n,iblk),alidfn(:,:,n,iblk),  &
                                 fswsfcn,           fswintn,             &
                                 fswthrun,          Iswabsn)
            endif

      !-----------------------------------------------------------------
      ! Atmosphere boundary layer calculation; compute coefficients
      ! for sensible and latent heat fluxes.
      !
      ! NOTE: The wind stress is computed here for later use.
      !-----------------------------------------------------------------

            if (trim(atmbndy) == 'constant') then
               call atmo_boundary_const(nx_block,      ny_block,        &
                                        'ice',          icells,         &
                                        indxi,          indxj,          &
                                        uatm(:,:,iblk), vatm(:,:,iblk), &
                                        wind(:,:,iblk), rhoa(:,:,iblk), &
                                        strairxn,       strairyn,       &
                                        lhcoef,         shcoef)
            else ! default
               call atmo_boundary_layer(nx_block,       ny_block,       &
                                        'ice',          icells,         &
                                        indxi,          indxj,          &
                                        trcrn(:,:,1,n,iblk),            &
                                        potT(:,:,iblk),                 &
                                        uatm(:,:,iblk), vatm(:,:,iblk), &
                                        wind(:,:,iblk), zlvl(:,:,iblk), &
                                        Qa  (:,:,iblk), rhoa(:,:,iblk), &
                                        strairxn,       strairyn,       &
                                        Trefn,          Qrefn,          &
                                        worka,          workb,          &
                                        lhcoef,         shcoef)
            endif ! atmbndy

      !-----------------------------------------------------------------
      ! Vertical thermodynamics: Heat conduction, growth and melting.
      !----------------------------------------------------------------- 

            il1 = ilyr1(n)
            il2 = ilyrn(n)
            sl1 = slyr1(n)
            sl2 = slyrn(n)

            call thermo_vertical                                       &
                            (nx_block,            ny_block,            &
                             dt,                  icells,              &
                             indxi,               indxj,               &
                             aicen(:,:,n,iblk),   trcrn(:,:,1,n,iblk), &
                             vicen(:,:,n,iblk),   vsnon(:,:,n,iblk),   &
                             eicen  (:,:,il1:il2,iblk),                &
                             esnon  (:,:,sl1:sl2,iblk),                &
                             flw    (:,:,iblk),   potT (:,:,iblk),     &
                             Qa     (:,:,iblk),   rhoa (:,:,iblk),     &
                             fsnow  (:,:,iblk),                        &
                             fbot,                Tbot,                &
                             lhcoef,              shcoef,              &
                             fswsfcn,             fswintn,             &
                             fswthrun,            Iswabsn,             &
                             fsensn,              flatn,               &
                             fswabsn,             flwoutn,             &
                             evapn,               freshn,              &
                             fsaltn,              fhocnn,              &
                             meltt   (:,:,iblk),  meltb   (:,:,iblk),  &
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
                            fsensn,             flatn,                &
                            fswabsn,            flwoutn,              &
                            evapn,                                    &
                            Trefn,              Qrefn,                &
                            freshn,             fsaltn,               &
                            fhocnn,             fswthrun,             &
                            alvdr   (:,:,iblk), alidr     (:,:,iblk), &
                            alvdf   (:,:,iblk), alidf     (:,:,iblk), &
                            strairxT(:,:,iblk), strairyT  (:,:,iblk), &
                            fsens   (:,:,iblk), flat      (:,:,iblk), &
                            fswabs  (:,:,iblk), flwout    (:,:,iblk), &
                            evap    (:,:,iblk),                       &
                            Tref    (:,:,iblk), Qref      (:,:,iblk), &
                            fresh   (:,:,iblk), fresh_hist(:,:,iblk), &
                            fsalt   (:,:,iblk), fsalt_hist(:,:,iblk), &
                            fhocn   (:,:,iblk), fhocn_hist(:,:,iblk), &
                            fswthru (:,:,iblk), fswthru_hist(:,:,iblk))

         enddo                  ! ncat

      !-----------------------------------------------------------------
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

         if (oceanmixed_ice) then
            do j = jlo, jhi
            do i = ilo, ihi
               if (hmix(i,j,iblk) > puny) then
                  sst(i,j,iblk) = sst(i,j,iblk) &
                       + (fhocn(i,j,iblk) + fswthru(i,j,iblk))*dt &
                       / (cprho*hmix(i,j,iblk))
               endif
            enddo
            enddo
         endif

!      call ice_timer_stop(timer_tmp)  ! temporary timer

      !-----------------------------------------------------------------
      ! Divide fluxes by ice area for the coupler, which assumes fluxes
      ! are per unit ice area.
      !-----------------------------------------------------------------

         call scale_fluxes (nx_block,            ny_block,           &
                            nghost,              tmask   (:,:,iblk), &
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

      enddo                      ! iblk

      call ice_timer_stop(timer_column) ! column physics
      call ice_timer_stop(timer_thermo) ! thermodynamics

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
!         hicen           ! ice thickness (m)

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
!      call ice_timer_start(timer_tmp)  ! temporary timer

      call init_flux_ocn        ! initialize ocean fluxes sent to coupler
      
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
            fresh_hist(i,j,iblk) = fresh_hist(i,j,iblk)  &
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
            call aggregate_area (nx_block, ny_block, &
                           aicen(:,:,:,iblk), &
                           aice(:,:,iblk),     aice0(:,:,iblk),  &
                           l_stop, &
                           istop,    jstop)

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

            if (icells > 0) &
            call linear_itd (nx_block, ny_block,       &
                             icells, indxi, indxj,     &
                             nghost,   trcr_depend,    &
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
                           tmask    (:,:,  iblk), dt,       &
                           aicen    (:,:,:,iblk),           &
                           trcrn    (:,:,:,:,iblk),         &
                           vicen    (:,:,:,iblk),           &
                           eicen    (:,:,:,iblk),           &
                           aice0    (:,:,  iblk),           &
                           aice     (:,:,  iblk),           &
                           frzmlt   (:,:,  iblk),           &
                           frazil   (:,:,  iblk),           &
                           frz_onset(:,:,  iblk), yday,     &
                           Tf       (:,:,  iblk), l_stop,   &
                           istop, jstop)

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
                            nghost,   dt,           &
                            fresh     (:,:,  iblk), &
                            fsalt     (:,:,  iblk), &    
                            fhocn     (:,:,  iblk), &
                            fresh_hist(:,:,  iblk), &
                            fsalt_hist(:,:,  iblk), &
                            fhocn_hist(:,:,  iblk), &
                            rside     (:,:,  iblk), &
                            meltl     (:,:,  iblk), &
                            aicen     (:,:,:,iblk), &
                            vicen     (:,:,:,iblk), &
                            vsnon     (:,:,:,iblk), &
                            eicen     (:,:,:,iblk), &
                            esnon     (:,:,:,iblk) )

         call ice_timer_stop(timer_thermo) ! thermodynamics

      !-----------------------------------------------------------------
      ! For the special case of a single category, adjust the area and
      ! volume (assuming that half the volume change decreases the
      ! thickness, and the other half decreases the area).  
      !-----------------------------------------------------------------

!NOTE - this does not work - hicen_init is not defined - ECH

!         if (ncat==1) &
!              call reduce_area (nx_block, ny_block,     &
!                                nghost,                 &
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
                           nghost,               dt,                   &
                           aicen   (:,:,:,iblk), trcrn (:,:,:,:,iblk), &
                           vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
                           eicen   (:,:,:,iblk), esnon (:,:,  :,iblk), &
                           aice0   (:,:,  iblk), aice      (:,:,iblk), &
                           trcr_depend,                                &
                           fresh   (:,:,  iblk), fresh_hist(:,:,iblk), &
                           fsalt   (:,:,  iblk), fsalt_hist(:,:,iblk), &
                           fhocn   (:,:,  iblk), fhocn_hist(:,:,iblk), &
                           l_stop,                                     &
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

      !-----------------------------------------------------------------
      ! Aggregate the updated state variables. 
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

!      call ice_timer_stop(timer_tmp)  ! temporary timer
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
                         fresh   (:,:,iblk),   fresh_hist(:,:,iblk),     &
                         fhocn   (:,:,iblk),   fhocn_hist(:,:,iblk))      

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

      enddo                     ! iblk

      call ice_timer_stop(timer_ridge)

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk), iblk)

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

         call cleanup_itd (nx_block,             ny_block,             &
                           nghost,               dt,                   &
                           aicen   (:,:,:,iblk), trcrn (:,:,:,:,iblk), &
                           vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
                           eicen   (:,:,:,iblk), esnon (:,:,  :,iblk), &
                           aice0   (:,:,  iblk), aice      (:,:,iblk), &
                           trcr_depend,                                &
                           fresh   (:,:,  iblk), fresh_hist(:,:,iblk), &
                           fsalt   (:,:,  iblk), fsalt_hist(:,:,iblk), &
                           fhocn   (:,:,  iblk), fhocn_hist(:,:,iblk), &
                           l_stop,                                     &
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

      !-----------------------------------------------------------------
      ! Aggregate the updated state variables. 
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

      enddo

      call ice_timer_stop(timer_column)

      call scale_hist_fluxes    ! to match coupler fluxes

      end subroutine step_dynamics

!=======================================================================

      end module CICE_RunMod

!=======================================================================
