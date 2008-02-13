!=======================================================================
!
!BOP
!
! !MODULE: ice_step_mod - thermodynamic and dynamics step routines
!
! !DESCRIPTION:
!
!  Contains CICE thermodynamic and dynamics step routines
!
! !REVISION HISTORY:
!  SVN:$Id: 
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2006 ECH: moved exit timeLoop to prevent execution of unnecessary timestep
! 2006 ECH: Streamlined for efficiency 
! 2006 ECH: Converted to free source form (F90)
! 2007 MV : Moved thermodynamic and dynamics step routines to ice_step_mod
!
! !INTERFACE:
!
      module ice_step_mod
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
#if (defined CCSM) || (defined SEQ_MCT)
      use ice_prescribed_mod
#endif

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: step_therm1, step_therm2, step_dynamics, &
                step_rad1, step_rad2
!
!EOP
!
!=======================================================================

      contains

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
         flatn       , & ! surface downward latent heat       (W/m^2)
         fswabsn     , & ! shortwave absorbed by ice          (W/m^2)
         flwoutn     , & ! upward LW at surface               (W/m^2)
         evapn       , & ! flux of vapor, atmos to ice   (kg m-2 s-1)
         freshn      , & ! flux of water, ice to ocean     (kg/m^2/s)
         fsaltn      , & ! flux of salt, ice to ocean      (kg/m^2/s)
         fhocnn      , & ! fbot corrected for leftover energy (W/m^2)
         strairxn    , & ! air/ice zonal  strss,              (N/m^2)
         strairyn    , & ! air/ice merdnl strss,              (N/m^2)
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

      call ice_timer_start(timer_column)  ! column physics
      call ice_timer_start(timer_thermo)  ! thermodynamics

      call init_history_therm    ! initialize thermo history variables
      call init_flux_ocn        ! initialize ocean fluxes sent to coupler

      if (oceanmixed_ice) &
           call ocean_mixed_layer (dt)   ! ocean surface fluxes and sst

!      call ice_timer_start(timer_tmp)  ! temporary timer

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
                                 nghost,             dt,                 &
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

      !-----------------------------------------------------------------
      ! Atmosphere boundary layer calculation; compute coefficients
      ! for sensible and latent heat fluxes.
      !
      ! NOTE: The wind stress is computed here for later use if 
      !       calc_strair = .true.   Otherwise, the wind stress
      !       components are set to the data values.
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

            if (.not.(calc_strair)) then
               strairxn(:,:) = strax(:,:,iblk)
               strairyn(:,:) = stray(:,:,iblk)
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

            call thermo_vertical                                       &
                            (nx_block,            ny_block,            &
                             dt,                  icells,              &
                             indxi,               indxj,               &
                             aicen(:,:,n,iblk),                        &
                             trcrn(:,:,nt_Tsfc,n,iblk),                &
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
                             fsensn,              flatn,               &
                             fswabsn,             flwoutn,             &
                             evapn,               freshn,              &
                             fsaltn,              fhocnn,              &
                             meltt   (:,:,iblk),  melts(:,:,iblk),     &
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

            call compute_ponds(nx_block, ny_block, nghost,              &
                               meltt_tmp, melts_tmp,                    &
                               frain(:,:,iblk),                         &
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
                            fsensn,             flatn,                &
                            fswabsn,            flwoutn,              &
                            evapn,                                    &
                            Trefn,              Qrefn,                &
                            freshn,             fsaltn,               &
                            fhocnn,             fswthrun(:,:,n,iblk), &
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

         if (prescribed_ice) then

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

         endif

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

      enddo              ! iblk

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

      enddo              ! iblk

      call ice_timer_stop(timer_column)

      end subroutine step_dynamics

!=======================================================================
!BOP
!
! !ROUTINE: step_rad1 - step pre-thermo radiation
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: David A. Bailey, NCAR
!
! !INTERFACE:

      subroutine step_rad1 (dt)
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

      ! snow variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsn             ! snow horizontal fraction
      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) :: &
         rhosnwn     , & ! snow density (kg/m3)
         rsnwn           ! snow grain radius (micro-meters)

      ! pond variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fpn         , & ! pond fraction
         hpn             ! pond depth (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         scale_factor
 
      real (kind=dbl_kind) :: netsw, netsw_old, ar

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      l_stop = .false.

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Compute cosine of solar zenith angle.
      ! This is used by the delta-Eddington shortwave module.
      ! Albedos are aggregated in merge_fluxes only for cells w/ coszen > 0.
      ! For basic shortwave, simply set coszen to a constant between 0 and 1.
      !-----------------------------------------------------------------

         if (trim(shortwave) == 'dEdd') then ! delta Eddington

            scale_factor(:,:) = c1

            do j = jlo, jhi
            do i = ilo, ihi
               if (aice(i,j,iblk) > c0) then
               netsw = swvdr(i,j,iblk)*(c1 - alvdr(i,j,iblk)) &
                     + swvdf(i,j,iblk)*(c1 - alvdf(i,j,iblk)) &
                     + swidr(i,j,iblk)*(c1 - alidr(i,j,iblk)) &
                     + swidf(i,j,iblk)*(c1 - alidf(i,j,iblk))
               netsw_old = c0
               do n=1, ncat
                  netsw_old = netsw_old + (fswsfcn(i,j,n,iblk) &
                            + fswintn(i,j,n,iblk) &
                            + fswthrun(i,j,n,iblk)) * aicen(i,j,n,iblk)
               enddo
               ar = c1 / aice(i,j,iblk)
               netsw_old = netsw_old * ar
               if (netsw_old > c0) then
                  scale_factor(i,j) = netsw / netsw_old
!                 if (my_task == 0 .and. i == 19 .and. j == 29) then
!                    print *,'coszen,netsw,netsw_old,scale_factor', &
!                         my_task, i, j, iblk, istep, &
!                         coszen(i,j,iblk),netsw,netsw_old,scale_factor(i,j)
!                    print *,'aicen,aice',aicen(i,j,n,iblk),aice(i,j,iblk)
!                 endif
               endif
               endif
            enddo               ! i
            enddo               ! j

         endif

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

            il1 = ilyr1(n)
            il2 = ilyrn(n)
            sl1 = slyr1(n)
            sl2 = slyrn(n)

            if (trim(shortwave) == 'dEdd') then   ! delta Eddington
              
               do ij=1,icells
                  i = indxi(ij)
                  j = indxj(ij)
                  fswsfcn(i,j,n,iblk) = scale_factor(i,j)*fswsfcn(i,j,n,iblk)
                  fswintn(i,j,n,iblk) = scale_factor(i,j)*fswintn(i,j,n,iblk)
                  fswthrun(i,j,n,iblk) = scale_factor(i,j)*fswthrun(i,j,n,iblk)
                  Sswabsn(i,j,sl1:sl2,iblk) = scale_factor(i,j)*Sswabsn(i,j,sl1:sl2,iblk)
                  Iswabsn(i,j,il1:il2,iblk) = scale_factor(i,j)*Iswabsn(i,j,il1:il2,iblk)
                  if (scale_factor(i,j) > 1000._dbl_kind) then
!                    print *,'fswsfcn,fswintn,fswthrun', &
!                       fswsfcn(i,j,n,iblk),fswintn(i,j,n,iblk),fswthrun(i,j,n,iblk)
!                    print *,'Sswabsn', Sswabsn(i,j,sl1:sl2,iblk)
!                    print *,'Iswabsn', Iswabsn(i,j,il1:il2,iblk)
                  endif
               enddo

            else

               Sswabsn(i,j,sl1:sl2,iblk) = c0

               call absorbed_solar  (nx_block,   ny_block,               &
                               icells,                                   &
                               indxi,      indxj,                        &
                               aicen(:,:,n,iblk),                        &
                               vicen(:,:,n,iblk),  vsnon(:,:,n,iblk),    &
                               swvdr(:,:,iblk),    swvdf(:,:,iblk),      &
                               swidr(:,:,iblk),     swidf(:,:,iblk),     &
                               alvdrni(:,:,n,iblk), alvdfni(:,:,n,iblk), &
                               alidrni(:,:,n,iblk), alidfni(:,:,n,iblk), &
                               alvdrns(:,:,n,iblk), alvdfns(:,:,n,iblk), &
                               alidrns(:,:,n,iblk), alidfns(:,:,n,iblk), &
                               fswsfcn(:,:,n,iblk), fswintn(:,:,n,iblk), &
                               fswthrun(:,:,n,iblk),                     &
                               Iswabsn(:,:,il1:il2,iblk))

            endif
         enddo                  ! ncat

      enddo                      ! iblk

      end subroutine step_rad1

!=======================================================================
!BOP
!
! !ROUTINE: step_rad2 - step pre-coupler radiation
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: David A. Bailey, NCAR
!
! !INTERFACE:

      subroutine step_rad2 (dt)
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

      ! snow variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsn             ! snow horizontal fraction
      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) :: &
         rhosnwn     , & ! snow density (kg/m3)
         rsnwn           ! snow grain radius (micro-meters)

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

      l_stop = .false.

      alvdr(:,:,:) = c0
      alvdf(:,:,:) = c0
      alidr(:,:,:) = c0
      alidf(:,:,:) = c0
      Sswabsn(:,:,:,:) = c0

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

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
                                 coszen(:,:,iblk), dt)

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

            il1 = ilyr1(n)
            il2 = ilyrn(n)
            sl1 = slyr1(n)
            sl2 = slyrn(n)

            if (trim(shortwave) == 'dEdd') then   ! delta Eddington

      ! note that rhoswn, rsnw, fp, hp and Sswabs ARE NOT dimensioned with ncat
      ! BPB 19 Dec 2006

               ! set snow properties
               call shortwave_dEdd_set_snow(nx_block, ny_block,           &
                                 icells,                                  &
                                 indxi,               indxj,              &
                                 aicen(:,:,n,iblk),   vsnon(:,:,n,iblk),  &
                                 trcrn(:,:,nt_Tsfc,n,iblk), fsn,          &
                                 rhosnwn,             rsnwn)


               if (.not. tr_pond) then

               ! set pond properties
               call shortwave_dEdd_set_pond(nx_block, ny_block,            &
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

               call shortwave_dEdd(nx_block,        ny_block,            &
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
                                 fswsfcn(:,:,n,iblk),fswintn(:,:,n,iblk),&
                                 fswthrun(:,:,n,iblk),                   &
                                 Sswabsn(:,:,sl1:sl2,iblk),              &
                                 Iswabsn(:,:,il1:il2,iblk))

! Special case of night to day

               do ij=1,icells
                  i = indxi(ij)
                  j = indxj(ij)
                  fswsfcn(i,j,n,iblk) = max(p01, fswsfcn(i,j,n,iblk))
               enddo

            else

               call compute_albedos (nx_block,   ny_block, &
                               icells,               &
                               indxi,      indxj,    &
                               aicen(:,:,n,iblk), vicen(:,:,n,iblk),    &
                               vsnon(:,:,n,iblk),                       &
                               trcrn(:,:,nt_Tsfc,n,iblk),               &
                               alvdrni(:,:,n,iblk),alidrni(:,:,n,iblk), &
                               alvdfni(:,:,n,iblk),alidfni(:,:,n,iblk), &
                               alvdrns(:,:,n,iblk),alidrns(:,:,n,iblk), &
                               alvdfns(:,:,n,iblk),alidfns(:,:,n,iblk), &
                               alvdrn(:,:,n,iblk),alidrn(:,:,n,iblk),   &
                               alvdfn(:,:,n,iblk),alidfn(:,:,n,iblk),   &
                               apondn(:,:,n,iblk),hpondn(:,:,n,iblk))

            endif


         ! Aggregate albedos for coupler

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               alvdf(i,j,iblk) = alvdf(i,j,iblk) &
                  + alvdfn(i,j,n,iblk)*aicen(i,j,n,iblk)
               alidf(i,j,iblk) = alidf(i,j,iblk) &
                  + alidfn(i,j,n,iblk)*aicen(i,j,n,iblk)
               alvdr(i,j,iblk) = alvdr(i,j,iblk) &
                  + alvdrn(i,j,n,iblk)*aicen(i,j,n,iblk)
               alidr(i,j,iblk) = alidr(i,j,iblk) &
                  + alidrn(i,j,n,iblk)*aicen(i,j,n,iblk)

            enddo

         enddo                  ! ncat

      !-----------------------------------------------------------------
      ! Divide fluxes by ice area for the coupler, which assumes fluxes
      ! are per unit ice area.
      !-----------------------------------------------------------------

         call scale_fluxes (nx_block,            ny_block,           &
                            nghost,              tmask   (:,:,iblk), &
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

      enddo                      ! iblk

      call scale_hist_fluxes     ! to match coupler fluxes
      
      end subroutine step_rad2

!=======================================================================

      end module ice_step_mod

!=======================================================================
