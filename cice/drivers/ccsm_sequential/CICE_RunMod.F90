!=======================================================================
!
!BOP
!
! !MODULE: CICE_RunMod
!
! !DESCRIPTION:
!
!  Contains CICE component driver routines particular to the specific
!  driver.
!
! !REVISION HISTORY:
!  SVN:$Id: CICE_RunMod.F90 138 2008-07-08 20:39:37Z eclare $
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2008 ECH: created module by moving subroutines from drivers/cice4/
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
!      use ice_exit
      use ice_fileunits
      use ice_flux
      use ice_forcing
      use ice_grid
      use ice_history
      use ice_restart
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

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: step_therm1, coupling_prep
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
      use ice_work, only: worka, workb
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

            if (.not.(calc_strair)) then
               ! Set to data values (on T points)
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

!echmod
      !-----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !-----------------------------------------------------------------
            scale_factor(i,j,iblk) = &
                       swvdr(i,j,iblk)*(c1 - alvdr_gbm(i,j,iblk)) &
                     + swvdf(i,j,iblk)*(c1 - alvdf_gbm(i,j,iblk)) &
                     + swidr(i,j,iblk)*(c1 - alidr_gbm(i,j,iblk)) &
                     + swidf(i,j,iblk)*(c1 - alidf_gbm(i,j,iblk))
!echmod

         enddo
         enddo

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

      enddo

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (scale_factor,     halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      call ice_timer_stop(timer_column)

      end subroutine coupling_prep

!=======================================================================

      end module CICE_RunMod

!=======================================================================
