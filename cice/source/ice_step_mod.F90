!=======================================================================
!
!BOP
!
! !MODULE: ice_step_mod
!
! !DESCRIPTION:
!
!  Contains CICE component driver routines common to all drivers.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2008 ECH: created module by moving subroutines from drivers/cice4/
!
! !INTERFACE:
!
      module ice_step_mod
!
! !USES:
!
!      use ice_age
      use ice_atmo
      use ice_calendar
      use ice_communicate
      use ice_diagnostics
      use ice_domain
      use ice_dyn_evp
      use ice_dyn_eap
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
      use ice_meltpond_cesm, only: compute_ponds_cesm, compute_ponds_simple
      use ice_meltpond_lvl, only: compute_ponds_lvl
      use ice_meltpond_topo, only: compute_ponds_topo
      use ice_restart_meltpond_lvl, only: ffracn, dhsn, rfracmin, rfracmax
      use ice_ocean
      use ice_orbital
      use ice_shortwave
      use ice_state
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use ice_zsalinity, only: first_ice

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: step_therm1, step_therm2, step_dynamics, &
                prep_radiation, step_radiation, post_thermo
!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: prep_radiation
!
! !DESCRIPTION:
!
! Scales radiation fields computed on the previous time step.
!
! !REVISION HISTORY:
!
! authors: David A. Bailey, NCAR
!
! !INTERFACE:

      subroutine prep_radiation (dt, iblk)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         k           , & ! vertical index       
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n               ! thickness category index

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      real (kind=dbl_kind) :: netsw 

      type (block) :: &
         this_block      ! block information for current block

      call ice_timer_start(timer_sw)      ! shortwave

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Compute netsw scaling factor (new netsw / old netsw)
      !-----------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi
            if (aice(i,j,iblk) > c0 .and. scale_factor(i,j,iblk) > puny) then
               netsw = swvdr(i,j,iblk)*(c1 - alvdr_gbm(i,j,iblk)) &
                     + swvdf(i,j,iblk)*(c1 - alvdf_gbm(i,j,iblk)) &
                     + swidr(i,j,iblk)*(c1 - alidr_gbm(i,j,iblk)) &
                     + swidf(i,j,iblk)*(c1 - alidf_gbm(i,j,iblk))
               scale_factor(i,j,iblk) = netsw / scale_factor(i,j,iblk)
            else
               scale_factor(i,j,iblk) = c1
            endif
            fswfac(i,j,iblk) = scale_factor(i,j,iblk) ! for history 
         enddo               ! i
         enddo               ! j

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
      ! Scale absorbed solar radiation for change in net shortwave
      !-----------------------------------------------------------------

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               fswsfcn(i,j,n,iblk)  = scale_factor(i,j,iblk)*fswsfcn (i,j,n,iblk)
               fswintn(i,j,n,iblk)  = scale_factor(i,j,iblk)*fswintn (i,j,n,iblk)
               fswthrun(i,j,n,iblk) = scale_factor(i,j,iblk)*fswthrun(i,j,n,iblk)
               do k = 1,nilyr+1
                  fswthruln(i,j,k,n,iblk) &
                                    = scale_factor(i,j,iblk)*fswthruln(i,j,k,n,iblk)
               enddo       !k

               Sswabsn(i,j,:,n,iblk) = &
                       scale_factor(i,j,iblk)*Sswabsn(i,j,:,n,iblk)
               Iswabsn(i,j,:,n,iblk) = &
                       scale_factor(i,j,iblk)*Iswabsn(i,j,:,n,iblk)
            enddo
         enddo                  ! ncat

      call ice_timer_stop(timer_sw)     ! shortwave

      end subroutine prep_radiation

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

      subroutine step_therm1 (dt, iblk)
!
! !USES:
!
      use ice_aerosol
      use ice_age
      use ice_firstyear
      use ice_work, only: worka, workb
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n               ! thickness category index

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

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         meltsn      , & ! snow melt in category n (m)
         vsnon_init  , & ! for aerosol mass budget
         rfrac           ! water fraction retained for melt ponds

      real (kind=dbl_kind) :: &
         raice       , & ! 1/aice
         pond            ! water retained in ponds (m)

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      l_stop = .false.

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
                                 ntrcr,              dt,                 &
                                 aice  (:,:,  iblk), frzmlt(:,:,  iblk), &
                                 vicen (:,:,:,iblk), vsnon (:,:,:,iblk), &
                                 trcrn (:,:,1:ntrcr,:,iblk),             &
                                 sst   (:,:,  iblk), Tf    (:,:,  iblk), &
                                 strocnxT(:,:,iblk), strocnyT(:,:,iblk), &
                                 Tbot,               fbot,               &
                                 rside (:,:,  iblk) )

         do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify cells with nonzero ice area
      !-----------------------------------------------------------------
            melttn(:,:,n,iblk)  = c0
            meltbn(:,:,n,iblk)  = c0
            congeln(:,:,n,iblk) = c0
            snoicen(:,:,n,iblk) = c0
            dsnown(:,:,n,iblk) = c0 
!            Tsf_icen(:,:,n,iblk) = c0
           
            icells = 0
            do j = jlo, jhi
            do i = ilo, ihi
#ifdef oned
               if (i==4.and.j==4) then
#else
               if (aicen(i,j,n,iblk) > puny) then
#endif
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif
            enddo               ! i
            enddo               ! j

            if ((calc_Tsfc .or. calc_strair) .and. icells > 0) then 

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
                                    trcrn(:,:,nt_Tsfc,n,iblk),      &
                                    potT(:,:,iblk), Qa  (:,:,iblk), &
                                    worka,          workb,          &
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
               ! Set to data values (on T points)
               strairxn(:,:) = strax(:,:,iblk)
               strairyn(:,:) = stray(:,:,iblk)
#else
               ! NEMO wind stress is supplied on u grid, multipied 
               ! by ice concentration and set directly in evp, so
               ! strairxT/yT = 0. Zero u-components here for safety.
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
            if (tr_FY) then
               call update_FYarea (nx_block, ny_block,      &
                                   dt, icells,              &
                                   indxi, indxj,            &
                                   lmask_n(:,:,iblk),       &
                                   lmask_s(:,:,iblk),       &
                                   trcrn(:,:,nt_FY,n,iblk))
            endif

      !-----------------------------------------------------------------
      ! Vertical thermodynamics: Heat conduction, growth and melting.
      !----------------------------------------------------------------- 

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

            vsnon_init(:,:) = vsnon(:,:,n,iblk)

            call thermo_vertical(nx_block,           ny_block,            &
                                dt,                  icells,              &
                                indxi,               indxj,               &
                                aicen(:,:,n,iblk),                        &
                                trcrn(:,:,:,n,iblk),                      &
                                vicen(:,:,n,iblk),   vsnon(:,:,n,iblk),   &
                                flw    (:,:,iblk),   potT (:,:,iblk),     &
                                Qa     (:,:,iblk),   rhoa (:,:,iblk),     &
                                fsnow  (:,:,iblk),   fpond (:,:,iblk),    &
                                fbot,                Tbot,                &
                                sss  (:,:,iblk),                          &
                                lhcoef,              shcoef,              &
                                fswsfcn(:,:,n,iblk), fswintn(:,:,n,iblk), &
                                fswthrun(:,:,n,iblk),                     &
                                Sswabsn(:,:,:,n,iblk),                    &
                                Iswabsn(:,:,:,n,iblk),                    &
                                fsurfn(:,:,n,iblk),  fcondtopn(:,:,n,iblk),&
                                fsensn,              flatn(:,:,n,iblk),   &
                                fswabsn,             flwoutn,             &
                                evapn,               freshn,              &
                                fsaltn,              fhocnn,              &
                                melttn(:,:,n,iblk),  meltsn,              &
                                meltbn(:,:,n,iblk),                       &
                                congeln(:,:,n,iblk), snoicen(:,:,n,iblk), &
                                mlt_onset(:,:,iblk), frz_onset(:,:,iblk), &
                                yday,                l_stop,              &
                                istop,               jstop,               &
                                dsnown(:,:,n,iblk),                       &
                                fsicen(:,:,n,iblk)) !,  Tsf_icen(:,:,n,iblk))

         if (l_stop) then
            write (nu_diag,*) 'istep1, my_task, iblk =', &
                               istep1, my_task, iblk
            write (nu_diag,*) 'category n = ', n
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
      ! Aerosol update
      !-----------------------------------------------------------------
         if (tr_aero .and. icells > 0) then

               call update_aerosol (nx_block, ny_block,                  &
                                    dt, icells,                          &
                                    indxi, indxj,                        &
                                    melttn(:,:,n,iblk),  meltsn,         &
                                    meltbn(:,:,n,iblk),                  &
                                    congeln(:,:,n,iblk),                 &
                                    snoicen(:,:,n,iblk),                 &
                                    fsnow(:,:,iblk),                     &
                                    trcrn(:,:,:,n,iblk),                 &
                                    aicen_init(:,:,n,iblk),              &
                                    vicen_init(:,:,n,iblk),              &
                                    vsnon_init(:,:),                     &
                                    vicen(:,:,n,iblk),                   &
                                    vsnon(:,:,n,iblk),                   &
                                    aicen(:,:,n,iblk),                   &
                                    faero_atm(:,:,:,iblk),               &
                                    faero_ocn(:,:,:,iblk))
         endif

      !-----------------------------------------------------------------
      ! Melt ponds
      ! If using tr_pond_cesm, the full calculation is performed here.
      ! If using tr_pond_topo, the rest of the calculation is done after
      ! the surface fluxes are merged, below.
      !-----------------------------------------------------------------

         if (tr_pond) then
            call ice_timer_start(timer_ponds)

            if (tr_pond_cesm) then
!               rfrac(:,:) = 0.15_dbl_kind + 0.7_dbl_kind * aicen(:,:,n,iblk)
               rfrac(:,:) = rfracmin + (rfracmax-rfracmin) * aicen(:,:,n,iblk) 
              call compute_ponds_cesm(nx_block, ny_block,                      &
                                       ilo, ihi, jlo, jhi,                      &
                                       rfrac, melttn(:,:,n,iblk), meltsn, frain(:,:,iblk),  &
                                       aicen (:,:,n,iblk), vicen (:,:,n,iblk),  &
                                       vsnon (:,:,n,iblk), trcrn (:,:,:,n,iblk))

            elseif (tr_pond_lvl) then
               rfrac(:,:) = rfracmin + (rfracmax-rfracmin) * aicen(:,:,n,iblk)
               call compute_ponds_lvl(nx_block, ny_block,                      &
                                      ilo, ihi, jlo, jhi,                      &
                                      rfrac,                                   &
                                      melttn(:,:,n,iblk), meltsn,              &
                                      frain (:,:,iblk),   Tair  (:,:,iblk),    &
                                      fsurfn(:,:,n,iblk),                      &
                                      dhsn  (:,:,n,iblk), ffracn(:,:,n,iblk),  &
                                      aicen (:,:,n,iblk), vicen (:,:,n,iblk),  &
                                      vsnon (:,:,n,iblk),                      &
                                      trcrn (:,:,:,n,iblk))

            elseif (tr_pond_topo .and. icells > 0) then
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

                  ! collect liquid water in ponds
                  ! assume salt still runs off

                  rfrac(i,j) = rfracmin + (rfracmax-rfracmin) * aicen(i,j,n,iblk)
                  pond = rfrac(i,j)/rhofresh * (melttn(i,j,n,iblk)*rhoi &
                       +                        meltsn(i,j       )*rhos &
                       +                        frain (i,j,iblk)*dt)

                  ! if pond does not exist, create new pond over full ice area
                  ! otherwise increase pond depth without changing pond area
                  if (trcrn(i,j,nt_apnd,n,iblk) < puny) then
                      trcrn(i,j,nt_hpnd,n,iblk) = c0
                      trcrn(i,j,nt_apnd,n,iblk) = c1
                  endif 
                  trcrn(i,j,nt_hpnd,n,iblk) = (pond &
                  + trcrn(i,j,nt_hpnd,n,iblk)*trcrn(i,j,nt_apnd,n,iblk)) &
                                            / trcrn(i,j,nt_apnd,n,iblk)
                  fpond(i,j,iblk) = fpond(i,j,iblk) &
                                  + pond * aicen(i,j,n,iblk) ! m
               enddo
            endif

            call ice_timer_stop(timer_ponds)
         endif

         if (tr_pond .and. trim(shortwave) /= 'dEdd') then
            call ice_timer_start(timer_ponds)

            rfrac(:,:) = c1

            call compute_ponds_simple(nx_block, ny_block,                      &
                                   ilo, ihi, jlo, jhi,                      &
                                   rfrac, melttn(:,:,n,iblk), meltsn, frain(:,:,iblk),  &
                                   aicen (:,:,n,iblk), vicen (:,:,n,iblk),  &
                                   vsnon (:,:,n,iblk), trcrn (:,:,:,n,iblk))

            call ice_timer_stop(timer_ponds)
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
                            fhocn   (:,:,iblk), fswthru   (:,:,iblk), &
                            melttn  (:,:,n,iblk), meltsn,             &
                            meltbn(:,:,n,iblk), congeln(:,:,n,iblk),  &
                            snoicen(:,:,n,iblk),                      &
                            meltt   (:,:,iblk),  melts   (:,:,iblk),  &
                            meltb   (:,:,iblk),                       &
                            congel  (:,:,iblk),  snoice  (:,:,iblk))

         enddo                  ! ncat

      !-----------------------------------------------------------------
      ! Calculate ponds from the topographic scheme
      !-----------------------------------------------------------------
            if (tr_pond_topo) then
               call ice_timer_start(timer_ponds)
               call compute_ponds_topo(nx_block, ny_block,                &
                                    ilo, ihi, jlo, jhi,                   &
                                    aice (:,:,  iblk), aicen(:,:,:,iblk), &
                                    vice (:,:,  iblk), vicen(:,:,:,iblk), &
                                    vsno (:,:,  iblk), vsnon(:,:,:,iblk), &
                                    trcrn(:,:,:,:,iblk),                  &
                                    potT(:,:,  iblk),  meltt(:,:,iblk),   &
                                    fsurf(:,:,iblk),   fpond(:,:,iblk))
               call ice_timer_stop(timer_ponds)
            endif

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

      subroutine step_therm2 (dt, iblk)
!
! !USES:
!
      use ice_therm_mushy, only: add_new_ice_mushy
      use ice_therm_bl99, only: add_new_ice_bl99
      use ice_therm_oned, only: diagnose_itd
      use ice_zbgc_public, only: ocean_bio
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, ij

      integer (kind=int_kind) :: &
         icells          ! number of ice/ocean cells 

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for ice/ocean cells

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts

      l_stop = .false.
      
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

      !-----------------------------------------------------------------
      ! Compute fractional ice area in each grid cell.
      !-----------------------------------------------------------------
         call aggregate_area (nx_block,          ny_block, &
                              aicen(:,:,:,iblk),           &
                              aice (:,:,  iblk), aice0(:,:,iblk))

         if (kitd == 1) then
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

            call linear_itd (nx_block, ny_block,             &
                             icells,   indxi, indxj,         &
                             ntrcr,    trcr_depend(1:ntrcr), &
                             aicen_init(:,:,:,iblk),         &
                             vicen_init(:,:,:,iblk),         &
                             aicen     (:,:,:,iblk),         &
                             trcrn     (:,:,1:ntrcr,:,iblk), & 
                             vicen     (:,:,:,iblk),         &
                             vsnon     (:,:,:,iblk),         &
                             aice      (:,:,  iblk),         &
                             aice0     (:,:,  iblk),         &
                             fpond     (:,:,  iblk),         &
                             l_stop,                         &
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

            endif ! icells

         endif  ! kitd = 1

         call ice_timer_stop(timer_catconv)    ! category conversions

      !-----------------------------------------------------------------
      ! Add frazil ice growing in leads.
      !-----------------------------------------------------------------

         ! identify ice-ocean cells
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
            
         if (ktherm == 2) then
         call add_new_ice_mushy (nx_block,              ny_block, &
                           ntrcr,                 icells,   &
                           indxi,                 indxj,    &
                           tmask     (:,:,  iblk), dt,      &
                           aicen     (:,:,:,iblk),          &
                           trcrn     (:,:,1:ntrcr,:,iblk),  &
                           vicen     (:,:,:,iblk),          &
                           aice0     (:,:,  iblk),          &
                           aice      (:,:,  iblk),          &
                           frzmlt    (:,:,  iblk),          &
                           frazil    (:,:,  iblk),          &
                           frz_onset (:,:,  iblk), yday,    &
                           update_ocn_f,                    &
                           fresh     (:,:,  iblk),          &
                           fsalt     (:,:,  iblk),          &
                           Tf        (:,:,  iblk),          &
                           sss       (:,:,  iblk),          &
                           phi_init, dSin0_frazil,          &
                           l_stop,                          &
                           istop                 , jstop)

         else

         call add_new_ice_bl99 (nx_block,              ny_block, &
                           ntrcr,                 icells,   &
                           indxi,                 indxj,    &
                           tmask     (:,:,  iblk), dt,      &
                           aicen     (:,:,:,iblk),          &
                           trcrn     (:,:,1:ntrcr,:,iblk),  &
                           vicen     (:,:,:,iblk),          &
                           aice0     (:,:,  iblk),          &
                           aice      (:,:,  iblk),          &
                           frzmlt    (:,:,  iblk),          &
                           frazil    (:,:,  iblk),          &
                           frz_onset (:,:,  iblk), yday,    &
                           update_ocn_f,                    &
                           fresh     (:,:,  iblk),          &
                           fsalt     (:,:,  iblk),          &
                           Tf        (:,:,  iblk),          &
                           sss       (:,:,  iblk),          &
                           salinz    (:,:,:,iblk), l_stop,  &
                           istop                 , jstop,   &
                           fsice(:,:,iblk), &
                           flux_bio(:,:,:,iblk),nbltrcr, &
                           ocean_bio(:,:,:,iblk))
         endif

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

#ifndef oned
         call lateral_melt (nx_block, ny_block,     &
                            ilo, ihi, jlo, jhi,     &
                            dt,                     &
                            fpond     (:,:,  iblk), &
                            fresh     (:,:,  iblk), &
                            fsalt     (:,:,  iblk), &    
                            fhocn     (:,:,  iblk), &
                            faero_ocn (:,:,:,iblk), &
                            rside     (:,:,  iblk), &
                            meltl     (:,:,  iblk), &
                            aicen     (:,:,:,iblk), &
                            vicen     (:,:,:,iblk), &
                            vsnon     (:,:,:,iblk), &
                            trcrn     (:,:,:,:,iblk),&
                            fsice     (:,:,  iblk), &
                            flux_bio(:,:,:,iblk),nbltrcr)

      !-----------------------------------------------------------------
      ! For the special case of a single category, adjust the area and
      ! volume (assuming that half the volume change decreases the
      ! thickness, and the other half decreases the area).  
      !-----------------------------------------------------------------

!echmod: test this
         if (ncat==1) &
             call reduce_area (nx_block, ny_block,     &
                               ilo, ihi, jlo, jhi,     &
                               tmask     (:,:,  iblk), &
                               aicen     (:,:,1,iblk), &
                               vicen     (:,:,1,iblk), &
                               aicen_init(:,:,1,iblk), &
                               vicen_init(:,:,1,iblk))
#endif
         
      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

         call cleanup_itd (nx_block,             ny_block,             &
                           ilo, ihi,             jlo, jhi,             &
                           dt,                   ntrcr,                &
                           aicen   (:,:,:,iblk),                       &
                           trcrn (:,:,1:ntrcr,:,iblk),                 &
                           vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
                           aice0   (:,:,  iblk), aice      (:,:,iblk), &
                           trcr_depend(1:ntrcr), fpond     (:,:,iblk), &
                           fresh   (:,:,  iblk), fsalt     (:,:,iblk), &
                           fhocn   (:,:,  iblk),                       &
                           faero_ocn(:,:,:,iblk),tr_aero,              &
                           tr_pond_topo,         heat_capacity,        &
                           nbltrcr ,             first_ice(:,:,iblk),  &
                           fsice(:,:,  iblk),    flux_bio(:,:,:,iblk), &
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
            call abort_ice ('ice: ITD cleanup error in step_therm2')
         endif

         !call diagnose_itd(nx_block, ny_block, aicen(:,:,:,1), vicen(:,:,:,1), vsnon(:,:,:,1), trcrn(:,:,:,:,1))

      end subroutine step_therm2

!=======================================================================
!BOP
!
! !ROUTINE: post_thermo: finalize thermo updates
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: Elizabeth Hunke, LANL
!
! !INTERFACE:

      subroutine post_thermo
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: & 
         iblk        , & ! block index 
         i,j             ! horizontal indices

      !-------------------------------------------------------------------
      ! Ghost cell updates for state variables.
      !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)
      call bound_state (aicen, trcrn, &
                        vicen, vsnon)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Aggregate the updated state variables (includes ghost cells). 
      !----------------------------------------------------------------- 
 
         call aggregate (nx_block,          ny_block,             &
                         aicen(:,:,:,iblk),                       &
                         trcrn(:,:,1:ntrcr,:,iblk),               &
                         vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
                         aice (:,:,  iblk),                       &
                         trcr (:,:,1:ntrcr,  iblk),               &
                         vice (:,:,  iblk), vsno (:,:,    iblk),  &
                         aice0(:,:,  iblk), tmask(:,:,    iblk),  &
                         ntrcr, trcr_depend(1:ntrcr)) 

      !-----------------------------------------------------------------
      ! Compute thermodynamic area and volume tendencies.
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            daidtt(i,j,iblk) = (aice(i,j,iblk) - daidtt(i,j,iblk)) / dt
            dvidtt(i,j,iblk) = (vice(i,j,iblk) - dvidtt(i,j,iblk)) / dt
         enddo
         enddo

      enddo ! iblk
      !$OMP END PARALLEL DO

      end subroutine post_thermo

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

      subroutine step_dynamics (dt, ndtd)

      use ice_state, only: nt_qsno, trcrn, vsnon, aicen
      use ice_domain_size, only: nslyr
      use ice_calendar, only: istep
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         ndtd    ! number of dynamics subcycles
!
!EOP
!
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: & 
         iblk        , & ! block index 
         i,j         , & ! horizontal indices
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      call init_history_dyn     ! initialize dynamic history variables

      !-----------------------------------------------------------------
      ! Elastic-viscous-plastic ice dynamics
      !-----------------------------------------------------------------

      if (kdyn == 1) call evp (dt)
      if (kdyn == 2) call eap (dt)

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

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call step_ridge (dt, ndtd, iblk)
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_stop(timer_ridge)

      !-------------------------------------------------------------------
      ! Ghost cell updates for state variables.
      !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)
      call bound_state (aicen, trcrn, &
                        vicen, vsnon)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Aggregate the updated state variables (includes ghost cells). 
      !----------------------------------------------------------------- 
 
         call aggregate (nx_block,          ny_block,             &
                         aicen(:,:,:,iblk),                       &
                         trcrn(:,:,1:ntrcr,:,iblk),               &
                         vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
                         aice (:,:,  iblk),                       &
                         trcr (:,:,1:ntrcr,  iblk),               &
                         vice (:,:,  iblk), vsno (:,:,    iblk),  &
                         aice0(:,:,  iblk), tmask(:,:,    iblk),  &
                         ntrcr, trcr_depend(1:ntrcr)) 

      !-----------------------------------------------------------------
      ! Compute dynamic area and volume tendencies.
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo,jhi
         do i = ilo,ihi
            dvidtd(i,j,iblk) = (vice(i,j,iblk) - dvidtd(i,j,iblk)) /dt
            daidtd(i,j,iblk) = (aice(i,j,iblk) - daidtd(i,j,iblk)) /dt
         enddo
         enddo

      enddo
      !$OMP END PARALLEL DO

      call ice_timer_stop(timer_column)

      end subroutine step_dynamics

!=======================================================================
!BOP
!
! !ROUTINE: step_ridge
!
! !DESCRIPTION:
!
! Computes sea ice mechanical deformation
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !INTERFACE:

      subroutine step_ridge (dt, ndtd, iblk)

! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         ndtd, & ! number of dynamics subcycles
         iblk            ! block index
!
!EOP
!
      real (kind=dbl_kind) :: &
         dtt      ! thermo time step

      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: & 
         i,j         , & ! horizontal indices
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts

         l_stop = .false.

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
                         dt,                   ndtd,                     &
                         ntrcr,                icells,                   &
                         indxi,                indxj,                    &
                         rdg_conv(:,:,  iblk), rdg_shear (:,:,  iblk),   &
                         aicen   (:,:,:,iblk),                           &
                         trcrn     (:,:,1:ntrcr,:,iblk),                 &
                         vicen   (:,:,:,iblk), vsnon     (:,:,:,iblk),   &
                         aice0   (:,:,  iblk),                           &
                         trcr_depend(1:ntrcr), l_stop,                   &
                         istop,                jstop,                    &   
                         dardg1dt(:,:,iblk),   dardg2dt  (:,:,iblk),     &
                         dvirdgdt(:,:,iblk),   opening   (:,:,iblk),     &
                         fpond   (:,:,iblk),                             &
                         fresh   (:,:,iblk),   fhocn     (:,:,iblk),     &
                         fsicen(:,:,:,iblk),   faero_ocn(:,:,:,iblk),    &
                         aparticn(:,:,:,iblk), krdgn     (:,:,:,iblk),   &
                         aredistn(:,:,:,iblk), vredistn  (:,:,:,iblk),   &
                         dardg1ndt(:,:,:,iblk),dardg2ndt (:,:,:,iblk),   &
                         dvirdgndt(:,:,:,iblk),                          &
                         araftn   (:,:,:,iblk),vraftn   (:,:,:,iblk))

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
    
      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

         dtt = dt * ndtd  ! for proper averaging over thermo timestep
         call cleanup_itd (nx_block,             ny_block,             &
                           ilo, ihi,             jlo, jhi,             &
                           dtt,                  ntrcr,                &
                           aicen   (:,:,:,iblk),                       &
                           trcrn (:,:,1:ntrcr,:,iblk),                 &
                           vicen   (:,:,:,iblk), vsnon (:,:,  :,iblk), &
                           aice0   (:,:,  iblk), aice      (:,:,iblk), &
                           trcr_depend(1:ntrcr), fpond     (:,:,iblk), &
                           fresh   (:,:,  iblk), fsalt     (:,:,iblk), &
                           fhocn   (:,:,  iblk),                       &
                           faero_ocn(:,:,:,iblk),tr_aero,              &
                           tr_pond_topo,         heat_capacity,        &
                           nbltrcr ,             first_ice(:,:,iblk),  &
                           fsice(:,:,  iblk),    flux_bio(:,:,:,iblk), &
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
            call abort_ice ('ice: ITD cleanup error in step_ridge')
         endif

      end subroutine step_ridge

!=======================================================================
!BOP
!
! !ROUTINE: step_radiation
!
! !DESCRIPTION:
!
! Computes radiation fields
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          David Bailey, NCAR
!          Elizabeth C. Hunke, LANL
!
! !INTERFACE:

      subroutine step_radiation (dt, iblk)

! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk            ! block index
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         n               ! thickness category index

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      call ice_timer_start(timer_sw)      ! shortwave

      ! Initialize
      do n = 1, ncat
      do j = 1, ny_block
      do i = 1, nx_block
         alvdrn(i,j,n,iblk) = c0
         alidrn(i,j,n,iblk) = c0
         alvdfn(i,j,n,iblk) = c0
         alidfn(i,j,n,iblk) = c0
         fswsfcn(i,j,n,iblk) = c0
         fswintn(i,j,n,iblk) = c0
         fswthrun(i,j,n,iblk) = c0
      enddo   ! i
      enddo   ! j
      enddo   ! ncat
      fswthruln(:,:,:,:,iblk) = c0
      Iswabsn(:,:,:,:,iblk) = c0
      Sswabsn(:,:,:,:,iblk) = c0

      if (trim(shortwave) == 'dEdd') then ! delta Eddington

         call run_dEdd (iblk)

      else  ! .not. dEdd

         call shortwave_ccsm3(iblk,                                   &
                              nx_block,     ny_block,                 &
                              aicen(:,:,:,iblk), vicen(:,:,:,iblk),   &
                              vsnon(:,:,:,iblk),                      &
                              trcrn(:,:,nt_Tsfc,:,iblk),              &
                              swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                              swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                              alvdrn(:,:,:,iblk),alidrn(:,:,:,iblk),  &
                              alvdfn(:,:,:,iblk),alidfn(:,:,:,iblk),  &
                              fswsfcn(:,:,:,iblk),fswintn(:,:,:,iblk),&
                              fswthrun(:,:,:,iblk),                   &
                              fswthruln(:,:,:,:,iblk),                &
                              Iswabsn(:,:,:,:,iblk),                  &
                              Sswabsn(:,:,:,:,iblk),                  &
                              albicen(:,:,:,iblk),albsnon(:,:,:,iblk),&
                              coszen(:,:,iblk))

      endif   ! shortwave

      call ice_timer_stop(timer_sw)     ! shortwave

      end subroutine step_radiation

!=======================================================================

      end module ice_step_mod

!=======================================================================
