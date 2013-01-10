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
! 2008 ECH: moved ESMF code to its own driver
!
! !INTERFACE:
!
      module CICE_RunMod
!
! !USES:
!
      use ice_aerosol
      use ice_age
      use ice_algae
      use ice_atmo
      use ice_brine
      use ice_calendar
      use ice_communicate
      use ice_diagnostics
      use ice_domain
      use ice_dyn_eap
      use ice_dyn_evp
      use ice_exit
      use ice_fileunits
      use ice_firstyear
      use ice_flux
      use ice_forcing
      use ice_grid
      use ice_history
      use ice_restart
      use ice_init
      use ice_itd
      use ice_kinds_mod
      use ice_lvl
      use ice_mechred
      use ice_meltpond_cesm
      use ice_meltpond_lvl
      use ice_meltpond_topo
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
      use ice_zbgc
      use ice_zsalinity
      use ice_domain_size, only: nltrcr

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: CICE_Run, ice_step
!
!EOP
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

   !--------------------------------------------------------------------
   ! timestep loop
   !--------------------------------------------------------------------

      timeLoop: do

         call ice_step

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date

         call calendar(time)    ! at the end of the timestep

         if (stop_now >= 1) exit timeLoop

#ifndef coupled
         call ice_timer_start(timer_couple)  ! atm/ocn coupling
         call get_forcing_atmo     ! atmospheric forcing from data
         call get_forcing_ocn(dt)  ! ocean forcing from data
!         if (tr_aero) call faero_data        ! aerosols
         if (tr_aero) call faero_default     ! aerosols
         if (read_Sin)  call get_ice_salinity  !update ice salinity from file
         if (tr_bgc_NO .OR. tr_bgc_Sil) call get_forcing_bgc    ! biogeochemistry
         call ice_timer_stop(timer_couple)   ! atm/ocn coupling
#endif

         call init_flux_atm     ! initialize atmosphere fluxes sent to coupler
         call init_flux_ocn     ! initialize ocean fluxes sent to coupler

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

      use ice_state, only: nt_qsno, trcrn
      use ice_domain_size, only: nslyr
      use ice_calendar, only: istep

      integer (kind=int_kind) :: &
         iblk        , & ! block index 
         k               ! dynamics supercycling index

      !-----------------------------------------------------------------
      ! restoring on grid boundaries
      !-----------------------------------------------------------------

         if (restore_ice) call ice_HaloRestore

      !-----------------------------------------------------------------
      ! initialize diagnostics
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics/history
         call init_mass_diags   ! diagnostics per timestep
         call init_history_therm
         call init_history_bgc
         call ice_timer_stop(timer_diags)   ! diagnostics/history

         call ice_timer_start(timer_column)  ! column physics
         call ice_timer_start(timer_thermo)  ! thermodynamics

         do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Scale radiation fields
      !-----------------------------------------------------------------

            if (calc_Tsfc) call prep_radiation (dt, iblk)

      !-----------------------------------------------------------------
      ! thermodynamics
      !-----------------------------------------------------------------
            
            call step_therm1     (dt, iblk) ! vertical thermodynamics
            call biogeochemistry (dt, iblk) ! biogeochemistry
            call step_therm2     (dt, iblk) ! ice thickness distribution thermo

         enddo ! iblk

         call post_thermo               ! finalize thermo update

         call ice_timer_stop(timer_thermo) ! thermodynamics
         call ice_timer_stop(timer_column) ! column physics

      !-----------------------------------------------------------------
      ! dynamics, transport, ridging
      !-----------------------------------------------------------------

         do k = 1, ndtd
            call step_dynamics (dt_dyn, ndtd)
         enddo

      !-----------------------------------------------------------------
      ! albedo, shortwave radiation
      !-----------------------------------------------------------------

         call ice_timer_start(timer_column)  ! column physics
         call ice_timer_start(timer_thermo)  ! thermodynamics

         do iblk = 1, nblocks

            call step_radiation (dt, iblk)

      !-----------------------------------------------------------------
      ! get ready for coupling and the next time step
      !-----------------------------------------------------------------

            call coupling_prep (iblk)

         enddo ! iblk

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (scale_factor,     halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_timer_stop(timer_bound)

         call ice_timer_stop(timer_thermo) ! thermodynamics
         call ice_timer_stop(timer_column) ! column physics

      !-----------------------------------------------------------------
      ! write data
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics
         if (mod(istep,diagfreq) == 0) call runtime_diags(dt) ! log file
         call ice_timer_stop(timer_diags)   ! diagnostics

         call ice_timer_start(timer_diags_bgc)
         if (mod(istep,diagfreq) == 0) then
            if (tr_bgc_S)                 call S_diags   (dt)
            if (solve_bgc .OR. tr_bgc_NO) call bgc_diags (dt)
         endif
         call ice_timer_stop(timer_diags_bgc)

         call ice_timer_start(timer_hist)   ! history
         call accum_hist (dt)               ! history file
         call ice_timer_stop(timer_hist)    ! history

         call ice_timer_start(timer_readwrite)  ! reading/writing
         if (write_restart == 1) then
            if (restart_ext) then
               call dumpfile_ext ! core variables for restarting
            else
               call dumpfile     ! core variables for restarting
            endif
            if (tr_iage)      call write_restart_age
            if (tr_FY)        call write_restart_FY
            if (tr_lvl)       call write_restart_lvl
            if (tr_pond_cesm) call write_restart_pond_cesm
            if (tr_pond_lvl)  call write_restart_pond_lvl
            if (tr_pond_topo) call write_restart_pond_topo
            if (nltrcr > 0)   call write_restart_bgc  ! for layer model
            if (tr_bgc_S)     call write_restart_S    ! for layer model
            if (kdyn == 2)    call write_restart_eap
         endif

         call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine ice_step
    
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

      subroutine coupling_prep (iblk)!

! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: & 
         iblk            ! block index 
!
!EOP
!
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: & 
         n           , & ! thickness category index
         i,j         , & ! horizontal indices
         k           , & ! tracer index
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      real (kind=dbl_kind) :: cszn ! counter for history averaging

      !-----------------------------------------------------------------
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

         if (oceanmixed_ice) &
         call ocean_mixed_layer (dt,iblk) ! ocean surface fluxes and sst

      !-----------------------------------------------------------------
      ! Aggregate albedos
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = c0
            alidf(i,j,iblk) = c0
            alvdr(i,j,iblk) = c0
            alidr(i,j,iblk) = c0

            albice(i,j,iblk) = c0
            albsno(i,j,iblk) = c0
            albpnd(i,j,iblk) = c0
            apeff_ai(i,j,iblk) = c0

            ! for history averaging
            cszn = c0
            if (coszen(i,j,iblk) > puny) cszn = c1
            do n = 1, nstreams
               albcnt(i,j,iblk,n) = albcnt(i,j,iblk,n) + cszn
            enddo
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

            if (coszen(i,j,iblk) > puny) then ! sun above horizon
            albice(i,j,iblk) = albice(i,j,iblk) &
               + albicen(i,j,n,iblk)*aicen(i,j,n,iblk)
            albsno(i,j,iblk) = albsno(i,j,iblk) &
               + albsnon(i,j,n,iblk)*aicen(i,j,n,iblk)
            albpnd(i,j,iblk) = albpnd(i,j,iblk) &
               + albpndn(i,j,n,iblk)*aicen(i,j,n,iblk)
            endif

            apeff_ai(i,j,iblk) = apeff_ai(i,j,iblk) &       ! for history
               + apeffn(i,j,n,iblk)*aicen(i,j,n,iblk)
         enddo
         enddo
         enddo

         do j = 1, ny_block
         do i = 1, nx_block

      !-----------------------------------------------------------------
      ! reduce fresh by fpond for coupling
      !-----------------------------------------------------------------

            fpond(i,j,iblk) = fpond(i,j,iblk) * rhofresh/dt
            fresh(i,j,iblk) = fresh(i,j,iblk) - fpond(i,j,iblk)

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

            alvdf_gbm  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_gbm  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_gbm  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_gbm  (i,j,iblk) = alidr  (i,j,iblk)
            fresh_gbm  (i,j,iblk) = fresh  (i,j,iblk)
            fsalt_gbm  (i,j,iblk) = fsalt  (i,j,iblk)
            fhocn_gbm  (i,j,iblk) = fhocn  (i,j,iblk)
            fswthru_gbm(i,j,iblk) = fswthru(i,j,iblk)
            fsice_gbm  (i,j,iblk) = fsice  (i,j,iblk)
            fsice_g_gbm  (i,j,iblk) = fsice_g  (i,j,iblk)

            do k = 1,nbltrcr
              flux_bio_gbm  (i,j,k,iblk) = flux_bio (i,j,k,iblk)
              flux_bio_g_gbm  (i,j,k,iblk) = flux_bio_g  (i,j,k,iblk)
            enddo

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
                            tmask    (:,:,iblk), nbltrcr,            &
                            aice     (:,:,iblk), Tf      (:,:,iblk), &
                            Tair     (:,:,iblk), Qa      (:,:,iblk), &
                            strairxT (:,:,iblk), strairyT(:,:,iblk), &
                            fsens    (:,:,iblk), flat    (:,:,iblk), &
                            fswabs   (:,:,iblk), flwout  (:,:,iblk), &
                            evap     (:,:,iblk),                     &
                            Tref     (:,:,iblk), Qref    (:,:,iblk), &
                            fresh    (:,:,iblk), fsalt   (:,:,iblk), &
                            fhocn    (:,:,iblk), fswthru (:,:,iblk), &
                            faero_ocn(:,:,:,iblk),                   &
                            alvdr    (:,:,iblk), alidr   (:,:,iblk), &
                            alvdf    (:,:,iblk), alidf   (:,:,iblk), &
                            fsice    (:,:,iblk), fsice_g (:,:,iblk), &
                            flux_bio(:,:,:,iblk),flux_bio_g(:,:,:,iblk))
 
!echmod - comment this out for efficiency, if .not. calc_Tsfc
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
!echmod

      end subroutine coupling_prep

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
