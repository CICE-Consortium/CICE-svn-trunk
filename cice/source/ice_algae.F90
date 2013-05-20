
!=======================================================================
!
!BOP
!
! !MODULE: ice_algae - Biogeochemistry
!
! !DESCRIPTION:
!
! Compute biogeochemistry vertically resolved.
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
! Vertical Transport of biological scalars is split between implicit diffusion
! (tridiag_solver) and advection.  A nondimensional z coordinate is used which 
! ranges from 0 (ice surface or snow/ice interface)  to 1 (ice bottom) and has 
! the number of layers nblyr while the ice has nilyr. Changes in ice thickness are 
! included as advective fluxes.  Initially, the only scalar modeled is nitrate (NO).  
! Reaction terms will follow original ice_bgc.F90.   
!  
!
!
! !INTERFACE:
!
      module ice_algae
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: nblyr, nilyr, nblyr_hist, max_blocks, nbltrcr
      use ice_blocks, only: nx_block, ny_block
      use ice_fileunits, only: nu_diag, nu_restart_bgc, nu_rst_pointer, nu_dump_bgc, flush_fileunit
      use ice_read_write, only: ice_open, ice_read, ice_write
      use ice_timers, only: timer_bgc2, timer_bgc3, ice_timer_start, ice_timer_stop
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_state, only: vicen, vice, trcr, ntrcr, ntraceb, nt_bgc_am_sk, nt_bgc_c_sk, &
           nt_bgc_chl_sk, nt_bgc_DMS_sk, nt_bgc_DMSPd_sk, nt_bgc_DMSPp_sk, nt_bgc_N, &
           nt_bgc_N_sk, nt_bgc_Nit_sk, nt_bgc_NO, nt_bgc_PON, nt_bgc_Sil, nt_bgc_Sil_sk, &
           nt_bgc_C, nt_bgc_chl, nt_bgc_DMS, nt_bgc_DMSPd, nt_bgc_DMSPp, nt_bgc_NH, &
           nt_bgc_NO, nt_bgc_PON, nt_bgc_Sil, nt_bgc_Nit_sk, nt_bgc_Sil_sk, nt_bgc_Nit_sk, &
           nt_bgc_Sil_sk
      use ice_zbgc_public, only: sil, nit, algalN, amm, dms, dmsp, R_C2N, R_chl2N, tr_bgc_C, &
           tr_bgc_chl, tr_bgc_DMSPd, tr_bgc_DMSPp, tr_bgc_N, tr_bgc_NH, tr_bgc_NO, tr_bgc_PON, &
           tr_bgc_Sil, tr_bgc_DMSPd, tr_bgc_DMSPp, tr_bgc_N, tr_bgc_NH, tr_bgc_NO, tr_bgc_PON, &
           tr_bgc_Sil, tr_bgc_N, tr_bgc_N_sk, tr_bgc_NO, tr_bgc_PON, tr_bgc_Sil, grid_o, &
           grid_o_t, fr_resp, nlt_bgc_DMS, nlt_bgc_DMSPd, nlt_bgc_N, nlt_bgc_NH, nlt_bgc_NO, &
           nlt_bgc_PON, nlt_bgc_Sil, tr_bgc_DMS, Dm, initbio_frac, min_salin, nlt_bgc_C, &
           nlt_bgc_chl, nlt_bgc_DMS, nlt_bgc_DMSPd, nlt_bgc_DMSPp, nlt_bgc_N, nlt_bgc_NH, &
           nlt_bgc_NO, nlt_bgc_PON, nlt_bgc_Sil, thin, nit_data_type, nit_file, &
           restore_bgc, sil_data_type, sil_file, zfswin, grownp, remap_layers_bgc_plus
!
!EOP
!
      implicit none

      private
      public :: get_forcing_bgc, bgc_diags, write_restart_bgc, &
                algal_dynamics, tracer_transport, read_restart_bgc

      real (kind=dbl_kind), parameter, private :: &
         R_Si2N  = 1.5_dbl_kind , & ! algal Si to N (mole/mole)
         R_S2N   = 0.03_dbl_kind, & ! algal S to N (mole/mole)
         zphimin  = 0.03_dbl_kind    ! minimum porosity for bgc only

!=======================================================================

      contains


!=======================================================================
!BOP
!
! !IROUTINE: get_forcing_bgc
!
! !INTERFACE:
!
      subroutine get_forcing_bgc
!
! !DESCRIPTION:
!
! Read and interpolate annual climatologies of silicate and nitrate.
! Restore model quantities to data if desired.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_calendar, only: dt, istep, mday, month, sec
      use ice_domain, only: nblocks
      use ice_forcing, only: trestore, trest, ocn_data_dir, &
          read_clim_data, interpolate_data, interp_coeff_monthly
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
          i, j, iblk  , & ! horizontal indices
          ixm,ixp     , & ! record numbers for neighboring months
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth        ! middle day of month

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
          nitdat      , & ! data value toward which nitrate is restored
          sildat          ! data value toward which silicate is restored

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks) :: &
         nit_data, & ! field values at 2 temporal data points
         sil_data

      logical (kind=log_kind) :: readm

      if (trim(nit_data_type) == 'clim'.or. &
          trim(sil_data_type) == 'clim') then

         nit_file = trim(ocn_data_dir)//'nitrate_WOA2005_surface_monthly' ! gx1 only
         sil_file = trim(ocn_data_dir)//'silicate_WOA2005_surface_monthly' ! gx1 only

         if (my_task == master_task .and. istep == 1) then
         if (trim(sil_data_type)=='clim') then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'silicate data interpolated to timestep:'
            write (nu_diag,*) trim(sil_file)
         endif
         if (trim(nit_data_type)=='clim') then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'nitrate data interpolated to timestep:'
            write (nu_diag,*) trim(nit_file)
            if (restore_bgc) write (nu_diag,*) &
              'bgc restoring timescale (days) =', trestore
         endif
         endif                     ! my_task, istep

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

         midmonth = 15          ! data is given on 15th of every month
!!!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

         ! Compute record numbers for surrounding months
         maxrec = 12
         ixm  = mod(month+maxrec-2,maxrec) + 1
         ixp  = mod(month,         maxrec) + 1
         if (mday >= midmonth) ixm = 99 ! other two points will be used
         if (mday <  midmonth) ixp = 99

         ! Determine whether interpolation will use values 1:2 or 2:3
         ! recslot = 2 means we use values 1:2, with the current value (2)
         !  in the second slot
         ! recslot = 1 means we use values 2:3, with the current value (2)
         !  in the first slot
         recslot = 1            ! latter half of month
         if (mday < midmonth) recslot = 2 ! first half of month

         ! Find interpolation coefficients
         call interp_coeff_monthly (recslot)

         readm = .false.
         if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

      endif   ! sil/nit_data_type


    !-------------------------------------------------------------------
    ! Read two monthly silicate values and interpolate.
    ! Restore toward interpolated value.
    !-------------------------------------------------------------------

      if (trim(nit_data_type)=='clim') then
         call read_clim_data (readm, 0, ixm, month, ixp, &
                              sil_file, sil_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (sil_data, sildat)

         if (restore_bgc) then
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sil(i,j,iblk) = sil(i,j,iblk)  &
                         + (sildat(i,j,iblk)-sil(i,j,iblk))*dt/trest
            enddo
            enddo
         enddo
         endif
      endif
    !-------------------------------------------------------------------
    ! Read two monthly nitrate values and interpolate.
    ! Restore toward interpolated value.
    !-------------------------------------------------------------------

      if (trim(nit_data_type)=='clim') then
         call read_clim_data (readm, 0, ixm, month, ixp, &
                              nit_file, nit_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (nit_data, nitdat)

         if (restore_bgc) then
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               nit(i,j,iblk) = nit(i,j,iblk)  &
                         + (nitdat(i,j,iblk)-nit(i,j,iblk))*dt/trest
            enddo
            enddo
         enddo
         endif
      endif

      end subroutine get_forcing_bgc

!=======================================================================
!BOP
!
! !ROUTINE: algal_dynamics
!
! !DESCRIPTION:
!
! Do biogeochemistry
!
! The philosophy in preliminary runs will be to build on
! the relatively well developed landfast geocycling mechanisms
! with adjustments as required in order to match sea ice data.
! Major starting points are thus Lavoie/Jin/Arrigo et al.
! N/C/Si and S cycling are all included.
! The latter is based almost entirely on the recent Stefels 07 review,
! Organism elemental ratios are rounded from Jin(Deal).
! An IARC/LANL EPSCoR proposal focusses on DMS.
! SE takes responsbility for development of the attached sulfur cycle.
! No models are currently available so that here the start point becomes
! the set of discussions included in observational papers.
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine algal_dynamics (nx_block, ny_block,   &
                                 icells,               &
                                 indxi,    indxj,      &
                                 hmix,     siln,       &
                                 nitn,     ammn,       &
                                 dmspn,    dmsn,       &
                                 flux_bion,            &
                                 meltb, congel,        &
                                 fswthrul, trcrn)
!
! !USES:
!
      use ice_calendar, only: dt
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         hmix   , & ! mixed layer depth
         meltb  , & ! bottom ice melt
         congel , & ! bottom ice growth 
         fswthrul    ! shortwave passing through ice to ocean

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         siln , & ! ocean silicate (mmol/m^3)
         ammn , & ! ocean ammonia/um (mmol/m^3)
         nitn , & ! ocean nitrate (mmol/m^3)
         dmspn, & ! ocean dmsp (mmol/m^3)
         dmsn     ! ocean dms  (mmol/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntraceb), &
         intent(inout) :: &
         flux_bion ! bgc fluxes to ocean
!
!  local parameters
!  sorted by layer -skeletal versus below ice supply
!
      real (kind=dbl_kind), parameter :: &
         sk_l       = 0.03_dbl_kind, &  ! skeletal layer thickness (m)
         pr_l       = 10.0_dbl_kind, &  ! product layer thickness (m) 
         R_C2N      = 7.0_dbl_kind , &  ! algal C to N (mole/mole)
         R_chl2N    = c1           , &  ! algal chlorophyll to N (mg/millimole)
         R_Si2N     = c1           , &  ! algal Si to N (mole/mole)
         R_S2N      = p1           , &  ! algal S to N (mole/mole)
         chlabs     = 0.03_dbl_kind, &  ! chlorophyll absorption (1/m(mg/m^3)) 
         mu_max     = 0.5_dbl_kind , &  ! cold, Eppley-like growth constant (1/day) 
         alpha      = p1           , &  ! P vs I slope in limit term (1/day/(W/m^2))
         P_max      = p1           , &  ! max production in limit term (1/day)
         K_Nit      = c1           , &  ! nitrate half saturation (mmol/m^3) 
         K_Am       = c1           , &  ! ammonia/um half saturation (mmol/m^3) 
         K_Sil      = c3           , &  ! silicon half saturation (mmol/m^3) 
         PV         = 1.e-6_dbl_kind,&  ! piston velocity for interface (m/s)
         N_congel   = c1           , &  ! required to trap initial conc (mmol/m^3)
         f_graze    = p1           , &  ! fraction grazed set constant per Arrigo
         f_graze_a  = 0.5_dbl_kind , &  ! fraction of grazing assimilated
         f_graze_s  = 0.5_dbl_kind , &  ! fraction of grazing spilled or slopped
         f_assim_e  = 0.5_dbl_kind , &  ! fraction of assimilation excreted         
         f_excrt_2S = c1           , &  ! excretion is efficient initially
         y_sk_DMS   = c1           , &  ! and conversion given high yield
         T_sk_conv=10800.0_dbl_kind, &  ! plus its fast (3 hours in seconds)
         T_sk_ox =864000.0_dbl_kind, &  ! DMS in turn oxidizes slowly (seconds)
         f_melt_2S  = p1           , &  ! melt by contrast entails sinkage losses
         y_pr_DMS   = c1           , &  ! but we begin again with unit yield
         T_pr_conv=10800.0_dbl_kind, &  ! and a similar conversion (seconds)
         T_pr_ox =864000.0_dbl_kind     ! plus round final time scale (seconds)

!
!  local variables
!
      integer (kind=int_kind) :: i, j, ij

      real (kind=dbl_kind) :: &
         N_sk_a     , &  ! skl algal nitrogen concentration on area (mmol/m^2)
         N_sk_v     , &  ! skl algal nitrogen concentration on volume (mmol/m^3) 
         C_sk_a     , &  ! skl algal carbon concentration on area (mmol/m^2)
         C_sk_v     , &  ! skl algal carbon concentration on volume (mmol/m^3)
         chl_sk_a   , &  ! skl algal chlorophyll concentration on area (mg/m^2)
         chl_sk_v   , &  ! skl algal chlorophyll concentration on volume (mg/m^3)
         Nit_sk_a   , &  ! skl nitrate concentration on area (mmol/m^2)
         Nit_sk_v   , &  ! skl nitrate concentration on volume (mmol/m^3) 
         Am_sk_a    , &  ! skl ammonia/um concentration on area (mmol/m^2)
         Am_sk_v    , &  ! skl ammonia/um concentration on volume (mmol/m^3) 
         Sil_sk_a   , &  ! skl silicon concentration on area (mmol/m^2)
         Sil_sk_v   , &  ! skl silicon concentration on volume (mmol/m^3) 
         DMSPp_sk_a , &  ! skl DMSPp concentration on area (mmol/m^2)
         DMSPp_sk_v , &  ! skl DMSPp concentration on volume (mmol/m^3)
         DMSPd_sk_a , &  ! skl DMSPd concentration on area (mmol/m^2)
         DMSPd_sk_v , &  ! skl DMSPd concentration on volume (mmol/m^3)
         DMS_sk_a   , &  ! skl DMS concentration on area (mmol/m^2)
         DMS_sk_v        ! skl DMS concentration on volume (mmol/m^3)

      real (kind=dbl_kind) :: &
         so_l     , &  ! source layer thickness (m)  
         op_dep   , &  ! attenuation exponent (optical depth)
         Iavg          ! attenuated Fswthrul (W/m^2)

      real (kind=dbl_kind) :: &
         L_lim    , &  ! overall light limitation 
         Nit_lim  , &  ! overall nitrate limitation
         Am_lim   , &  ! overall ammonia/um limitation
         N_lim    , &  ! overall nitrogen species limitation
         Sil_lim  , &  ! overall silicon limitation
         growmax_N, &  ! maximum growth rate in N currency (mol/m^2/s)
         grow_N   , &  ! true growth rate in N currency (mmol/m^2/s)
         potU_Nit , &  ! potential nitrate uptake (mmol/m^2/s)
         potU_Am  , &  ! potential ammonia/um uptake (mmol/m^2/s)
         potU_Sil , &  ! potential silicon uptake (mmol/m^2/s)
         U_Nit    , &  ! actual nitrate uptake (mmol/m^2/s)
         U_Am     , &  ! actual ammonia/um uptake (mmol/m^2/s)
         U_Sil         ! actual silicon uptake (mmol/m^2/s)

      real (kind=dbl_kind) :: &
         F_Nit    , &  ! flux out of the sk layer (mmol/m^2/s)
         F_Am     , &  ! flux out of the sk layer (mmol/m^2/s)
         F_Sil    , &  ! flux out of the sk layer (mmol/m^2/s)
         F_DMSP   , &  ! flux out of the sk layer (mmol/m^2/s)
         F_DMS    , &  ! flux out of the sk layer (mmol/m^2/s)
         cong_alg , &  ! algae captured with new ice growth (mmol/m^2)
         f_melt        ! vertical melt fraction of skeletal layer

!  source terms underscore s, removal underscore r

      real (kind=dbl_kind) :: &
         N_sk_s_p     , &  ! algal nitrogen photosynthesis (mmol/m^2)
         N_sk_s_c     , &  ! algal nitrogen via congel (mmol/m^2)
         N_sk_r_g     , &  ! algal nitrogen losses to grazing (mmol/m^2)
         N_sk_r_m     , &  ! algal nitrogen losses to melting (mmol/m^2) 
         N_sk_s       , &  ! net algal nitrogen sources (mmol/m^2)
         N_sk_r       , &  ! net algal nitrogen removal (mmol/m^2)
         C_sk_s       , &  ! net algal carbon sources (mmol/m^2)
         C_sk_r       , &  ! net algal carbon removal (mmol/m^2
         chl_sk_s     , &  ! net algal chlorophyll sources (mmol/m^2)
         chl_sk_r     , &  ! net algal chlorophyll removal (mmol/m^2)
         Nit_sk_f     , &  ! skl nitrate flux (mmol/m^2)
         Nit_sk_r_p   , &  ! skl nitrate uptake by algae (mmol/m^2)
         Nit_sk_s     , &  ! net skl nitrate sources (mmol/m^2)
         Nit_sk_r     , &  ! net skl nitrate removal (mmol/m^2)
         Am_sk_s_e    , &  ! skl ammonia/um source from excretion (mmol/m^2)
         Am_sk_f      , &  ! skl ammonia/um flux (mmol/m^2)
         Am_sk_r_p    , &  ! skl ammonia/um uptake by algae (mmol/m^2)
         Am_sk_s      , &  ! net skl ammonia/um sources (mmol/m^2)
         Am_sk_r      , &  ! net skl ammonia/um removal (mmol/m^2)
         Sil_sk_f     , &  ! skl silicon flux (mmol/m^2)
         Sil_sk_r_p   , &  ! skl silicon uptake by algae (mmol/m^2)
         Sil_sk_s     , &  ! net skl silicon sources (mmol/m^2)
         Sil_sk_r     , &  ! net skl silicon removal (mmol/m^2)
         DMSPp_sk_r_g , &  ! algal DMSP losses to grazing (mmol/m^2)
         DMSPp_sk_r_m , &  ! algal DMSP losses to melting (mmol/m^2)
         DMSPp_sk_s   , &  ! net algal DMSP sources (mmol/m^2)
         DMSPp_sk_r   , &  ! net algal DMSP removal (mmol/m^2)
         DMSPd_sk_s_s , &  ! skl dissolved DMSP from spillage (mmol/m^2)
         DMSPd_sk_s_e , &  ! skl dissolved DMSP from excretion (mmol/m^2)
         DMSPd_sk_f   , &  ! skl dissolved DMSP flux (mmmol/m^2)
         DMSPd_sk_r_c , &  ! skl dissolved DMSP conversion (mmol/m^2)
         DMSPd_sk_s   , &  ! net skl dissolved DMSP sources (mmol/m^2)
         DMSPd_sk_r   , &  ! net skl dissolved DMSP removal (mmol/m^2)
         DMS_sk_s_c   , &  ! skl DMS source via conversion (mmol/m^2)
         DMS_sk_f     , &  ! skl DMS flux (mmol/m^2)
         DMS_sk_r_o   , &  ! skl DMS losses due to oxidation (mmol/m^2)
         DMS_sk_s     , &  ! net skl DMS sources (mmol/m^2)
         DMS_sk_r     , &  ! net skl DMS removal (mmol/m^2)
         DMSP_pr_s_m  , &  ! prod dissolved DMSP from melting (mmol/m^2)
         DMSP_pr_f    , &  ! prod dissolved DMSP flux (mmol/m^2)
         DMSP_pr_r_c  , &  ! prod dissolved DMSP conversion (mmol/m^2)
         DMSP_pr_s    , &  ! net prod dissolved DMSP sources (mmol/m^2)
         DMSP_pr_r    , &  ! net prod dissolved DMSP removal (mmol/m^2)
         DMS_pr_s_c   , &  ! prod DMS source via conversion (mmol/m^2)
         DMS_pr_f     , &  ! prod DMS flux (mmol/m^2)
         DMS_pr_r_o   , &  ! prod DMS losses due to oxidation (mmol/m^2)
         DMS_pr_s     , &  ! net prod DMS sources (mmol/m^2)
         DMS_pr_r          ! net prod DMS removal (mmol/m^2)

!
!EOP
!
!  begin building biogeochemistry terms

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         so_l = hmix(i,j)

!  map to local variable names, provide conversions to volume data

         N_sk_a    = trcrn(i,j,nt_bgc_N_sk)
         N_sk_v    = N_sk_a / sk_l
         C_sk_a    = trcrn(i,j,nt_bgc_C_sk)
         C_sk_v    = C_sk_a / sk_l
         chl_sk_a  = trcrn(i,j,nt_bgc_chl_sk)
         chl_sk_v  = chl_sk_a / sk_l
         Nit_sk_a  = trcrn(i,j,nt_bgc_Nit_sk)
         Nit_sk_v  = Nit_sk_a / sk_l
         Am_sk_a   = trcrn(i,j,nt_bgc_Am_sk)
         Am_sk_v   = Am_sk_a / sk_l
         Sil_sk_a  = trcrn(i,j,nt_bgc_Sil_sk)
         Sil_sk_v  = Sil_sk_a / sk_l
         DMSPp_sk_a = trcrn(i,j,nt_bgc_DMSPp_sk)
         DMSPp_sk_v = DMSPp_sk_a / sk_l
         DMSPd_sk_a = trcrn(i,j,nt_bgc_DMSPd_sk)
         DMSPd_sk_v = DMSPd_sk_a / sk_l
         DMS_sk_a   = trcrn(i,j,nt_bgc_DMS_sk)
         DMS_sk_v   = DMS_sk_a / sk_l

!  geographic term construction for skeletal layer

!  radiation factors

         chl_sk_a = R_chl2N * N_sk_a 
         op_dep   = chlabs * chl_sk_a
         if (op_dep > puny) then
           Iavg   = fswthrul(i,j) * (c1 - exp(-op_dep)) / op_dep
         else
           Iavg   = c0
         endif

!  compute growth rate and related -note time conversion 1/(86400s/d)

         L_lim    = min(alpha*Iavg,P_max) / P_max
         Nit_lim  = Nit_sk_v/(Nit_sk_v + K_Nit)
         Am_lim   = Am_sk_v/(Am_sk_v + K_Am) 
         Nit_lim  = min(Nit_lim, c1 - Am_lim)
         N_lim    = Am_lim + Nit_lim
         Sil_lim  = Sil_sk_v/(Sil_sk_v + K_Sil)
         growmax_N= mu_max * N_sk_a / secday
         grow_N   = min(L_lim,N_lim,Sil_lim) * growmax_N
         potU_Nit = Nit_lim          * growmax_N
         potU_Am  = Am_lim           * growmax_N 
         potU_Sil = Sil_lim * R_Si2N * growmax_N
         U_Am     = min(grow_N, potU_Am)
         U_Nit    = grow_N - U_Am
         U_Sil    = min(R_Si2N * grow_N, potU_Sil)  

!  Nutrient exchange between skeletal layer and ocean (mmol/m^2/s).
!  Also sulfur compounds move into a below-ice product bin.
!  In cice, all fluxes are positive downward.
!  Piston velocity already has units m per second.

         F_Nit   = PV * (Nit_sk_v - nitn(i,j))
         F_Am    = PV * (Am_sk_v  - ammn(i,j))
         F_Sil   = PV * (Sil_sk_v - siln(i,j))
         F_DMSP  = PV * (DMSPd_sk_v - dmspn(i,j))
         F_DMS   = PV * (DMS_sk_v - dmsn(i,j))

!  congelation additions, set to round value from Arrigo initially

         cong_alg = congel(i,j) * N_congel 

!  loss fractions -note, grazing is a fixed portion of growth

         f_melt = min(c1,max(meltb(i,j)/sk_l,c0))

!  Since the framework remains incomplete at this point,
!  it is assumed as a starting expedient that 
!  DMSP loss to melting results in 10% conversion to DMS
!  which is then given a ten day removal constant.
!  Grazing losses are channeled into rough spillage and assimilation
!  then following ammonia there is some recycling.

         N_sk_s_p = grow_N * dt
         N_sk_s_c = cong_alg
         N_sk_r_g = f_graze * N_sk_s_p
         N_sk_r_m = f_melt * N_sk_a
         N_sk_s = N_sk_s_p + N_sk_s_c
         N_sk_r = N_sk_r_g + N_sk_r_m

         C_sk_s = R_C2N * N_sk_s
         C_sk_r = R_C2N * N_sk_r

         chl_sk_s = R_chl2N * N_sk_s
         chl_sk_r = R_chl2N * N_sk_r

         Nit_sk_f   =-F_Nit * dt 
         Nit_sk_r_p = U_Nit * dt
         Nit_sk_s = Nit_sk_f
         Nit_sk_r = Nit_sk_r_p

         Am_sk_s_e = f_assim_e * f_graze_a * N_sk_r_g
         Am_sk_f   =-F_Am * dt 
         Am_sk_r_p = U_Am * dt
         Am_sk_s = Am_sk_s_e + Am_sk_f
         Am_sk_r = Am_sk_r_p

         Sil_sk_f   =-F_Sil * dt 
         Sil_sk_r_p = U_Sil * dt
         Sil_sk_s = Sil_sk_f
         Sil_sk_r = Sil_sk_r_p

         DMSPp_sk_r_g = R_S2N * N_sk_r_g
         DMSPp_sk_r_m = R_S2N * N_sk_r_m
         DMSPp_sk_s = R_S2N * N_sk_s
         DMSPp_sk_r = DMSPp_sk_r_g + DMSPp_sk_r_m

         DMSPd_sk_s_s = f_graze_s * R_S2N * N_sk_r_g
         DMSPd_sk_s_e = f_excrt_2S * f_assim_e * f_graze_a * R_S2N * N_sk_r_g
         DMSPd_sk_f   =-F_DMSP * dt
         DMSPd_sk_r_c = (c1/T_sk_conv) * (DMSPd_sk_v*sk_l) * dt
         DMSPd_sk_s = DMSPd_sk_s_s + DMSPd_sk_s_e + DMSPd_sk_f
         DMSPd_sk_r = DMSPd_sk_r_c     

         DMS_sk_s_c = y_sk_DMS * DMSPd_sk_r_c
         DMS_sk_f   =-F_DMS * dt
         DMS_sk_r_o = (c1/T_sk_ox) * (DMS_sk_v*sk_l) * dt
         DMS_sk_s = DMS_sk_s_c + DMS_sk_f
         DMS_sk_r = DMS_sk_r_o

         DMSP_pr_s_m = f_melt_2S * DMSPp_sk_r_m
         DMSP_pr_f   = F_DMSP * dt
         DMSP_pr_r_c = (c1/T_pr_conv) * (dmspn(i,j)*pr_l) * dt 
         DMSP_pr_s = DMSP_pr_s_m + DMSP_pr_f
         DMSP_pr_r = DMSP_pr_r_c

         DMS_pr_s_c = y_pr_DMS * DMSP_pr_r_c
         DMS_pr_f   = F_DMS * dt 
         DMS_pr_r_o = (c1/T_pr_ox) * (dmsn(i,j)*pr_l) * dt
         DMS_pr_s = DMS_pr_s_c + DMS_pr_f
         DMS_pr_r = DMS_pr_r_o

! place back onto grid

         trcrn(i,j,nt_bgc_N_sk)    = N_sk_a    + N_sk_s    - N_sk_r
         trcrn(i,j,nt_bgc_C_sk)    = C_sk_a    + C_sk_s    - C_sk_r
         trcrn(i,j,nt_bgc_chl_sk)  = chl_sk_a  + chl_sk_s  - chl_sk_r
         trcrn(i,j,nt_bgc_Nit_sk)  = Nit_sk_a  + Nit_sk_s  - Nit_sk_r
         trcrn(i,j,nt_bgc_Am_sk)   = Am_sk_a   + Am_sk_s   - Am_sk_r
         trcrn(i,j,nt_bgc_Sil_sk)  = Sil_sk_a  + Sil_sk_s  - Sil_sk_r
         trcrn(i,j,nt_bgc_DMSPp_sk)= DMSPp_sk_a+ DMSPp_sk_s- DMSPp_sk_r
         trcrn(i,j,nt_bgc_DMSPd_sk)= DMSPd_sk_a+ DMSPd_sk_s- DMSPd_sk_r
         trcrn(i,j,nt_bgc_DMS_sk)  = DMS_sk_a  + DMS_sk_s  - DMS_sk_r
         nitn (i,j)                = nitn (i,j)+ F_Nit *dt/so_l
         ammn (i,j)                = ammn (i,j)+ F_Am  *dt/so_l
         siln (i,j)                = siln (i,j)+ F_Sil *dt/so_l
         dmspn(i,j)                = dmspn(i,j)+(DMSP_pr_s - DMSP_pr_r)/pr_l
         dmsn (i,j)                = dmsn(i,j) +(DMS_pr_s  - DMS_pr_r )/pr_l 

         trcrn(i,j,nt_bgc_N_sk)    = max(trcrn(i,j,nt_bgc_N_sk)   , c0)
         trcrn(i,j,nt_bgc_C_sk)    = max(trcrn(i,j,nt_bgc_C_sk)   , c0)
         trcrn(i,j,nt_bgc_chl_sk)  = max(trcrn(i,j,nt_bgc_chl_sk) , c0)
         trcrn(i,j,nt_bgc_Nit_sk)  = max(trcrn(i,j,nt_bgc_Nit_sk) , c0)
         trcrn(i,j,nt_bgc_Am_sk)   = max(trcrn(i,j,nt_bgc_Am_sk)  , c0)
         trcrn(i,j,nt_bgc_Sil_sk)  = max(trcrn(i,j,nt_bgc_Sil_sk) , c0)
         trcrn(i,j,nt_bgc_DMSPp_sk)= max(trcrn(i,j,nt_bgc_DMSPp_sk), c0)
         trcrn(i,j,nt_bgc_DMSPd_sk)= max(trcrn(i,j,nt_bgc_DMSPd_sk), c0)
         trcrn(i,j,nt_bgc_DMS_sk)  = max(trcrn(i,j,nt_bgc_DMS_sk) , c0)
         nitn(i,j) = max(nitn(i,j), c0)
         ammn(i,j) = max(ammn(i,j), c0)
         siln(i,j) = max(siln(i,j), c0)
         dmspn(i,j)= max(dmspn(i,j),c0)
         dmsn(i,j) = max(dmsn(i,j), c0) 

         flux_bion(i,j,nlt_bgc_NO   ) = flux_bion(i,j,nlt_bgc_NO   ) + F_Nit
         flux_bion(i,j,nlt_bgc_NH   ) = flux_bion(i,j,nlt_bgc_NH   ) + F_Am
         flux_bion(i,j,nlt_bgc_Sil  ) = flux_bion(i,j,nlt_bgc_Sil  ) + F_Sil
         flux_bion(i,j,nlt_bgc_DMSPp) = flux_bion(i,j,nlt_bgc_DMSPp) + F_DMSP
         flux_bion(i,j,nlt_bgc_DMS  ) = flux_bion(i,j,nlt_bgc_DMS  ) + F_DMS

      enddo

      end subroutine algal_dynamics

!=======================================================================
!BOP
!
! !ROUTINE: algal_biochemistry
!
! !DESCRIPTION:
!
! Do biogeochemistry  without including ice-ocean fluxes, melt/growth losses,
!   or interior ice transport processes
!
! After Elliott:
!
! The philosophy in preliminary runs will be to build on
! the relatively well developed landfast geocycling mechanisms
! with adjustments as required in order to match sea ice data.
! Major starting points are thus Lavoie/Jin/Arrigo et al.
! N/C/Si and S cycling are all included.
! The latter is based almost entirely on the recent Stefels 07 review,
! Organism elemental ratios are rounded from Jin(Deal).
! An IARC/LANL EPSCoR proposal focusses on DMS.
! SE takes responsbility for development of the attached sulfur cycle.
! No models are currently available so that here the start point becomes
! the set of discussions included in observational papers.
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine algal_biochemistry (nx_block, ny_block,   &
                                 icells,               &
                                 indxi,    indxj,      &
                               !  hmix,     siln,       &
                               !  nitn,     ammn,       &
                               !  dmspn,    dmsn,       &
                               !  fsiln, fnitn, fammn,  &
                               !  fdmspn, fdmsn,        &
                               !  meltb, congel,        &
                                 fswthrul, trcrn)
!
! !USES:
!
      use ice_calendar, only: dt
!      use ice_zbgc_public, only: R_C2N, R_chl2N
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
       !  hmix   , & ! mixed layer depth
       !  meltb  , & ! bottom ice melt
       !  congel , & ! bottom ice growth 
         fswthrul    ! shortwave passing through ice to ocean

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn

  !    real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
  !       siln , & ! ocean silicate (mmol/m^3)
  !       ammn , & ! ocean ammonia/um (mmol/m^3)
  !       nitn , & ! ocean nitrate (mmol/m^3)
  !       dmspn, & ! ocean dmsp (mmol/m^3)
  !       dmsn,  & ! ocean dms  (mmol/m^3)
  !       fsiln, & ! ocean silicate flux (mmol/m^2/s)
  !       fammn, & ! ocean ammonia flux (mmol/m^2/s)
  !       fnitn, & ! ocean nitrate flux (mmol/m^2/s)
  !       fdmspn,& ! ocean dmsp flux (mmol/m^2/s)
  !       fdmsn    ! ocean dms flux (mmmol/m^2/s)
!
!  local parameters
!  sorted by layer -skeletal versus below ice supply
!
!
!  Parameters shared by both skeletal layer and zbgc models 

      real (kind=dbl_kind), parameter :: &
!njeffery (now in ice_zbgc_public with different values)
!         R_C2N      = 7.0_dbl_kind , &  ! algal C to N (mole/mole)
!njeffery (R_chl2N = 3)
!         R_chl2N    = c1           , &  ! algal chlorophyll to N (mg/millimole)
!njeffery (R_Si2N = 1.5)
!         R_Si2N     = c1           , &  ! algal Si to N (mole/mole)
!njeffery (R_S2N = 0.03)
!         R_S2N      = p1           , &  ! algal S to N (mole/mole)
         chlabs     = 0.03_dbl_kind, &  ! chlorophyll absorption (1/m(mg/m^3)) 
!njeffery  chlabs      = 9.e-4_dbl_kind, & ! chlorophyll absorption (1/(mg/m^3)) where the 3 cm 
         mu_max     = 0.5_dbl_kind , &  ! cold, Eppley-like growth constant (1/day)
         K_Nit      = c1           , &  ! nitrate half saturation (mmol/m^3) 
         K_Am       = c1           , &  ! ammonia/um half saturation (mmol/m^3) 
         K_Sil      = c3           , &  ! silicon half saturation (mmol/m^3) 
!njeffery K_Sil = 4.0_dbl_kind , & ! silicon half saturation (mmol/m^3) 
         f_graze    = p1           , &  ! fraction grazed set constant per Arrigo
         f_graze_s  = 0.5_dbl_kind , &  ! fraction of grazing spilled or slopped   
         f_excrt_2S = c1           , &  ! excretion is efficient initially
         y_sk_DMS   = c1           , &  ! and conversion given high yield
         T_sk_conv=10800.0_dbl_kind, &  ! plus its fast (3 hours in seconds)
         T_sk_ox =864000.0_dbl_kind     ! DMS in turn oxidizes slowly (seconds)
!njeffery         t_sk_conv  = 10.0_dbl_kind, &  ! at a Stefels rate (days)
!njeffery         t_sk_ox    = 10.0_dbl_kind     ! DMS in turn oxidizes slowly (days)

! Parameters in just skeletal layer
 
      real (kind=dbl_kind), parameter :: &
         sk_l       = 0.03_dbl_kind, &  ! skeletal layer thickness (m)
         pr_l       = 10.0_dbl_kind, &  ! product layer thickness (m)
         alpha      = p1           , &  ! P vs I slope in limit term (1/day/(W/m^2)) 
         P_max      = p1           , &  ! max production in limit term (1/day)
         f_graze_a  = 0.5_dbl_kind , &  ! fraction of grazing assimilated
         f_assim_e  = 0.5_dbl_kind      ! fraction of assimilation excreted      

! Parameters in just zbgc model

      real (kind=dbl_kind), parameter :: &
         T_bot        = -1.8_dbl_kind, & ! interface close to freezing (C)    
         T_max        = -1.8_dbl_kind, & !   !maximum growth at Tmax
         alpha2max    = 0.8_dbl_kind, &   !
         beta2max     = 0.018_dbl_kind,&  ! corresponding light inhibition (1/W/m^2)
         grow_Tdep    = 0.0633_dbl_kind,& ! and its T dependence (1/C)
         fr_graze_tot = 0.25_dbl_kind, & !combine fr_graze_a*fr_assim_e
         inhib        =-1.46_dbl_kind, &  ! inhibition by ammonia (per conc)      
         mort_pre   = 0.0208_dbl_kind, &  ! mort_pre'*exp(-mort_Tdep*T_max) = 0.022 = mort_pre_old
                                          !prefix to mortality, will remin (1/day)
         mort_Tdep  = 0.03_dbl_kind , &  ! and its T dependence (1/C) 
         fr_mort2min= c1           , &  ! plus fractionation to remin
         t_nitrif   = 67.0_dbl_kind     ! nitrification time scale (days)    


! Parameters for transport or product layer
 
 !     real (kind=dbl_kind), parameter :: &
      !   PV         = 1.e-6_dbl_kind,&  ! piston velocity for interface (m/s)
      !   N_congel   = c1           , &  ! required to trap initial conc (mmol/m^3)
      !   f_melt_2S  = p1           , &  ! melt by contrast entails sinkage losses
      !   y_pr_DMS   = c1           , &  ! but we begin again with unit yield
      !   T_pr_conv=10800.0_dbl_kind, &  ! and a similar conversion (seconds)
      !   T_pr_ox =864000.0_dbl_kind     ! plus round final time scale (seconds)
!njeffery         t_pr_conv  = 10.0_dbl_kind, &  ! and a similar conversion (days)
!njeffery         t_pr_ox    = 10.0_dbl_kind     ! plus round final time scale (days)

!
!
!  local variables
!
      integer (kind=int_kind) :: i, j, ij

      real (kind=dbl_kind) :: &
         N_sk_a     , &  ! skl algal nitrogen concentration on area (mmol/m^2)
         N_sk_v     , &  ! skl algal nitrogen concentration on volume (mmol/m^3) 
         C_sk_a     , &  ! skl algal carbon concentration on area (mmol/m^2)
         C_sk_v     , &  ! skl algal carbon concentration on volume (mmol/m^3)
         chl_sk_a   , &  ! skl algal chlorophyll concentration on area (mg/m^2)
         chl_sk_v   , &  ! skl algal chlorophyll concentration on volume (mg/m^3)
         Nit_sk_a   , &  ! skl nitrate concentration on area (mmol/m^2)
         Nit_sk_v   , &  ! skl nitrate concentration on volume (mmol/m^3) 
         Am_sk_a    , &  ! skl ammonia/um concentration on area (mmol/m^2)
         Am_sk_v    , &  ! skl ammonia/um concentration on volume (mmol/m^3) 
         Sil_sk_a   , &  ! skl silicon concentration on area (mmol/m^2)
         Sil_sk_v   , &  ! skl silicon concentration on volume (mmol/m^3) 
         DMSPp_sk_a , &  ! skl DMSPp concentration on area (mmol/m^2)
         DMSPp_sk_v , &  ! skl DMSPp concentration on volume (mmol/m^3)
         DMSPd_sk_a , &  ! skl DMSPd concentration on area (mmol/m^2)
         DMSPd_sk_v , &  ! skl DMSPd concentration on volume (mmol/m^3)
         DMS_sk_a   , &  ! skl DMS concentration on area (mmol/m^2)
         DMS_sk_v        ! skl DMS concentration on volume (mmol/m^3)

      real (kind=dbl_kind) :: &
         so_l     , &  ! source layer thickness (m)  
         op_dep   , &  ! attenuation exponent (optical depth)
         Iavg          ! attenuated Fswthrul (W/m^2)

      real (kind=dbl_kind) :: &
         L_lim    , &  ! overall light limitation 
         Nit_lim  , &  ! overall nitrate limitation
         Am_lim   , &  ! overall ammonia/um limitation
         N_lim    , &  ! overall nitrogen species limitation
         Sil_lim  , &  ! overall silicon limitation
         growmax_N, &  ! maximum growth rate in N currency (mol/m^2/s)
         grow_N   , &  ! true growth rate in N currency (mmol/m^2/s)
         potU_Nit , &  ! potential nitrate uptake (mmol/m^2/s)
         potU_Am  , &  ! potential ammonia/um uptake (mmol/m^2/s)
         potU_Sil , &  ! potential silicon uptake (mmol/m^2/s)
         U_Nit    , &  ! actual nitrate uptake (mmol/m^2/s)
         U_Am     , &  ! actual ammonia/um uptake (mmol/m^2/s)
         U_Sil         ! actual silicon uptake (mmol/m^2/s)

      real (kind=dbl_kind) :: &
         F_Nit    , &  ! flux out of the sk layer (mmol/m^2/s)
         F_Am     , &  ! flux out of the sk layer (mmol/m^2/s)
         F_Sil    , &  ! flux out of the sk layer (mmol/m^2/s)
         F_DMSP   , &  ! flux out of the sk layer (mmol/m^2/s)
         F_DMS    , &  ! flux out of the sk layer (mmol/m^2/s)
         cong_alg , &  ! algae captured with new ice growth (mmol/m^2)
         f_melt        ! vertical melt fraction of skeletal layer

!  source terms underscore s, removal underscore r

      real (kind=dbl_kind) :: &
         N_sk_s_p     , &  ! algal nitrogen photosynthesis (mmol/m^2)
         N_sk_s_c     , &  ! algal nitrogen via congel (mmol/m^2)
         N_sk_r_g     , &  ! algal nitrogen losses to grazing (mmol/m^2)
         N_sk_r_m     , &  ! algal nitrogen losses to melting (mmol/m^2) 
         N_sk_s       , &  ! net algal nitrogen sources (mmol/m^2)
         N_sk_r       , &  ! net algal nitrogen removal (mmol/m^2)
         C_sk_s       , &  ! net algal carbon sources (mmol/m^2)
         C_sk_r       , &  ! net algal carbon removal (mmol/m^2
         chl_sk_s     , &  ! net algal chlorophyll sources (mmol/m^2)
         chl_sk_r     , &  ! net algal chlorophyll removal (mmol/m^2)
         Nit_sk_f     , &  ! skl nitrate flux (mmol/m^2)
         Nit_sk_r_p   , &  ! skl nitrate uptake by algae (mmol/m^2)
         Nit_sk_s     , &  ! net skl nitrate sources (mmol/m^2)
         Nit_sk_r     , &  ! net skl nitrate removal (mmol/m^2)
         Am_sk_s_e    , &  ! skl ammonia/um source from excretion (mmol/m^2)
         Am_sk_f      , &  ! skl ammonia/um flux (mmol/m^2)
         Am_sk_r_p    , &  ! skl ammonia/um uptake by algae (mmol/m^2)
         Am_sk_s      , &  ! net skl ammonia/um sources (mmol/m^2)
         Am_sk_r      , &  ! net skl ammonia/um removal (mmol/m^2)
         Sil_sk_f     , &  ! skl silicon flux (mmol/m^2)
         Sil_sk_r_p   , &  ! skl silicon uptake by algae (mmol/m^2)
         Sil_sk_s     , &  ! net skl silicon sources (mmol/m^2)
         Sil_sk_r     , &  ! net skl silicon removal (mmol/m^2)
         DMSPp_sk_r_g , &  ! algal DMSP losses to grazing (mmol/m^2)
         DMSPp_sk_r_m , &  ! algal DMSP losses to melting (mmol/m^2)
         DMSPp_sk_s   , &  ! net algal DMSP sources (mmol/m^2)
         DMSPp_sk_r   , &  ! net algal DMSP removal (mmol/m^2)
         DMSPd_sk_s_s , &  ! skl dissolved DMSP from spillage (mmol/m^2)
         DMSPd_sk_s_e , &  ! skl dissolved DMSP from excretion (mmol/m^2)
         DMSPd_sk_f   , &  ! skl dissolved DMSP flux (mmmol/m^2)
         DMSPd_sk_r_c , &  ! skl dissolved DMSP conversion (mmol/m^2)
         DMSPd_sk_s   , &  ! net skl dissolved DMSP sources (mmol/m^2)
         DMSPd_sk_r   , &  ! net skl dissolved DMSP removal (mmol/m^2)
         DMS_sk_s_c   , &  ! skl DMS source via conversion (mmol/m^2)
         DMS_sk_f     , &  ! skl DMS flux (mmol/m^2)
         DMS_sk_r_o   , &  ! skl DMS losses due to oxidation (mmol/m^2)
         DMS_sk_s     , &  ! net skl DMS sources (mmol/m^2)
         DMS_sk_r     , &  ! net skl DMS removal (mmol/m^2)
         DMSP_pr_s_m  , &  ! prod dissolved DMSP from melting (mmol/m^2)
         DMSP_pr_f    , &  ! prod dissolved DMSP flux (mmol/m^2)
         DMSP_pr_r_c  , &  ! prod dissolved DMSP conversion (mmol/m^2)
         DMSP_pr_s    , &  ! net prod dissolved DMSP sources (mmol/m^2)
         DMSP_pr_r    , &  ! net prod dissolved DMSP removal (mmol/m^2)
         DMS_pr_s_c   , &  ! prod DMS source via conversion (mmol/m^2)
         DMS_pr_f     , &  ! prod DMS flux (mmol/m^2)
         DMS_pr_r_o   , &  ! prod DMS losses due to oxidation (mmol/m^2)
         DMS_pr_s     , &  ! net prod DMS sources (mmol/m^2)
         DMS_pr_r          ! net prod DMS removal (mmol/m^2)

!
!EOP
!
!  begin building biogeochemistry terms

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

     !    so_l = hmix(i,j)

!  map to local variable names, provide conversions to volume data

         N_sk_a    = trcrn(i,j,nt_bgc_N_sk)
         N_sk_v    = N_sk_a / sk_l
         C_sk_a    = trcrn(i,j,nt_bgc_C_sk)
         C_sk_v    = C_sk_a / sk_l
         chl_sk_a  = trcrn(i,j,nt_bgc_chl_sk)
         chl_sk_v  = chl_sk_a / sk_l
         Nit_sk_a  = trcrn(i,j,nt_bgc_Nit_sk)
         Nit_sk_v  = Nit_sk_a / sk_l
         Am_sk_a   = trcrn(i,j,nt_bgc_Am_sk)
         Am_sk_v   = Am_sk_a / sk_l
         Sil_sk_a  = trcrn(i,j,nt_bgc_Sil_sk)
         Sil_sk_v  = Sil_sk_a / sk_l
         DMSPp_sk_a = trcrn(i,j,nt_bgc_DMSPp_sk)
         DMSPp_sk_v = DMSPp_sk_a / sk_l
         DMSPd_sk_a = trcrn(i,j,nt_bgc_DMSPd_sk)
         DMSPd_sk_v = DMSPd_sk_a / sk_l
         DMS_sk_a   = trcrn(i,j,nt_bgc_DMS_sk)
         DMS_sk_v   = DMS_sk_a / sk_l

!  geographic term construction for skeletal layer

!  radiation factors

         chl_sk_a = R_chl2N * N_sk_a 
         op_dep   = chlabs * chl_sk_a
         if (op_dep > puny) then
           Iavg   = fswthrul(i,j) * (c1 - exp(-op_dep)) / op_dep
         else
           Iavg   = c0
         endif

!  compute growth rate and related -note time conversion 1/(86400s/d)

         L_lim    = min(alpha*Iavg,P_max) / P_max
         Nit_lim  = Nit_sk_v/(Nit_sk_v + K_Nit)
         Am_lim   = Am_sk_v/(Am_sk_v + K_Am) 
         Nit_lim  = min(Nit_lim, c1 - Am_lim)
         N_lim    = Am_lim + Nit_lim
         Sil_lim  = Sil_sk_v/(Sil_sk_v + K_Sil)
         growmax_N= mu_max * N_sk_a / secday
         grow_N   = min(L_lim,N_lim,Sil_lim) * growmax_N
         potU_Nit = Nit_lim          * growmax_N
         potU_Am  = Am_lim           * growmax_N 
         potU_Sil = Sil_lim * R_Si2N * growmax_N
         U_Am     = min(grow_N, potU_Am)
         U_Nit    = grow_N - U_Am
         U_Sil    = min(R_Si2N * grow_N, potU_Sil)  

!  Nutrient exchange between skeletal layer and ocean (mmol/m^2/s).
!  Also sulfur compounds move into a below-ice product bin.
!  In cice, all fluxes are positive downward.
!  Piston velocity already has units m per second.

   !      F_Nit   = PV * (Nit_sk_v - nitn(i,j))
   !      F_Am    = PV * (Am_sk_v  - ammn(i,j))
   !      F_Sil   = PV * (Sil_sk_v - siln(i,j))
   !      F_DMSP  = PV * (DMSPd_sk_v - dmspn(i,j))
   !      F_DMS   = PV * (DMS_sk_v - dmsn(i,j))

!  congelation additions, set to round value from Arrigo initially
           cong_alg = c0
!   !      cong_alg = congel(i,j) * N_congel 

!  loss fractions -note, grazing is a fixed portion of growth
         f_melt = c0
!         f_melt = min(c1,max(meltb(i,j)/sk_l,c0))
!
!  Since the framework remains incomplete at this point,
!  it is assumed as a starting expedient that 
!  DMSP loss to melting results in 10% conversion to DMS
!  which is then given a ten day removal constant.
!  Grazing losses are channeled into rough spillage and assimilation
!  then following ammonia there is some recycling.

         N_sk_s_p = grow_N * dt
         N_sk_s_c = cong_alg
         N_sk_r_g = f_graze * N_sk_s_p
         N_sk_r_m = c0 !f_melt * N_sk_a
         N_sk_s = N_sk_s_p + N_sk_s_c
         N_sk_r = N_sk_r_g + N_sk_r_m

         C_sk_s = R_C2N * N_sk_s
         C_sk_r = R_C2N * N_sk_r

         chl_sk_s = R_chl2N * N_sk_s
         chl_sk_r = R_chl2N * N_sk_r

         Nit_sk_f   =-F_Nit * dt 
         Nit_sk_r_p = U_Nit * dt
         Nit_sk_s = Nit_sk_f
         Nit_sk_r = Nit_sk_r_p

         Am_sk_s_e = f_assim_e * f_graze_a * N_sk_r_g
         Am_sk_f   =-F_Am * dt 
         Am_sk_r_p = U_Am * dt
         Am_sk_s = Am_sk_s_e + Am_sk_f
         Am_sk_r = Am_sk_r_p

         Sil_sk_f   =-F_Sil * dt 
         Sil_sk_r_p = U_Sil * dt
         Sil_sk_s = Sil_sk_f
         Sil_sk_r = Sil_sk_r_p

         DMSPp_sk_r_g = R_S2N * N_sk_r_g
         DMSPp_sk_r_m = R_S2N * N_sk_r_m
         DMSPp_sk_s = R_S2N * N_sk_s
         DMSPp_sk_r = DMSPp_sk_r_g + DMSPp_sk_r_m

         DMSPd_sk_s_s = f_graze_s * R_S2N * N_sk_r_g
         DMSPd_sk_s_e = f_excrt_2S * f_assim_e * f_graze_a * R_S2N * N_sk_r_g
         DMSPd_sk_f   =-F_DMSP * dt
         DMSPd_sk_r_c = (c1/T_sk_conv) * (DMSPd_sk_v*sk_l) * dt
         DMSPd_sk_s = DMSPd_sk_s_s + DMSPd_sk_s_e + DMSPd_sk_f
         DMSPd_sk_r = DMSPd_sk_r_c     

         DMS_sk_s_c = y_sk_DMS * DMSPd_sk_r_c
         DMS_sk_f   =-F_DMS * dt
         DMS_sk_r_o = (c1/T_sk_ox) * (DMS_sk_v*sk_l) * dt
         DMS_sk_s = DMS_sk_s_c + DMS_sk_f
         DMS_sk_r = DMS_sk_r_o

         DMSP_pr_s_m = c0 !f_melt_2S * DMSPp_sk_r_m
         DMSP_pr_f   = F_DMSP * dt
         DMSP_pr_r_c = c0 !(c1/T_pr_conv) * (dmspn(i,j)*pr_l) * dt 
         DMSP_pr_s = DMSP_pr_s_m + DMSP_pr_f
         DMSP_pr_r = DMSP_pr_r_c

         DMS_pr_s_c = c0 !y_pr_DMS * DMSP_pr_r_c
         DMS_pr_f   = F_DMS * dt 
         DMS_pr_r_o = c0 !(c1/T_pr_ox) * (dmsn(i,j)*pr_l) * dt
         DMS_pr_s = DMS_pr_s_c + DMS_pr_f
         DMS_pr_r = DMS_pr_r_o

! place back onto grid

         trcrn(i,j,nt_bgc_N_sk)    = N_sk_a    + N_sk_s    - N_sk_r
         trcrn(i,j,nt_bgc_C_sk)    = C_sk_a    + C_sk_s    - C_sk_r
         trcrn(i,j,nt_bgc_chl_sk)  = chl_sk_a  + chl_sk_s  - chl_sk_r
         trcrn(i,j,nt_bgc_Nit_sk)  = Nit_sk_a  + Nit_sk_s  - Nit_sk_r
         trcrn(i,j,nt_bgc_Am_sk)   = Am_sk_a   + Am_sk_s   - Am_sk_r
         trcrn(i,j,nt_bgc_Sil_sk)  = Sil_sk_a  + Sil_sk_s  - Sil_sk_r
         trcrn(i,j,nt_bgc_DMSPp_sk)= DMSPp_sk_a+ DMSPp_sk_s- DMSPp_sk_r
         trcrn(i,j,nt_bgc_DMSPd_sk)= DMSPd_sk_a+ DMSPd_sk_s- DMSPd_sk_r
         trcrn(i,j,nt_bgc_DMS_sk)  = DMS_sk_a  + DMS_sk_s  - DMS_sk_r
        ! nitn (i,j)                = nitn (i,j)+ F_Nit *dt/so_l
        ! ammn (i,j)                = ammn (i,j)+ F_Am  *dt/so_l
        ! siln (i,j)                = siln (i,j)+ F_Sil *dt/so_l
        ! dmspn(i,j)                = dmspn(i,j)+(DMSP_pr_s - DMSP_pr_r)/pr_l
        ! dmsn (i,j)                = dmsn(i,j) +(DMS_pr_s  - DMS_pr_r )/pr_l 

        ! fnitn(i,j)                = fnitn(i,j) + F_Nit  ! for pop coupling
        ! fammn(i,j)                = fammn(i,j) + F_Am   ! for pop coupling
        ! fsiln(i,j)                = fsiln(i,j) + F_Sil  ! for pop coupling
        ! fdmspn(i,j)               = fdmspn(i,j)+ F_DMSP ! for pop coupling
        ! fdmsn(i,j)                = fdmsn(i,j) + F_DMS  ! for pop coupling

         trcrn(i,j,nt_bgc_N_sk)    = max(trcrn(i,j,nt_bgc_N_sk)   , c0)
         trcrn(i,j,nt_bgc_C_sk)    = max(trcrn(i,j,nt_bgc_C_sk)   , c0)
         trcrn(i,j,nt_bgc_chl_sk)  = max(trcrn(i,j,nt_bgc_chl_sk) , c0)
         trcrn(i,j,nt_bgc_Nit_sk)  = max(trcrn(i,j,nt_bgc_Nit_sk) , c0)
         trcrn(i,j,nt_bgc_Am_sk)   = max(trcrn(i,j,nt_bgc_Am_sk)  , c0)
         trcrn(i,j,nt_bgc_Sil_sk)  = max(trcrn(i,j,nt_bgc_Sil_sk) , c0)
         trcrn(i,j,nt_bgc_DMSPp_sk)= max(trcrn(i,j,nt_bgc_DMSPp_sk), c0)
         trcrn(i,j,nt_bgc_DMSPd_sk)= max(trcrn(i,j,nt_bgc_DMSPd_sk), c0)
         trcrn(i,j,nt_bgc_DMS_sk)  = max(trcrn(i,j,nt_bgc_DMS_sk) , c0)
        ! nitn(i,j) = max(nitn(i,j), c0)
        ! ammn(i,j) = max(ammn(i,j), c0)
        ! siln(i,j) = max(siln(i,j), c0)
        ! dmspn(i,j)= max(dmspn(i,j),c0)
        ! dmsn(i,j) = max(dmsn(i,j), c0) 

      enddo

      end subroutine algal_biochemistry

!=======================================================================
!BOP
!
! !ROUTINE: tracer_transport
!
! !DESCRIPTION:
!
! Solve the scalar vertical diffusion equation implicitly using 
! tridiag_solver. Calculate the diffusivity from temperature and salinity.
! 
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice with 
! dynamic salinity or the height ratio == hinS/vicen*aicen, where hinS is the 
! height of the brine surface relative to the bottom of the ice.  This volume fraction
! may be > 1 in which case there is brine above the ice surface (meltponds). 
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine tracer_transport (nx_block, ny_block,    &
                                  icells,   dt,         &
                                  indxii,   indxjj,    &   
                                  aicen,  vicen,  & 
                                  hice_old,     nitn,  & 
                                  flux_bio, flux_bio_g, &
                                  ammn, siln,  dmspn,  dmsn,  & 
                                  algalNn, zphin, iphin, & 
                                  trcrn,  iDin,   sss, &
                                  fswthrul, &
                                  growN, upNOn, upNHn,  &
                                  dh_top, dh_bot, dh_bot_chl, &
                                  zfswin,  n_cat, &
                                  first_ice,  TLAT, TLON,      &
                                  hbri, hbri_old, &
                                  bgrid, igrid, cgrid)     
!
! !USES:
!
      use ice_therm_shared, only: solve_Sin
      use ice_calendar,    only: istep1, time
!      use ice_zbgc_public, only: remap_layers_bgc_plus, min_salin, thin, Dm
!
! !INPUT/OUTPUT PARAMETERS:                                
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells              ! number of true cells with aicen > 0
                              
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxii, indxjj ! compressed indices for icells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt            ! time step 

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid           ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid          ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid           !CICE vertical coordinate   

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         sss         , & ! ocean salinity (ppt)
         TLAT        ,&
         TLON        ,&
         hbri         , &   ! brine height  (m)
         hbri_old      , &   ! brine height  (m)
         hice_old
 

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         dh_top       , &  !  brine change in top and bottom for diagnostics (m)
         dh_bot       , &      !
         dh_bot_chl  ! (bottom) change in hinS

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         first_ice            ! for first category ice only .true. initialized values should be used

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr+1), &
         intent(in) :: &
         fswthrul         ! Short wave radiation at each ice layer (W/m^2)  

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr + 1), &
         intent(out) :: &
         zfswin           ! Short wave flux on bio grid (W/m^2)  1 is surface point 
                         !2:nblyr+1 interior grid points
       
      real (kind=dbl_kind), dimension (nx_block,ny_block,ntraceb), intent(inout) :: &
         flux_bio, &   !ocean tracer flux (mmol/m^2/s)
         flux_bio_g    !ocean tracer flux from gravity drainage (and molecular diffusion)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &   
!change to  inout when updating ocean fields
         nitn        , & ! ocean nitrate (mmol/m^3) 
         ammn        , & ! ocean ammonium (mmol/m^3) 
         siln        , & ! ocean silicate (mmol/m^3) 
         algalNn     , & ! ocean algal nitrogen (mmol/m^3)
         dmspn       , & ! ocean DMSP (mmol/m^3)
         dmsn            ! ocean DMS (mmol/m^3)  

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr), intent(inout) :: &
         growN       , & ! algal growth rate at each layer (mmol/m^3/s)
         upNOn        , & ! algal nitrate uptake rate at each layer (mmol/m^3/s)
         upNHn            ! algal growth rate at each layer (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(in) :: &
         zphin            ! Porosity of layers

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(in) :: &
         iphin       , &  ! Porosity on the igrid   
         iDin            !  Diffusivity on the igrid   (1/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn           ! for salinity (ppt * vicen) and NO (mmol/m^3 * vicen ) 

      integer (kind=int_kind), intent(in) :: &
         n_cat           ! category number (for use in test conservation)

      
!
!EOP
!

     
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k, m, mm     ! vertical biology layer index 
                         !
  
      real (kind=dbl_kind), dimension (icells,nblyr) :: &       
         tracer_loss_z   ! loss term (tracer that remains in ice during flushing)   

      real (kind=dbl_kind), dimension(icells) :: &
         hin         , & ! ice thickness (m)        
         hin_old     , & ! ice thickness before current melt/growth (m)
         hinS        , & ! brine thickness (m)
         hinS_old    , & ! brine thickness before current melt/growth (m)
         surface_S       ! salinity of ice above hin > hinS
  

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2) :: &
         zphin_N          ! zphi >= zphimin

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1) :: &
         iDin_N , &      ! (1/s) diffusivity for algae with no bottom flux
         iphin_N

      real (kind=dbl_kind), dimension(icells,nblyr_hist):: & 
         NOin        , & ! Local layer nitrate values (mmol/m^3)
         Nin         , & ! Local layer algal N values (mmol/m^3)
         Silin       , & ! Local layer silicate values (mmol/m^3)
         Cin         , & ! Local layer carbon values (mmol/m^3)
         chlin       , & ! Local layer chlorophyll values (mg/m^3)
         DMSPpin     , & ! Local layer DMSPp values (mmol/m^3)
         DMSPdin     , & ! Local layer DMSPd values (mmol/m^3)
         DMSin        , &! Local layer DMS values (mmol/m^3)
         NHin         , &! Local layer ammonium values (mmol/m^3)
         PONin           ! Local layer PON values (mmol/m^3)

      real (kind=dbl_kind), dimension(icells,nblyr_hist):: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix equation

      real (kind=dbl_kind), dimension(icells,nblyr):: &
         react_N  ,& ! sources and sinks for N bulk equation matrix 
         react_DMSPp ! sources and sinks for DMSPp bulk equation matrix

      real (kind=dbl_kind), dimension(icells,nblyr,ntraceb):: &
         react  , &   ! biological sources and sinks for equation matrix
         react_phi    ! react*zphi
 
      real (kind=dbl_kind), dimension(icells, nblyr_hist,ntraceb):: &
         in_init      , & ! Initial concentration from previous time-step
         biomat       , & ! Matrix output     
         dbio          ! biomat(ij,k,mm) - in_init(ij,k,mm)

      real (kind=dbl_kind), dimension(icells, nblyr, ntraceb):: &
         bulk_biomat      , & ! Final ice bulk concentration from previous time-step
         bulk_biomat_o        ! initial ice bulk concentration

      real (kind=dbl_kind), dimension(icells, ntraceb):: &
         bulk_top      , & ! top bulk concentration k =1 
         bulk_bot      , &  ! bot bulk concentration
         dC_dx_bot         ! bottom tracer gradient

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: &
         trtmp0       , &  ! temporary, remapped tracers     !need extra 
         trtmp            ! temporary, remapped tracers     ! for nblyr_hist = nblyr+2

      logical (kind=log_kind), dimension(icells,ntraceb) :: & 
         good_bio_numerics      !


! local parameters
                           
      !-----------------------------------------------------------------------------
      !  Data from Cottier et al, 1999, JGR for fitting zbgc physical parameters
      !-----------------------------------------------------------------------------
      integer, parameter :: &
         nt_zfswin = 1        ! for interpolation of short wave
    
    
   call ice_timer_start(timer_bgc2)  
 
  !------------------------------------------------------------------------ 
  !define gravity drainage diffusivity on bio grid 
  !if using tr_salinity. 
  !----------------------------------------------------------------------


  !-------------------------------------
  ! Initialize 
  !------------------------------------

        NOin(:,:) = c0
        Nin(:,:) = c0
        NHin(:,:) = c0
        Silin(:,:) = c0
        Cin(:,:) = c0
        chlin(:,:) = c0
        DMSPpin(:,:) = c0
        DMSPdin(:,:) = c0
        DMSin(:,:) = c0
        PONin(:,:) = c0
        iDin_N(:,:,:) = c0

        do k = 1, nblyr
           do ij = 1, icells   
              i = indxii(ij)
              j = indxjj(ij)  
              zphin_N(i,j,k+2) = max(zphimin,zphin(i,j,k+2))
              zphin_N(i,j,k+1) = max(zphimin,zphin(i,j,k+1))
              zphin_N(i,j,k) = max(zphimin,zphin(i,j,k))
              iphin_N(i,j,k) = max(zphimin,iphin(i,j,k))
              iphin_N(i,j,k+1) = max(zphimin,iphin(i,j,k+1))          
              iDin_N(i,j,k+1) = iphin_N(i,j,k+1)*Dm/hbri(i,j)**2
       

               if (tr_bgc_NO) then
                     if (.NOT. first_ice(i,j)) then  !brine concentration
                        NOin(ij,k+1) = trcrn(i,j,nt_bgc_NO+k-1)/zphin_N(i,j,k+1)

                        if (trcrn(i,j,nt_bgc_NO+k-1) > 1000.0_dbl_kind .AND. tr_bgc_N) then
                          write(nu_diag,*)'Bulk Nitrate solution error initially, first_ice(i,j):',first_ice(i,j)      
                          write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,mm:'&
                                           ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                            TLON(i,j)*rad_to_deg,istep1,mm
                          write(nu_diag,*)'h_bot_chl(i,j),dh_top(i,j)'&
                                           ,dh_bot_chl(i,j),dh_top(i,j) 
                          write(nu_diag,*)'ij,k,NOin(ij,k+1)'
                          write(nu_diag,*)ij,k,NOin(ij,k+1)
                          write(nu_diag,*)'k,trcrn(i,j,nt_bgc_NO+k-1),zphin_N(i,j,k+1)' 
                          write(nu_diag,*)k,trcrn(i,j,nt_bgc_NO+k-1),zphin_N(i,j,k+1)
                          write(nu_diag,*)'hinS(ij),hinS_old(ij)' 
                          write(nu_diag,*)hinS(ij),hinS_old(ij)
                          call abort_ice ('ice_bgc.F90: BGC error1')
                        endif

                    else
                       trcrn(i,j,nt_bgc_NO+k-1) = nitn(i,j)*initbio_frac   !initialize
                       NOin(ij,k+1) = trcrn(i,j,nt_bgc_NO+k-1)/zphin_N(i,j,k+1) 

                       if (trcrn(i,j,nt_bgc_NO+k-1) > 1000.0_dbl_kind .AND. tr_bgc_N) then
                         write(nu_diag,*)'Bulk Nitrate solution error initially, first_ice(i,j):',first_ice(i,j)      
                         write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,mm:'&
                                         ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                         TLON(i,j)*rad_to_deg,istep1,mm
                         write(nu_diag,*)'h_bot_chl(i,j),dh_top(i,j)'&
                                         ,dh_bot_chl(i,j),dh_top(i,j) 
                         write(nu_diag,*)'ij,k,NOin(ij,k+1)'
                         write(nu_diag,*)ij,k,NOin(ij,k+1)
                         write(nu_diag,*)'k,trcrn(i,j,nt_bgc_NO+k-1),zphin_N(i,j,k+1)' 
                         write(nu_diag,*)k,trcrn(i,j,nt_bgc_NO+k-1),zphin_N(i,j,k+1)
                         write(nu_diag,*)'hinS(ij),hinS_old(ij)' 
                         write(nu_diag,*)hinS(ij),hinS_old(ij)
                         call abort_ice ('ice_bgc.F90: BGC error1')
                       endif
   
                   endif
                endif
                if (tr_bgc_N) then
                     if (.NOT. first_ice(i,j) ) then 
                          Nin(ij,k+1) = trcrn(i,j,nt_bgc_N+k-1)/zphin_N(i,j,k+1)
                     else
                         trcrn(i,j,nt_bgc_N+k-1) = algalNn(i,j)  !*initbio_frac   
                         Nin(ij,k+1) = trcrn(i,j,nt_bgc_N+k-1)/zphin_N(i,j,k+1)  
                     endif      
                endif
                if (tr_bgc_NH) then
                     if (.NOT. first_ice(i,j) ) then 
                          NHin(ij,k+1) = trcrn(i,j,nt_bgc_NH+k-1)/zphin_N(i,j,k+1) 
                     else
                         trcrn(i,j,nt_bgc_NH+k-1) = ammn(i,j)*initbio_frac     
                         NHin(ij,k+1) = trcrn(i,j,nt_bgc_NH+k-1)/zphin_N(i,j,k+1)  
                     endif
                endif
                if (tr_bgc_Sil) then
                     if (.NOT. first_ice(i,j) ) then 
                          Silin(ij,k+1) = trcrn(i,j,nt_bgc_Sil+k-1)/zphin_N(i,j,k+1)  
                     else
                         trcrn(i,j,nt_bgc_Sil+k-1) = siln(i,j)*initbio_frac    
                         Silin(ij,k+1) = trcrn(i,j,nt_bgc_Sil+k-1)/zphin_N(i,j,k+1)  
                     endif
                endif
                if (tr_bgc_C) then
                     if (.NOT. first_ice(i,j) ) then 
                          Cin(ij,k+1) = trcrn(i,j,nt_bgc_C+k-1)/zphin_N(i,j,k+1)  
                     else
                         trcrn(i,j,nt_bgc_C+k-1) = R_C2N*algalNn(i,j) !*initbio_frac   
                         Cin(ij,k+1) = trcrn(i,j,nt_bgc_C+k-1)/zphin_N(i,j,k+1)  
                     endif
                endif
                if (tr_bgc_chl) then
                     if (.NOT. first_ice(i,j) ) then 
                          chlin(ij,k+1) = trcrn(i,j,nt_bgc_chl+k-1)/zphin_N(i,j,k+1) 
                     else
                         trcrn(i,j,nt_bgc_chl+k-1) = R_chl2N*algalNn(i,j)  !*initbio_frac    
                         chlin(ij,k+1) = trcrn(i,j,nt_bgc_chl+k-1)/zphin_N(i,j,k+1)  
                     endif 
                endif
                if (tr_bgc_DMSPp) then
                     if (.NOT. first_ice(i,j) ) then
                          DMSPpin(ij,k+1) = trcrn(i,j,nt_bgc_DMSPp+k-1)/zphin_N(i,j,k+1)  
                     else
                         trcrn(i,j,nt_bgc_DMSPp+k-1) = dmspn(i,j) !*initbio_frac   
                         DMSPpin(ij,k+1) = trcrn(i,j,nt_bgc_DMSPp+k-1)/zphin_N(i,j,k+1)  
                     endif
                endif
                if (tr_bgc_DMSPd) then
                     if (.NOT. first_ice(i,j) ) then 
                          DMSPdin(ij,k+1) = trcrn(i,j,nt_bgc_DMSPd+k-1)/zphin_N(i,j,k+1)  
                     else
                         trcrn(i,j,nt_bgc_DMSPd+k-1) = dmspn(i,j)*initbio_frac    
                         DMSPdin(ij,k+1) = trcrn(i,j,nt_bgc_DMSPd+k-1)/zphin_N(i,j,k+1)  
                     endif
                endif
                if (tr_bgc_DMS) then
                     if (.NOT. first_ice(i,j) ) then
                          DMSin(ij,k+1) = trcrn(i,j,nt_bgc_DMS+k-1)/zphin_N(i,j,k+1)   
                     else
                         trcrn(i,j,nt_bgc_DMS+k-1) = dmsn(i,j)*initbio_frac    
                         DMSin(ij,k+1) = trcrn(i,j,nt_bgc_DMS+k-1)/zphin_N(i,j,k+1)  
                     endif   
                endif
                if (tr_bgc_PON) then
                     if (.NOT. first_ice(i,j) ) then 
                          PONin(ij,k+1) = trcrn(i,j,nt_bgc_PON+k-1)/zphin_N(i,j,k+1)   
                     else
                         trcrn(i,j,nt_bgc_PON+k-1) =  nitn(i,j)*initbio_frac      
                         PONin(ij,k+1) = trcrn(i,j,nt_bgc_PON+k-1)/zphin_N(i,j,k+1)  
                     endif   
                endif
           
       enddo             !ij
       enddo          !k
  
      !-----------------------------------------------------------------
      !     boundary conditions
      !-----------------------------------------------------------------
       do ij = 1, icells
          
          i = indxii(ij)
          j = indxjj(ij)            
          
           surface_S(ij) = min_salin
           hinS(ij) = hbri(i,j)
           hin(ij) = vicen(i,j)/aicen(i,j)
           hin_old(ij) = hice_old(i,j)
           hinS_old(ij) = hbri_old(i,j)

           !---------------------------------
           ! diagnostic variable
           !----------------------------------
           
           !changed iphin_N to zphin_N

           NOin(ij,1)   = NOin(ij,2)
           Nin(ij,1)    = Nin(ij,2) 
           NHin(ij,1)   = NHin(ij,2)
           Silin(ij,1)  = Silin(ij,2)
           Cin(ij,1)    = Cin(ij,2)
           chlin(ij,1)  = chlin(ij,2) 
           DMSPpin(ij,1)= DMSPpin(ij,2)
           DMSPdin(ij,1)= DMSPdin(ij,2)
           DMSin(ij,1)  = DMSin(ij,2) 
           PONin(ij,1)  = PONin(ij,2)

           if (dh_top(i,j) < c0) then  ! .AND. sss(i,j) .GE. min_salin) then
                NOin(ij,1)   = min_salin/zphin_N(i,j,2)*nitn(i,j)/sss(i,j)
                Nin(ij,1)    = min_salin/zphin_N(i,j,2)*algalNn(i,j)/sss(i,j) 
                Cin(ij,1)    = R_C2N* Nin(ij,1) 
                chlin(ij,1)  = R_chl2N * Nin(ij,1)
                NHin(ij,1)   = min_salin/zphin_N(i,j,2)* ammn(i,j)/sss(i,j)  
                Silin(ij,1)  = min_salin/zphin_N(i,j,2)* siln(i,j)/sss(i,j) 
                DMSPpin(ij,1)= min_salin/zphin_N(i,j,2)*dmspn(i,j)/sss(i,j)
                DMSPdin(ij,1)= min_salin/zphin_N(i,j,2)*dmspn(i,j)/sss(i,j)
                DMSin(ij,1)  = min_salin/zphin_N(i,j,2)*dmsn(i,j)/sss(i,j)
                PONin(ij,1) =  NOin(ij,1) !min_salin/zphin_N(i,j,2)*algalNn(i,j)/sss(i,j)  !

            endif

        
            NOin(ij,nblyr_hist)   = nitn(i,j) 
            Nin(ij,nblyr_hist)    = algalNn(i,j)  
            NHin(ij,nblyr_hist)   = ammn(i,j) 
            Silin(ij,nblyr_hist)  = siln(i,j) 
            Cin(ij,nblyr_hist)    = R_C2N * Nin(ij,nblyr_hist)
            chlin(ij,nblyr_hist)  = R_chl2N * Nin(ij,nblyr_hist)
            DMSPpin(ij,nblyr_hist)= dmspn(i,j)
            DMSPdin(ij,nblyr_hist)= dmspn(i,j)
            DMSin(ij,nblyr_hist)  = dmsn(i,j)
            PONin(ij,nblyr_hist)  = nitn(i,j)   !algalNn(i,j) !c0


       enddo             !ij
   
      !-----------------------------------------------------------------
      ! Interpolate shortwave flux, fswthrul (defined at top to bottom with nilyr+1 
      !  evenly spaced  with spacing = (1/nilyr) to grid variable zfswin:
      !-----------------------------------------------------------------
      
        trtmp(:,:,:) = c0     
        do k = 1, nilyr
         do ij = 1, icells    
            i = indxii(ij)
            j = indxjj(ij) 
               
            ! contains cice values (fswthrul(i,j,1) is surface value)
             trtmp0(i,j,nt_zfswin+k-1) = fswthrul(i,j,k+1) 
          enddo
         enddo
        
        call remap_layers_bgc_plus (nx_block,ny_block,        &
                             indxii,   indxjj,           &
                             icells,                   &
                             ntrcr,                    &
                             nilyr,                    &
                             nt_zfswin,                  &
                             trtmp0,    trtmp,          &
                             0,        nblyr,          &
                             hin, hinS,         &
                             cgrid(2:nilyr+1),         &
                             bgrid(2:nblyr+1), surface_S )

       do k = 1,nblyr
       do ij = 1, icells  
           i = indxii(ij)
           j = indxjj(ij) 
          
           zfswin(i,j,k+1) = trtmp(i,j,nt_zfswin+k-1)
           zfswin(i,j,1) = fswthrul(i,j,1)
           
       enddo                    !ij
       enddo
  
      !-----------------------------------------------------------------
      ! Initialize Biology  or better, combine this with above
      !-----------------------------------------------------------------
 
      dbio(:,:,:) = c0

      do k = 1, nblyr_hist
         do ij = 1, icells    
            i = indxii(ij)
            j = indxjj(ij)!!
                           ! save initial values
           if (tr_bgc_NO)then             
               in_init(ij,k,nlt_bgc_NO) = NOin(ij,k)        
               biomat(ij,k,nlt_bgc_NO) = NOin(ij,k)     
            endif
            if (tr_bgc_NH) then
               in_init(ij,k,nlt_bgc_NH) = NHin(ij,k)      
               biomat(ij,k,nlt_bgc_NH) = NHin(ij,k)      
            endif
            if (tr_bgc_N) then
                in_init(ij,k,nlt_bgc_N)  = Nin(ij,k)      
                biomat(ij,k,nlt_bgc_N)  = Nin(ij,k)      
            endif
            if (tr_bgc_Sil) then
               in_init(ij,k,nlt_bgc_Sil)= Silin(ij,k)  
               biomat(ij,k,nlt_bgc_Sil)= Silin(ij,k)      
            endif
            if (tr_bgc_C) then
               in_init(ij,k,nlt_bgc_C)  = Cin(ij,k)  
               biomat(ij,k,nlt_bgc_C)  = Cin(ij,k)      
             endif
            if (tr_bgc_chl) then
                in_init(ij,k,nlt_bgc_chl)= chlin(ij,k)      
                biomat(ij,k,nlt_bgc_chl)= chlin(ij,k)      
            endif
            if (tr_bgc_DMSPp) then
                in_init(ij,k,nlt_bgc_DMSPp) = DMSPpin(ij,k)   
                biomat(ij,k,nlt_bgc_DMSPp) = DMSPpin(ij,k)      
            endif
            if (tr_bgc_DMSPd) then
                in_init(ij,k,nlt_bgc_DMSPd) = DMSPdin(ij,k) 
                biomat(ij,k,nlt_bgc_DMSPd) = DMSPdin(ij,k)      
            endif
            if (tr_bgc_DMS) then
               in_init(ij,k,nlt_bgc_DMS) = DMSin(ij,k)  
               biomat(ij,k,nlt_bgc_DMS) = DMSin(ij,k)           
            endif
            if (tr_bgc_PON) then
              in_init(ij,k,nlt_bgc_PON) = PONin(ij,k) 
              biomat(ij,k,nlt_bgc_PON) = PONin(ij,k) 
            endif
  
         enddo             !ij 
      enddo                !k
  

      !-----------------------------------------------------------------
      ! Compute biological reaction terms: react(1:nblyr,mm) 
      !        and algal growth rates: growN(i,j,1:nblyr)
      !        in_init(2:nblyr+1),     zfswin (2:nblyr+1), upNOn(1:nblyr)
      !----------------------------------------------------------------- 
         react(:,:,:) = c0   !reaction terms
         react_phi(:,:,:) = c0   !reaction terms
         growN(:,:,:) = c0   !algal growth rate
         upNOn (:,:,:) = c0   !algal NO uptake rate
         upNHn (:,:,:) = c0   !algal NH uptake rate

     call ice_timer_stop(timer_bgc2)

     call ice_timer_start(timer_bgc3)

    !  if (solve_bgc) then
    !    call zbiochemistry
    !  endif  

      !-----------------------------------------------------------------
      ! Compute elements of tridiagonal matrix for each tracer
      !
      ! Note, the Darcy velocity is defined as positive down
      !-----------------------------------------------------------------

      do mm = 1, ntraceb

       if (mm .NE. nlt_bgc_C .AND. mm .NE. nlt_bgc_chl ) then
    
        if (mm .NE. nlt_bgc_N .AND. mm .NE. nlt_bgc_DMSPp ) then !.AND. mm .NE. nlt_bgc_PON

            call get_matrix_elements_calc_bgc &
                                     (nx_block, ny_block,         &
                                      icells,                     &
                                      indxii,    indxjj,            &
                                      igrid, bgrid,              &
                                      dt,     zphin_N,               &
                                      iDin,   iphin_N,              &
                                      in_init(:,:,mm), hinS_old,   &
                                      hinS,                        &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs,              &
                                      react(:,:,mm),   &
                                      dh_top,dh_bot)
      

         !------------------------------------------------------------------------
         !       
         !  Use molecular diffusion only for DMSPp and N  and PON
         !------------------------------------------------------------------------

        else

            call get_matrix_elements_calc_bgc &
                                     (nx_block, ny_block,         &
                                      icells,                     &
                                      indxii,    indxjj,            &
                                      igrid, bgrid,              &
                                      dt,     zphin_N,               &
                                      iDin_N,   iphin_N,            &
                                      in_init(:,:,mm), hinS_old,   &
                                      hinS,                        &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs,              &
                                      react(:,:,mm), &
                                      dh_top,dh_bot_chl)

        endif
        
      !-----------------------------------------------------------------
      ! Solve tridiagonal matrix to obtain the new tracers
      !-----------------------------------------------------------------
         call tridiag_solver (icells,             &
                              nblyr_hist, sbdiag, &
                              diag,     spdiag,   &
                              rhs,      biomat(:,:,mm))

       endif  

       do ij = 1, icells  
         i = indxii(ij)
         j = indxjj(ij) 
         bulk_top(ij,mm) = biomat(ij,1,mm)*zphin_N(i,j,1)
         bulk_bot(ij,mm) = biomat(ij,nblyr_hist,mm)*zphin_N(i,j,nblyr+1)
         if (dh_bot(i,j) < c0) then
           bulk_bot(ij,mm) = biomat(ij,nblyr+1,mm)*zphin_N(i,j,nblyr+1)
         endif
         
         dC_dx_bot(ij,mm) = (biomat(ij,nblyr_hist,mm) - biomat(ij,nblyr+1,mm)) &
                           /(bgrid(nblyr_hist)-bgrid(nblyr+1))
       enddo  !ij
      enddo                !mm


      !-----------------------------------------------------------------
      ! Reload NOin from the matrix solution
      !-----------------------------------------------------------------
     


      do k = 2,nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells  
            
            i = indxii(ij)
            j = indxjj(ij) 
              do mm = 1, ntraceb
                 bulk_biomat(ij,k-1,mm) = biomat(ij,k,mm)*zphin_N(i,j,k)
                 bulk_biomat_o(ij,k-1,mm)= in_init(ij,k,mm)*zphin_N(i,j,k)
                 react_phi(ij,k-1,mm) = react(ij,k-1,mm)*zphin_N(i,j,k)
              enddo

              if (tr_bgc_NO) then
                 if ( bulk_biomat(ij,k-1,nlt_bgc_NO) > 1000.0_dbl_kind) then
               
                   write(nu_diag,*)'Bulk Nitrate solution error'      
                   write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,mm:'&
                                   ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                    TLON(i,j)*rad_to_deg,istep1,mm
                   write(nu_diag,*)'h_bot_chl(i,j),dh_top(i,j)'&
                                   ,dh_bot_chl(i,j),dh_top(i,j) 
                   write(nu_diag,*)'ij,k,nlt_bgc_NO'
                   write(nu_diag,*) ij,k,nlt_bgc_NO
                   write(nu_diag,*)'biomat(ij,mm,nlt_bgc_NO),mm = 1,nblyr+2'
                   write(nu_diag,*)(biomat(ij,mm,nlt_bgc_NO),mm = 1,nblyr+2)
                   write(nu_diag,*)'in_init(ij,k,nlt_bgc_NO),zphin_N(i,j,k),ntraceb,nbltrcr'
                   write(nu_diag,*) in_init(ij,k,nlt_bgc_NO),zphin_N(i,j,k),ntraceb,nbltrcr 
                   write(nu_diag,*)'iphin_N(i,j,mm),mm =1,nblyr+1' 
                   write(nu_diag,*)(iphin_N(i,j,mm),mm =1,nblyr+1)
                   write(nu_diag,*) 'zphin_N(i,j,mm),mm=1,nblyr+2'
                   write(nu_diag,*) (zphin_N(i,j,mm),mm=1,nblyr+2)

                   write(nu_diag,*)'hinS(ij),hinS_old(ij),hin_old(ij),hin(ij)' 
                   write(nu_diag,*)hinS(ij),hinS_old(ij),hin_old(ij),hin(ij)
                   call abort_ice ('ice_bgc.F90: BGC error')
                  
                  endif   
                 NOin(ij,k) = max(biomat(ij,k,nlt_bgc_NO),c0)
              endif
              if (tr_bgc_N)   Nin(ij,k) = max(biomat(ij,k,nlt_bgc_N),c0) 
              if (tr_bgc_NH)  NHin(ij,k) = max(biomat(ij,k,nlt_bgc_NH),c0)
              if (tr_bgc_Sil)  Silin(ij,k) = max(biomat(ij,k,nlt_bgc_Sil),c0)
              if (tr_bgc_DMSPp)  DMSPpin(ij,k) = max(biomat(ij,k,nlt_bgc_DMSPp),c0)
              if (tr_bgc_DMSPd)  DMSPdin(ij,k) = max(biomat(ij,k,nlt_bgc_DMSPd),c0) 
              if (tr_bgc_DMS)  DMSin(ij,k) = max(biomat(ij,k,nlt_bgc_DMS),c0)
              if (tr_bgc_PON)  PONin(ij,k) = max(biomat(ij,k,nlt_bgc_PON),c0)
         
         enddo       ! ij
      enddo          ! k


      !-----------------------------------------------------------------
      ! Define the fluxes into the ocean
      !-----------------------------------------------------------------        
        do  mm = 1, ntraceb
         if (mm == nlt_bgc_N .OR. mm == nlt_bgc_DMSPp) then

           do ij = 1, icells 
              i = indxii(ij)
              j = indxjj(ij)

 
              call  calc_bgc_fluxes &
                                     (hinS(ij), hinS_old(ij), bulk_biomat(ij,2:nblyr+1,mm), &
                                      bulk_biomat_o(ij,2:nblyr+1,mm), &
                                      bulk_top(ij,mm), bulk_bot(ij,mm), dh_top(i,j), dh_bot_chl(i,j), &
                                      iDin_N(i,j,nblyr+1), dC_dx_bot(ij,mm), dt, &
                                      igrid,good_bio_numerics(ij,mm),flux_bio(i,j,mm),flux_bio_g(i,j,mm),&
                                      react_phi(ij,:,mm),aicen(i,j)) 
     
              if (.NOT. good_bio_numerics(ij,mm)) then
                write(nu_diag,*)'1st Calc_bgc_fluxes: Category,i,j,ij,TLAT,TLON,'
                write(nu_diag,*)'istep1,mm, nlt_bgc_NO:'&
                 ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                TLON(i,j)*rad_to_deg,istep1,mm, nlt_bgc_NO
                write(nu_diag,*)'dh_bot_chl(i,j),dh_top(i,j)'&
                ,dh_bot_chl(i,j),dh_top(i,j)
                write(nu_diag,*)'dC_dx_bot(ij,mm),iDin_N(i,j,nblyr+1)'
                write(nu_diag,*)dC_dx_bot(ij,mm),iDin_N(i,j,nblyr+1)
                write(nu_diag,*)'biomat(ij,m,nlt_bgc_NO),m = 1,nblyr+2'
                write(nu_diag,*)(biomat(ij,m,nlt_bgc_NO),m = 1,nblyr+2)
                write(nu_diag,*)' in_init(ij,m,nlt_bgc_NO),m = 1,nblyr'
                write(nu_diag,*)(in_init(ij,m,nlt_bgc_NO),m = 1,nblyr)
                call abort_ice ('ice: calc_bgc_fluxes ')
              endif
              

           
           enddo  !ij

         else
         
            do ij = 1, icells 
              i = indxii(ij)
              j = indxjj(ij)
           

              call  calc_bgc_fluxes &
                                     (hinS(ij), hinS_old(ij), bulk_biomat(ij,:,mm), &
                                      bulk_biomat_o(ij,:,mm), &
                                      bulk_top(ij,mm), bulk_bot(ij,mm), dh_top(i,j), dh_bot(i,j), &
                                      iDin(i,j,nblyr+1),&
                                      dC_dx_bot(ij,mm), dt, &
                                      igrid,good_bio_numerics(ij,mm),flux_bio(i,j,mm),flux_bio_g(i,j,mm),&
                                      react_phi(ij,:,mm),aicen(i,j))            

              if (.NOT. good_bio_numerics(ij,mm)) then
                write(nu_diag,*)'2nd Calc_bgc_fluxes: Category,i,j,ij,TLAT,TLON'
                write(nu_diag,*)'istep1,mm, nlt_bgc_NO:'&
                 ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                TLON(i,j)*rad_to_deg,istep1,mm, nlt_bgc_NO
                write(nu_diag,*)'dh_bot(i,j),dh_top(i,j)'&
                ,dh_bot(i,j),dh_top(i,j)
                write(nu_diag,*)'dC_dx_bot(ij,mm),iDin(i,j,nblyr+1)'
                write(nu_diag,*)dC_dx_bot(ij,mm),iDin(i,j,nblyr+1)
                write(nu_diag,*)'biomat(ij,m,nlt_bgc_NO),m = 1,nblyr+2'
                write(nu_diag,*)(biomat(ij,m,nlt_bgc_NO),m = 1,nblyr+2)
                write(nu_diag,*)' in_init(ij,m,nlt_bgc_NO),m = 1,nblyr'
                write(nu_diag,*)(in_init(ij,m,nlt_bgc_NO),m = 1,nblyr)
                call abort_ice ('ice: calc_bgc_fluxes ')
              endif

            enddo  !ij


         endif 
        enddo   !mm

      !-----------------------------------------------------------------
      ! Update the ocean nutrient concentration
      !  Use flux_bio_g
      !-----------------------------------------------------------------
        
         !nitn(:,:) = c1*nitn(:,:) + c0        
         !ammn(:,:) = c1*ammn(:,:) + c0             
         !siln(:,:) = c1*siln(:,:) + c0             
         !algalNn(:,:) = c1*algalNn(:,:) + c0          
         !dmspn(:,:) = c1*dmspn(:,:) + c0            
         !dmsn(:,:) = c1*dmsn(:,:) + c0   


      !-----------------------------------------------------------------
      ! Update the tracer variable
      !-----------------------------------------------------------------
    

      do k = 1,nblyr                  !back to bulk quantity
         do ij = 1, icells 
            i = indxii(ij)
            j = indxjj(ij)
            if (hinS_old(ij) > thin .AND. hinS(ij) > thin) then
              if (tr_bgc_NO ) trcrn(i,j,nt_bgc_NO+k-1) =  NOin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_N)  trcrn(i,j,nt_bgc_N+k-1) =  Nin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_NH) trcrn(i,j,nt_bgc_NH+k-1) =  NHin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_C) trcrn(i,j,nt_bgc_C+k-1) =  R_C2N * Nin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_chl) trcrn(i,j,nt_bgc_chl+k-1) =  R_chl2N * Nin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_sil) trcrn(i,j,nt_bgc_Sil+k-1) =  Silin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_DMSPp) trcrn(i,j,nt_bgc_DMSPp+k-1) =  DMSPpin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_DMSPd) trcrn(i,j,nt_bgc_DMSPd+k-1) =  DMSPdin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_DMS) trcrn(i,j,nt_bgc_DMS+k-1) =  DMSin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_PON) trcrn(i,j,nt_bgc_PON+k-1) =  PONin(ij,k+1)*zphin_N(i,j,k+1)
            endif

         enddo           !ij
      enddo        !k   

      !-----------------------------------------------------------------
      !  couple back to ice_therm
      !-----------------------------------------------------------------
      
     
   
  !  endif   !ntraceb and solve_bgc



     call ice_timer_stop(timer_bgc3)

 
770 format (I6,D16.6)        
781 format (I6,I6,I6)
790 format (I6,I6)
791 format (f24.17)
792 format (2D16.6)
793 format (3D16.6)
794 format (4D15.5)
800 format (F10.4)

      end subroutine tracer_transport

!=======================================================================
!BOP
!
! !ROUTINE: algal_dyn
!
! !DESCRIPTION:
!
! Do biogeochemistry from subroutine algal_dynamics
! by Scott Elliott: updated to 175
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine algal_dyn (nx_block, ny_block,       &
                                 icells,              &
                                 indxi,    indxj,     &
                                 fswthrul, reactb,      &
                                 ltrcrn, ntr, growN,&
                                 upNOn, upNHn )
!
! !USES:
!
      use ice_calendar, only: dt
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells,  &            ! number of cells with aicen > puny
         ntr                   ! number of layer tracers

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny


      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fswthrul    ! average shortwave passing through to current level in  ice

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         growN,  &    !  algal growth rate   (mmol/m^3/s)
         upNOn ,  &    !  algal NO uptake  rate   (mmol/m^3/s)
         upNHn         !  algal NH uptake rate    (mmol/m^3/s)


      real (kind=dbl_kind), dimension(icells,ntr), &
         intent(out) :: &
         reactb      ! biological reaction terms

      real (kind=dbl_kind), dimension(icells,ntr), &
         intent(in) :: &
         ltrcrn       ! concentrations  in layer


!------------------------------------------------------
!
! Scott's new parameters
!------------------------------------------------------
      real (kind=dbl_kind), parameter :: & 
         sk_l       = 0.03_dbl_kind, &  ! skeletal layer thickness (m)
         pr_l       = 10.0_dbl_kind, &  ! product layer thickness (m) 
         T_bot      =-1.8_dbl_kind , &  ! interface close to freezing (C)

!----------- later will replace T_bot with layer temperature *****  
!         chlabs     = 0.03_dbl_kind, &  ! chlorophyll absorption (1/m(mg/m^3)) 
         chlabs      = 9.e-4_dbl_kind, & ! chlorophyll absorption (1/(mg/m^3)) where the 3 cm 
!                                         assumed skeletal layer thickness provides the conversion
!                                         of units
!njmod-------------   redefined so that coefficient out front is the maximum growth rate
         mu_max     = 0.5_dbl_kind, &  ! 0.5_dbl_kind, &  ! maximum growth rate (1/day)   !like Meibings     
         T_max      = -1.8_dbl_kind, & !   !maximum growth at Tmax
!njmod         grow_pre   = 1.68_dbl_kind, &  ! prefix to growth constant (1/day) 
                                        ! Note: this is so large because of the
                                        ! inhibition term 
         grow_Tdep  = 0.0633_dbl_kind,& ! and its T dependence (1/C)
        ! fr_resp    = 0.05_dbl_kind, &  ! respiration fraction  !moved to top of subroutine
         fr_graze   = p1           , &  ! A93 val for S, set to zero in Jin 06
         fr_graze_s = 0.5_dbl_kind , &  ! fraction of grazing spilled or slopped

!njmod------------------------ combine fr_graze_a*fr_assim_e    
!         fr_graze_a = 0.5_dbl_kind , &  ! fraction of grazing assimilated
!         fr_assim_e = 0.5_dbl_kind , &  ! fraction of assimilation excreted  
         fr_graze_tot = 0.25_dbl_kind, &    

         alpha2max  = 0.8_dbl_kind, &   !
         beta2max   = 0.018_dbl_kind,&  ! corresponding light inhibition (1/W/m^2)
         K_Nit      = c1           , &  ! nitrate half saturation (mmol/m^3) 
         K_Am       = c1           , &  ! ammonium half saturation (mmol/m^3) 
         K_Sil      = 4.0_dbl_kind , &  ! silicon half saturation (mmol/m^3) 
         inhib      =-1.46_dbl_kind, &  ! inhibition by ammonia (per conc)

!njmod------------------ pull out exp(mort_Tdep*T_max)  in the definition
!         mort_pre   = 0.022_dbl_kind, &  !mort_pre_old: prefix to mortality, will remin (1/day)     
         mort_pre   = 0.0208_dbl_kind, &  ! mort_pre'*exp(-mort_Tdep*T_max) = 0.022 = mort_pre_old
                                          !prefix to mortality, will remin (1/day)
         mort_Tdep  = 0.03_dbl_kind , &  ! and its T dependence (1/C) 
         fr_mort2min= c1           , &  ! plus fractionation to remin
         t_nitrif   = 67.0_dbl_kind     ! nitrification time scale (days)

      real (kind=dbl_kind), parameter :: &
         fr_excrt_2S= c1           , &  ! excretion is efficient initially
         y_sk_DMS   = c1           , &  ! and conversion given high yield
         t_sk_conv  = 10.0_dbl_kind, &  ! at a Stefels rate (days)
         t_sk_ox    = 10.0_dbl_kind     ! DMS in turn oxidizes slowly (days)

     ! real (kind=dbl_kind), parameter :: &
     !    chl_pr_v   = 0.1_dbl_kind , &  ! fixed nondiatom chlorophyll in ml (mg/m^3)
     !    R_chl2N_nd = 3.0_dbl_kind , &  ! shade adaptation below (mg/millimole)
     !    R_C2N_nd   = 7.0_dbl_kind , &  ! open water ratio (mole/mole)
     !    t_pr_dsrp  = 10.0_dbl_kind     ! disruption time scale (days)
         
     ! real (kind=dbl_kind), parameter :: &
     !    R_S2N_nd   = 0.03_dbl_kind, &  ! open water ratio nondiatoms (mole/mole)
     !    y_pr_DMS   = c1           , &  ! but we begin again with unit yield
     !    t_pr_conv  = 10.0_dbl_kind, &  ! and a similar conversion (days)
     !    t_pr_ox    = 10.0_dbl_kind     ! plus round final time scale (days)

!
!  local variables
!
      integer (kind=int_kind) :: i, j, ij

      real (kind=dbl_kind) :: &
         Nin        , &     ! algal nitrogen concentration on volume (mmol/m^3) 
         Cin        , &     ! algal carbon concentration on volume (mmol/m^3)
         chlin      , &     ! algal chlorophyll concentration on volume (mg/m^3)
         NOin       , &     ! nitrate concentration on volume (mmol/m^3) 
         NHin       , &     ! ammonia/um concentration on volume (mmol/m^3) 
         Silin      , &     ! silicon concentration on volume (mmol/m^3) 
         DMSPpin    , &     ! DMSPp concentration on volume (mmol/m^3)
         DMSPdin    , &     ! DMSPd concentration on volume (mmol/m^3)
         DMSin      , &     ! DMS concentration on volume (mmol/m^3)
         PONin              ! PON concentration on volume (mmol/m^3)

     
!-------------------------------------------------------------
!
! More of Scott's new local variables
!--------------------------------------------------------------

      real (kind=dbl_kind) :: &
         narrow     , &  ! spectral narrowing to PAR in some runs
         I_nar_MBJ  , &  ! narrowed intensity follows MBJ form (W/m^2) 
         op_dep     , &  ! bottom layer attenuation exponent (optical depth)
         Iavg_loc        ! bottom layer attenuated Fswthrul (W/m^2)

      real (kind=dbl_kind) :: &
         Sal_brine, &  ! brine salinity per MBJ (A93)(ppt weight)  
                       ! use this value from diffuse bio
         Sal_lim  , &  ! overall salinity limitation per MBJ (A93)
         L_lim    , &  ! overall light limitation 
         Nit_lim  , &  ! overall nitrate limitation
         Am_lim   , &  ! overall ammonium limitation
         N_lim    , &  ! overall nitrogen species limitation
         Sil_lim  , &  ! overall silicon limitation
         fr_Nit   , &  ! fraction of local ecological growth as nitrate
         fr_Am    , &  ! fraction of local ecological growth as ammonia
      !   mu_max   , &  ! maximum growth rate constant (1/s)
         growmax_N, &  ! maximum growth rate in N currency (mmol/m^2/s)
         grow_N_ecolim,   & ! add ecological limitations (mmol/m^2/s)
         U_Nit_ecolim,    & ! compute potential uptakes (mmol/m^2/s)
         U_Am_ecolim,     & ! compute potential uptakes (mmol/m^2/s)
         U_Sil_ecolim,    & ! compute potential uptakes (mmol/m^2/s)
         U_Nit_msslim,    & ! then numerical mass limits (mmol/m^2/s)
         U_Am_msslim,     & ! then numerical mass limits (mmol/m^2/s)
         U_Sil_msslim,    & ! then numerical mass limits (mmol/m^2/s)
         Umss2eco_Nit,    & ! and ratio of mass available to ecolimit
         Umss2eco_Am,     & ! and ratio of mass available to ecolimit
         Umss2eco_Sil,    & ! and ratio of mass available to ecolimit
         U_Nit_totlim,    & ! then select overall limit (mmol/m^2/s) 
         U_Am_totlim,     & ! then select overall limit (mmol/m^2/s) 
         U_Sil_totlim,    & ! then select overall limit (mmol/m^2/s) 
         grow_N_totNlim,  & ! sum over nitrogen species (mmol/m^2/s)
         grow_N_totSilim, & ! also account for silicate (mmol/m^2/s)
         grow_N   , &  ! true growth rate in N currency (mmol/m^2/s)
         U_Nit    , &  ! actual nitrate uptake (mmol/m^2/s)
         U_Am     , &  ! actual ammonium uptake (mmol/m^2/s)
         U_Sil         ! actual silicon uptake (mmol/m^2/s)

      real (kind=dbl_kind) :: &
         resp     , &  ! respiration (mmol/m^2/s)
         graze    , &  ! grazing (mmol/m^2/s)
         mort     , &  ! sum of mortality and excretion (mmol/m^2/s)
         nitrif        ! nitrification (mmol/m^2/s)


!  source terms underscore s, removal underscore r

!\\\ change units to mmol/m^3
      real (kind=dbl_kind) :: &
         N_s_p     , &  ! algal nitrogen photosynthesis (mmol/m^2)
         N_s_c     , &  ! algal nitrogen via congel (mmol/m^2)
         N_r_g     , &  ! algal nitrogen losses to grazing (mmol/m^2)
         N_r_r     , &  ! algal nitrogen losses to respiration (mmol/m^2)
         N_r_mo    , &  ! algal nitrogen losses to mortality (mmol/m^2)
         N_s       , &  ! net algal nitrogen sources (mmol/m^2)
         N_r       , &  ! net algal nitrogen removal (mmol/m^2)
         C_s       , &  ! net algal carbon sources (mmol/m^2)
         C_r       , &  ! net algal carbon removal (mmol/m^2)
         chl_s     , &  ! net algal chlorophyll sources (mmol/m^2)
         chl_r     , &  ! net algal chlorophyll removal (mmol/m^2)
         NO_s_n   , &  ! skl nitrate from nitrification (mmol/m^2)
         NO_s_r   , &  ! skl nitrate from respiration (mmol/m^2)
         NO_r_p   , &  ! skl nitrate uptake by algae (mmol/m^2)
         NO_s     , &  ! net skl nitrate sources (mmol/m^2)
         NO_r     , &  ! net skl nitrate removal (mmol/m^2)
         NH_s_e    , &  ! skl ammonium source from excretion (mmol/m^2)
         NH_s_r    , &  ! skl ammonium source from respiration (mmol/m^2)
         NH_s_mo   , &  ! skl ammonium source from mort/remin (mmol/m^2) 
         NH_r_p    , &  ! skl ammonium uptake by algae (mmol/m^2)
         NH_r_n    , &  ! skl ammonium removal to nitrification (mmol/m^2)
         NH_s      , &  ! net skl ammonium sources (mmol/m^2)
         NH_r      , &  ! net skl ammonium removal (mmol/m^2)
         Sil_s_r   , &  ! skl silicon from respiration (mmol/m^2)
         Sil_r_p   , &  ! skl silicon uptake by algae (mmol/m^2)
         Sil_s     , &  ! net skl silicon sources (mmol/m^2)
         Sil_r          ! net skl silicon removal (mmol/m^2)

      real (kind=dbl_kind) :: &
         DMSPp_s_p , &  ! algal DMSP to match photosynthesis (mmol/m^2)
         DMSPp_s_c , &  ! algal DMSP to match congel (mmol/m^2)
         DMSPp_r_g , &  ! algal DMSP losses to grazing (mmol/m^2)
         DMSPp_r_r , &  ! algal DMSP losses to respiration (mmol/m^2)
         DMSPp_r_mo, &  ! algal DMSP losses to mortality (mmol/m^2)
         DMSPp_s   , &  ! net algal DMSP sources (mmol/m^2)
         DMSPp_r   , &  ! net algal DMSP removal (mmol/m^2)
         DMSPd_s_s , &  ! skl dissolved DMSP from grazing spillage (mmol/m^2)
         DMSPd_s_e , &  ! skl dissolved DMSP from zooplankton excretion (mmol/m^2)
         DMSPd_s_mo, &  ! skl dissolved DMSP from MBJ algal mortexc (mmol/m^2)
         DMSPd_r_c , &  ! skl dissolved DMSP conversion (mmol/m^2)
         DMSPd_s   , &  ! net skl dissolved DMSP sources (mmol/m^2)
         DMSPd_r   , &  ! net skl dissolved DMSP removal (mmol/m^2)
         DMS_s_c   , &  ! skl DMS source via conversion (mmol/m^2)
         DMS_r_o   , &  ! skl DMS losses due to oxidation (mmol/m^2)
         DMS_s     , &  ! net skl DMS sources (mmol/m^2)
         DMS_r     , &  ! net skl DMS removal (mmol/m^2)
         PON_s_z   , &  ! PON source as zooplankton (mmol/m^3)
         PON_s_d   , &  ! PON source as detritus (mmol/m^3)
         PON_s     , &  ! net  PON sources (mmol/m^3)
         PON_r          ! net  PON removal (mmol/m^3)

      real (kind=dbl_kind) :: &
         DMSP_pr_s_nd , &  ! product layer dissolved DMSP from local bio (mmol/m^2)
         DMSP_pr_s_me , &  ! product layer dissolved DMSP from melting (mmol/m^2)
         DMSP_pr_r_c  , &  ! product layer dissolved DMSP conversion (mmol/m^2)
         DMSP_pr_s    , &  ! net product dissolved DMSP sources (mmol/m^2)
         DMSP_pr_r    , &  ! net product dissolved DMSP removal (mmol/m^2)
         DMS_pr_s_c   , &  ! product layer DMS source via conversion (mmol/m^2)
         DMS_pr_r_o   , &  ! product layer DMS losses due to oxidation (mmol/m^2)
         DMS_pr_s     , &  ! net product DMS sources (mmol/m^2)
         DMS_pr_r          ! net product DMS removal (mmol/m^2)


      logical (kind=log_kind) :: &   
         write_flag           ! set to true at each timestep        

!
!EOP
!
!  begin building biogeochemistry terms

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

      write_flag = .true.

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

!  map to local variable names, provide conversions to volume data
         Nin = c0
         Cin = c0
         chlin = c0
         NOin = c0
         NHin = c0
         Silin = c0
         DMSPpin = c0
         DMSPdin = c0
         DMSin = c0
         PONin = c0

         if (tr_bgc_N)         Nin       = ltrcrn(ij,nlt_bgc_N)
        ! if (tr_bgc_C)         Cin       = ltrcrn(ij,nlt_bgc_C)
        ! if (tr_bgc_chl) then
        !      chlin     = ltrcrn(ij,nlt_bgc_chl)
        ! else
              chlin = R_chl2N * Nin   !on volume
        ! endif 
         if (tr_bgc_NO)        NOin      = ltrcrn(ij,nlt_bgc_NO)
         if (tr_bgc_NH)        NHin      = ltrcrn(ij,nlt_bgc_NH)
         if (tr_bgc_Sil)       Silin     = ltrcrn(ij,nlt_bgc_Sil)
        ! if (tr_bgc_DMSPp)     DMSPpin   = ltrcrn(ij,nlt_bgc_DMSPp)
         if (tr_bgc_DMSPd)     DMSPdin   = ltrcrn(ij,nlt_bgc_DMSPd)
         if (tr_bgc_DMS)       DMSin     = ltrcrn(ij,nlt_bgc_DMS)
         if (tr_bgc_PON)       PONin     = ltrcrn(ij,nlt_bgc_PON)


!  radiation factors

        I_nar_MBJ = fswthrul(i,j)      

        op_dep   = chlabs * chlin   !! now chlabs is in 1/(mg/m^3)

!njmod-------------  this condition does not effect the L_lim term so 
! remove it.
!        if (op_dep > 0.01) then

!  Okhotsk maxima causes a reevaluation.
!  The new concept is, late algae at the bottom of the bottom strongly attenuated.
!  Since k about 0.03 1/m(mg/m3), efold at about 30 mg/m3
!  values of order hundreds will be shut down...
!  Newest algae experience super aborption above because they sit low.
!  More than perhaps two efolds and light falls below half value.


           Iavg_loc   = I_nar_MBJ * exp(-op_dep)
!         else
!           Iavg_loc   = I_nar_MBJ   
!         endif



!  compute growth rate and related -note time conversion 1/(86400s/d)
!  MBJ limitation forms maintained below
!  Not clear that the nitrogen term remains below unity, discuss and alter
!  Can be enforced through Gregg minimization as applied in early LANL versions
!  Also maintain Arrigo salinity restriction though bottom layers close to FP
!  This may be useful as the single ice bottom bin gives way to NJ brine convection
!  It appears MBJ inserted option for an L05 ice freeze/melt rate limit
!  but then left this turned off -here we mimic
!        Sal_brine =-3.9921 - 22.7*T_bot - 1.0015*T_bot**2 - 0.02*T_bot**3
!        Sal_lim   = 1.1e-2 + 3.012e-2*S_brine + 1.0342e-3*S_brine*S_brine &
!                   -4.6033e-5*S_brine**3 + 4.926e-7*S_brine**4 - 1.659e-9*S_brine**5
       Sal_lim   = c1 ! close to unity in any case, start simple


!njmod--------------  the inihibition term means that L_lim no longer peaks at 1
!       but we want all limitation terms to be [0,1]
!          
!       L_lim     = (c1 - exp(-alpha2max*Iavg_loc)) * exp(-beta2max*Iavg_loc)      

        L_lim     = (c1 - exp(-alpha2max*Iavg_loc))   ! no inhibition

!-----------------------------------------------------------------------
!njmod---------------  Again the limitation term, done this way, does not produce
!       N_lim  in [0,1].  Also, N_lim has local unphysical maxima.   
!    
!        Nit_lim   = (NOin/(NOin + K_Nit)) * exp(inhib*NHin)
!
 
        Nit_lim =   NOin/(NOin + K_Nit)
        Am_lim = c0
        if (tr_bgc_NH) then
            Am_lim    = NHin/(NHin + K_Am)
            N_lim     = (Nit_lim + Am_lim)/c2  !bounded by [0,1]
        else
            N_lim = Nit_lim
        endif
        Sil_lim = c1
        if (tr_bgc_Sil)  Sil_lim   = Silin/(Silin + K_Sil)  
    
!njmod------------ don't need ifthen
!  use inhibition here

         fr_Am    = c0  !p5
         if (tr_bgc_NH) then 
          fr_Am = p5
          if (Nit_lim*exp(inhib*NHin) + Am_lim > c0 ) &
            fr_Am  = Am_lim/(Nit_lim*exp(inhib*NHin) + Am_lim)   
         endif
        
         fr_Nit    = c1 - fr_Am 


!  Growth and uptake computed within the bottom layer 
!  Note here per A93 discussions and MBJ model, salinity is a universal restriction
!  Comparison with available column nutrients inserted but in tests had no effect
!  Primary production reverts to SE form, see MBJ below and be careful

!------- change units from Scott  to mmol/m^3/s

!      mu_max      = (grow_pre / secday)  * exp(grow_Tdep*T_bot) * Sal_lim
!                    ! pretty large.  something like 1.5 /day
! njmod ------------------------  No salinity limitation or light limitation
!              mu_max is already the maximum growth rated for T_bot <= T_max
!              change T_bot to grid T

       growmax_N   = mu_max / secday * exp(grow_Tdep * (T_bot - T_max))* Nin 
!-------------------------------------------------------------------------
!
!  Uptake should not exceed current values
!

        grow_N_ecolim = min(L_lim,N_lim,Sil_lim) * growmax_N
        U_Sil_ecolim = grow_N_ecolim* R_si2N
        if (tr_bgc_Sil) U_Sil_ecolim = min(grow_N_ecolim * R_Si2N, Silin/dt)
        U_Nit_ecolim = min(grow_N_ecolim * fr_Nit, NOin/dt)  
        U_Am_ecolim  = min(grow_N_ecolim * fr_Am , NHin/dt)    
 
        grow_N_ecolim = min(U_Sil_ecolim/R_Si2N,U_Nit_ecolim + U_Am_ecolim)

        if (tr_bgc_NH) then
          fr_Am = p5
          if (grow_N_ecolim > c0) fr_Am = min(U_Am_ecolim/grow_N_ecolim,c1)
        endif

        fr_Nit = c1-fr_Am

        U_Nit  = fr_Nit * grow_N_ecolim
        U_Am   = fr_Am  * grow_N_ecolim
        U_Sil  = R_Si2N * grow_N_ecolim


        grow_N_totSilim = U_Sil / R_Si2N  !don't need
        grow_N = grow_N_ecolim   

        
        resp        = fr_resp* grow_N
        graze       = fr_graze*grow_N

!-------------- now mort_pre is the maximum mortality rate for T_bot <= T_max)
!njmod  
        mort        = mort_pre * exp(mort_Tdep*(T_bot-T_max)) * Nin / secday 
        nitrif      = c0 !(NHin / t_nitrif) / secday

  
!------ history variable --------------

        growN(i,j) = grow_N 
        if (Nin > c0) growN(i,j) = grow_N/Nin  ! specific growth rate (per s)
         upNOn(i,j) = U_Nit
         upNHn(i,j) = U_Am

!  Since the framework remains incomplete at this point,
!  it is assumed as a starting expedient that 
!  DMSP loss to melting results in 10% conversion to DMS
!  which is then given a ten day removal constant.
!  Grazing losses are channeled into rough spillage and assimilation
!  then following ammonia there is some recycling.
!
!  Algal reaction term
!       N_react = (grow_N*(c1 -fr_resp - fr_graze) - mort)*dt

        N_s_p = grow_N * dt                  
        N_s_c = c0                           ! fall seeding off to start
        N_r_g = graze * dt                   ! primarily for S, MBJ leaves out
        N_r_r = resp * dt
        N_r_mo= mort * dt
        N_s = N_s_p
        N_r = N_r_g + N_r_r + N_r_mo 

!  Don't need these now
!       C_react = R_C2N * N_react 
! and
!       chl_react = R_chl2N * N_react
!
!        C_s = R_C2N * N_s
!        C_r = R_C2N * N_r
!
!         chl_s = R_chl2N * N_s
!         chl_r = R_chl2N * N_r

! Nitrate reaction term
!        NO_react = (nitrif - fr_Nit*grow_N)*dt

         NO_s_n = nitrif * dt                 
         NO_s_r = c0                         
         NO_r_p = U_Nit * dt                
         NO_s = NO_s_n + NO_s_r
         NO_r = NO_r_p

! Ammonium reaction term
!        NH_react = (-nitrif - (c1-fr_Nit - fr_resp - &
!           fr_graze*fr_graze_tot)*grow_N + mort*fr_mort2min)*dt  

!        NH_s_r = fr_Am * NH_r_r           ! likely frac 1 but could follow J06

         NH_s_r = c1 * N_r_r 
!njmod--------------  fr_assim_e*fr_graze_a always appear together,  redefine
!       fr_graze_tot = fr_assim_e*fr_graze_a
              
!         NH_s_e = fr_assim_e * fr_graze_a * N_r_g
         NH_s_e = fr_graze_tot * N_r_g
         NH_s_mo= fr_mort2min * N_r_mo
         NH_r_p = U_Am * dt                
         NH_r_n = nitrif *dt                  
         NH_s = NH_s_r + NH_s_e + NH_s_mo
         NH_r = NH_r_p + NH_r_n 

! Silica  reaction term
!        Sil_react = - R_Si2N * grow_N * dt

 !       Sil_s_r = R_Si2N * N_r_r          ! L05 say no but could follow J06
         Sil_s_r = c0                          
         Sil_r_p = U_Sil * dt
         Sil_s = Sil_s_r
         Sil_r = Sil_r_p 

!  Sulfur cycle begins here
!  Grazing losses are channeled into rough spillage and assimilation
!  then onward and the MBJ mortality channel is included
!  It is assumed as a starting expedient that 
!  DMSP loss to melting gives partial conversion to DMS in product layer
!  which then undergoes Stefels removal.

! DMSPp reaction term: Don't need
!        DMSPp_react = R_S2N * N_react
!
!         DMSPp_s_p = R_S2N * N_s_p
!         DMSPp_s_c = R_S2N * N_s_c 
!         DMSPp_r_g = R_S2N * N_r_g
!         DMSPp_r_r = R_S2N * N_r_r
!         DMSPp_r_mo= R_S2N * N_r_mo
!         DMSPp_s = DMSPp_s_p + DMSPp_s_c
!         DMSPp_r = DMSPp_r_g + DMSPp_r_r+DMSPp_r_mo

! DMSPd  reaction term
!        DMSPd_react = R_S2N*((fr_graze_s+fr_excrt_2S*fr_graze_tot)*fr_graze*grow_N + &
!                  fr_mort2min*mort)*dt - [\DMSPd]/t_sk_conv*dt

         DMSPd_s_s = fr_graze_s * R_S2N * N_r_g
!         DMSPd_s_e = fr_excrt_2S * fr_assim_e * fr_graze_a * R_S2N * N_r_g
         DMSPd_s_e = fr_excrt_2S * fr_graze_tot * R_S2N * N_r_g
         DMSPd_s_mo= fr_mort2min * R_S2N * N_r_mo
         DMSPd_r_c = (c1/t_sk_conv) * (c1/secday)  * (DMSPdin) * dt
         DMSPd_s = DMSPd_s_s + DMSPd_s_e + DMSPd_s_mo 
         DMSPd_r = DMSPd_r_c
!
! DMS reaction term
!       DMS_react = ([\DMSPd]*y_sk_DMS/t_sk_conv - c1/t_sk_ox *[\DMS])*dt


 
         DMS_s_c = y_sk_DMS * DMSPd_r_c 
         DMS_r_o = (c1/t_sk_ox) * (c1/secday) * (DMSin) * dt 
         DMS_s = DMS_s_c  
         DMS_r = DMS_r_o 

!  for mixed layer sulfur chemistry, fluxes kept separate for ice area weighting
!   units are different here, but not sure if they need to be changed
!   no fluxes into the product layer here
!
       !  DMSP_pr_s_nd= (c1/t_pr_dsrp) * (c1/secday) * &
       !                (chl_pr_v*pr_l) * (R_S2N_nd/R_chl2N_nd) * dt 
       !  DMSP_pr_s_me= fr_melt_2S * DMSPp_sk_r_me
       !  DMSP_pr_r_c = (c1/t_pr_conv)*(c1/secday) * (dmsp(i,j,iblk)*pr_l) * dt 
      !   DMSP_pr_f   = F_DMSP * dt
       !  DMSP_pr_s = DMSP_pr_s_nd                ! + DMSP_pr_s_me + DMSP_pr_f
       !  DMSP_pr_r = DMSP_pr_r_c

       !  DMS_pr_s_c = y_pr_DMS * DMSP_pr_r_c
       !  DMS_pr_r_o = (c1/t_pr_ox) * (c1/secday) * (dms(i,j,iblk)*pr_l) * dt
       !  DMS_pr_f   = F_DMS * dt 
       !  DMS_pr_s = DMS_pr_s_c                   ! + DMS_pr_f
       !  DMS_pr_r = DMS_pr_r_o



!\\\\\\ Added PON  ! needs to be updated  
!        currently using PON to shadow nitrate without reactions
        ! PON_s_z = N_r_g * f_graze_a 
        ! PON_s_d = N_r_g * (c1 - f_graze_a)
        ! PON_s   = PON_s_z + PON_s_d
        ! PON_r   = N_r_g * f_graze_a * f_assim_e
                  


! Define reaction terms

         if (tr_bgc_N) reactb(ij,nlt_bgc_N)    =  N_s    - N_r 
      !   if (tr_bgc_C) reactb(ij,nlt_bgc_C)    =  C_s    - C_r
      !   if (tr_bgc_chl)reactb(ij,nlt_bgc_chl)  =  chl_s  - chl_r
         if (tr_bgc_NO)reactb(ij,nlt_bgc_NO)   =  NO_s  - NO_r
         if (tr_bgc_NH)reactb(ij,nlt_bgc_NH)   =  NH_s   - NH_r 

         if (tr_bgc_Sil)reactb(ij,nlt_bgc_Sil)  =  Sil_s  - Sil_r

       !  if (tr_bgc_DMSPp)reactb(ij,nlt_bgc_DMSPp)=  DMSPp_s- DMSPp_r
         if (tr_bgc_DMSPd)reactb(ij,nlt_bgc_DMSPd)=  DMSPd_s- DMSPd_r

         if (tr_bgc_DMS)reactb(ij,nlt_bgc_DMS)  =  DMS_s  - DMS_r

         if (tr_bgc_PON)reactb(ij,nlt_bgc_PON)  =  c0 !  
! - (reactb(ij,nlt_bgc_N) + reactb(ij,nlt_bgc_NO) + reactb(ij,nlt_bgc_NH))


        ! if (growN(i,j) > 1.0e-6_dbl_kind .AND. write_flag) then
        !        write(nu_diag, *) 'reactb(ij,nlt_bgc_N), ij:',reactb(ij,nlt_bgc_N),ij
        !        write(nu_diag, *) 'growN(i,j), i,j:',growN(i,j),i,j
        !        write(nu_diag, *) 'fswthrul(ij):',fswthrul(ij)
        !        write(nu_diag, *) 'Nin,NOin,NHin,Silin:',Nin,NOin,NHin,Silin
        !        write_flag = .false.
        ! endif



      enddo


      end subroutine algal_dyn
!
!=======================================================================
!BOP
!
! !ROUTINE: get_matrix_elements for bgc equation - compute tridiagonal matrix elements
!
! !DESCRIPTION:
!
! Compute terms in tridiagonal matrix that will be solved to find
!  the bgc vertical profile
!
! !REVISION HISTORY:
!
! authors     Nicole Jeffery, LANL
!
!
! November 2008 by N. Jeffery, modified for bgc layer model
!
! !INTERFACE:
!
      subroutine get_matrix_elements_calc_bgc &
                                     (nx_block, ny_block,         &
                                      icells,                     &
                                      indxi,   indxj,             &
                                      igrid, bgrid,              &
                                      dt,     zphin_N,               &
                                      iDin,    iphin_N,             &
                                      NOin_init, hin_old,         &
                                      hin,                        &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs,              &
                                      reactb, dht,dhb)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj  ! compressed indices 

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(in) :: &
         zphin_N             ! Porosity of layers

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hin         , & ! ice thickness (m)                  (m)
         hin_old      ! ice thickness before growth/melt (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(in) :: &
         iphin_N           ! Porosity of layers on interface

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dht, dhb           !

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(in) :: &
         iDin              ! Diffusivity on the igrid (m^2/s)

      real (kind=dbl_kind), dimension (nblyr_hist), intent(in) :: &
         bgrid            ! biology nondimensional grid layer points 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid         ! biology grid interface points 

      real (kind=dbl_kind), dimension (icells,nblyr), intent(in) :: &
         reactb          ! Reaction terms from  algal_dyn

      real (kind=dbl_kind), dimension (icells,nblyr_hist), &
         intent(in) :: &
         NOin_init         ! bgc concentration at beginning of time step


      real (kind=dbl_kind), dimension (icells,nblyr_hist), &
         intent(out) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

 
!
!EOP
!

      real (kind=dbl_kind), dimension (icells,nblyr) :: &
         vel              ! advective velocity times dt (m)

     

      real (kind=dbl_kind), dimension (icells,nblyr+1) :: &
         ivel   

   
      real (kind=dbl_kind), dimension (nblyr_hist) :: &
         bgrid_temp   ! biology nondimensional grid layer points 


      real (kind=dbl_kind) :: &
         Crank           , &! Crank-Nicolson parameter  where
                         !C = 1 is fully implicit advection
                         !and C = 0 is fully explicit advection
         rhs_sp, rhs_sb,&! components of the rhs tri-diagonal matrix:
                         ! sp is super-diagonal, sb is sub-diagonal 
         rhs_diag        ! and diag is diagonal
 
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij,           & ! horizontal indices, combine i and j loops
         k, ks, ki, kr   ! vertical indices and row counters

      real (kind=dbl_kind), dimension (icells):: &
         dvel_dx        ! derivative of the velocity wrt x
      
      
      !-----------------------------------------------------------------
      ! Initialize matrix elements.
      !-----------------------------------------------------------------

      do k = 1, nblyr_hist

         if (k == 1) then
           do ij = 1, icells
             sbdiag(ij,k) = c0
             diag  (ij,k) = c1
             spdiag(ij,k) = c0
             rhs   (ij,k) = NOin_init(ij,k)  
           enddo
          else
           do ij = 1, icells
             sbdiag(ij,k) = c0
             diag  (ij,k) = c1
             spdiag(ij,k) = c0
             rhs   (ij,k) = NOin_init(ij,k)  
           enddo
         endif
         bgrid_temp(k) = bgrid(k)
      enddo

            
      Crank =  c1           ! C = c1 is fully implicit advection
                          ! C = c0 is fully explicit advection
      !-----------------------------------------------------------------
      ! Compute matrix elements
      !-----------------------------------------------------------------



802 format (5d15.5)


       

       do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (hin_old(ij) > thin) then

           bgrid_temp(1) =  c0 - grid_o_t/hin_old(ij)        
           bgrid_temp(nblyr_hist) = c1 + grid_o/hin_old(ij) 
           
             dvel_dx(ij) = (dhb(i,j)-dht(i,j))/hin_old(ij) 

             do k = 1, nblyr+1
               ivel(ij,k) =  (igrid(k)*dhb(i,j) - (igrid(k)-c1)*(dht(i,j)))/hin_old(ij)  
               if (k < nblyr+1) then
                  vel(ij,k) = (bgrid_temp(k+1)*(dhb(i,j)) - &
                          (bgrid_temp(k+1) - c1)* (dht(i,j)) )/hin_old(ij)  
               endif
             enddo    !k


           do k = 2, nblyr+1


           !=====================================================================================
           ! Let's model advection term at boundaries a little differently:
           ! wd(zphi*N)/dx = d(w*zphi*N)/dx - dw/dx * zphi * N
           !====================================================================================


             if (k == 2) then
               
                 if (dht(i,j) > c0) then   

                    sbdiag(ij,2) =  c0 
                                    
                                   
    

                    diag(ij,2)   = c1 + Crank * dvel_dx(ij) + dt/zphin_N(i,j,k)/(igrid(k) - &
                              igrid(k-1))* (iDin(i,j,k)/(bgrid_temp(k+1)- &
                              bgrid_temp(k))) - Crank * ((ivel(ij,2)*iphin_N(i,j,2)) *&
                               (bgrid_temp(3)-igrid(2))/(bgrid_temp(3)-bgrid_temp(2))- (ivel(ij,1)*zphin_N(i,j,2)))/&
                              (igrid(2)-igrid(1))/zphin_N(i,j,2)        



                    spdiag(ij,2) = - Crank *(ivel(ij,2)*iphin_N(i,j,2))/&
                             zphin_N(i,j,k)/(bgrid_temp(k+1) - &   
                             bgrid_temp(2)) * (igrid(2) - bgrid_temp(2))/(igrid(2)-igrid(1)) -  &
                             dt *iDin(i,j,k)/zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1))/(bgrid_temp(k+1) - bgrid_temp(k))
                 
                    rhs_sb      = c0
                    rhs_sp      =  (c1-Crank)*(ivel(ij,2)*iphin_N(i,j,2))/zphin_N(i,j,k)/(bgrid_temp(3) - &
                             bgrid_temp(2)) * (igrid(2) - bgrid_temp(2))/(igrid(2)-igrid(1))
                  

                    rhs_diag    = - (c1-Crank) *dvel_dx(ij)* NOin_init(ij,k) + (c1-Crank) *&
                                ((ivel(ij,2)*iphin_N(i,j,2)) *&
                               (bgrid_temp(3)-igrid(2))/(bgrid_temp(3)-bgrid_temp(2))- (ivel(ij,1)*iphin_N(i,j,1)))/&
                              (igrid(2)-igrid(1))/zphin_N(i,j,2)* NOin_init(ij,k) + &
                                 NOin_init(ij,k) + reactb(ij,k-1) 
                 else   
                 

                    sbdiag(ij,2) = Crank *(ivel(ij,1)*zphin_N(i,j,2))/&
                           zphin_N(i,j,2)/(bgrid_temp(k+1) - igrid(1))
 
                    diag(ij,2)   = c1+ Crank * dvel_dx(ij) + dt/zphin_N(i,j,k)/(igrid(k) - &
                              igrid(k-1))* (iDin(i,j,k)/(bgrid_temp(k+1)- &
                              bgrid_temp(k))) - Crank* (ivel(ij,2)*iphin_N(i,j,2))/zphin_N(i,j,2)* &
                              (bgrid_temp(3)-igrid(2))/(bgrid_temp(3)-bgrid_temp(2))/(igrid(2)-igrid(1))


                    spdiag(ij,2) = - Crank *(ivel(ij,2)*iphin_N(i,j,2))/zphin_N(i,j,k)/(bgrid_temp(3) - &
                             bgrid_temp(2)) * (igrid(2) - bgrid_temp(2))/(igrid(2)-igrid(1)) -  &
                             dt *iDin(i,j,k)/zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1))/(bgrid_temp(k+1) - bgrid_temp(k)) 

                     rhs_sb      =  -(c1 - Crank)*(ivel(ij,1)*zphin_N(i,j,2))/(igrid(2) - &
                              igrid(1))/zphin_N(i,j,2) 

                     rhs_sp      =  (c1-Crank)*(ivel(ij,2)*iphin_N(i,j,2))/zphin_N(i,j,k)/(bgrid_temp(3) - &
                             bgrid_temp(2)) * (igrid(2) - bgrid_temp(2))/(igrid(2)-igrid(1))

                     rhs_diag    = -(c1-Crank) * dvel_dx(ij)* NOin_init(ij,k) + (c1-Crank)*&
                               (ivel(ij,2)*iphin_N(i,j,2))/zphin_N(i,j,2)* &
                              (bgrid_temp(3)-igrid(2))/(bgrid_temp(3)-bgrid_temp(2))/(igrid(2)-igrid(1))* &
                                NOin_init(ij,k) + NOin_init(ij,k) + reactb(ij,k-1) 

              
                  endif

 

             elseif (k == nblyr+1) then

                if (dhb(i,j) < 0) then    !bottom melt

                   sbdiag(ij,nblyr+1) = Crank *(ivel(ij,nblyr)*iphin_N(i,j,nblyr))/&
                          zphin_N(i,j,k)/(igrid(nblyr+1) - &
                          igrid(nblyr))*(bgrid_temp(nblyr+1)-igrid(nblyr))/(bgrid_temp(nblyr+1)-bgrid_temp(nblyr)) -  & 
                          dt *iDin(i,j,k-1)/zphin_N(i,j,k)/(igrid(k) - &
                          igrid(k-1))/(bgrid_temp(k) - bgrid_temp(k-1)) 

                   diag(ij,nblyr+1)   = c1 + Crank* dvel_dx(ij) + dt/zphin_N(i,j,k)/(igrid(k) - &
                              igrid(k-1))* (iDin(i,j,k)/(bgrid_temp(k+1)- &    
                              igrid(k)) + iDin(i,j,k-1)/(bgrid_temp(k)-  &  
                              bgrid_temp(k-1))) - Crank * ((ivel(ij,nblyr+1)*zphin_N(i,j,nblyr+1)) - &
                              (ivel(ij,nblyr)*iphin_N(i,j,nblyr))*(igrid(nblyr)-bgrid_temp(nblyr))/&
                              (bgrid_temp(nblyr+1)-bgrid_temp(nblyr)))/zphin_N(i,j,nblyr+1)/&
                             (igrid(nblyr+1)-igrid(nblyr))


                   spdiag(ij,nblyr+1) = - dt *iDin(i,j,k)/zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1))/(bgrid_temp(k+1) - igrid(k))


                   rhs_sb      =  -(c1 - Crank)*(ivel(ij,nblyr)*iphin_N(i,j,nblyr))/&
                          zphin_N(i,j,k)/(igrid(nblyr+1) - &
                          igrid(nblyr))*(bgrid_temp(nblyr+1)-igrid(nblyr))/(bgrid_temp(nblyr+1)-bgrid_temp(nblyr))

                   rhs_sp      = c0
                 
                   rhs_diag    = -(c1-Crank)* dvel_dx(ij)*NOin_init(ij,k) + (c1-Crank)* &
                             ((ivel(ij,nblyr+1)*iphin_N(i,j,nblyr+1)) - &
                              (ivel(ij,nblyr)*iphin_N(i,j,nblyr))*(igrid(nblyr)-bgrid_temp(nblyr))/&
                              (bgrid_temp(nblyr+1)-bgrid_temp(nblyr)))/zphin_N(i,j,nblyr+1)/&
                             (igrid(nblyr+1)-igrid(nblyr))*NOin_init(ij,k) + &
                             NOin_init(ij,k) + reactb(ij,k-1)

                else

                    sbdiag(ij,nblyr+1) = Crank *(ivel(ij,nblyr)*iphin_N(i,j,nblyr))/&
                         zphin_N(i,j,k)/(igrid(k) - &
                         igrid(k-1))*(bgrid_temp(nblyr+1)-igrid(nblyr))/(bgrid_temp(nblyr+1)-bgrid_temp(nblyr)) -  & 
                          dt *iDin(i,j,k-1)/zphin_N(i,j,k)/(igrid(k) - &
                          igrid(k-1))/(bgrid_temp(k) - bgrid_temp(k-1))

                   diag(ij,nblyr+1)   = c1 + Crank* dvel_dx(ij) + dt/zphin_N(i,j,k)/(igrid(k) - &
                              igrid(k-1))* (iDin(i,j,k)/(bgrid_temp(k+1)- &
                              igrid(k)) + iDin(i,j,k-1)/(bgrid_temp(k)-  & 
                              bgrid_temp(k-1))) + Crank * (ivel(ij,nblyr)*iphin_N(i,j,nblyr))/&
                              zphin_N(i,j,nblyr+1)* &
                              (igrid(nblyr)-bgrid_temp(nblyr))/(bgrid_temp(nblyr+1)-bgrid_temp(nblyr))/&
                              (igrid(nblyr+1)-igrid(nblyr))

                  spdiag(ij,nblyr+1) = - Crank *(ivel(ij,nblyr+1)*iphin_N(i,j,nblyr+1))/&
                             zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1)) -  &
                             dt *iDin(i,j,k)/zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1))/(bgrid_temp(k+1) - igrid(k))



                 rhs_sb      =  -(c1 - Crank)*(ivel(ij,nblyr)*iphin_N(i,j,nblyr))/&
                           zphin_N(i,j,k)/(igrid(k) - &
                          igrid(k-1))*(bgrid_temp(nblyr+1)-igrid(nblyr))/(bgrid_temp(nblyr+1)-bgrid_temp(nblyr))

                 rhs_sp      =  (c1-Crank)*(ivel(ij,nblyr+1)*iphin_N(i,j,nblyr+1))/&
                             zphin_N(i,j,k)/(igrid(k) -  igrid(k-1))
                 
                 rhs_diag    =  -(c1-Crank)* dvel_dx(ij)* NOin_init(ij,k) - (c1-Crank)* &
                               (ivel(ij,nblyr)*iphin_N(i,j,nblyr))/zphin_N(i,j,nblyr+1)* &
                              (igrid(nblyr)-bgrid_temp(nblyr))/(bgrid_temp(nblyr+1)-bgrid_temp(nblyr))/&
                              (igrid(nblyr+1)-igrid(nblyr))* NOin_init(ij,k) + &
                              NOin_init(ij,k) + reactb(ij,k-1)
               endif

             else    
 

                 sbdiag(ij,k) = Crank *(vel(ij,k-1)*zphin_N(i,j,k-1))/&
                           zphin_N(i,j,k)/(bgrid_temp(k+1) - bgrid_temp(k-1)) -  & 
                           dt *iDin(i,j,k-1)/zphin_N(i,j,k)/(igrid(k) - &
                           igrid(k-1))/(bgrid_temp(k) - bgrid_temp(k-1))
  
                 diag(ij,k)   = c1 + dt/zphin_N(i,j,k)/(igrid(k) - &
                              igrid(k-1))* (iDin(i,j,k)/(bgrid_temp(k+1)- &
                              bgrid_temp(k)) + iDin(i,j,k-1)/(bgrid_temp(k)-  &
                              bgrid_temp(k-1)))

                 spdiag(ij,k) = - Crank *(vel(ij,k-1)*zphin_N(i,j,k+1))/zphin_N(i,j,k)/(bgrid_temp(k+1) - &
                             bgrid_temp(k-1)) - dt *iDin(i,j,k)/zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1))/(bgrid_temp(k+1) - bgrid_temp(k)) 


                 rhs_sb      =  -(c1 - Crank)*(vel(ij,k-1)*zphin_N(i,j,k-1))/&
                           zphin_N(i,j,k)/(bgrid_temp(k+1)-bgrid_temp(k-1))

                 rhs_sp      =  (c1-Crank)*(vel(ij,k-1)*zphin_N(i,j,k+1))/&
                            zphin_N(i,j,k)/(bgrid_temp(k+1)- bgrid_temp(k-1)) 

                 rhs_diag    = NOin_init(ij,k) + reactb(ij,k-1)  

            endif 

            rhs(ij,k)   = rhs_sb * NOin_init(ij,k-1) + rhs_diag &
                           + rhs_sp * NOin_init(ij,k+1)



           enddo               !k
      endif                    !hin_old > thin
      enddo                    ! ij


!790 format (I6,I6)
!781 format (I6,I6,I6)
!792 format (2D15.2)
!793 format (3D15.2)
!794 format (4D15.2)
!795 format (2D15.2,F10.4)
!800 format (F10.4)

      end subroutine get_matrix_elements_calc_bgc

!=======================================================================
!BOP
!
! !ROUTINE: tridiag_solver - tridiagonal matrix solver
!
! !DESCRIPTION:
!
! Tridiagonal matrix solver-- for salinity
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine tridiag_solver (icells,    &
                                 nmat,     sbdiag,   &
                                 diag,     spdiag,   &
                                 rhs,      xout)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
          icells      ! number of cells with aicen > puny


      integer (kind=int_kind), intent(in) :: &
         nmat            ! matrix dimension

      real (kind=dbl_kind), dimension (icells,nmat), &
           intent(in) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      real (kind=dbl_kind), dimension (icells,nmat), &
           intent(inout) :: &
         xout            ! solution vector
!
!EOP
!
      integer (kind=int_kind) :: &
         ij   , & ! horizontal index, combines i and j loops
         k               ! row counter

      real (kind=dbl_kind), dimension (icells) :: &
         wbeta           ! temporary matrix variable

      real (kind=dbl_kind), dimension(icells,nmat):: &
         wgamma          ! temporary matrix variable

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         wbeta(ij) = diag(ij,1)
         xout(ij,1) = rhs(ij,1) / wbeta(ij)
      enddo                     ! ij

      do k = 2, nmat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            wgamma(ij,k) = spdiag(ij,k-1) / wbeta(ij)
            wbeta(ij) = diag(ij,k) - sbdiag(ij,k)*wgamma(ij,k)
            xout(ij,k) = (rhs(ij,k) - sbdiag(ij,k)*xout(ij,k-1)) &
                         / wbeta(ij)
         enddo                  ! ij
      enddo                     ! k

      do k = nmat-1, 1, -1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            xout(ij,k) = xout(ij,k) - wgamma(ij,k+1)*xout(ij,k+1)
         enddo                  ! ij
      enddo                     ! k

      end subroutine tridiag_solver

!=======================================================================
!BOP
!
! !IROUTINE: bgc_diags - writes max,min,global sums to standard out
!
! !INTERFACE:
!
      subroutine bgc_diags (dt)
!
! !DESCRIPTION:
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW
!          Nicole Jeffery, LANL
!
! !USES:
!
      use ice_broadcast, only: broadcast_scalar
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, pbloc, plat, plon
      use ice_domain_size, only: ncat, nltrcr
      use ice_grid, only: lmask_n, lmask_s, tarean, tareas, grid_type
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, k, n, nn, ii,jj, iblk

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         pN_sk, pC_sk, pchl_sk, pNit_sk, pAm_sk, pSil_sk, &
         pDMSPp_sk, pDMSPd_sk, pDMS_sk, &
         pNit_ac, pAm_ac, pSil_ac, pDMSP_ac, pDMS_ac

      ! vertical  fields of category 1 at diagnostic points for bgc layer model
      real (kind=dbl_kind), dimension(npnt,nblyr) :: &
         pNO,pN, pSil, pPON, pgrowN

      real (kind=dbl_kind), dimension(npnt,nblyr+1) :: &
         pzfswin      

      !-----------------------------------------------------------------
      ! biogeochemical state of the ice
      !-----------------------------------------------------------------

      if (print_points) then

      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               iblk = pbloc(n)

               if (tr_bgc_N_sk) then   !bgc

                 pN_sk(n) = trcr(i,j,nt_bgc_N_sk,iblk)         ! biogeochemistry 
                 pC_sk(n) = trcr(i,j,nt_bgc_C_sk,iblk)         ! biogeochemistry 
                 pchl_sk(n) = trcr(i,j,nt_bgc_chl_sk,iblk)     ! biogeochemistry 
                 pNit_sk(n) = trcr(i,j,nt_bgc_Nit_sk,iblk)     ! biogeochemistry 
                 pAm_sk(n) = trcr(i,j,nt_bgc_Am_sk,iblk)       ! biogeochemistry 
                 pSil_sk(n) = trcr(i,j,nt_bgc_Sil_sk,iblk)     ! biogeochemistry 
                 pDMSPp_sk(n) = trcr(i,j,nt_bgc_DMSPp_sk,iblk) ! biogeochemistry 
                 pDMSPd_sk(n) = trcr(i,j,nt_bgc_DMSPd_sk,iblk) ! biogeochemistry 
                 pDMS_sk(n) = trcr(i,j,nt_bgc_DMS_sk,iblk)     ! biogeochemistry 
                 pNit_ac(n) = nit(i,j,iblk)

                endif !tr_bgc_N_sk              
               
                
               if (nltrcr > 0) then   ! layer bgc
               ! for layer model bgc

                 do k = 1, nblyr+1
                   pzfswin(n,k) = c0
                   do nn = 1,ncat
                     pzfswin(n,k) = pzfswin(n,k) +  zfswin(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                   enddo !nn
                   if (vice(i,j,iblk) > c0) then
                     pzfswin(n,k) = pzfswin(n,k)/vice(i,j,iblk)
                   endif !vice
                 enddo  !k
                 do k = 1,nblyr
                       pNO(n,k) = c0
                       pSil(n,k) = c0
                       pN(n,k) = c0
                       pPON(n,k) = c0
                       pgrowN(n,k) =  growNp(i,j,k,iblk)

                   if (vice(i,j,iblk) > c0) pgrowN(n,k) = pgrowN(n,k)/vice(i,j,iblk)                   
                   if (tr_bgc_NO) pNO(n,k) =  trcr(i,j,nt_bgc_NO+k-1,iblk)                   
                   if (tr_bgc_Sil)pSil(n,k)=  trcr(i,j,nt_bgc_Sil+k-1,iblk)     
                   if (tr_bgc_N)  pN(n,k)=  trcr(i,j,nt_bgc_N+k-1,iblk)    
                   if (tr_bgc_PON)  pPON(n,k)=  trcr(i,j,nt_bgc_PON+k-1,iblk)
                  
                 enddo   !k
               endif !(nltrcr > 0)

            endif                 ! my_task = pmloc

            if (tr_bgc_N_sk) then                  ! bgc
             call broadcast_scalar(pN_sk    (n), pmloc(n))             
             call broadcast_scalar(pC_sk    (n), pmloc(n))             
             call broadcast_scalar(pchl_sk  (n), pmloc(n))             
             call broadcast_scalar(pNit_sk  (n), pmloc(n))             
             call broadcast_scalar(pAm_sk   (n), pmloc(n))             
             call broadcast_scalar(pSil_sk  (n), pmloc(n))             
             call broadcast_scalar(pDMSPp_sk(n), pmloc(n))             
             call broadcast_scalar(pDMSPd_sk(n), pmloc(n))             
             call broadcast_scalar(pDMS_sk  (n), pmloc(n))             
             call broadcast_scalar(pNit_ac  (n), pmloc(n))             
             call broadcast_scalar(pAm_ac   (n), pmloc(n))             
             call broadcast_scalar(pSil_ac  (n), pmloc(n))             
             call broadcast_scalar(pDMSP_ac (n), pmloc(n))             
             call broadcast_scalar(pDMS_ac  (n), pmloc(n))             
            endif   !tr_bgc_sk

           if (nltrcr > 0) then                   !layer bgc
             do k = 1,nblyr+1
              call broadcast_scalar(pzfswin (n,k), pmloc(n))
             enddo !k
             do k = 1,nblyr
              call broadcast_scalar(pNO (n,k), pmloc(n))
              call broadcast_scalar(pSil (n,k), pmloc(n))
              call broadcast_scalar(pN (n,k), pmloc(n))
              call broadcast_scalar(pPON (n,k), pmloc(n))
              call broadcast_scalar(pgrowN (n,k), pmloc(n))
             enddo  !k
           endif   !nltrcr

          
         enddo                  ! npnt

      endif                     ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

     if (my_task == master_task) then

       call flush_fileunit(nu_diag)

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

       if (print_points) then

        write(nu_diag,*) '                         '
        write(nu_diag,902) '       Lat, Long         ',plat(1),plon(1), &
                                                       plat(2),plon(2)
        write(nu_diag,903) '  my_task, iblk, i, j     ', &
                              pmloc(1),pbloc(1),piloc(1),pjloc(1), &
                              pmloc(2),pbloc(2),piloc(2),pjloc(2)

      if (tr_bgc_N_sk) then  
        write(nu_diag,*) '----------bgc----------'
        write(nu_diag,900) 'nitrogen, skeletal     = ',pN_sk(1),pN_sk(2)
        write(nu_diag,900) 'carbon, skeletal       = ',pC_sk(1),pC_sk(2)
        write(nu_diag,900) 'chlorophyll, skeletal  = ',pchl_sk(1),pchl_sk(2)
        write(nu_diag,900) 'nitrate, skeletal      = ',pNit_sk(1),pNit_sk(2)
        write(nu_diag,900) 'ammonia/um, skeletal   = ',pAm_sk(1),pAm_sk(2)
        write(nu_diag,900) 'silicon, skeletal      = ',pSil_sk(1),pSil_sk(2)
        write(nu_diag,900) 'DMSPp, skeletal        = ',pDMSPp_sk(1),pDMSPp_sk(2)
        write(nu_diag,900) 'DMSPd, skeletal        = ',pDMSPd_sk(1),pDMSPd_sk(2)
        write(nu_diag,900) 'DMS, skeletal          = ',pDMS_sk(1),pDMS_sk(2)
        write(nu_diag,900) 'nitrate, mixed layer   = ',pNit_ac(1),pNit_ac(2)
        write(nu_diag,900) 'ammonia/um, mixed layer= ',pAm_ac(1),pAm_ac(2)
        write(nu_diag,900) 'silicon, mixed layer   = ',pSil_ac(1),pSil_ac(2)
        write(nu_diag,900) 'DMSP, mixed layer      = ',pDMSP_ac(1),pDMSP_ac(2)
        write(nu_diag,900) 'DMS, mixed layer       = ',pDMS_ac(1),pDMS_ac(2)
       endif

       if (nltrcr > 0) then   
        write(nu_diag,*) '                         '
        write(nu_diag,*) ' Top down bgc Layer Model'
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'NO(1) Bulk NO3   ','NO(2) Bulk NO3 '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pNO(n,k),n=1,2), k = 1,nblyr)              
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'PON(1) NO3 tracer   ','PON(2) NO3 tracer '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pPON(n,k),n=1,2), k = 1,nblyr)              
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'growN(1) specific growth  ','growN(2) specific growth '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pgrowN(n,k),n=1,2), k = 1,nblyr) 
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'zfswin(1) PAR  ','zfswin(2) PAR '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pzfswin(n,k),n=1,2), k = 1,nblyr+1)              
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'N(1) Bulk algalN   ','N(2) Bulk algalN '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pN(n,k),n=1,2), k = 1,nblyr)              
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'Sil(1) Bulk Silicate   ','Sil(2) Bulk Silicate '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pSil(n,k),n=1,2), k = 1,nblyr)  

        endif                   ! nltrcr
       endif                    ! print_points
      endif                     ! my_task = master_task

  802 format (f24.17,2x,f24.17)
  803 format (a25,2x,a25)
  900 format (a25,2x,f24.17,2x,f24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
  903 format (a25,5x,i4,1x,i4,1x,i4,1x,i4,7x,i4,1x,i4,1x,i4,1x,i4)

      end subroutine bgc_diags


!=======================================================================
!BOP
!
! !ROUTINE: calc_bgc_fluxes -- of Lax-Wendroff numerics for thick ice
!                          
!
! !DESCRIPTION:
!
!  Written in flux conservative form 
!                d(h zphi c)/dt = d(v zphi c + D*h dc/dx)/dx + zphi hR  
!                where v = (x-1) dzt/dt - x dzb/dt (which includes flushing),
!                D includes the gravity drainage eddy term plus molecular = (De + zphi Dm)/h^2, 
!                and R is the algal chemistry = react/dt
!               
! !REVISION HISTORY:
!
! authors     Nicole Jeffery, LANL
!
!
!
! !INTERFACE:
!
      subroutine calc_bgc_fluxes &
                                     (hin, hin_old, Cin, Cin_old, Cin_top, Cin_bot, dzt, dzb, Dm_bot, &
                                      dC_dx_bot, dt, &
                                      igrid,good_numerics,flux_o_tot,flux_o_g,react,aicen) 
!
! !USES:
!
! 
!
! !INPUT/OUTPUT PARAMETERS:
!
  
   


      real (kind=dbl_kind), dimension(nblyr), intent(in) :: &
         Cin           ,&  ! bulk c = zphi*c
         Cin_old       , & ! old bulk c = zphi_old * c_old
         react             ! reaction term

      real (kind=dbl_kind), intent(in) :: &
         Cin_top,        &  !bulk top conc.
         Cin_bot,        &
         aicen,  &
         hin_old   , &       ! brine thickness (m) 
         hin,      &       ! new brine thickness (m)
         dt, &        ! time-step 
         dzt, &
         dzb, &
         dC_dx_bot, &      !brine trace conc gradient
         Dm_bot          ! bottom molecular diffusive flux

      !---------------------------------------------------
      !  These fluxes are +ive into ocean (different from below)
      !---------------------------------------------------
      real (kind=dbl_kind), intent(inout) :: &
         flux_o_tot, flux_o_g
                          ! tracer flux, gravity+molecular drainage flux ,
                          ! and boundary flux to ocean (kg/m^2/s)                                  
      real (kind=dbl_kind), dimension (nblyr + 1), intent(in) :: &
         igrid            ! biology nondimensional grid interface points 

      logical (kind=log_kind), intent(inout) :: &   
         good_numerics          ! true if conservation satisfied within error

!
!
!EOP
!
     integer (kind=int_kind) :: &
         k, m, k_t, k_b ! vertical biology layer index 

    
     real (kind=dbl_kind) :: &
         sum_old      , &  !
         sum_new      , &  !
         sum_true     , &  !
         sum_true2     , &  !
         sum_react     , & !
         dh           , &  ! hin-hin_old
         dsum         , &  ! sum_new - sum_true
         dsum2        , &  ! sum_true2 - sum_new
         htrue        , &  !
         !-------------------------------------------------------------------
         ! Net fluxes (flux, fluxg, fluxb, fluxf) are defined positive into the ice
         ! Bottom fluxes (eg fluxg_b) are positive down into the ocean
         ! Top fluxes (eg fluxg_t) are positive down into the ice
         !-------------------------------------------------------------------
         fluxg, fluxg_b, fluxg_t, & !gravity drainage net  flux 
         fluxb, fluxb_b, fluxb_t, & !boundary motion net flux, 
         flux_o_b, & ! positive into ocean
         flux   , & ! total flux just gravity and boundary
         diff_dt   !uses flux. in units of kg/m^2/s
        

       


     real (kind=log_kind), parameter :: &
         accuracy = 2.0e-5   !1.0e-5 kg/m^2/s difference between boundary fluxes 
                             
       
         good_numerics = .true.
         dh = hin - hin_old
         htrue = hin_old + dzb - dzt
         sum_old = c0
         sum_new = c0
         sum_true = c0
         sum_true2 = c0
         sum_react = c0

         do k = 1, nblyr
          sum_old = sum_old + Cin_old(k)*hin_old*(igrid(k+1)-igrid(k))  !hin_old
          sum_new = sum_new + Cin(k)*hin*(igrid(k+1)-igrid(k))  
          sum_true = sum_true + Cin(k)*htrue*(igrid(k+1)-igrid(k))  
          sum_react = sum_react + react(k)*hin_old*(igrid(k+1)-igrid(k))
          if (Cin(k) > c0) then
            sum_true2 = sum_true2 +Cin(k)*hin*(igrid(k+1)-igrid(k))
          endif 
         enddo
         dsum = sum_new - sum_true     !dh correction for algae
         dsum2 = sum_true2 - sum_new   !negative correction

         fluxb_b = dzb*Cin_bot  !dzb < 0 is meltb > 0 or downward flushing
         fluxb_t = dzt*Cin_top  !dzt < 0 is meltwater added to surface
         fluxb = fluxb_b - fluxb_t  + dsum + dsum2! salt flux into ice
        
         fluxg_b = Dm_bot*hin_old*dC_dx_bot*dt
         fluxg_t = c0
         fluxg = fluxg_b - fluxg_t

         
         flux = fluxb + fluxg  !tot salt flux into the ice
     
     !-------------------------------------
     !  Ocean flux: positive into the ocean
     !-------------------------------------    

         flux_o_b = -fluxb/dt   !all boundary additions come from ocean
         flux_o_g = -fluxg/dt
         flux_o_tot = -flux/dt

         diff_dt = (sum_true2 - sum_old -sum_react- flux)/dt*aicen   
  
     !if (abs(diff_dt) > accuracy ) then
     !      good_numerics = .false.
     !      write(nu_diag,*) 'Conservation failure in bgc,diff_dt,accuracy:',diff_dt, accuracy
     !      write(nu_diag,*) 'Gravity flux*dt,Boundary flux*dt, total flux*dt:'&
     !                        ,fluxg,fluxb,flux
     !      write(nu_diag,*) 'Gravity flux* dt top, bottom:',fluxg_t,fluxg_b
     !      write(nu_diag,*) 'Boundary flux*dt top, bottom:',fluxb_t,fluxb_b
     !      write(nu_diag,*) 'dzt, dzb:', dzt, dzb
     !      write(nu_diag,*) 'hin, hin_old, dh:',hin, hin_old, dh
     !      write(nu_diag,*) 'Boundary values:, Cin_top,Cin_bot',Cin_top,Cin_bot
     !      write(nu_diag,*) 'Near Boundary values:, Cin(1),Cin(nblyr)',Cin(1),Cin(nblyr)
     !      write(nu_diag,*) 'Total initial tracer',sum_old
     !      write(nu_diag,*) 'Total final  tracer',sum_true2
     !      write(nu_diag,*) 'Total reactions',sum_react
     !      write(nu_diag,*) 'sum_new and sum_true',sum_new,sum_true
     !      write(nu_diag,*) 'fluxes into ocean (kg/m^2/s):'
     !      write(nu_diag,*) 'Gravity flux ',flux_o_g
     !      write(nu_diag,*) 'Boundary flux ',flux_o_b
     !      write(nu_diag,*) 'Total flux_o_tot',flux_o_tot
     !      write(nu_diag,*) 'aicen:',aicen
     !endif


     end subroutine calc_bgc_fluxes

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_bgc - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_bgc(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for a bgc restart
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size, only: ncat
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_domain, only: nblocks
      use ice_restart, only: lenstr, restart_dir, restart_file
      use ice_state, only: trcrn
!      use ice_zbgc_public, only:tr_bgc_N, tr_bgc_NO,tr_bgc_NH, &
!                  tr_bgc_Sil, tr_bgc_DMSPp, tr_bgc_DMSPd,  &
!                  tr_bgc_PON, nit, amm, sil, dmsp, dms, algalN
                  
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         iyear = nyr + year_init - 1
         imonth = month
         iday = mday
         
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.bgc.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_bgc,filename,0)

      if (my_task == master_task) then
        write(nu_dump_bgc) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
        write(nu_diag,*) 'BGC Restart written ',istep1,time,time_forc
      endif

      diag = .true.
      !--------------------------
      !bgc_stuff
      !---------------------------
      do n = 1, ncat
         do k = 1,nblyr
         if (tr_bgc_N) call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_N+k-1,n,:),'ruf8',diag)

         enddo
         do k = 1,nblyr
      
         if (tr_bgc_NO) call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_NO+k-1,n,:),'ruf8',diag)

         enddo
         do k = 1,nblyr
         if (tr_bgc_NH) call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_NH+k-1,n,:),'ruf8',diag)

         enddo
         do k = 1,nblyr
         if (tr_bgc_Sil) call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Sil+k-1,n,:),'ruf8',diag)

         enddo
         do k = 1,nblyr
         if (tr_bgc_DMSPp) call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPp+k-1,n,:),'ruf8',diag)

         enddo
         do k = 1,nblyr
         if (tr_bgc_DMSPd) call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPd+k-1,n,:),'ruf8',diag)

         enddo
         do k = 1,nblyr
         if (tr_bgc_PON) call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_PON+k-1,n,:),'ruf8',diag)

         enddo
      enddo
      
      !---------------------------------------------------      
      ! ocean values
      !-----------------------------------------------------
      if (tr_bgc_N) call ice_write(nu_dump_bgc,0,algalN,'ruf8',diag)
      if (tr_bgc_NO) call ice_write(nu_dump_bgc,0,nit,'ruf8',diag)
      if (tr_bgc_NH) call ice_write(nu_dump_bgc,0,amm,'ruf8',diag)
      if (tr_bgc_Sil) call ice_write(nu_dump_bgc,0,sil,'ruf8',diag)
      if (tr_bgc_DMSPp) call ice_write(nu_dump_bgc,0,dmsp,'ruf8',diag)
      if (tr_bgc_DMSPd) call ice_write(nu_dump_bgc,0,dms,'ruf8',diag)

      if (my_task == master_task) close(nu_dump_bgc)

      end subroutine write_restart_bgc

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_bgc - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_bgc(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for a bgc restart
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size, only: ncat
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init, &
                              istep0
      use ice_domain, only: nblocks
      use ice_restart, only: lenstr, restart_dir, restart_file, pointer_file, runtype
      use ice_state, only: trcrn
      use ice_exit, only: abort_ice
!      use ice_zbgc_public, only: tr_bgc_N, tr_bgc_NO,tr_bgc_NH, tr_bgc_C, tr_bgc_chl,&
!                  tr_bgc_Sil, tr_bgc_DMSPp, tr_bgc_DMSPd, tr_bgc_PON,  &
!                  nit, amm, sil, dmsp, dms, &
!                  algalN, R_C2N, R_chl2N, R_S2N
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: &
         filename, filename0, string1, string2

      logical (kind=log_kind) :: &
         diag, hit_eof

      if (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         read(nu_rst_pointer,'(a)') filename0
         filename = trim(filename0)
         close(nu_rst_pointer)

         ! reconstruct path/file
         n = index(filename0,trim(restart_file))
         if (n == 0) call abort_ice('bgc restart: filename discrepancy')
         string1 = trim(filename0(1:n-1))
         string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
         write(filename,'(a,a,a,a)') &
            string1(1:lenstr(string1)), &
            restart_file(1:lenstr(restart_file)),'.bgc', &
            string2(1:lenstr(string2))
      endif ! master_task

      call ice_open(nu_restart_bgc,filename,0)

      if (my_task == master_task) then
        read(nu_restart_bgc) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
        write(nu_diag,*) 'BGC Restart read at istep=',istep0,time,time_forc
      endif

      diag = .true.

      !-----------------------------------------------------------------
    
      do n = 1, ncat
         if (tr_bgc_N) then
           if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' N for each bgc layer'
           do k = 1,nblyr
             call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_N+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
                if (tr_bgc_C .OR. tr_bgc_chl) then
                  do iblk = 1, max_blocks
                  do i = 1,nx_block
                  do j = 1,ny_block
                  if (tr_bgc_C) trcrn(i,j,nt_bgc_C+k-1,n,iblk) = trcrn(i,j,nt_bgc_N+k-1,n,iblk)*R_C2N
                  if (tr_bgc_chl) trcrn(i,j,nt_bgc_chl+k-1,n,iblk) = trcrn(i,j,nt_bgc_N+k-1,n,iblk)*R_chl2N
                  enddo  !j
                  enddo  !i
                  enddo  !iblk
                endif
          enddo
         endif
         
         if (tr_bgc_NO) then
           if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' NO for each bgc layer'
           do k = 1,nblyr
             call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_NO+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
           enddo
         endif
         if (tr_bgc_NH) then
           if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' NH for each bgc layer'
           do k = 1,nblyr
             call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_NH+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
           enddo
         endif
         if (tr_bgc_Sil) then
           if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' Sil for each bgc layer'
           do k = 1,nblyr
            call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Sil+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
           enddo
          endif
         if (tr_bgc_DMSPp) then
           if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' DMSPp for each bgc layer'
           do k = 1,nblyr
             call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPp+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
           enddo
         endif
         if (tr_bgc_DMSPd) then
           if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' DMSPd for each bgc layer'
           do k = 1,nblyr
              call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPd+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
           enddo
         endif
         if (tr_bgc_PON) then
           if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' PON for each bgc layer'
           do k = 1,nblyr
              call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_PON+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
           enddo
         endif

      enddo
     
      if (my_task == master_task) &
         write(nu_diag,*) 'ocean bgc fields:  algalN, nit, amm, sil, dmsp, dms'
      
      
       if (tr_bgc_N)call ice_read(nu_restart_bgc,0,algalN,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
       if (tr_bgc_NO)call ice_read(nu_restart_bgc,0,nit,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
       if (tr_bgc_NH)call ice_read(nu_restart_bgc,0,amm,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
       if (tr_bgc_Sil)call ice_read(nu_restart_bgc,0,sil,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
       if (tr_bgc_DMSPp)call ice_read(nu_restart_bgc,0,dmsp,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
       if (tr_bgc_DMSPd)call ice_read(nu_restart_bgc,0,dms,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
     
       if (my_task == master_task) close(nu_restart_bgc)

      end subroutine read_restart_bgc

!=======================================================================

      end module ice_algae

!=======================================================================
