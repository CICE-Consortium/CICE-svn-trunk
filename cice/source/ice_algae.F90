!=======================================================================
!
! Compute biogeochemistry vertically resolved.
!
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
      module ice_algae

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: nblyr, nilyr, nblyr_hist, max_blocks, nbltrcr
      use ice_blocks, only: nx_block, ny_block
      use ice_fileunits, only: nu_diag, nu_restart_bgc, nu_rst_pointer, &
          nu_dump_bgc, flush_fileunit
      use ice_read_write, only: ice_open, ice_read, ice_write
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_zbgc_public
      use ice_state, only: vicen, vice, trcr, ntrcr, ntraceb, nt_bgc_am_sk, &
          nt_bgc_c_sk, nt_bgc_chl_sk, nt_bgc_DMS_sk, nt_bgc_DMSPd_sk, &
          nt_bgc_DMSPp_sk, nt_bgc_N, nt_bgc_N_sk, nt_bgc_Nit_sk, nt_bgc_NO, &
          nt_bgc_PON, nt_bgc_Sil, nt_bgc_Sil_sk,nt_bgc_C, nt_bgc_chl, &
          nt_bgc_DMS, nt_bgc_DMSPd, nt_bgc_DMSPp, nt_bgc_NH, nt_bgc_NO, &
          nt_bgc_PON, nt_bgc_Sil, nt_bgc_Nit_sk, nt_bgc_Sil_sk, nt_bgc_Nit_sk, &
          nt_bgc_Sil_sk, nt_bgc_S

      implicit none

      private
      public :: get_forcing_bgc, bgc_diags, write_restart_bgc, &
                algal_dyn, read_restart_bgc, &
                z_biogeochemistry, skl_biogeochemistry

      real (kind=dbl_kind), parameter, private :: &
         R_Si2N   = 1.5_dbl_kind,  & ! algal Si to N (mole/mole)
         zphimin  = 0.01_dbl_kind    ! minimum porosity for bgc only

!=======================================================================

      contains

!=======================================================================
!
! Read and interpolate annual climatologies of silicate and nitrate.
! Restore model quantities to data if desired.
!
! author: Elizabeth C. Hunke, LANL

      subroutine get_forcing_bgc

      use ice_calendar, only: dt, istep, mday, month, sec
      use ice_domain, only: nblocks
      use ice_flux, only:  hmix, sss
      use ice_forcing, only: trestore, trest, ocn_data_dir, &
          read_clim_data, interpolate_data, interp_coeff_monthly

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

         nit_file = trim(bgc_data_dir)//'nitrate_WOA2005_surface_monthly' ! gx1 only
         sil_file = trim(bgc_data_dir)//'silicate_WOA2005_surface_monthly' ! gx1 only

         if (my_task == master_task .and. istep == 1) then
         if (trim(sil_data_type)=='clim' .AND. (tr_bgc_Sil_sk .OR. tr_bgc_Sil)) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'silicate data interpolated to timestep:'
            write (nu_diag,*) trim(sil_file)
         endif
         if (trim(nit_data_type)=='clim' .AND. (tr_bgc_Nit_sk .OR. tr_bgc_NO)) then
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
         if (mday >= midmonth) ixm = -99 ! other two points will be used
         if (mday <  midmonth) ixp = -99

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

      if (trim(sil_data_type)=='clim'  .AND. (tr_bgc_Sil_sk .OR. tr_bgc_Sil)) then
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

      if (trim(nit_data_type)=='clim' .AND. (tr_bgc_Nit_sk .OR. tr_bgc_NO)) then
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

      elseif (trim(nit_data_type) == 'sss'  .AND. (tr_bgc_NO .OR. tr_bgc_Nit_sk)) then
          
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  nit(i,j,iblk) =  sss(i,j,iblk)        
               enddo
               enddo
            enddo
                 
      endif

      end subroutine get_forcing_bgc

!=======================================================================
!
! After algal_dynamics but separates the skeletal layer 
! biochemistry from the physics 
! 
      subroutine skl_biogeochemistry (nx_block, ny_block, &
                                 icells,  n_cat, dt,             &
                                 indxi,    indxj,         &
                                 flux_bio,   ocean_bio,   &
                                 hmix, aicen,         &
                                 meltb,    congel,        &
                                 fswthrul, first_ice, &
                                 trcrn, Iavgn, &
                                 grow_Cn, fN_partn)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells, n_cat         ! number of cells with aicen > puny 
                               !category

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt         ! time step 

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         first_ice  ! initialized values should be used

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         hmix   , & ! mixed layer depth
         aicen  , & ! ice area 
         meltb  , & ! bottom ice melt
         congel , & ! bottom ice growth 
         fswthrul   ! shortwave passing through ice to ocean

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn
   
    !-------------------------
    ! history variables
    !----------------------------   	

     real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         Iavgn,   & ! attenuated Fswthru  (W/m^2)
         grow_Cn, & ! specific growth (1/s)  
         fN_partn   ! particulate N flux  (mmol/m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbltrcr), intent(inout) :: &
         flux_bio,& ! ocean tracer flux (mmol/m^2/s) positive into ocean
         ocean_bio  ! ocean tracer concentration (mmol/m^3)

      !  local variables

      integer (kind=int_kind) :: i, j, ij, nn

      real (kind=dbl_kind), dimension(icells,nbltrcr):: &
         react      , & ! biological sources and sinks for equation matrix
         in_init    , & ! Initial concentration from previous time-step
         congel_alg     ! (mg N/m^2/s) congelation flux contribution to ice algae 
                        ! (used as initialization)

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         growN      ,&  !  algal growth rate    (mmol/m^3/s)
         upNOn      ,&  !  algal NO uptake rate (mmol/m^3/s)
         upNHn          !  algal NH uptake rate (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nbltrcr):: &
         flux_bio_temp  !  tracer flux to ocean (mmol/m^2/s)

      real (kind=dbl_kind), parameter :: &
         PVc = 1.e-6_dbl_kind           , & ! type 'constant' piston velocity for interface (m/s) 
         PV_scale_growth = p5           , & ! scale factor in Jin code PV during ice growth
         PV_scale_melt = p05            , & ! scale factor in Jin code PV during ice melt
         MJ1 = 9.667e-9_dbl_kind        , & ! (m/s) coefficients in Jin2008
         MJ2 = 38.8_dbl_kind            , & ! (1) from:4.49e-4_dbl_kind*secday   
         MJ3 = 1.04e7_dbl_kind          , & ! 1/(m/s) from: 1.39e-3_dbl_kind*secday^2  
         PV_frac_max = 0.9_dbl_kind         ! Maximum Piston velocity is 90% of skeletal layer/dt

      real (kind=dbl_kind), dimension(icells) :: &
         PVt       , & ! type 'Jin2006' piston velocity (m/s) 
         ice_growth, & ! Jin2006 definition: either congel rate or bottom melt rate  (m/s)
         f_meltn       ! vertical melt fraction of skeletal layer in dt

      real (kind=dbl_kind), dimension(nbltrcr) :: &
         PVn    , & ! 1 for tracers that flow with the brine and 0 otherwise
         fmeltn     ! 0 for tracers that flow and  1 for tracers that cling

      real (kind=dbl_kind):: &
         init_temp  ! tracer Concentration after biochemistry (mmol/m^3) 

  !-------------------------------------
  ! Initialize 
  !------------------------------------
          
  in_init(:,:) = c0
  bgc_flux_type = 'Jin2006'    ! or 'constant'
  congel_alg(:,:) = c0
  PVn(:) = c1
  PVt(:) = c0
  fmeltn(:) = c0
  f_meltn(:) = c0
  react(:,:) = c0
  fmeltn(nlt_bgc_N) = c1

  do nn = 1,nbltrcr          
       if (bgc_tracer_type(nn) < p5) PVn(nn) = c0
  enddo 

 !----------------------------------------------------------------------
 ! 'Jin2006':
 ! 1. congel/melt dependent piston velocity (PV) for both
 ! growth and melt
 ! 2. If congel > melt use 'congel' and melt > congel use 'melt'
 ! 3. For algal N, PV for ice growth only provides a seeding concentration 
 ! 4. Melt affects nutrients and algae in the same manner through PV(melt)
 !-----------------------------------------------------------------------
 if (trim(bgc_flux_type) == 'Jin2006') then  

   do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
 
         if (first_ice(i,j)) then     
          trcrn(i,j,nt_bgc_N_sk)                          = ocean_bio(i,j,nlt_bgc_N)*sk_l/phi_sk
          if (tr_bgc_Nit_sk) trcrn(i,j,nt_bgc_Nit_sk)     = ocean_bio(i,j,nlt_bgc_NO)*sk_l/phi_sk
          if (tr_bgc_Am_sk)  trcrn(i,j,nt_bgc_Am_sk)      = ocean_bio(i,j,nlt_bgc_NH)*sk_l/phi_sk
          if (tr_bgc_Sil_sk) trcrn(i,j,nt_bgc_Sil_sk)     = ocean_bio(i,j,nlt_bgc_Sil)*sk_l/phi_sk
          if (tr_bgc_C_sk)   trcrn(i,j,nt_bgc_C_sk)       = ocean_bio(i,j,nlt_bgc_C)*sk_l/phi_sk 
          if (tr_bgc_chl_sk) trcrn(i,j,nt_bgc_chl_sk)     = ocean_bio(i,j,nlt_bgc_chl)*sk_l/phi_sk
          if (tr_bgc_DMSPp_sk) trcrn(i,j,nt_bgc_DMSPp_sk) = ocean_bio(i,j,nlt_bgc_DMSPp)*sk_l/phi_sk 
          if (tr_bgc_DMSPd_sk) trcrn(i,j,nt_bgc_DMSPd_sk) = ocean_bio(i,j,nlt_bgc_DMSPd)*sk_l/phi_sk
          if (tr_bgc_DMS_sk) trcrn(i,j,nt_bgc_DMS_sk)     = ocean_bio(i,j,nlt_bgc_DMS)*sk_l/phi_sk
         endif

         in_init(ij,nlt_bgc_N)                           = trcrn(i,j,nt_bgc_N_sk)/sk_l
         if (tr_bgc_Nit_sk) in_init(ij,nlt_bgc_NO)       = trcrn(i,j,nt_bgc_Nit_sk)/sk_l
         if (tr_bgc_Am_sk)  in_init(ij,nlt_bgc_NH)       = trcrn(i,j,nt_bgc_Am_sk)/sk_l
         if (tr_bgc_Sil_sk) in_init(ij,nlt_bgc_Sil)      = trcrn(i,j,nt_bgc_Sil_sk)/sk_l
         if (tr_bgc_C_sk)   in_init(ij,nlt_bgc_C)        = trcrn(i,j,nt_bgc_C_sk)/sk_l
         if (tr_bgc_chl_sk) in_init(ij,nlt_bgc_chl)      = trcrn(i,j,nt_bgc_chl_sk)/sk_l
         if (tr_bgc_DMSPp_sk) in_init(ij,nlt_bgc_DMSPp)  = trcrn(i,j,nt_bgc_DMSPp_sk)/sk_l
         if (tr_bgc_DMSPd_sk) in_init(ij,nlt_bgc_DMSPd)  = trcrn(i,j,nt_bgc_DMSPd_sk)/sk_l
         if (tr_bgc_DMS_sk) in_init(ij,nlt_bgc_DMS)      = trcrn(i,j,nt_bgc_DMS_sk)/sk_l

         do nn = 1,nbltrcr
           if (in_init(ij,nn) < c0) then
             write(nu_diag,*)'Before: initial sk_bgc < 0, ij,nn,nbltrcr,in_init(ij,nn)', &
                    ij,nn,nbltrcr,in_init(ij,nn)
             call abort_ice ('ice_bgc.F90: BGC error1')
           endif
         enddo
         ice_growth(ij) = (congel(i,j)-meltb(i,j))/dt
         if (ice_growth(ij) > c0) then         ! ice_growth(ij) = congel(i,j)/dt
             PVt(ij) =-min(abs(PV_scale_growth*(MJ1 + MJ2*ice_growth(ij) - MJ3*ice_growth(ij)**2)), &
                      PV_frac_max*sk_l/dt)  
         else                                   !ice_growth(ij) = -meltb(i,j)/dt
             PVt(ij) = min(abs(PV_scale_melt*(MJ2*ice_growth(ij)-MJ3*ice_growth(ij)**2)), &
                      PV_frac_max*sk_l/dt)
         endif
   
       !---------------------------------------------------------
       ! Algae melt like nutrients
       !---------------------------------------------------------        
        if  (ice_growth(ij) < c0) then   !melt:  flux from ice to ocean
 
             f_meltn(ij) =  PVt(ij)*in_init(ij,nlt_bgc_N)   !for algae only

        elseif (ice_growth(ij) > c0 .AND. in_init(ij,nlt_bgc_N) < ocean_bio(i,j,nlt_bgc_N)/phi_sk) then

       !---------------------------------------------------------
       ! Growth only contributes to seeding from ocean 
       !---------------------------------------------------------
           congel_alg(ij,nlt_bgc_N) = (ocean_bio(i,j,nlt_bgc_N)/phi_sk - in_init(ij,nlt_bgc_N))*sk_l/dt
     
        endif  !PVt > c0       
    enddo   !ij

 !----------------------------------------------------------------------
 ! 'constant':
 ! 1. Constant PV  for congel > melt
 ! 2. For algae, PV for ice growth only provides a seeding concentration 
 ! 3. Melt loss (f_meltn) affects algae only and is proportional to melt
 !-----------------------------------------------------------------------
  else   !PV_type = 'constant'


    do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (congel(i,j) > meltb(i,j)) PVt(ij) =-PVc  

         if (first_ice(i,j) ) then     
          trcrn(i,j,nt_bgc_N_sk)                          = ocean_bio(i,j,nlt_bgc_N)*sk_l/phi_sk
          if (tr_bgc_Nit_sk)   trcrn(i,j,nt_bgc_Nit_sk)   = ocean_bio(i,j,nlt_bgc_NO)*sk_l/phi_sk
          if (tr_bgc_Am_sk)    trcrn(i,j,nt_bgc_Am_sk)    = ocean_bio(i,j,nlt_bgc_NH)*sk_l/phi_sk
          if (tr_bgc_Sil_sk)   trcrn(i,j,nt_bgc_Sil_sk)   = ocean_bio(i,j,nlt_bgc_Sil)*sk_l/phi_sk
          if (tr_bgc_C_sk)     trcrn(i,j,nt_bgc_C_sk)     = ocean_bio(i,j,nlt_bgc_C)*sk_l/phi_sk 
          if (tr_bgc_chl_sk)   trcrn(i,j,nt_bgc_chl_sk)   = ocean_bio(i,j,nlt_bgc_chl)*sk_l/phi_sk
          if (tr_bgc_DMSPp_sk) trcrn(i,j,nt_bgc_DMSPp_sk) = ocean_bio(i,j,nlt_bgc_DMSPp)*sk_l/phi_sk 
          if (tr_bgc_DMSPd_sk) trcrn(i,j,nt_bgc_DMSPd_sk) = ocean_bio(i,j,nlt_bgc_DMSPd)*sk_l/phi_sk
          if (tr_bgc_DMS_sk)   trcrn(i,j,nt_bgc_DMS_sk)   = ocean_bio(i,j,nlt_bgc_DMS)*sk_l/phi_sk
         endif

         in_init(ij,nlt_bgc_N)                          = trcrn(i,j,nt_bgc_N_sk)/sk_l
         if (tr_bgc_Nit_sk)   in_init(ij,nlt_bgc_NO)    = trcrn(i,j,nt_bgc_Nit_sk)/sk_l
         if (tr_bgc_Am_sk)    in_init(ij,nlt_bgc_NH)    = trcrn(i,j,nt_bgc_Am_sk)/sk_l
         if (tr_bgc_Sil_sk)   in_init(ij,nlt_bgc_Sil)   = trcrn(i,j,nt_bgc_Sil_sk)/sk_l
         if (tr_bgc_C_sk)     in_init(ij,nlt_bgc_C)     = trcrn(i,j,nt_bgc_C_sk)/sk_l
         if (tr_bgc_chl_sk)   in_init(ij,nlt_bgc_chl)   = trcrn(i,j,nt_bgc_chl_sk)/sk_l
         if (tr_bgc_DMSPp_sk) in_init(ij,nlt_bgc_DMSPp) = trcrn(i,j,nt_bgc_DMSPp_sk)/sk_l
         if (tr_bgc_DMSPd_sk) in_init(ij,nlt_bgc_DMSPd) = trcrn(i,j,nt_bgc_DMSPd_sk)/sk_l
         if (tr_bgc_DMS_sk)   in_init(ij,nlt_bgc_DMS)   = trcrn(i,j,nt_bgc_DMS_sk)/sk_l
        
        do nn = 1,nbltrcr
           if (in_init(ij,nn) < c0) then
             write(nu_diag,*)'Before: initial sk_bgc < 0, ij,nn,nbltrcr,in_init(ij,nn)', &
                    ij,nn,nbltrcr,in_init(ij,nn)
             call abort_ice ('ice_bgc.F90: BGC error1')
           endif
        enddo

        if (congel(i,j) .GE. meltb(i,j) .AND. in_init(ij,nlt_bgc_N) < ocean_bio(i,j,nlt_bgc_N)/phi_sk) then
             congel_alg(ij,nlt_bgc_N) = (ocean_bio(i,j,nlt_bgc_N)/phi_sk - in_init(ij,nlt_bgc_N))*sk_l/dt
        elseif (meltb(i,j) > congel(i,j)) then
             f_meltn(ij) = min(in_init(ij,nlt_bgc_N), min(c1,meltb(i,j)/sk_l)*in_init(ij,nlt_bgc_N))*sk_l/dt
        endif
      enddo    !ij

 endif

!  begin building biogeochemistry terms

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            
       call algal_dyn &
                        (nx_block, ny_block,        &
                         icells,                    &
                         indxi,    indxj,           &
                         fswthrul, react,           & 
                         in_init, nbltrcr,          &
                         grow_Cn, upNOn,            &
                         upNHn,                     &
                         tr_bgc_N_sk,tr_bgc_Nit_sk, &
                         tr_bgc_Am_sk,tr_bgc_Sil_sk,&
                         tr_bgc_C_sk,tr_bgc_chl_sk, &
                         tr_bgc_DMSPp_sk,           & 
                         tr_bgc_DMSPd_sk,           &
                         tr_bgc_DMS_sk)

!  Compute new tracer concencentrations
  
     do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
        
          do nn = 1,nbltrcr

           !--------------------------------------------------------------------
           ! if PVt(ij) > 0, ie melt, then ocean_bio term drops out (MJ2006)
           ! Combine boundary fluxes
           !-------------------------------------------------------------------        
           
           PVn(nn) = SIGN(PVn(nn),PVt(ij))
           init_temp = max(c0,in_init(ij,nn) + react(ij,nn))
           flux_bio_temp(nn) = (PVn(nn)*PVt(ij)*init_temp -PVn(nn)*min(c0,PVt(ij))*ocean_bio(i,j,nn))+ &
                                + f_meltn(ij)*fmeltn(nn) - congel_alg(ij,nn)

           if (init_temp - flux_bio_temp(nn)*dt/sk_l < c0) then
               flux_bio_temp(nn) = init_temp*sk_l/dt*(c1-puny)
           endif

           in_init(ij,nn) = init_temp - flux_bio_temp(nn)*dt/sk_l
           flux_bio(i,j,nn) = flux_bio(i,j,nn) +  flux_bio_temp(nn)*phi_sk  

           !--------------------------------------------------
           ! Uncomment to update ocean concentration
           ! Currently not coupled with ocean biogeochemistry
           !----------------------------------------------------
           !ocean_bio(i,j,nn) = ocean_bio(i,j,nn) + flux_bio(i,j,nn)/hmix(i,j)*aicen(i,j)

           if (in_init(ij,nn) < c0) then
                write(nu_diag,*)'after algal fluxes sk_bgc < 0, ij,nn,in_init(ij,nn),flux_bio(i,j,nn)', &
                    ij,nn,in_init(ij,nn),flux_bio(i,j,nn)
                write(nu_diag,*) 'init_temp,flux_bio_temp(nn),f_meltn(ij), congel_alg(ij,nn),PVt(ij),PVn(nn):'
                write(nu_diag,*) init_temp,flux_bio_temp(nn),f_meltn(ij), congel_alg(ij,nn),PVt(ij),PVn(nn)
                write(nu_diag,*) 'congel(i,j),meltb(i,j):'
                write(nu_diag,*) congel(i,j),meltb(i,j)
                call abort_ice ('ice_bgc.F90: BGC error3')
           endif
         
          enddo  !nbltrcr

           !------------------------
           ! combine reactions 
           ! and back to (m^2)
           !-------------------------

         trcrn(i,j,nt_bgc_N_sk) = in_init(ij,nlt_bgc_N)*sk_l
         if (tr_bgc_Nit_sk)trcrn(i,j,nt_bgc_Nit_sk) =  sk_l*in_init(ij,nlt_bgc_NO)
         if (tr_bgc_Am_sk)trcrn(i,j,nt_bgc_Am_sk) = sk_l*in_init(ij,nlt_bgc_NH)
         if (tr_bgc_Sil_sk)trcrn(i,j,nt_bgc_Sil_sk) = sk_l*in_init(ij,nlt_bgc_Sil)
         if (tr_bgc_C_sk) trcrn(i,j,nt_bgc_C_sk) =   trcrn(i,j,nt_bgc_N_sk)*R_C2N
         if (tr_bgc_chl_sk)  trcrn(i,j,nt_bgc_chl_sk) =  trcrn(i,j,nt_bgc_N_sk)*R_chl2N
         if (tr_bgc_DMSPp_sk) trcrn(i,j,nt_bgc_DMSPp_sk) = sk_l*in_init(ij,nlt_bgc_DMSPp)
         if (tr_bgc_DMSPd_sk) trcrn(i,j,nt_bgc_DMSPd_sk) = sk_l*in_init(ij,nlt_bgc_DMSPd)
         if (tr_bgc_DMS_sk)trcrn(i,j,nt_bgc_DMS_sk) = sk_l*in_init(ij,nlt_bgc_DMS)
        
    enddo !icells

    end subroutine skl_biogeochemistry

!=======================================================================
!
! Solve the scalar vertical diffusion equation implicitly using 
! tridiag_solver. Calculate the diffusivity from temperature and salinity.
! 
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice with 
! dynamic salinity or the height ratio == hinS/vicen*aicen, where hinS is the 
! height of the brine surface relative to the bottom of the ice.  This volume fraction
! may be > 1 in which case there is brine above the ice surface (meltponds). 
! 
      subroutine z_biogeochemistry (nx_block, ny_block,      &
                                  icells, n_cat, dt,         &
                                  indxii,   indxjj,          &   
                                  aicen,  vicen,             & 
                                  hice_old,     nitn,        & 
                                  flux_bio, flux_bio_g,      &
                                  ammn, siln,  dmspn,  dmsn, & 
                                  algalNn, zphin, iphin,     & 
                                  trcrn,  iDin,   sss,       &
                                  fswthrul,                  &
                                  growN, upNOn, upNHn,       &
                                  dh_top, dh_bot,            &
                                  dh_top_chl, dh_bot_chl,    &
                                  zfswin,                    &
                                  first_ice,  TLAT, TLON,    &
                                  hbri, hbri_old, darcy_V,   &
                                  darcy_V_chl, bgrid, igrid, &
                                  cgrid, flood_val,zphi_min)     

      use ice_calendar,    only: istep1, time
      use ice_therm_shared, only: solve_Sin, ktherm
      use ice_state, only: aice, nt_sice      

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells, n_cat         ! number of true cells with aicen > 0
                              
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxii, indxjj ! compressed indices for icells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt             ! time step 

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bgrid          ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid          ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid          ! CICE vertical coordinate   

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         sss         , & ! ocean salinity (ppt)
         TLAT        , &
         TLON        , &
         hbri        , & ! brine height  (m)
         hbri_old    , & ! brine height  (m)
         hice_old    , & ! ice height (m)
         darcy_V     , & ! darcy velocity
         darcy_V_chl , & ! darcy velocity for algae
         zphi_min    , & ! surface porosity
         dh_top      , & ! change in brine top (m)
         dh_bot      , & ! change in brine bottom (m)
         dh_bot_chl  , & ! change in brine bottom (m) felt by algae
         dh_top_chl      ! change in brine top (m) felt by algae

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         first_ice   , &  ! .true. for first ice
         flood_val        ! .true. if ocean floods surface    

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr+1), &
         intent(in) :: &
         fswthrul         ! Short wave radiation at each ice layer (W/m^2)  

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr + 1), &
         intent(out) :: &
         zfswin           ! Short wave flux on bio grid (W/m^2)  1 is surface point 
                          ! 2:nblyr+1 interior grid points
       
      real (kind=dbl_kind), dimension (nx_block,ny_block,nbltrcr), intent(inout) :: &
         flux_bio, &   ! total ocean tracer flux (mmol/m^2/s)
         flux_bio_g    ! ocean tracer flux from gravity drainage &
                       ! and molecular diffusion (mmol/m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &   
         !change to  inout when updating ocean fields
         nitn        , & ! ocean nitrate (mmol/m^3) 
         ammn        , & ! ocean ammonium (mmol/m^3) 
         siln        , & ! ocean silicate (mmol/m^3) 
         algalNn     , & ! ocean algal nitrogen (mmol/m^3)
         dmspn       , & ! ocean DMSP (mmol/m^3)
         dmsn            ! ocean DMS (mmol/m^3)  

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr), intent(inout) :: &
         growN       , & ! algal growth rate          (mmol/m^3/s)
         upNOn       , & ! algal nitrate uptake rate  (mmol/m^3/s)
         upNHn           ! algal ammonium uptake rate (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(in) :: &
         zphin            ! Porosity on the bgrid

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(in) :: &
         iphin       , &  ! Porosity on the igrid   
         iDin             ! Diffusivity/h on the igrid (1/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn           ! tracers: for salinity (ppt * vicen) and NO (mmol/m^3 * vicen)    

      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k, m, mm        ! vertical biology layer index 

      real (kind=dbl_kind), dimension(icells) :: &
         hin         , & ! ice thickness (m)        
         hin_old     , & ! ice thickness before current melt/growth (m)
         hinS        , & ! brine thickness (m)
         hinS_old    , & ! brine thickness before current melt/growth (m)
         surface_S       ! salinity of ice above hin > hinS
   
      real (kind=dbl_kind) :: &
         dh_add      , & ! Ice algae added throughout ice during growth (m) 
                         ! max(0,dh_bot-dh_bot_chl)
         Dm_bot      , & ! Bottom diffusivity/h (1/s)
         dh_bot_all      ! dh_bot or dh_bot_chl   

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2) :: &
         zphin_N         ! porosity for tracer model has minimum 
                         ! zphin_N >= zphimin

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1) :: &
         iDin_N      , & ! (1/s) diffusivity/h for algae (molecular only)
         iphin_N         ! tracer porosity on the igrid

      real (kind=dbl_kind), dimension(icells,nblyr_hist):: & 
         NOin        , & ! Local layer nitrate values (mmol/m^3)
         Nin         , & ! Local layer algal N values (mmol/m^3)
         Silin       , & ! Local layer silicate values (mmol/m^3)
         Cin         , & ! Local layer carbon values (mmol/m^3)
         chlin       , & ! Local layer chlorophyll values (mg/m^3)
         DMSPpin     , & ! Local layer DMSPp values (mmol/m^3)
         DMSPdin     , & ! Local layer DMSPd values (mmol/m^3)
         DMSin       , & ! Local layer DMS values (mmol/m^3)
         NHin        , & ! Local layer ammonium values (mmol/m^3)
         PONin           ! Local layer PON values (mmol/m^3)

      real (kind=dbl_kind), dimension(icells,nblyr_hist):: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix equation

      real (kind=dbl_kind), dimension(icells,nblyr,nbltrcr):: &
         react       , & ! biological sources and sinks for equation matrix
         react_phi       ! react*zphin_N
 
      real (kind=dbl_kind), dimension(icells, nblyr_hist,nbltrcr):: &
         in_init     , & ! Initial concentration from previous time-step
         biomat          ! Matrix output    

      real (kind=dbl_kind), dimension(icells, nblyr, nbltrcr):: &
         bulk_biomat , & ! Final ice bulk concentration from previous time-step
         bulk_biomat_o   ! initial ice bulk concentration

      real (kind=dbl_kind):: &
         bulk_top    , & ! top bulk concentration  
         bulk_bot    , & ! bot bulk concentration
         dC_dx_bot       ! bottom tracer gradient

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: &
         trtmp0      , & ! temporary, remapped tracers
         trtmp           ! temporary, remapped tracers

      real (kind=dbl_kind), dimension(icells):: &
         thin_check      ! c0 if thin and c1 if not

      logical (kind=log_kind), dimension(icells,nbltrcr) :: & 
         good_bio_numerics ! .false. if tracer conservation fails     
   
      logical (kind=log_kind):: & 
         check_conserve   ! .true. if check tracer conservation

      ! local parameters
         
      integer, parameter :: &
         nt_zfswin = 1    ! for interpolation of short wave to bgrid
 
  !------------------------------------------------------------------------ 
  ! define gravity drainage diffusivity on bio grid 
  ! if using tr_salinity. 
  !----------------------------------------------------------------------
  !-------------------------------------
  ! Initialize 
  !------------------------------------
        thin_check(:) = c0
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
              zphin_N(i,j,k)   = max(zphimin,zphin(i,j,k))
              iphin_N(i,j,k)   = max(zphimin,iphin(i,j,k))
              iphin_N(i,j,k+1) = max(zphimin,iphin(i,j,k+1))
              zphin_N(i,j,2) = zphi_min(i,j)
              zphin_N(i,j,1) = zphi_min(i,j)
              iphin_N(i,j,1) = zphi_min(i,j)

              iDin_N(i,j,k+1) = iphin_N(i,j,k+1)*Dm/hbri(i,j)**2 
   
                if (first_ice(i,j)) then                             !initialize
                       trcrn(i,j,nt_bgc_NO+k-1) = nitn(i,j)*initbio_frac  
                       NOin(ij,k+1) = trcrn(i,j,nt_bgc_NO+k-1)/zphin_N(i,j,k+1) 

                       if (trcrn(i,j,nt_bgc_NO+k-1) > 1000.0_dbl_kind) then
                          write(nu_diag,*)'Bulk Nitrate solution error initially, first_ice(i,j):',first_ice(i,j)      
                          write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,mm:'&
                                         ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                         TLON(i,j)*rad_to_deg,istep1,mm
                          write(nu_diag,*)'dh_bot_chl(i,j),dh_top(i,j)'&
                                         ,dh_bot_chl(i,j),dh_top(i,j) 
                          write(nu_diag,*)'ij,k,NOin(ij,k+1)'
                          write(nu_diag,*)ij,k,NOin(ij,k+1)
                          write(nu_diag,*)'k,trcrn(i,j,nt_bgc_NO+k-1),zphin_N(i,j,k+1)' 
                          write(nu_diag,*)k,trcrn(i,j,nt_bgc_NO+k-1),zphin_N(i,j,k+1)
                        if (ktherm == 1 ) then
                          write(nu_diag,*)'trcrn(i,j,nt_bgc_S+nblyr-1),trcrn(i,j,nt_bgc_S+nblyr-2)'
                          write(nu_diag,*)trcrn(i,j,nt_bgc_S+nblyr-1),trcrn(i,j,nt_bgc_S+nblyr-2)
                        elseif (ktherm == 2) then
                          write(nu_diag,*)'trcrn(i,j,nt_sice+nilyr-1),trcrn(i,j,nt_sice+nilyr-2)'
                          write(nu_diag,*)trcrn(i,j,nt_sice+nilyr-1),trcrn(i,j,nt_sice+nilyr-2)
                        endif

                          write(nu_diag,*)'hinS(ij),hinS_old(ij)' 
                          write(nu_diag,*)hinS(ij),hinS_old(ij)
                          call abort_ice ('ice_bgc.F90: BGC error1')
                       endif
   
                else
                       NOin(ij,k+1) = trcrn(i,j,nt_bgc_NO+k-1)/zphin_N(i,j,k+1)

                       if (trcrn(i,j,nt_bgc_NO+k-1) > 1000.0_dbl_kind .AND. tr_bgc_N) then
                          write(nu_diag,*)'Bulk Nitrate solution error initially, first_ice(i,j):',first_ice(i,j)      
                          write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,mm:'&
                                           ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                            TLON(i,j)*rad_to_deg,istep1,mm
                          write(nu_diag,*)'dh_bot_chl(i,j),dh_top(i,j)'&
                                           ,dh_bot_chl(i,j),dh_top(i,j) 
                          write(nu_diag,*)'ij,k,NOin(ij,k+1)'
                          write(nu_diag,*)ij,k,NOin(ij,k+1)
                          write(nu_diag,*)'k,trcrn(i,j,nt_bgc_NO+k-1),zphin_N(i,j,k+1)' 
                          write(nu_diag,*)k,trcrn(i,j,nt_bgc_NO+k-1),zphin_N(i,j,k+1)
                        if (ktherm == 1 ) then
                          write(nu_diag,*)'trcrn(i,j,nt_bgc_S+nblyr-1),trcrn(i,j,nt_bgc_S+nblyr-2)'
                          write(nu_diag,*)trcrn(i,j,nt_bgc_S+nblyr-1),trcrn(i,j,nt_bgc_S+nblyr-2)
                        elseif (ktherm == 2) then
                          write(nu_diag,*)'trcrn(i,j,nt_sice+nilyr-1),trcrn(i,j,nt_sice+nilyr-2)'
                          write(nu_diag,*)trcrn(i,j,nt_sice+nilyr-1),trcrn(i,j,nt_sice+nilyr-2)
                        endif
                          write(nu_diag,*)'hinS(ij),hinS_old(ij)'
                          write(nu_diag,*)hinS(ij),hinS_old(ij)
                          call abort_ice ('ice_bgc.F90: BGC error1')
                       endif
                     
                endif
                if (tr_bgc_N) then
                     if (first_ice(i,j)) then 
                         trcrn(i,j,nt_bgc_N+k-1) = algalNn(i,j) 
                         Nin(ij,k+1) = trcrn(i,j,nt_bgc_N+k-1)/zphin_N(i,j,k+1) 
                     else 
                          Nin(ij,k+1) = trcrn(i,j,nt_bgc_N+k-1)/zphin_N(i,j,k+1)
                     endif      
                endif
                if (tr_bgc_NH) then
                     if (first_ice(i,j)) then 
                         trcrn(i,j,nt_bgc_NH+k-1) = ammn(i,j)*initbio_frac     
                         NHin(ij,k+1) = trcrn(i,j,nt_bgc_NH+k-1)/zphin_N(i,j,k+1) 
                     else
                         NHin(ij,k+1) = trcrn(i,j,nt_bgc_NH+k-1)/zphin_N(i,j,k+1) 
                     endif
                endif
                if (tr_bgc_Sil) then
                     if (first_ice(i,j)) then 
                         trcrn(i,j,nt_bgc_Sil+k-1) = siln(i,j)*initbio_frac    
                         Silin(ij,k+1) = trcrn(i,j,nt_bgc_Sil+k-1)/zphin_N(i,j,k+1)
                     else
                         Silin(ij,k+1) = trcrn(i,j,nt_bgc_Sil+k-1)/zphin_N(i,j,k+1)   
                     endif
                endif
                if (tr_bgc_C) then
                     if (first_ice(i,j)) then 
                         trcrn(i,j,nt_bgc_C+k-1) = R_C2N*algalNn(i,j)
                         Cin(ij,k+1) = trcrn(i,j,nt_bgc_C+k-1)/zphin_N(i,j,k+1)  
                     else
                         Cin(ij,k+1) = trcrn(i,j,nt_bgc_C+k-1)/zphin_N(i,j,k+1)  
                     endif
                endif
                if (tr_bgc_chl) then
                     if (first_ice(i,j)) then 
                         trcrn(i,j,nt_bgc_chl+k-1) = R_chl2N*algalNn(i,j) 
                         chlin(ij,k+1) = trcrn(i,j,nt_bgc_chl+k-1)/zphin_N(i,j,k+1) 
                     else 
                         chlin(ij,k+1) = trcrn(i,j,nt_bgc_chl+k-1)/zphin_N(i,j,k+1) 
                     endif 
                endif
                if (tr_bgc_DMSPp) then
                     if (first_ice(i,j)) then
                         trcrn(i,j,nt_bgc_DMSPp+k-1) = dmspn(i,j)
                         DMSPpin(ij,k+1) = trcrn(i,j,nt_bgc_DMSPp+k-1)/zphin_N(i,j,k+1) 
                     else
                         DMSPpin(ij,k+1) = trcrn(i,j,nt_bgc_DMSPp+k-1)/zphin_N(i,j,k+1)    
                     endif
                endif
                if (tr_bgc_DMSPd) then
                     if (first_ice(i,j)) then
                         trcrn(i,j,nt_bgc_DMSPd+k-1) = dmspn(i,j)*initbio_frac    
                         DMSPdin(ij,k+1) = trcrn(i,j,nt_bgc_DMSPd+k-1)/zphin_N(i,j,k+1) 
                     else
                         DMSPdin(ij,k+1) = trcrn(i,j,nt_bgc_DMSPd+k-1)/zphin_N(i,j,k+1)   
                     endif
                endif
                if (tr_bgc_DMS) then
                     if (first_ice(i,j)) then
                         trcrn(i,j,nt_bgc_DMS+k-1) = dmsn(i,j)*initbio_frac    
                         DMSin(ij,k+1) = trcrn(i,j,nt_bgc_DMS+k-1)/zphin_N(i,j,k+1)  
                     else
                         DMSin(ij,k+1) = trcrn(i,j,nt_bgc_DMS+k-1)/zphin_N(i,j,k+1)    
                     endif   
                endif
                if (tr_bgc_PON) then
                     if (first_ice(i,j)) then
                         trcrn(i,j,nt_bgc_PON+k-1) =  nitn(i,j)*initbio_frac      
                         PONin(ij,k+1) = trcrn(i,j,nt_bgc_PON+k-1)/zphin_N(i,j,k+1)  
                     else
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
           !---------------------------------
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

           if (dh_top(i,j) +darcy_V(i,j)/zphi_min(i,j)*dt < c0 .AND. .NOT. flood_val(i,j)) then 
                NOin(ij,1)   = min_bgc*nitn(i,j)/zphi_min(i,j) 
                NHin(ij,1)   = min_bgc*ammn(i,j)/zphi_min(i,j)
                Silin(ij,1)  = min_bgc*siln(i,j)/zphi_min(i,j)
                DMSPdin(ij,1)= min_bgc*dmspn(i,j)/zphi_min(i,j)
                DMSin(ij,1)  = min_bgc*dmsn(i,j)/zphi_min(i,j)
                PONin(ij,1) =  NOin(ij,1) 
           elseif (dh_top(i,j) + darcy_V(i,j)/zphi_min(i,j)*dt < c0 .AND. flood_val(i,j)) then              
                NOin(ij,1)   = nitn(i,j)/zphi_min(i,j)
                NHin(ij,1)   = ammn(i,j)/zphi_min(i,j)
                Silin(ij,1)  = siln(i,j)/zphi_min(i,j)
                DMSPdin(ij,1)= dmspn(i,j)/zphi_min(i,j)
                DMSin(ij,1)  = dmsn(i,j)/zphi_min(i,j)
                PONin(ij,1) =  NOin(ij,1) 
           endif
           if (dh_top_chl(i,j) +darcy_V_chl(i,j)/zphi_min(i,j)*dt < c0 .AND. .NOT. flood_val(i,j)) then
                Nin(ij,1)    = min_bgc*algalNn(i,j)/zphi_min(i,j)
                Cin(ij,1)    = R_C2N* Nin(ij,1) 
                chlin(ij,1)  = R_chl2N * Nin(ij,1)
                DMSPpin(ij,1)= min_bgc*dmspn(i,j)/zphi_min(i,j)
           elseif (dh_top_chl(i,j) + darcy_V_chl(i,j)/zphi_min(i,j)*dt < c0 .AND. flood_val(i,j)) then  
                Nin(ij,1)    = algalNn(i,j)/zphi_min(i,j)
                Cin(ij,1)    = R_C2N* Nin(ij,1) 
                chlin(ij,1)  = R_chl2N * Nin(ij,1)
                DMSPpin(ij,1)= dmspn(i,j)/zphi_min(i,j)
           endif

           if ((dh_bot(i,j) + darcy_V(i,j)*dt > c0)) then 
               NOin(ij,nblyr_hist)   = nitn(i,j) 
               NHin(ij,nblyr_hist)   = ammn(i,j) 
               Silin(ij,nblyr_hist)  = siln(i,j) 
               DMSPdin(ij,nblyr_hist)= dmspn(i,j)
               DMSin(ij,nblyr_hist)  = dmsn(i,j)
               PONin(ij,nblyr_hist)  = nitn(i,j)  
            else  
               NOin(ij,nblyr_hist)   = NOin(ij,nblyr+1)
               NHin(ij,nblyr_hist)   = NHin(ij,nblyr+1)
               Silin(ij,nblyr_hist)  = Silin(ij,nblyr+1)
               DMSPdin(ij,nblyr_hist)= DMSPdin(ij,nblyr+1)
               DMSin(ij,nblyr_hist)  =  DMSin(ij,nblyr+1)
               PONin(ij,nblyr_hist)  = PONin(ij,nblyr+1) 
            endif
           if ((dh_bot_chl(i,j) + darcy_V_chl(i,j)*dt > c0)) then  
               Nin(ij,nblyr_hist)    = algalNn(i,j)  
               Cin(ij,nblyr_hist)    = R_C2N * Nin(ij,nblyr_hist)
               chlin(ij,nblyr_hist)  = R_chl2N * Nin(ij,nblyr_hist)
               DMSPpin(ij,nblyr_hist)= dmspn(i,j)
            else  
               Nin(ij,nblyr_hist)    = Nin(ij,nblyr+1)
               Cin(ij,nblyr_hist)    = R_C2N * Nin(ij,nblyr_hist)
               chlin(ij,nblyr_hist)  = R_chl2N * Nin(ij,nblyr_hist)
               DMSPpin(ij,nblyr_hist)= DMSPpin(ij,nblyr+1)
            endif
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
        
        call remap_layers_bgc_plus (nx_block,ny_block, &
                             indxii,   indxjj,         &
                             icells,                   &
                             ntrcr,                    &
                             nilyr,                    &
                             nt_zfswin,                &
                             trtmp0,    trtmp,         &
                             0,        nblyr,          &
                             hin, hinS,                &
                             cgrid(2:nilyr+1),         &
                             bgrid(2:nblyr+1), surface_S)

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
 
      do k = 1, nblyr_hist
         do ij = 1, icells    
            i = indxii(ij)
            j = indxjj(ij)!!
                           ! save initial values
            in_init(ij,k,nlt_bgc_NO) = NOin(ij,k)        
            biomat(ij,k,nlt_bgc_NO)  = NOin(ij,k)

            if (tr_bgc_NH) then
               in_init(ij,k,nlt_bgc_NH) = NHin(ij,k)      
               biomat(ij,k,nlt_bgc_NH)  = NHin(ij,k)      
            endif
            if (tr_bgc_N) then
               in_init(ij,k,nlt_bgc_N)  = Nin(ij,k)      
               biomat(ij,k,nlt_bgc_N)   = Nin(ij,k)      
            endif
            if (tr_bgc_Sil) then
               in_init(ij,k,nlt_bgc_Sil)= Silin(ij,k)  
               biomat(ij,k,nlt_bgc_Sil) = Silin(ij,k)      
            endif
            if (tr_bgc_C) then
               in_init(ij,k,nlt_bgc_C)  = Cin(ij,k)  
               biomat(ij,k,nlt_bgc_C)   = Cin(ij,k)      
             endif
            if (tr_bgc_chl) then
               in_init(ij,k,nlt_bgc_chl)= chlin(ij,k)      
               biomat(ij,k,nlt_bgc_chl) = chlin(ij,k)      
            endif
            if (tr_bgc_DMSPp) then
                in_init(ij,k,nlt_bgc_DMSPp)= DMSPpin(ij,k)   
                biomat(ij,k,nlt_bgc_DMSPp) = DMSPpin(ij,k)      
            endif
            if (tr_bgc_DMSPd) then
                in_init(ij,k,nlt_bgc_DMSPd)= DMSPdin(ij,k) 
                biomat(ij,k,nlt_bgc_DMSPd) = DMSPdin(ij,k)      
            endif
            if (tr_bgc_DMS) then
               in_init(ij,k,nlt_bgc_DMS) = DMSin(ij,k)  
               biomat(ij,k,nlt_bgc_DMS)  = DMSin(ij,k)           
            endif
            if (tr_bgc_PON) then
              in_init(ij,k,nlt_bgc_PON) = PONin(ij,k) 
              biomat(ij,k,nlt_bgc_PON)  = PONin(ij,k) 
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

     if (solve_zbgc) then

       do k = 2, nblyr+1   
       
         call algal_dyn &
                        (nx_block, ny_block,           &
                         icells,                       &
                         indxii,    indxjj,            &
                         zfswin(:,:,k), react(:,k-1,:),& 
                         in_init(:,k,:),ntraceb,       &
                         growN(:,:,k-1),               &
                         upNOn(:,:,k-1),               &
                         upNHn(:,:,k-1),               &
                         tr_bgc_N,tr_bgc_NO,           &
                         tr_bgc_NH,tr_bgc_Sil,         &
                         tr_bgc_C,tr_bgc_chl,          &
                         tr_bgc_DMSPp, tr_bgc_DMSPd,   & 
                         tr_bgc_DMS, tr_bgc_PON)
          
          
         do ij = 1, icells    
            i = indxii(ij)
            j = indxjj(ij)
            dh_add = max(c0,dh_bot(i,j)-dh_bot_chl(i,j))
            react(ij,k-1, nlt_bgc_N) = react(ij,k-1,nlt_bgc_N) + &
                       dh_add*algalNn(i,j)/zphin_N(i,j,k)/ &
                       hinS(ij)*(igrid(k)-igrid(k-1))      
         enddo     !ij
       enddo       !k

      endif        !solve_zbgc

      !-----------------------------------------------------------------
      ! Compute elements of tridiagonal matrix for each tracer
      !-----------------------------------------------------------------
      do mm = 1, nbltrcr

       if (bgc_tracer_type(mm) > p5) then  !c1
          
            call get_matrix_elements_calc_bgc &
                                     (nx_block, ny_block,        &
                                      icells,                    &
                                      indxii,    indxjj,         &
                                      igrid, bgrid,              &
                                      dt,     zphin_N,           &
                                      iDin,   iphin_N,           &
                                      in_init(:,:,mm), hinS_old, &
                                      hinS,                      &
                                      sbdiag,   diag,            &
                                      spdiag,   rhs, darcy_V,    &
                                      react(:,:,mm),             &
                                      dh_top,dh_bot,             &
                                      thin_check,mm)
      
         !------------------------------------------------------------------------
         !
         !  Use molecular diffusion only for DMSPp and N
         !------------------------------------------------------------------------

        elseif (bgc_tracer_type(mm) < p5) then  !c0

            call get_matrix_elements_calc_bgc &
                                     (nx_block, ny_block,         &
                                      icells,                     &
                                      indxii,    indxjj,          &
                                      igrid, bgrid,               &
                                      dt,     zphin_N,            &
                                      iDin_N,   iphin_N,          &
                                      in_init(:,:,mm), hinS_old,  &
                                      hinS,                       &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs,  darcy_V_chl,&
                                      react(:,:,mm)            ,  &
                                      dh_top_chl,dh_bot_chl,      &
                                      thin_check,mm)

        endif
        
      !-----------------------------------------------------------------
      ! Solve tridiagonal matrix to obtain the new tracers
      !-----------------------------------------------------------------
         call tridiag_solver (icells,                  &
                              nblyr_hist, sbdiag,      &
                              diag,     spdiag,        &
                              rhs,      biomat(:,:,mm))

        do ij = 1, icells  
          i = indxii(ij)
          j = indxjj(ij) 
         
          bulk_top = biomat(ij,1,mm)*zphin_N(i,j,1)
          bulk_bot = biomat(ij,nblyr_hist,mm)*zphin_N(i,j,nblyr+1)
          if (dh_bot(i,j) < c0) then
            bulk_bot = biomat(ij,nblyr+1,mm)*zphin_N(i,j,nblyr+1)
          endif
          dC_dx_bot = (biomat(ij,nblyr_hist,mm) - biomat(ij,nblyr+1,mm)) &
                           /(bgrid(nblyr_hist)-bgrid(nblyr+1))
          Dm_bot = bgc_tracer_type(mm)*iDin(i,j,nblyr+1) + &
                   (c1-bgc_tracer_type(mm))*iDin_N(i,j,nblyr+1)
          dh_bot_all = bgc_tracer_type(mm)*dh_bot(i,j) + &
                (c1-bgc_tracer_type(mm))*dh_bot_chl(i,j)
          flux_bio_g(i,j,mm) = flux_bio_g(i,j,mm) -thin_check(ij)*Dm_bot*hinS_old(ij)*dC_dx_bot
          flux_bio(i,j,mm) = flux_bio(i,j,mm) + thin_check(ij)*(flux_bio_g(i,j,mm)+&
                    (dh_bot_all*bulk_bot-dh_top(i,j)*bulk_top)/dt)   
          if (hinS_old(ij) < thin)&
                Call thin_ice_flux(hinS(ij),hinS_old(ij),zphin_N(i,j,:),in_init(ij,:,mm), &
                            flux_bio(i,j,mm), igrid, dt)
          !-------------------------------
          ! update ocean concentration
          ! Not coupled to ocean biogeochemistry
          !-------------------------------------
          !ocean_bio(i,j,mm) = ocean_bio(i,j,mm) + flux_bio(i,j,mm)/hmix(i,j)*aicen(i,j)   
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
              do mm = 1, nbltrcr
                 bulk_biomat(ij,k-1,mm)  = biomat(ij,k,mm)*zphin_N(i,j,k)
                 bulk_biomat_o(ij,k-1,mm)= in_init(ij,k,mm)*zphin_N(i,j,k)
                 react_phi(ij,k-1,mm)    = react(ij,k-1,mm)*zphin_N(i,j,k)
              enddo
              if ( bulk_biomat(ij,k-1,nlt_bgc_NO) > 1000.0_dbl_kind) then              
                   write(nu_diag,*)'Bulk Nitrate solution error'      
                   write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,mm:'&
                                   ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                    TLON(i,j)*rad_to_deg,istep1,mm
                   write(nu_diag,*)'dh_bot_chl(i,j),dh_top(i,j)'&
                                   ,dh_bot_chl(i,j),dh_top(i,j) 
                   write(nu_diag,*)'ij,k,nlt_bgc_NO'
                   write(nu_diag,*) ij,k,nlt_bgc_NO
                   write(nu_diag,*)'biomat(ij,mm,nlt_bgc_NO),mm = 1,nblyr+2'
                   write(nu_diag,*)(biomat(ij,mm,nlt_bgc_NO),mm = 1,nblyr+2)
                   write(nu_diag,*)'in_init(ij,k,nlt_bgc_NO),zphin_N(i,j,k),nbltrcr,nbltrcr'
                   write(nu_diag,*) in_init(ij,k,nlt_bgc_NO),zphin_N(i,j,k),nbltrcr,nbltrcr 
                   write(nu_diag,*)'iphin_N(i,j,mm),mm =1,nblyr+1' 
                   write(nu_diag,*)(iphin_N(i,j,mm),mm =1,nblyr+1)
                   write(nu_diag,*) 'zphin_N(i,j,mm),mm=1,nblyr+2'
                   write(nu_diag,*) (zphin_N(i,j,mm),mm=1,nblyr+2)
                        if (ktherm == 1 ) then
                          write(nu_diag,*)'trcrn(i,j,nt_bgc_S-1+mm),mm=1,nblyr+2'
                          write(nu_diag,*) ( trcrn(i,j,nt_bgc_S-1+mm),mm=1,nblyr+2)
                        elseif (ktherm == 2) then
                          write(nu_diag,*)'trcrn(i,j,nt_sice-1+mm),mm=1,nilyr+2'
                          write(nu_diag,*) ( trcrn(i,j,nt_sice-1+mm),mm=1,nilyr+2)
                        endif
                   write(nu_diag,*) 'Sin: trcrn(i,j,nt_bgc_S-1+mm),mm=1,nblyr+2'
                   write(nu_diag,*) ( trcrn(i,j,nt_bgc_S-1+mm),mm=1,nblyr+2)
                   write(nu_diag,*)'hinS(ij),hinS_old(ij),hin_old(ij),hin(ij)' 
                   write(nu_diag,*)hinS(ij),hinS_old(ij),hin_old(ij),hin(ij)
                   call abort_ice ('ice_bgc.F90: BGC error')
              endif   
              NOin(ij,k) = max(biomat(ij,k,nlt_bgc_NO),c0)
              if (tr_bgc_N)    Nin(ij,k)     = max(biomat(ij,k,nlt_bgc_N),c0) 
              if (tr_bgc_NH)   NHin(ij,k)    = max(biomat(ij,k,nlt_bgc_NH),c0)
              if (tr_bgc_Sil)  Silin(ij,k)   = max(biomat(ij,k,nlt_bgc_Sil),c0)
              if (tr_bgc_DMSPp)DMSPpin(ij,k) = max(biomat(ij,k,nlt_bgc_DMSPp),c0)
              if (tr_bgc_DMSPd)DMSPdin(ij,k) = max(biomat(ij,k,nlt_bgc_DMSPd),c0) 
              if (tr_bgc_DMS)  DMSin(ij,k)   = max(biomat(ij,k,nlt_bgc_DMS),c0)
              if (tr_bgc_PON)  PONin(ij,k)   = max(biomat(ij,k,nlt_bgc_PON),c0)   
         enddo       ! ij
      enddo          ! k

      !-----------------------------------------------------------------
      ! Check flux conservation of tracers
      !-----------------------------------------------------------------
      check_conserve = .false.
      good_bio_numerics(:,:) = .true.
      if (check_conserve) then        
        do  mm = 1, nbltrcr
           do ij = 1, icells 
              i = indxii(ij)
              j = indxjj(ij)

              if (thin_check(ij) > c0) call check_conserve_bgc &
                                     (hinS(ij), hinS_old(ij), bulk_biomat(ij,1:nblyr,mm), &
                                      bulk_biomat_o(ij,1:nblyr,mm), dt, &
                                      igrid,good_bio_numerics(ij,mm),flux_bio(i,j,mm),&
                                      react_phi(ij,:,mm),aicen(i,j)) 
     
              if (.NOT. good_bio_numerics(ij,mm)) then
                write(nu_diag,*)'1st Calc_bgc_fluxes: Category,i,j,ij,TLAT,TLON,'
                write(nu_diag,*)'istep1,mm, nlt_bgc_NO:'&
                 ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                TLON(i,j)*rad_to_deg,istep1,mm, nlt_bgc_NO
                write(nu_diag,*)'hinS(ij),hinS_old(ij),aicen(i,j)'
                write(nu_diag,*) hinS(ij),hinS_old(ij),aicen(i,j)
                write(nu_diag,*)'dh_bot_chl(i,j),dh_top(i,j)'&
                ,dh_bot_chl(i,j),dh_top(i,j)
                write(nu_diag,*)'iDin_N(i,j,nblyr+1)'
                write(nu_diag,*)iDin_N(i,j,nblyr+1)
                write(nu_diag,*)'flux_bio(i,j,mm)'
                write(nu_diag,*) flux_bio(i,j,mm)
                write(nu_diag,*)'react_phi(ij,m,mm),m = 1,nblyr'
                write(nu_diag,*)(react_phi(ij,m,mm),m = 1,nblyr)
                write(nu_diag,*)'bulk_biomat(ij,m,mm),m = 1,nblyr'
                write(nu_diag,*)(bulk_biomat(ij,m,mm),m = 1,nblyr)
                write(nu_diag,*)'bulk_biomat_o(ij,m,mm),m = 1,nblyr'
                write(nu_diag,*)(bulk_biomat_o(ij,m,mm),m = 1,nblyr)
                call abort_ice ('ice: calc_bgc_fluxes ')
              endif
                         
           enddo  !ij
        enddo   !mm
      endif     !check_conserve
      
      !-----------------------------------------------------------------
      ! Update the tracer variable
      !-----------------------------------------------------------------
    
       do k = 1,nblyr                  !back to bulk quantity
         do ij = 1, icells 
            i = indxii(ij)
            j = indxjj(ij)
            if (hinS_old(ij) > thin) then
              trcrn(i,j,nt_bgc_NO+k-1)                     =  NOin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_N)    trcrn(i,j,nt_bgc_N+k-1)     =  Nin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_NH)   trcrn(i,j,nt_bgc_NH+k-1)    =  NHin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_C)    trcrn(i,j,nt_bgc_C+k-1)     =  R_C2N * Nin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_chl)  trcrn(i,j,nt_bgc_chl+k-1)   =  R_chl2N * Nin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_sil)  trcrn(i,j,nt_bgc_Sil+k-1)   =  Silin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_DMSPp)trcrn(i,j,nt_bgc_DMSPp+k-1) =  DMSPpin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_DMSPd)trcrn(i,j,nt_bgc_DMSPd+k-1) =  DMSPdin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_DMS)  trcrn(i,j,nt_bgc_DMS+k-1)   =  DMSin(ij,k+1)*zphin_N(i,j,k+1)
              if (tr_bgc_PON)  trcrn(i,j,nt_bgc_PON+k-1)   =  PONin(ij,k+1)*zphin_N(i,j,k+1)
            endif
         enddo           !ij
      enddo        !k   
   
770 format (I6,D16.6)        
781 format (I6,I6,I6)
790 format (I6,I6)
791 format (f24.17)
792 format (2D16.6)
793 format (3D16.6)
794 format (4D15.5)
800 format (F10.4)

      end subroutine z_biogeochemistry

!=======================================================================
!
! Do biogeochemistry from subroutine algal_dynamics
! by Scott Elliott: updated to 175
! 
      subroutine algal_dyn (nx_block, ny_block,       &
                                 icells,              &
                                 indxi,    indxj,     &
                                 fswthrul, reactb,    &
                                 ltrcrn, ntr, growN,  &
                                 upNOn, upNHn,        &
                                 tr_bio_N,tr_bio_NO,  &
                                 tr_bio_NH,tr_bio_Sil,&
                                 tr_bio_C,tr_bio_chl, &
                                 tr_bio_DMSPp,        &
                                 tr_bio_DMSPd,        &
                                 tr_bio_DMS, PON_flag)

      use ice_calendar, only: dt

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of cells with aicen > puny
         ntr                   ! number of layer tracers

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fswthrul   ! average shortwave passing through current ice layer (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         growN,  &  !  algal specific growth rate   (/s)
         upNOn,  &  !  algal NO uptake  rate   (mmol/m^3/s)
         upNHn      !  algal NH uptake rate    (mmol/m^3/s)

      real (kind=dbl_kind), dimension(icells,ntr), &
         intent(inout) :: &
         reactb     ! biological reaction terms (mmol/m3)

      real (kind=dbl_kind), dimension(icells,ntr), &
         intent(in) :: &
         ltrcrn     ! concentrations  in layer

      logical (kind=log_kind), intent(in):: &
         tr_bio_N, tr_bio_NO  , &   ! tracer flags for z_bgc or skl_bgc: 
         tr_bio_NH, tr_bio_Sil, &   ! algal nitrogen, nitrate, ammonium, silicate
         tr_bio_C, tr_bio_chl , &   ! algal carbon, algal chlorophyll, DMSP particulate
         tr_bio_DMSPp, tr_bio_DMSPd, & !DMSP dissolved, DMS
         tr_bio_DMS

      logical (kind=log_kind), intent(in), optional:: &
         PON_flag ! if true, add an extra tracer
                  ! currently for z_bgc only.

      !  local variables

      logical (kind=log_kind):: &
         tr_bio_PON ! carries a non-reactive nitrate tracer for testing physics parametrizations 

      real (kind=dbl_kind), parameter :: & 
         T_bot      = -1.8_dbl_kind , & ! interface close to freezing (C)
         chlabs     = 9.e-4_dbl_kind, & ! chlorophyll absorption (1/(mg/m^3)) where a 3 cm 
                                        ! skeletal layer thickness is assumed
         mu_max     = 1.5_dbl_kind,   & ! (1/day) maximum growth rate 'Jin2006'     
         T_max      = -1.8_dbl_kind,  & ! (C) maximum growth at Tmax
         op_dep_min = 0.1          ,  & ! 0.01 Light attenuates for optical depths exceeding
                                        ! op_dep_min (small effect unless chla > 100 mg/m^2)
         grow_Tdep  = 0.0633_dbl_kind,& ! (1/C)and its T dependence
         fr_graze   = p1           ,  & ! A93 val for S, set to zero in Jin 06
         fr_graze_s = 0.5_dbl_kind ,  & ! fraction of grazing spilled or slopped
         fr_graze_a = 0.5_dbl_kind ,  & ! fraction of grazing assimilated
         fr_graze_e = 0.5_dbl_kind ,  & ! fraction of assimilation excreted  
         alpha2max  = 0.8_dbl_kind,   & ! light limitation (1/(W/m^2))
         beta2max   = 0.018_dbl_kind, & ! corresponding light inhibition (1/W/m^2)
         K_Nit      = c1           ,  & ! nitrate half saturation (mmol/m^3) 
         K_Am       = c1           ,  & ! ammonium half saturation (mmol/m^3) 
         K_Sil      = 4.0_dbl_kind ,  & ! silicon half saturation (mmol/m^3)
         mort_pre   = 0.0208_dbl_kind,& ! (1/day) prefix to mortality
         mort_Tdep  = 0.03_dbl_kind , & ! (1/C) T dependence of mortality
         fr_mort2min= c1            , & ! fractionation to remin
         t_nitrif   = 67.0_dbl_kind , & ! (days) nitrification time scale
         max_loss   = 0.9               ! restrict uptake to 90% of remaining value 

      real (kind=dbl_kind), parameter :: &
         fr_excrt_2S= c1           , &  ! excretion is efficient initially
         y_sk_DMS   = c1           , &  ! and conversion given high yield
         t_sk_conv  = 10.0_dbl_kind, &  ! at a Stefels rate (days)
         t_sk_ox    = 10.0_dbl_kind     ! DMS in turn oxidizes slowly (days)

     ! real (kind=dbl_kind), parameter :: &
     !    pr_l       = 10.0_dbl_kind, & ! product layer thickness (m) 
     !    chl_pr_v   = 0.1_dbl_kind , & ! fixed nondiatom chlorophyll in ml (mg/m^3)
     !    R_chl2N_nd = 3.0_dbl_kind , & ! shade adaptation below (mg/millimole)
     !    R_C2N_nd   = 7.0_dbl_kind , & ! open water ratio (mole/mole)
     !    t_pr_dsrp  = 10.0_dbl_kind    ! disruption time scale (days)
         
     ! real (kind=dbl_kind), parameter :: &
     !    R_S2N_nd   = 0.03_dbl_kind, & ! open water ratio nondiatoms (mole/mole)
     !    y_pr_DMS   = c1           , & ! but we begin again with unit yield
     !    t_pr_conv  = 10.0_dbl_kind, & ! and a similar conversion (days)
     !    t_pr_ox    = 10.0_dbl_kind    ! plus round final time scale (days)

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

      real (kind=dbl_kind) :: &
         op_dep     , &  ! bottom layer attenuation exponent (optical depth)
         Iavg_loc        ! bottom layer attenuated Fswthrul (W/m^2)

      real (kind=dbl_kind) :: &
         L_lim    , &  ! overall light limitation 
         Nit_lim  , &  ! overall nitrate limitation
         Am_lim   , &  ! overall ammonium limitation
         N_lim    , &  ! overall nitrogen species limitation
         Sil_lim  , &  ! overall silicon limitation
         fr_Nit   , &  ! fraction of local ecological growth as nitrate
         fr_Am    , &  ! fraction of local ecological growth as ammonia
         growmax_N, &  ! maximum growth rate in N currency (mmol/m^3/s)
         grow_N   , &  ! true growth rate in N currency (mmol/m^3/s)
         potU_Nit , &  ! potential nitrate uptake (mmol/m^3/s)
         potU_Am  , &  ! potential ammonium uptake (mmol/m^3/s)
         U_Nit    , &  ! actual nitrate uptake (mmol/m^3/s)
         U_Am     , &  ! actual ammonium uptake (mmol/m^3/s)
         U_Sil         ! actual silicon uptake (mmol/m^3/s)

      real (kind=dbl_kind) :: &
         resp     , &  ! respiration (mmol/m^3/s)
         graze    , &  ! grazing (mmol/m^3/s)
         mort     , &  ! sum of mortality and excretion (mmol/m^3/s)
         nitrif        ! nitrification (mmol/m^3/s)

!  source terms underscore s, removal underscore r

      real (kind=dbl_kind) :: &
         N_s_p     , &  ! algal nitrogen photosynthesis (mmol/m^3)
         N_r_g     , &  ! algal nitrogen losses to grazing (mmol/m^3)
         N_r_r     , &  ! algal nitrogen losses to respiration (mmol/m^3)
         N_r_mo    , &  ! algal nitrogen losses to mortality (mmol/m^3)
         N_s       , &  ! net algal nitrogen sources (mmol/m^3)
         N_r       , &  ! net algal nitrogen removal (mmol/m^3)
         C_s       , &  ! net algal carbon sources (mmol/m^3)
         C_r       , &  ! net algal carbon removal (mmol/m^3)
         NO_s_n    , &  ! skl nitrate from nitrification (mmol/m^3)
         NO_s_r    , &  ! skl nitrate from respiration (mmol/m^3)
         NO_r_p    , &  ! skl nitrate uptake by algae (mmol/m^3)
         NO_s      , &  ! net skl nitrate sources (mmol/m^3)
         NO_r      , &  ! net skl nitrate removal (mmol/m^3)
         NH_s_e    , &  ! skl ammonium source from excretion (mmol/m^3)
         NH_s_r    , &  ! skl ammonium source from respiration (mmol/m^3)
         NH_s_mo   , &  ! skl ammonium source from mort/remin (mmol/m^3) 
         NH_r_p    , &  ! skl ammonium uptake by algae (mmol/m^3)
         NH_r_n    , &  ! skl ammonium removal to nitrification (mmol/m^3)
         NH_s      , &  ! net skl ammonium sources (mmol/m^3)
         NH_r      , &  ! net skl ammonium removal (mmol/m^3)
         Sil_s_r   , &  ! skl silicon from respiration (mmol/m^3)
         Sil_r_p   , &  ! skl silicon uptake by algae (mmol/m^3)
         Sil_s     , &  ! net skl silicon sources (mmol/m^3)
         Sil_r          ! net skl silicon removal (mmol/m^3)

      real (kind=dbl_kind) :: &
         DMSPd_s_s , &  ! skl dissolved DMSP from grazing spillage (mmol/m^3)
         DMSPd_s_e , &  ! skl dissolved DMSP from zooplankton excretion (mmol/m^3)
         DMSPd_s_mo, &  ! skl dissolved DMSP from MBJ algal mortexc (mmol/m^3)
         DMSPd_r_c , &  ! skl dissolved DMSP conversion (mmol/m^3)
         DMSPd_s   , &  ! net skl dissolved DMSP sources (mmol/m^3)
         DMSPd_r   , &  ! net skl dissolved DMSP removal (mmol/m^3)
         DMS_s_c   , &  ! skl DMS source via conversion (mmol/m^3)
         DMS_r_o   , &  ! skl DMS losses due to oxidation (mmol/m^3)
         DMS_s     , &  ! net skl DMS sources (mmol/m^3)
         DMS_r     , &  ! net skl DMS removal (mmol/m^3)
         PON_s_z   , &  ! PON source as zooplankton (mmol/m^3)
         PON_s_d   , &  ! PON source as detritus (mmol/m^3)
         PON_s     , &  ! net  PON sources (mmol/m^3)
         PON_r          ! net  PON removal (mmol/m^3)

     ! real (kind=dbl_kind) :: &
     !    DMSP_pr_s_nd , &  ! product layer dissolved DMSP from local bio (mmol/m^3)
     !    DMSP_pr_s_me , &  ! product layer dissolved DMSP from melting (mmol/m^3)
     !    DMSP_pr_r_c  , &  ! product layer dissolved DMSP conversion (mmol/m^3)
     !    DMSP_pr_s    , &  ! net product dissolved DMSP sources (mmol/m^3)
     !    DMSP_pr_r    , &  ! net product dissolved DMSP removal (mmol/m^3)
     !    DMS_pr_s_c   , &  ! product layer DMS source via conversion (mmol/m^3)
     !    DMS_pr_r_o   , &  ! product layer DMS losses due to oxidation (mmol/m^3)
     !    DMS_pr_s     , &  ! net product DMS sources (mmol/m^3)
     !    DMS_pr_r          ! net product DMS removal (mmol/m^3)

      logical (kind=log_kind) :: &   
         write_flag           ! set to true at each timestep        

!  begin building biogeochemistry terms

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

      write_flag = .true.
      
      if (present(PON_flag)) then
         tr_bio_PON = PON_flag
      else
         tr_bio_PON = .false.
      endif

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

         NOin     = ltrcrn(ij,nlt_bgc_NO)
         if (tr_bio_N)         Nin      = ltrcrn(ij,nlt_bgc_N)
         if (tr_bio_C)         Cin      = ltrcrn(ij,nlt_bgc_C)
         if (tr_bio_NH)        NHin     = ltrcrn(ij,nlt_bgc_NH)
         if (tr_bio_Sil)       Silin    = ltrcrn(ij,nlt_bgc_Sil)
         if (tr_bio_DMSPp)     DMSPpin  = ltrcrn(ij,nlt_bgc_DMSPp)
         if (tr_bio_DMSPd)     DMSPdin  = ltrcrn(ij,nlt_bgc_DMSPd)
         if (tr_bio_DMS)       DMSin    = ltrcrn(ij,nlt_bgc_DMS)
         if (tr_bio_PON)       PONin    = ltrcrn(ij,nlt_bgc_PON)
         chlin = R_chl2N * Nin  

!----------------------
! Light limitation
!----------------------- 

        op_dep   = chlabs * chlin   

!  Okhotsk maxima causes a reevaluation.
!  The new concept is, late algae at the bottom of the bottom strongly attenuated.
!  Since k about 0.03 1/m(mg/m3), efold at about 30 mg/m3
!  values of order hundreds will be shut down...
!  Newest algae experience super aborption above because they sit low.
!  More than perhaps two efolds and light falls below half value.

         if (op_dep > op_dep_min) then
           Iavg_loc   = fswthrul(i,j)*  (c1 - exp(-op_dep)) / op_dep
         else
           Iavg_loc   = fswthrul(i,j)
         endif

!  With light inhibition
!           
!       L_lim     = (c1 - exp(-alpha2max*Iavg_loc)) * exp(-beta2max*Iavg_loc)      
!
!  Without light inhibition
!           
        L_lim     = (c1 - exp(-alpha2max*Iavg_loc)) 

!--------------------------- 
! Nutrient limitation
!----------------------------

        Nit_lim =   NOin/(NOin + K_Nit)
        Am_lim = c0
        if (tr_bio_NH) then
            Am_lim    = NHin/(NHin + K_Am)
            N_lim     = min(c1,Nit_lim + Am_lim)  
        else
            N_lim = Nit_lim
        endif
        Sil_lim = c1
        if (tr_bio_Sil)  Sil_lim   = Silin/(Silin + K_Sil)
 
!  Growth and uptake computed within the bottom layer 
!  Note here per A93 discussions and MBJ model, salinity is a universal restriction
!  Comparison with available column nutrients inserted but in tests had no effect
!  Primary production reverts to SE form, see MBJ below and be careful

       growmax_N   = mu_max / secday * exp(grow_Tdep * (T_bot - T_max))* Nin 
       grow_N   = min(L_lim,N_lim,Sil_lim) * growmax_N
       potU_Nit = Nit_lim          * growmax_N
       potU_Am  = Am_lim           * growmax_N 
       U_Am     = min(grow_N, potU_Am)
       U_Nit    = grow_N - U_Am
       U_Sil    = R_Si2N * grow_N

       if (tr_bio_Sil) U_Sil = min(U_Sil, max_loss * Silin/dt)
       U_Nit = min(U_Nit, max_loss * NOin/dt)  
       U_Am  = min(U_Am, max_loss * NHin/dt)    
 
       grow_N = min(U_Sil/R_Si2N,U_Nit + U_Am)
       fr_Am = c0
       if (tr_bio_NH) then
          fr_Am = p5
          if (grow_N > c0) fr_Am = min(U_Am/grow_N,c1)
       endif
       fr_Nit = c1-fr_Am
       U_Nit  = fr_Nit * grow_N
       U_Am   = fr_Am  * grow_N
       U_Sil  = R_Si2N * grow_N
    
       resp        = fr_resp* grow_N
       graze       = fr_graze*grow_N
       mort        = mort_pre * exp(mort_Tdep*(T_bot-T_max)) * Nin / secday 
       nitrif      = c0 ! (NHin / t_nitrif) / secday
 
! history variables

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
        N_r_g = graze * dt 
        N_r_r = resp * dt
        N_r_mo= mort * dt
        N_s = N_s_p
        N_r = N_r_g + N_r_r + N_r_mo 

!  Carbon chemistry 
!       C_react = R_C2N * N_react 

        C_s = R_C2N * N_s
        C_r = R_C2N * N_r

! Nitrate reaction term
!       NO_react = (nitrif - fr_Nit*grow_N)*dt

        NO_s_n = nitrif * dt               
        NO_r_p = U_Nit * dt                
        NO_s = NO_s_n
        NO_r = NO_r_p

! Ammonium reaction term
!        NH_react = (-nitrif - (c1-fr_Nit - fr_resp - &
!           fr_graze*fr_graze_e*fr_graze_a)*grow_N + mort*fr_mort2min)*dt  

        NH_s_r =  N_r_r 
        NH_s_e = fr_graze_e*fr_graze_a* N_r_g
        NH_s_mo= fr_mort2min * N_r_mo
        NH_r_p = U_Am * dt                
        NH_r_n = nitrif *dt                  
        NH_s = NH_s_r + NH_s_e + NH_s_mo
        NH_r = NH_r_p + NH_r_n 

! Silica  reaction term
!        Sil_react = - R_Si2N * grow_N * dt
     
        Sil_r_p = U_Sil * dt
        Sil_s = c0
        Sil_r = Sil_r_p 

!  Sulfur cycle begins here
!  Grazing losses are channeled into rough spillage and assimilation
!  then onward and the MBJ mortality channel is included
!  It is assumed as a starting expedient that 
!  DMSP loss to melting gives partial conversion to DMS in product layer
!  which then undergoes Stefels removal.


! DMSPd  reaction term
!        DMSPd_react = R_S2N*((fr_graze_s+fr_excrt_2S*fr_graze_e*fr_graze_a)*fr_graze*grow_N + &
!                  fr_mort2min*mort)*dt - [\DMSPd]/t_sk_conv*dt

        DMSPd_s_s = fr_graze_s * R_S2N * N_r_g
        DMSPd_s_e = fr_excrt_2S * fr_graze_e * fr_graze_a * R_S2N * N_r_g
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
       !  DMSP_pr_f   = F_DMSP * dt
       !  DMSP_pr_s = DMSP_pr_s_nd                ! + DMSP_pr_s_me + DMSP_pr_f
       !  DMSP_pr_r = DMSP_pr_r_c

       !  DMS_pr_s_c = y_pr_DMS * DMSP_pr_r_c
       !  DMS_pr_r_o = (c1/t_pr_ox) * (c1/secday) * (dms(i,j,iblk)*pr_l) * dt
       !  DMS_pr_f   = F_DMS * dt 
       !  DMS_pr_s = DMS_pr_s_c                   ! + DMS_pr_f
       !  DMS_pr_r = DMS_pr_r_o



! Add PON 
!        currently using PON to shadow nitrate without reactions
!        PON_react = N_r_g *( c1- fr_graze_e*fr_graze_a)

         PON_s_z = N_r_g * fr_graze_a 
         PON_s_d = N_r_g * (c1 - fr_graze_a)
         PON_s   = PON_s_z + PON_s_d
         PON_r   = N_r_g * fr_graze_a* fr_graze_e
                 
! Define reaction terms

         reactb(ij,nlt_bgc_NO)                    = NO_s  - NO_r
         if (tr_bio_N) reactb(ij,nlt_bgc_N)       = N_s   - N_r 
         if (tr_bio_C) reactb(ij,nlt_bgc_C)       = C_s    - C_r
         if (tr_bio_NH)reactb(ij,nlt_bgc_NH)      = NH_s   - NH_r
         if (tr_bio_Sil)reactb(ij,nlt_bgc_Sil)    = Sil_s  - Sil_r
         if (tr_bio_DMSPd)reactb(ij,nlt_bgc_DMSPd)= DMSPd_s- DMSPd_r
         if (tr_bio_DMS)reactb(ij,nlt_bgc_DMS)    = DMS_s  - DMS_r
         if (tr_bio_PON)reactb(ij,nlt_bgc_PON)    =  c0 !PON_s - PON_r

         !if (growN(i,j) > c0 .AND. write_flag) then
         !       write(nu_diag, *) 'reactb(ij,nlt_bgc_N), ij:',reactb(ij,nlt_bgc_N),ij
         !       write(nu_diag, *) 'growN(i,j), i,j:',growN(i,j),i,j
         !       write(nu_diag, *) 'fswthrul(i,j):',fswthrul(i,j)
         !       write(nu_diag, *) 'Nin,NOin,NHin,Silin:',Nin,NOin,NHin,Silin
         !       write_flag = .false.
         !endif

      enddo


      end subroutine algal_dyn
!
!=======================================================================
!
! Compute terms in tridiagonal matrix that will be solved to find
!  the bgc vertical profile
!
! November 2008 by N. Jeffery, modified for bgc layer model
!
      subroutine get_matrix_elements_calc_bgc &
                                     (nx_block, ny_block,         &
                                      icells,                     &
                                      indxi,   indxj,             &
                                      igrid, bgrid,               &
                                      dt,     zphin_N,            &
                                      iDin,    iphin_N,           &
                                      NOin_init, hin_old,         &
                                      hin,                        &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs,  darcy_V,    &
                                      reactb,  dht,dhb,           &
                                      thin_check,mm)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells ,mm            ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj  ! compressed indices 

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(in) :: &
         zphin_N         ! Porosity of layers where zphin >= zphimin

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hin         , & ! brine thickness (m)                  (m)
         hin_old         ! brine thickness before growth/melt (m)

      real (kind=dbl_kind), dimension (icells), intent(inout) :: &
         thin_check      ! c0 if thin and c1 if > thin       

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(in) :: &
         iphin_N     , & ! Porosity with min condition on igrid
         iDin            ! Diffusivity on the igrid (m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dht         , & ! Change in brine top (m)
         dhb         , & ! Change in brine bottom (m)
         darcy_V         ! Darcy velocity (m/s)

      real (kind=dbl_kind), dimension (nblyr_hist), intent(in) :: &
         bgrid           ! biology nondimensional grid layer points 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid           ! biology interface grid points 

      real (kind=dbl_kind), dimension (icells,nblyr), intent(in) :: &
         reactb          ! Reaction terms from  algal_dyn (mmol/m^3)

      real (kind=dbl_kind), dimension (icells,nblyr_hist), &
         intent(in) :: &
         NOin_init       ! bgc concentration (mmol/m^3) at start of dt 

      real (kind=dbl_kind), dimension (icells,nblyr_hist), &
         intent(out) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      ! local variables

      real (kind=dbl_kind), dimension (icells,nblyr) :: &
         vel              ! advective velocity times dt (m)

      real (kind=dbl_kind), dimension (icells,nblyr+1) :: &
         ivel             ! advective velocity*dt on igrid (m)
   
      real (kind=dbl_kind), dimension (icells,nblyr_hist) :: &
         bgrid_temp       ! biology nondimensional grid layer points 

      real (kind=dbl_kind) :: &
         rhs_sp, rhs_sb,& ! components of the rhs tri-diagonal matrix:
                          ! sp is super-diagonal, sb is sub-diagonal 
         rhs_diag         ! and diag is diagonal (mmol/m^3)
 
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij,           & ! horizontal indices, combine i and j loops
         k, ks, ki, kr   ! vertical indices and row counters

      real (kind=dbl_kind), dimension (icells):: &
         dhf             ! darcy_V*dt/hin_old
            
      !-----------------------------------------------------------------
      ! Initialize matrix elements.
      !-----------------------------------------------------------------

      sbdiag(:,:)= c0
      diag(:,:) = c1
      spdiag(:,:) = c0
      rhs(:,:) = c0

      do k = 1, nblyr_hist
           do ij = 1, icells
             rhs   (ij,k) = NOin_init(ij,k)  
             bgrid_temp(ij,k) = bgrid(k)
           enddo
      enddo
     
      !-----------------------------------------------------------------
      ! Compute matrix elements
      !-----------------------------------------------------------------

       do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)  
         dhf(ij) = darcy_V(i,j)*dt/hin_old(ij)
         if (hin_old(ij) > thin .OR. thin_check(ij) > p5) then
           thin_check(ij) = c1
           bgrid_temp(ij,1) =  c0 - grid_o_t/hin_old(ij)
           bgrid_temp(ij,nblyr_hist) = c1 + grid_o/hin_old(ij)
           do k = 1, nblyr+1
               ivel(ij,k) =  (igrid(k)*dhb(i,j) - (igrid(k)-c1)*(dht(i,j)))/hin_old(ij) 
               if (k < nblyr+1) then
                  vel(ij,k) = (bgrid_temp(ij,k+1)*(dhb(i,j)) - &
                          (bgrid_temp(ij,k+1) - c1)* (dht(i,j)) )/hin_old(ij) 
               endif
           enddo    !k
           do k = 2, nblyr+1

             if (k == 2) then               
                 if (dht(i,j) + darcy_V(i,j)*dt/iphin_N(i,j,1) > c0) then   

                    sbdiag(ij,2) =  c0 
                 
                    diag(ij,2)   = c1 + dt/zphin_N(i,j,k)/(igrid(k) - igrid(k-1))* &
                             (iDin(i,j,k)/(bgrid_temp(ij,k+1)- bgrid_temp(ij,k))) + &
                                (ivel(ij,1)*(zphin_N(i,j,2)-iphin_N(i,j,1)/c2)+ dhf(ij)/c2)/&
                              (zphin_N(i,j,2)*igrid(2))


                    spdiag(ij,2) = -(ivel(ij,1)*zphin_N(i,j,1) + dhf(ij))/c2/&
                             zphin_N(i,j,k)/igrid(2)  -  &
                             dt *iDin(i,j,k)/zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1))/(bgrid_temp(ij,k+1) - bgrid_temp(ij,k))
                 else   
                    sbdiag(ij,2) =(ivel(ij,1)*zphin_N(i,j,1)+ dhf(ij))/&
                           zphin_N(i,j,2)/(igrid(2) - bgrid_temp(ij,1))
 
                    diag(ij,2)   = c1 + dt/zphin_N(i,j,k)/(igrid(k) - &
                              igrid(k-1))* (iDin(i,j,k)/(bgrid_temp(ij,k+1)- &
                              bgrid_temp(ij,k))) - &
                              (ivel(ij,1)*zphin_N(i,j,1)+ dhf(ij))/c2/zphin_N(i,j,2)/ &
                              (igrid(2)-bgrid_temp(ij,1))


                    spdiag(ij,2) = - (ivel(ij,1)*zphin_N(i,j,1) + dhf(ij))/c2/ &
                            zphin_N(i,j,2)/ (igrid(2)-bgrid_temp(ij,1)) - &
                             dt *iDin(i,j,k)/zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1))/(bgrid_temp(ij,k+1) - bgrid_temp(ij,k)) 
                  endif
             elseif (k == nblyr+1) then
                if (dhb(i,j) +darcy_V(i,j)*dt .LE. 0) then    !bottom melt /iphin_N(i,j,nblyr+1)

                   sbdiag(ij,nblyr+1) = (vel(ij,nblyr)*(iphin_N(i,j,nblyr)) + dhf(ij))/c2/&
                          zphin_N(i,j,k)/(igrid(nblyr+1) - igrid(nblyr)) - & 
                          dt *iDin(i,j,k-1)/zphin_N(i,j,k)/(igrid(k) - &
                          igrid(k-1))/(bgrid_temp(ij,k) - bgrid_temp(ij,k-1)) 

                   diag(ij,nblyr+1)   = c1 +  dt/zphin_N(i,j,k)/(igrid(k) - &
                              igrid(k-1))* (iDin(i,j,k)/(bgrid_temp(ij,k+1)- &    
                              igrid(k)) + iDin(i,j,k-1)/(bgrid_temp(ij,k)-  &  
                              bgrid_temp(ij,k-1))) - &
                              (vel(ij,nblyr)*(iphin_N(i,j,nblyr+1) - zphin_N(i,j,nblyr)/c2) + &
                              dhf(ij)/c2)/&
                             (igrid(nblyr+1)-igrid(nblyr))/zphin_N(i,j,nblyr+1)

                   spdiag(ij,nblyr+1) = - dt *iDin(i,j,k)/zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1))/(bgrid_temp(ij,k+1) - igrid(k))

                else
                 sbdiag(ij,nblyr+1) = (vel(ij,nblyr)*iphin_N(i,j,nblyr) + dhf(ij))/c2/&
                         zphin_N(i,j,k)/(bgrid_temp(ij,nblyr+2)-igrid(nblyr)) -  & 
                          dt *iDin(i,j,k-1)/zphin_N(i,j,k)/(igrid(k) - &
                          igrid(k-1))/(bgrid_temp(ij,k) - bgrid_temp(ij,k-1))

                 diag(ij,nblyr+1)   = c1 +  dt/zphin_N(i,j,k)/(igrid(k) - &
                              igrid(k-1))* (iDin(i,j,k)/(bgrid_temp(ij,k+1)- &
                              igrid(k)) + iDin(i,j,k-1)/(bgrid_temp(ij,k)-  & 
                              bgrid_temp(ij,k-1))) + &
                              (vel(ij,nblyr)*iphin_N(i,j,nblyr) + dhf(ij))/c2/&
                              zphin_N(i,j,nblyr+1)/(bgrid_temp(ij,nblyr+2)-igrid(nblyr))

                 spdiag(ij,nblyr+1) = - (vel(ij,nblyr) + dhf(ij)/iphin_N(i,j,nblyr+1))/&
                             zphin_N(i,j,k)/(bgrid_temp(ij,nblyr+2)-igrid(nblyr)) -  &
                             dt *iDin(i,j,k)/zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1))/(bgrid_temp(ij,k+1) - igrid(k))
               endif
             else    
                 sbdiag(ij,k) = (vel(ij,k-1)*iphin_N(i,j,k-1)+ dhf(ij))/c2/&
                           zphin_N(i,j,k)/(igrid(k) - igrid(k-1)) -  & 
                           dt *iDin(i,j,k-1)/zphin_N(i,j,k)/(igrid(k) - &
                           igrid(k-1))/(bgrid_temp(ij,k) - bgrid_temp(ij,k-1))
  
                 diag(ij,k)   = c1 + dt/zphin_N(i,j,k)/(igrid(k) - &
                              igrid(k-1))* (iDin(i,j,k)/(bgrid_temp(ij,k+1)- &
                              bgrid_temp(ij,k)) + iDin(i,j,k-1)/(bgrid_temp(ij,k)-  &
                              bgrid_temp(ij,k-1))) - vel(ij,k-1)*(iphin_N(i,j,k)-iphin_N(i,j,k-1))/&
                              c2/(igrid(k)-igrid(k-1))/zphin_N(i,j,k)

                 spdiag(ij,k) = - (vel(ij,k-1)*iphin_N(i,j,k) + dhf(ij))/c2/zphin_N(i,j,k)/&
                             (igrid(k)-igrid(k-1)) - dt *iDin(i,j,k)/zphin_N(i,j,k)/(igrid(k) - &
                             igrid(k-1))/(bgrid_temp(ij,k+1) - bgrid_temp(ij,k)) 
            endif 
            rhs(ij,k)   = NOin_init(ij,k) + reactb(ij,k-1) 

           enddo               !k
      endif                    !hin_old > thin
      enddo                    ! ij

802 format (5d15.5)

      end subroutine get_matrix_elements_calc_bgc

!=======================================================================
!
! Tridiagonal matrix solver-- for salinity
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
      subroutine tridiag_solver (icells,    &
                                 nmat,     sbdiag,   &
                                 diag,     spdiag,   &
                                 rhs,      xout)

      integer (kind=int_kind), intent(in) :: &
          icells         ! number of cells with aicen > puny

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

      ! local variables     

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
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW
!          Nicole Jeffery, LANL

      subroutine bgc_diags (dt)

      use ice_broadcast, only: broadcast_scalar
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, pbloc, &
                                plat, plon
      use ice_domain_size, only: ncat, nltrcr
      use ice_grid, only: lmask_n, lmask_s, tarean, tareas, grid_type
      use ice_state, only: aice, aicen 

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, nn, ii,jj, iblk

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         pN_sk, pNit_sk, pAm_sk, pSil_sk, &
         pDMSPp_sk, pDMSPd_sk, pDMS_sk, pN_ac, &
         pNit_ac, pAm_ac, pSil_ac, pDMSP_ac, pDMS_ac, pN_tot, &
         pflux_NO, pflux_N, pflux_Sil, pflux_NH

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
               pNit_ac(n)  = ocean_bio(i,j,nlt_bgc_NO,iblk)   !nit(i,j,iblk)
               pAm_ac(n)   = ocean_bio(i,j,nlt_bgc_NH,iblk)   !amm(i,j,iblk)
               pSil_ac(n)  = ocean_bio(i,j,nlt_bgc_Sil,iblk)  !sil(i,j,iblk)
               pDMSP_ac(n) = ocean_bio(i,j,nlt_bgc_DMSPp,iblk)!dmsp(i,j,iblk)
               pDMS_ac(n)  = ocean_bio(i,j,nlt_bgc_DMS,iblk)  !dms(i,j,iblk)
               pN_ac(n)    = ocean_bio(i,j,nlt_bgc_N,iblk)    !algalN(i,j,iblk)             
               !----------------------------------------------------------------------
               ! fluxes in mmol/m^2/d
               ! concentrations are bulk in mmol/m^3
               !---------------------------------------------------------------------
               if (solve_skl_bgc) then
                 pNit_sk(n) = c0
                 pAm_sk(n) = c0
                 pSil_sk(n) = c0
                 pDMSPp_sk(n) = c0
                 pDMSPd_sk(n) = c0
                 pDMS_sk(n) = c0
               
                 pN_sk(n) = trcr(i,j,nt_bgc_N_sk,iblk)*phi_sk/sk_l                
                 pflux_N(n) = flux_bio(i,j,nlt_bgc_N,iblk)*mps_to_cmpdy/c100 
                 if (tr_bgc_Nit_sk) then
                    pNit_sk(n)  = trcr(i,j,nt_bgc_Nit_sk,iblk)*phi_sk/sk_l   
                    pflux_NO(n) = flux_bio(i,j,nlt_bgc_NO,iblk)*mps_to_cmpdy/c100 
                 endif
                 if (tr_bgc_Am_sk) then
                    pAm_sk(n)   = trcr(i,j,nt_bgc_Am_sk,iblk)*phi_sk/sk_l        
                    pflux_NH(n) = flux_bio(i,j,nlt_bgc_NH,iblk) *mps_to_cmpdy/c100 
                 endif
                 if (tr_bgc_Sil_sk) then
                    pSil_sk(n)   = trcr(i,j,nt_bgc_Sil_sk,iblk)*phi_sk/sk_l       
                    pflux_Sil(n) = flux_bio(i,j,nlt_bgc_Sil,iblk) *mps_to_cmpdy/c100 
                 endif
                 if (tr_bgc_DMSPp_sk) pDMSPp_sk(n) = trcr(i,j,nt_bgc_DMSPp_sk,iblk)*phi_sk/sk_l  
                 if (tr_bgc_DMSPd_sk) pDMSPd_sk(n) = trcr(i,j,nt_bgc_DMSPd_sk,iblk)*phi_sk/sk_l 
                 if (tr_bgc_DMS_sk)   pDMS_sk(n)   = trcr(i,j,nt_bgc_DMS_sk,iblk)*phi_sk/sk_l  

               elseif (tr_bgc_NO) then   ! zbgc
                 pN_tot(n) = c0
                 pflux_NO(n)                =   flux_bio(i,j,nlt_bgc_NO,iblk)
                 if (tr_bgc_Sil)pflux_Sil(n)=   flux_bio(i,j,nlt_bgc_Sil,iblk)
                 if (tr_bgc_N)  pflux_N(n)  =   flux_bio(i,j,nlt_bgc_N,iblk)
                 if (tr_bgc_NH) pflux_NH(n) =   flux_bio(i,j,nlt_bgc_NH,iblk)

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
                   pNO(n,k)                  =  trcr(i,j,nt_bgc_NO+k-1,iblk)                   
                   if (tr_bgc_Sil) pSil(n,k) =  trcr(i,j,nt_bgc_Sil+k-1,iblk)     
                   if (tr_bgc_N)   pN(n,k)   =  trcr(i,j,nt_bgc_N+k-1,iblk)    
                   if (tr_bgc_PON) pPON(n,k) =  trcr(i,j,nt_bgc_PON+k-1,iblk)
                  
                   if (aice(i,j,iblk) > c0) pN_tot(n) = pN_tot(n) + &
                          pN(n,k)/real(nblyr,kind=dbl_kind)*vice(i,j,iblk)/aice(i,j,iblk)
                 enddo   !k
                endif
             endif                 ! my_task = pmloc
 
             call broadcast_scalar(pN_ac    (n), pmloc(n))     
             call broadcast_scalar(pNit_ac  (n), pmloc(n))             
             call broadcast_scalar(pAm_ac   (n), pmloc(n))             
             call broadcast_scalar(pSil_ac  (n), pmloc(n))             
             call broadcast_scalar(pDMSP_ac (n), pmloc(n))             
             call broadcast_scalar(pDMS_ac  (n), pmloc(n))
             call broadcast_scalar(pflux_N  (n), pmloc(n))     
             call broadcast_scalar(pflux_NO (n), pmloc(n))             
             call broadcast_scalar(pflux_NH (n), pmloc(n))             
             call broadcast_scalar(pflux_Sil(n), pmloc(n))

            if (solve_skl_bgc) then              ! skl_bgc
             call broadcast_scalar(pN_sk    (n), pmloc(n))            
             call broadcast_scalar(pNit_sk  (n), pmloc(n))             
             call broadcast_scalar(pAm_sk   (n), pmloc(n))             
             call broadcast_scalar(pSil_sk  (n), pmloc(n))             
             call broadcast_scalar(pDMSPp_sk(n), pmloc(n))             
             call broadcast_scalar(pDMSPd_sk(n), pmloc(n))             
             call broadcast_scalar(pDMS_sk  (n), pmloc(n))        
            endif   !tr_bgc_sk

           if (tr_bgc_NO) then                   !  z_bgc
              call broadcast_scalar(pN_tot  (n), pmloc(n))
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
            endif
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
       write(nu_diag,*) '--ocean mixed layer----'
       write(nu_diag,900) 'algal N(mmol/m^3)      = ',pN_ac(1),pN_ac(2)
       write(nu_diag,900) 'nitrate(mmol/m^3)      = ',pNit_ac(1),pNit_ac(2)
       write(nu_diag,900) 'ammonia/um(mmol/m^3)   = ',pAm_ac(1),pAm_ac(2)
       write(nu_diag,900) 'silicon(mmol/m^3)      = ',pSil_ac(1),pSil_ac(2)
       if (tr_bgc_DMS_sk .OR. tr_bgc_DMS) then
         write(nu_diag,900) 'DMSP(mmol/m^3)         = ',pDMSP_ac(1),pDMSP_ac(2)
         write(nu_diag,900) 'DMS(mmol/m^3)          = ',pDMS_ac(1),pDMS_ac(2)
       endif
       write(nu_diag,*) '                         '
       write(nu_diag,*) '--ice-ocean fluxes----'
       write(nu_diag,900) 'algalN flux(mmol/m^2/d)= ',pflux_N(1),pflux_N(2)
       write(nu_diag,900) 'nit. flux(mmol/m^2/d)  = ',pflux_NO(1),pflux_NO(2)
       write(nu_diag,900) 'amm. flux(mmol/m^2/d)  = ',pflux_NH(1),pflux_NH(2)
       write(nu_diag,900) 'sil. flux(mmol/m^2/d)  = ',pflux_Sil(1),pflux_Sil(2)
       write(nu_diag,*) '                         '
      if (solve_skl_bgc) then  
        write(nu_diag,*) '-Bulk Bottom bgc-------'
        write(nu_diag,900) 'nitrogen, (mmol/m^3)   = ',pN_sk(1),pN_sk(2)
        write(nu_diag,900) 'nitrate, (mmol/m^3)    = ',pNit_sk(1),pNit_sk(2)
        write(nu_diag,900) 'ammonia/um, (mmol/m^3) = ',pAm_sk(1),pAm_sk(2)
        write(nu_diag,900) 'silicon, (mmol/m^3)    = ',pSil_sk(1),pSil_sk(2)
        if (tr_bgc_DMS_sk) then
          write(nu_diag,900) 'DMSPp, (mmol/m^3)      = ',pDMSPp_sk(1),pDMSPp_sk(2)
          write(nu_diag,900) 'DMSPd, (mmol/m^3)      = ',pDMSPd_sk(1),pDMSPd_sk(2)
          write(nu_diag,900) 'DMS, (mmol/m^3)        = ',pDMS_sk(1),pDMS_sk(2)
        endif
        write(nu_diag,*) '                         '
      elseif (tr_bgc_NO) then
        write(nu_diag,*) '------ zbgc  Model------' 
        if (tr_bgc_N) then
          write(nu_diag,900) 'tot bulk N(mmol/m^2) = ',pN_tot(1),pN_tot(2)
          write(nu_diag,*) '                         '
        endif
        write(nu_diag,803) 'NO(1) Bulk NO3   ','NO(2) Bulk NO3 '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pNO(n,k),n=1,2), k = 1,nblyr)              
        write(nu_diag,*) '                           '
        if (tr_bgc_PON) then
          write(nu_diag,803) 'PON(1) NO3 tracer   ','PON(2) NO3 tracer '
          write(nu_diag,*) '---------------------------------------------------'
          write(nu_diag,802) ((pPON(n,k),n=1,2), k = 1,nblyr)              
          write(nu_diag,*) '                         '
        endif
        if (solve_zbgc) then
          write(nu_diag,803) 'growN(1) specific growth  ','growN(2) specific growth '
          write(nu_diag,*) '---------------------------------------------------'
          write(nu_diag,802) ((pgrowN(n,k),n=1,2), k = 1,nblyr) 
          write(nu_diag,*) '                         '
          write(nu_diag,803) 'zfswin(1) PAR  ','zfswin(2) PAR '
          write(nu_diag,*) '---------------------------------------------------'
          write(nu_diag,802) ((pzfswin(n,k),n=1,2), k = 1,nblyr+1)              
          write(nu_diag,*) '                         '
        endif
        if (tr_bgc_N) then
          write(nu_diag,803) 'N(1) Bulk algalN   ','N(2) Bulk algalN '
          write(nu_diag,*) '---------------------------------------------------'
          write(nu_diag,802) ((pN(n,k),n=1,2), k = 1,nblyr)              
          write(nu_diag,*) '                         '
        endif
        if (tr_bgc_Sil) then
          write(nu_diag,803) 'Sil(1) Bulk Silicate   ','Sil(2) Bulk Silicate '
          write(nu_diag,*) '---------------------------------------------------'
          write(nu_diag,802) ((pSil(n,k),n=1,2), k = 1,nblyr)  
        endif
      endif                   ! tr_bgc_NO
      endif                   ! print_points
      endif                   ! my_task = master_task 

  802 format (f24.17,2x,f24.17)
  803 format (a25,2x,a25)
  900 format (a25,2x,f24.17,2x,f24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
  903 format (a25,5x,i4,1x,i4,1x,i4,1x,i4,7x,i4,1x,i4,1x,i4,1x,i4)

      end subroutine bgc_diags


!=======================================================================
!
!  Written in flux conservative form 
!           d(h zphi c)/dt = d(v zphi c + D*h dc/dx)/dx + zphi hR  
!           where v = (x-1) dzt/dt - x dzb/dt (which includes flushing),
!           D includes the gravity drainage eddy term plus molecular = (De + zphi Dm)/h^2, 
!           and R is the algal chemistry = react/dt
!
! authors     Nicole Jeffery, LANL

      subroutine check_conserve_bgc &
                                     (hin, hin_old, Cin, Cin_old, dt, &
                                      igrid,good_numerics,flux_o_tot, &
                                      react,aicen) 

      real (kind=dbl_kind), dimension(nblyr), intent(in) :: &
         Cin           , & ! bulk c = zphi*c
         Cin_old       , & ! old bulk c = zphi_old * c_old
         react             ! reaction term

      real (kind=dbl_kind), intent(in) :: &
         aicen         , & ! ice area
         hin_old       , & ! old brine thickness (m) 
         hin           , & ! new brine thickness (m)
         dt                ! time-step

      real (kind=dbl_kind), intent(in) :: &
         flux_o_tot       ! tracer flux, gravity+molecular drainage flux,
                          ! and boundary flux to ocean (kg/m^2/s)  
                          ! positive into the ocean  

      real (kind=dbl_kind), dimension (nblyr + 1), intent(in) :: &
         igrid            ! biology nondimensional grid interface points 

      logical (kind=log_kind), intent(inout) :: &   
         good_numerics    ! true if conservation satisfied within error

     ! local variables

     integer (kind=int_kind) :: &
         k         ! vertical biology layer index 
   
     real (kind=dbl_kind) :: &
         sum_old      , &  ! total initial ice tracer conc. (mmol/m^2)
         sum_new      , &  ! total final ice tracer conc. (mmol/m^2)
         sum_true     , &  !
         sum_react    , & !
         dh           , &  ! hin-hin_old
         dsum         , &  ! sum_true - sum_new
         htrue        , &  !
         diff_dt           !difference in units of (kg/m^2/s)
       
     real (kind=log_kind), parameter :: &
         accuracy = 2.0e-5 ! (kg/m^2/s)
                             
         good_numerics = .true.
         dh = hin - hin_old
         sum_old = c0
         sum_new = c0
         sum_true = c0
         sum_react = c0

         do k = 1, nblyr
          sum_old = sum_old + Cin_old(k)*hin_old*(igrid(k+1)-igrid(k))  !hin_old
          sum_new = sum_new + Cin(k)*hin*(igrid(k+1)-igrid(k)) 
          sum_react = sum_react + react(k)*hin_old*(igrid(k+1)-igrid(k))
          if (Cin(k) > c0) then
            sum_true = sum_true +Cin(k)*hin*(igrid(k+1)-igrid(k))
          endif 
         enddo
         dsum = sum_true - sum_new   !negative correction
     !-------------------------------------
     !  Ocean flux: positive into the ocean
     !-------------------------------------    
         diff_dt = ((sum_true - sum_old -sum_react)/dt + flux_o_tot)*aicen   
  
     if (abs(diff_dt) > accuracy ) then
         good_numerics = .false.
         write(nu_diag,*) 'Conservation failure in bgc,diff_dt,accuracy:',diff_dt, accuracy
         write(nu_diag,*) 'Total initial tracer',sum_old
         write(nu_diag,*) 'Total final1  tracer',sum_true
         write(nu_diag,*) 'Total final2 (with possible negative values)',sum_new
         write(nu_diag,*) 'final1-final2,dsum',dsum
         write(nu_diag,*) 'Total reactions',sum_react
         write(nu_diag,*) 'Total flux_o_tot',flux_o_tot
         write(nu_diag,*) 'aicen:',aicen
     endif

     end subroutine check_conserve_bgc

!=======================================================================
!
! Find ice-ocean flux when ice is thin and internal dynamics/reactions are
!             assumed to be zero
!
! authors     Nicole Jeffery, LANL

      subroutine thin_ice_flux (hin, hin_old, phin, Cin, flux_o_tot,igrid,dt) 

      real (kind=dbl_kind), dimension(nblyr), intent(in) :: &
         phin

      real (kind=dbl_kind), dimension(nblyr), intent(in) :: &
         Cin              ! initial concentration

      real (kind=dbl_kind), intent(in) :: &
         hin_old   , &     ! brine thickness (m) 
         hin       , &     ! new brine thickness (m)
         dt                ! time step

      real (kind=dbl_kind), intent(inout) :: &
         flux_o_tot       ! tracer flux, gravity+molecular drainage flux ,
                          ! and boundary flux to ocean (kg/m^2/s)  
                          ! positive into the ocean  

      real (kind=dbl_kind), dimension (nblyr + 1), intent(in) :: &
         igrid            ! biology nondimensional grid interface points 

     ! local variables

     integer (kind=int_kind) :: &
         k         ! vertical biology layer index
   
     real (kind=dbl_kind) :: &
         sum_bio        !
   
         sum_bio = c0
         do k = 1, nblyr
          sum_bio = sum_bio + Cin(k)*phin(k)*(igrid(k+1)-igrid(k))
         enddo

         flux_o_tot = flux_o_tot -(hin-hin_old)*sum_bio/dt

     end subroutine thin_ice_flux

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
! Dumps all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_bgc(filename_spec)

      use ice_domain_size, only: ncat
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_domain, only: nblocks
      use ice_restart, only: lenstr, restart_dir, restart_file
      use ice_state, only: trcrn

      character(len=char_len_long), intent(in), optional :: filename_spec

      ! local variables

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
      if (tr_bgc_NO) then
        if (tr_bgc_N) then
         do k = 1,nblyr
              call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_N+k-1,n,:),'ruf8',diag)
         enddo
        endif
        do k = 1,nblyr
            call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_NO+k-1,n,:),'ruf8',diag)
        enddo
        if (tr_bgc_NH) then
         do k = 1,nblyr 
            call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_NH+k-1,n,:),'ruf8',diag)
         enddo
        endif
         if (tr_bgc_Sil) then
         do k = 1,nblyr 
            call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Sil+k-1,n,:),'ruf8',diag)
         enddo
         endif
         if (tr_bgc_DMSPp) then
         do k = 1,nblyr 
           call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPp+k-1,n,:),'ruf8',diag)
         enddo
         endif
         if (tr_bgc_DMSPd) then
         do k = 1,nblyr 
            call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPd+k-1,n,:),'ruf8',diag)
         enddo
         endif
         if (tr_bgc_PON) then
         do k = 1,nblyr
            call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_PON+k-1,n,:),'ruf8',diag)
         enddo
         endif
        elseif (solve_skl_bgc) then
         call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_N_sk,n,:),'ruf8',diag)
         if (tr_bgc_C_sk)  call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_C_sk,n,:),'ruf8',diag)
         if (tr_bgc_chl_sk)call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_chl_sk,n,:),'ruf8',diag)
         if (tr_bgc_Nit_sk)call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Nit_sk,n,:),'ruf8',diag)
         if (tr_bgc_Am_sk) call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Am_sk,n,:),'ruf8',diag)
         if (tr_bgc_Sil_sk)call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Sil_sk,n,:),'ruf8',diag)
         if (tr_bgc_DMSPp_sk)call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPp_sk,n,:),'ruf8',diag)
         if (tr_bgc_DMSPd_sk)call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPd_sk,n,:),'ruf8',diag)
         if (tr_bgc_DMS_sk) call ice_write(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMS_sk,n,:),'ruf8',diag)
        endif  ! tr_bgc_NO
      enddo
      
      !---------------------------------------------------      
      ! ocean values
      !-----------------------------------------------------
      if (tr_bgc_N .OR. tr_bgc_N_sk) call ice_write(nu_dump_bgc,0,algalN,'ruf8',diag)
      if (tr_bgc_NO .OR. tr_bgc_Nit_sk) call ice_write(nu_dump_bgc,0,nit,'ruf8',diag)
      if (tr_bgc_NH .OR. tr_bgc_Am_sk) call ice_write(nu_dump_bgc,0,amm,'ruf8',diag)
      if (tr_bgc_Sil .OR. tr_bgc_Sil_sk) call ice_write(nu_dump_bgc,0,sil,'ruf8',diag)
      if (tr_bgc_DMSPp .OR. tr_bgc_DMSPp_sk) call ice_write(nu_dump_bgc,0,dmsp,'ruf8',diag)
      if (tr_bgc_DMS .OR. tr_bgc_DMS_sk) call ice_write(nu_dump_bgc,0,dms,'ruf8',diag)

      if (my_task == master_task) close(nu_dump_bgc)

      end subroutine write_restart_bgc

!=======================================================================
!
! Reads all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_bgc(filename_spec)

      use ice_domain_size, only: ncat
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init, &
                              istep0
      use ice_domain, only: nblocks
      use ice_restart, only: lenstr, restart_dir, restart_file, pointer_file, runtype
      use ice_state, only: trcrn
      use ice_exit, only: abort_ice

      character(len=char_len_long), intent(in), optional :: filename_spec

      ! local variables

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
      do n = 1, ncat
       if (tr_bgc_NO) then
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
         
         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' NO for each bgc layer'
          do k = 1,nblyr
             call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_NO+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
         enddo
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

       elseif (solve_skl_bgc) then
         call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_N_sk,n,:),'ruf8',diag)
         if (tr_bgc_C_sk) call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_C_sk,n,:),'ruf8',diag)
         if (tr_bgc_chl_sk) call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_chl_sk,n,:),'ruf8',diag)
         if (tr_bgc_Nit_sk) call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Nit_sk,n,:),'ruf8',diag)
         if (tr_bgc_Am_sk) call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Am_sk,n,:),'ruf8',diag)
         if (tr_bgc_Sil_sk)call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Sil_sk,n,:),'ruf8',diag)
         if(tr_bgc_DMSPp_sk)call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPp_sk,n,:),'ruf8',diag)
         if (tr_bgc_DMSPd_sk)call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPd_sk,n,:),'ruf8',diag)
         if (tr_bgc_DMS_sk) call ice_read(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMS_sk,n,:),'ruf8',diag)
       endif  !tr_bgc_NO
      enddo
     
      if (my_task == master_task) &
        write(nu_diag,*) 'ocean bgc fields:  algalN, nit, amm, sil, dmsp, dms'
      if (tr_bgc_N .OR. tr_bgc_N_sk)call ice_read(nu_restart_bgc,0,algalN,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      if (tr_bgc_NO .OR. tr_bgc_Nit_sk) call ice_read(nu_restart_bgc,0,nit,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      if (tr_bgc_NH .OR. tr_bgc_Am_sk) call ice_read(nu_restart_bgc,0,amm,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      if (tr_bgc_Sil .OR. tr_bgc_Sil_sk) call ice_read(nu_restart_bgc,0,sil,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      if (tr_bgc_DMSPp .OR. tr_bgc_DMSPp_sk) call ice_read(nu_restart_bgc,0,dmsp,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      if (tr_bgc_DMS .OR. tr_bgc_DMS_sk)call ice_read(nu_restart_bgc,0,dms,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
       if (my_task == master_task) close(nu_restart_bgc)

      end subroutine read_restart_bgc

!=======================================================================

      end module ice_algae

!=======================================================================
