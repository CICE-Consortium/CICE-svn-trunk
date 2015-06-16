!  SVN:$Id$
!=======================================================================
!
! Compute biogeochemistry in the skeletal layer
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_algae

      use ice_kinds_mod
      use ice_domain_size, only: nblyr, nilyr, max_blocks,max_algae, &
                     max_aero, max_doc, max_dic, max_don, max_fe
      use ice_blocks, only: nx_block, ny_block
      use ice_fileunits, only: nu_diag, nu_restart_bgc, nu_rst_pointer, &
          nu_dump_bgc, flush_fileunit
      use ice_read_write, only: ice_open, ice_read, ice_write
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_zbgc_shared 
      use ice_state, only: ntrcr, nt_bgc_S

      implicit none

      private 
      public :: get_forcing_bgc, bgc_diags, write_restart_bgc, &
                algal_dyn, read_restart_bgc, &
                skl_biogeochemistry, z_biogeochemistry, &
                get_atm_bgc, update_snow_bgc, fzaero_data, &
                init_bgc_data
      save

      integer (kind=int_kind) :: &
           bgcrecnum = 0   ! old record number (save between steps)

!=======================================================================

      contains

!=======================================================================
!
! Read and interpolate annual climatologies of silicate and nitrate.
! Restore model quantities to data if desired.
!
! author: Elizabeth C. Hunke, LANL

      subroutine get_forcing_bgc

      use ice_calendar, only: dt, istep, mday, month, sec, yday
      use ice_constants, only: c0, field_loc_center, field_type_scalar, &
                               secday
      use ice_domain, only: nblocks
      use ice_flux, only: sss
      use ice_forcing, only: trestore, trest, fyear, &
          read_clim_data_nc, interpolate_data, &
          interp_coeff_monthly, interp_coeff,  &
          read_data_nc_point, c1intp, c2intp

      integer (kind=int_kind) :: &
          i, j, k,iblk, & ! horizontal indices
          ixm,ixp, ixx, & ! record numbers for neighboring months
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth    , & ! middle day of month
          recnum      , & ! record number
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          sec_day     , & !  fix time to noon
          ks              ! bgc tracer index (bio_index_o)

      character (char_len_long) :: & 
         met_file,   &    ! netcdf filename
         fieldname        ! field name in netcdf file

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
          nitdat      , & ! data value toward which nitrate is restored
          sildat          ! data value toward which silicate is restored

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), save :: &
         nit_data, & ! field values at 2 temporal data points
         sil_data
          
      real (kind=dbl_kind), dimension(2), save :: &
         sil_data_p      , &  ! field values at 2 temporal data points
         nit_data_p           ! field values at 2 temporal data points

      real (kind=dbl_kind) :: &
          sec1hr              ! number of seconds in 1 hour

      logical (kind=log_kind) :: readm, read1

      if (.not. trim(nit_data_type)=='ISPOL' .AND. .not. trim(sil_data_type)=='ISPOL') then 
      if (trim(nit_data_type) == 'clim'.or. &
          trim(sil_data_type) == 'clim') then

         nit_file = trim(bgc_data_dir)//'nitrate_climatologyWOA_gx1v6f.nc'
                             !'nitrate_WOA2005_surface_monthly'  ! gx1 only
         sil_file = trim(bgc_data_dir)//'silicate_climatologyWOA_gx1v6f.nc'
                             !'silicate_WOA2005_surface_monthly' ! gx1 only

         if (my_task == master_task .and. istep == 1) then
         if (trim(sil_data_type)=='clim' .AND. tr_bgc_Sil) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'silicate data interpolated to timestep:'
            write (nu_diag,*) trim(sil_file)
         endif
         if (trim(nit_data_type)=='clim' .AND. tr_bgc_Nit) then
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

      endif   ! 'clim prep' sil/nit_data_type

    !-------------------------------------------------------------------
    ! Read two monthly silicate values and interpolate.
    ! Restore toward interpolated value.
    !-------------------------------------------------------------------
    
      if (trim(sil_data_type)=='clim'  .AND. tr_bgc_Sil) then
        ! call read_clim_data (readm, 0, ixm, month, ixp, &
        !                      sil_file,  sil_data, &
        !                      field_loc_center, field_type_scalar)
         fieldname = 'silicate'
         call read_clim_data_nc (readm, 0, ixm, month, ixp, &
                              sil_file, fieldname, sil_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (sil_data, sildat)
         
         if (istep == 1 .or. .NOT. restore_bgc) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sil(i,j,iblk) = sildat(i,j,iblk)
               ks = 2*max_algae + max_doc + 3 + max_dic
               ocean_bio_all(i,j,ks,iblk) = sil(i,j,iblk)                       !Sil
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         elseif (restore_bgc) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sil(i,j,iblk) = sil(i,j,iblk)  &
                         + (sildat(i,j,iblk)-sil(i,j,iblk))*dt/trest
               ks = 2*max_algae + max_doc + 3 + max_dic
               ocean_bio_all(i,j,ks,iblk) = sil(i,j,iblk)                       !Sil
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         endif  !restore
        elseif (tr_bgc_Sil) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sil(i,j,iblk) = 25.0_dbl_kind
               ks = 2*max_algae + max_doc + 3 + max_dic
               ocean_bio_all(i,j,ks,iblk) = sil(i,j,iblk)                       !Sil
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

      endif  !tr_bgc_Sil
    !-------------------------------------------------------------------
    ! Read two monthly nitrate values and interpolate.
    ! Restore toward interpolated value.
    !-------------------------------------------------------------------

      if (trim(nit_data_type)=='clim' .AND. tr_bgc_Nit) then 
        ! call read_clim_data (readm, 0, ixm, month, ixp, &
        !                      nit_file, nit_data, &
        !                      field_loc_center, field_type_scalar)
         fieldname = 'nitrate'
         call read_clim_data_nc (readm, 0, ixm, month, ixp, &
                              nit_file, fieldname, nit_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (nit_data, nitdat)

         if (istep == 1 .or. .NOT. restore_bgc) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               nit(i,j,iblk) = nitdat(i,j,iblk)
               ks = max_algae + 1
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit
               ks =  2*max_algae + max_doc + 7 + max_dic
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         elseif (restore_bgc ) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               nit(i,j,iblk) = nit(i,j,iblk)  &
                         + (nitdat(i,j,iblk)-nit(i,j,iblk))*dt/trest     
               ks = max_algae + 1
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit
               ks =  2*max_algae + max_doc + 7 + max_dic
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
        endif  !restore_bgc

      elseif (trim(nit_data_type) == 'sss'  .AND.  tr_bgc_Nit) then 
          
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               nit(i,j,iblk) =  sss(i,j,iblk)      
               ks = max_algae + 1
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit 
               ks =  2*max_algae + max_doc + 7 + max_dic
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON      
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

      elseif (tr_bgc_Nit) then 
       
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               nit(i,j,iblk) = 12.0_dbl_kind
               ks = max_algae + 1
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit 
               ks =  2*max_algae + max_doc + 7 + max_dic
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON      
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

      endif   !tr_bgc_Nit

    !-------------------------------------------------------------------
    ! Data from Papdimitrious et al., 2007, Limnol. Oceanogr. 
    ! and WOA at 68oS, 304.5oE : 
    ! daily data located at the end of the 24-hour period. 
    !-------------------------------------------------------------------

      elseif (trim(nit_data_type) == 'ISPOL' .or. trim(sil_data_type) == 'ISPOL') then 

         nit_file = trim(bgc_data_dir)//'nutrients_daily_ISPOL_WOA_field3.nc'
         sil_file = trim(bgc_data_dir)//'nutrients_daily_ISPOL_WOA_field3.nc' 

         if (my_task == master_task .and. istep == 1) then
         if (tr_bgc_Sil) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'silicate data interpolated to timestep:'
            write (nu_diag,*) trim(sil_file)
         endif
         if (tr_bgc_Nit) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'nitrate data interpolated to timestep:'
            write (nu_diag,*) trim(nit_file)
            if (restore_bgc) write (nu_diag,*) &
              'bgc restoring timescale (days) =', trestore
         endif
         endif                     ! my_task, istep

        dataloc = 2                          ! data located at end of interval
        sec1hr = secday                      ! seconds in day
        maxrec = 365                         ! 

        ! current record number
        recnum = int(yday)   

        ! Compute record numbers for surrounding data (2 on each side)
        ixm = mod(recnum+maxrec-2,maxrec) + 1
        ixx = mod(recnum-1,       maxrec) + 1
       
        recslot = 2
        ixp = -99
        call interp_coeff (recnum, recslot, sec1hr, dataloc)

        read1 = .false.
        if (istep==1 .or. bgcrecnum .ne. recnum) read1 = .true.
 
 
        if (tr_bgc_Sil) then
          met_file = sil_file
          fieldname= 'silicate' 

          call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, sil_data_p, &
                    field_loc_center, field_type_scalar)
      
          sil(:,:,:) =  c1intp * sil_data_p(1) &
                  + c2intp * sil_data_p(2)
         endif
         if (tr_bgc_Nit) then
           met_file = nit_file
           fieldname= 'nitrate' 

           call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, nit_data_p, &
                    field_loc_center, field_type_scalar)
      
           nit(:,:,:) =  c1intp * nit_data_p(1) &
                  + c2intp * nit_data_p(2)
         endif
         
       !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               ks = 2*max_algae + max_doc + 3 + max_dic
               ocean_bio_all(i,j,ks,iblk) = sil(i,j,iblk)                       !Sil  
               ks = max_algae + 1
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !nit
               ks =  2*max_algae + max_doc + 7 + max_dic
               ocean_bio_all(i,j,ks,iblk) = nit(i,j,iblk)                       !PON    
            enddo
            enddo
         enddo
       !$OMP END PARALLEL DO

       ! Save record number for next time step
       bgcrecnum = recnum

      endif

      end subroutine get_forcing_bgc

!=======================================================================
!
! skeletal layer biochemistry
! 
      subroutine skl_biogeochemistry (nx_block, ny_block,  &
                                      icells,   dt,        &
                                      indxi,    indxj,     &
                                      nbtrcr,   n_algae,   &
                                      flux_bio, ocean_bio, &
                                      hmix,     aicen,     &
                                      meltb,    congel,    &
                                      fswthru,  first_ice, &
                                      trcrn,    upNOn,     &
                                      upNHn,    grow_alg_skl, &
                                      nbtrcr_sw, trcrn_sw,  hin)

      use ice_constants, only: p5, p05, p1, c1, c0, puny, c10
      use ice_state, only: nt_bgc_N

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of cells with aicen > puny 
         nbtrcr , n_algae  , & ! number of bgc tracers and number algae
         nbtrcr_sw             ! nilyr+nslyr+2 for chlorophyhll

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt         ! time step 

      real (kind=dbl_kind), dimension(nx_block*ny_block), intent(in) :: &
         hin        ! ice thickness (m)

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         first_ice  ! initialized values should be used

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         hmix   , & ! mixed layer depth
         aicen  , & ! ice area 
         meltb  , & ! bottom ice melt
         congel , & ! bottom ice growth 
         fswthru    ! shortwave passing through ice to ocean


      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn      !bulk concentration per m^3
       
      real (kind=dbl_kind), dimension(nx_block,ny_block,nbtrcr_sw), &
         intent(inout) :: &
         trcrn_sw   ! total bulk chlorophyll mg chl/m^3 
                    ! nonzero in bottom 3 cm  (scaled to the bottom 
                    ! of shortwave grid)
       
      ! history variables
  
      real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr), intent(inout) :: &
         flux_bio,& ! ocean tracer flux (mmol/m^2/s) positive into ocean
         ocean_bio  ! ocean tracer concentration (mmol/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,n_algae), intent(out) :: &
         grow_alg_skl, & ! tot algal growth rate (mmol/m^3/s)  
         upNOn       , & !  algal NO uptake rate (mmol/m^3/s) 
         upNHn           !  algal NH uptake rate (mmol/m^3/s) 

      !  local variables

      integer (kind=int_kind) :: i, j, ij, nn, mm

      real (kind=dbl_kind), dimension(icells,nbtrcr):: &
         react  , & ! biological sources and sinks (mmol/m^3)
         cinit  , & ! initial brine concentration*sk_l (mmol/m^2)
         cinit_v, & ! initial brine concentration (mmol/m^3)
         congel_alg, & ! congelation flux contribution to ice algae (mmol/m^2 s) 
                       ! (used as initialization)
         f_meltn       ! vertical melt fraction of skeletal layer in dt

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         Zoo_skl    ! N losses from zooplankton/bacteria ... (mmol/m^3)

      real (kind=dbl_kind), dimension (nx_block*ny_block) :: &
         iTin      

      real (kind=dbl_kind), dimension (nbtrcr):: &
         flux_bio_temp, & ! tracer flux to ocean (mmol/m^2 s)
         PVflag       , & ! 1 for tracers that flow with the brine, 0 otherwise
         cling            ! 1 for tracers that cling, 0 otherwise

      real (kind=dbl_kind), parameter :: &
         PVc = 1.e-6_dbl_kind           , & ! type 'constant' piston velocity for interface (m/s) 
         PV_scale_growth = p5           , & ! scale factor in Jin code PV during ice growth
         PV_scale_melt = p05            , & ! scale factor in Jin code PV during ice melt
         growth_max = 1.85e-10_dbl_kind , & ! PVt function reaches maximum here.  (m/s)
         Tin_bot = -1.8_dbl_kind        , & ! temperature of the ice bottom (oC)
         MJ1 = 9.667e-9_dbl_kind        , & ! (m/s) coefficients in Jin2008
         MJ2 = 38.8_dbl_kind            , & ! (1) from:4.49e-4_dbl_kind*secday   
         MJ3 = 1.04e7_dbl_kind          , & ! 1/(m/s) from: 1.39e-3_dbl_kind*secday^2  
         PV_frac_max = 0.9_dbl_kind         ! Maximum Piston velocity is 90% of skeletal layer/dt

      real (kind=dbl_kind), dimension(icells) :: &
         PVt       , & ! type 'Jin2006' piston velocity (m/s) 
         ice_growth    ! Jin2006 definition: either congel rate or bottom melt rate  (m/s)

      real (kind=dbl_kind):: &
         grow_val  , & ! (m/x)
         rphi_sk   , & ! 1 / skeletal layer porosity
         cinit_tmp      ! temporary variable for concentration (mmol/m^2)

      logical (kind=log_kind), dimension(nx_block,ny_block) :: conserve_N

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         Nerror        ! change in total nitrogen from reactions
  !-------------------------------------
  ! Initialize 
  !------------------------------------
     conserve_N(:,:) = .true.
     Zoo_skl(:,:) = c0

     do nn = 1, nbtrcr 
        cinit     (:,nn) = c0
        cinit_v   (:,nn) = c0 
        congel_alg(:,nn) = c0
        f_meltn   (:,nn) = c0
        react     (:,nn) = c0
        PVflag(nn) = c1
        cling(nn) = c0
  !--------------------------------------------------
  ! only the dominant tracer_type affects behavior
  !  < 0 is purely mobile:  > 0 stationary behavior
  ! NOTE: retention times are not used in skl model
  !---------------------------------------------------
        if (bgc_tracer_type(nn) >= c0) then  
            PVflag(nn) = c0
            cling(nn) = c1
         endif
     enddo    !nn

     rphi_sk = c1/phi_sk
     PVt (:) = c0
     iTin(:) = Tin_bot

     do nn = 1, nbtrcr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
     do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
 
         ice_growth(ij) = (congel(i,j)-meltb(i,j))/dt
         if (first_ice(i,j)) then 
               trcrn(i,j,bio_index(nn))   = ocean_bio(i,j,nn)! * sk_l*rphi_sk
         endif ! first_ice
         cinit(ij,nn)     = trcrn(i,j,bio_index(nn)) * sk_l * rphi_sk
         cinit_v(ij,nn) = cinit(ij,nn)/sk_l
         if (cinit(ij,nn) < c0) then
               write(nu_diag,*)'initial sk_bgc < 0, ij,nn,nbtrcr,cinit(ij,nn)', &
                    ij,nn,nbtrcr,cinit(ij,nn)
               call abort_ice ('ice_bgc.F90: BGC error1')
         endif  
     enddo     ! ij
     enddo     ! nbtrcr

      !-------------------------------------------------------------------
      ! 'Jin2006':
      ! 1. congel/melt dependent piston velocity (PV) for growth and melt
      ! 2. If congel > melt use 'congel'; if melt > congel use 'melt'
      ! 3. For algal N, PV for ice growth only provides a seeding concentration 
      ! 4. Melt affects nutrients and algae in the same manner through PV(melt)
      !-------------------------------------------------------------------

     if (trim(bgc_flux_type) == 'Jin2006') then  
         
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
 
            if (ice_growth(ij) > c0) then  ! ice_growth(ij) = congel(i,j)/dt
               grow_val = min(ice_growth(ij),growth_max)  
               PVt(ij) = -min(abs(PV_scale_growth*(MJ1 + MJ2*grow_val &
                                                       - MJ3*grow_val**2)), &
                              PV_frac_max*sk_l/dt)  
            else                           ! ice_growth(ij) = -meltb(i,j)/dt
               PVt(ij) =  min(abs(PV_scale_melt  *(      MJ2*ice_growth(ij) &
                                                       - MJ3*ice_growth(ij)**2)), &
                              PV_frac_max*sk_l/dt)
            endif
            do nn = 1,nbtrcr  
            if (bgc_tracer_type(nn) >= c0) then
            if (ice_growth(ij) < c0) then ! flux from ice to ocean
               ! Algae and clinging tracers melt like nutrients              
               f_meltn(ij,nn) = PVt(ij)*cinit_v(ij,nn) ! for algae only
            elseif (ice_growth(ij) > c0 .AND. &
               cinit(ij,nn) < ocean_bio(i,j,nn)*sk_l/phi_sk) then
               ! Growth only contributes to seeding from ocean 
               congel_alg(ij,nn) = (ocean_bio(i,j,nn)*sk_l/phi_sk &
                                        - cinit(ij,nn))/dt
            endif ! PVt > c0  
            endif     
            enddo
         enddo    ! ij

      !----------------------------------------------------------------------
      ! 'constant':
      ! 1. Constant PV for congel > melt
      ! 2. For algae, PV for ice growth only provides a seeding concentration 
      ! 3. Melt loss (f_meltn) affects algae only and is proportional to melt
      !-----------------------------------------------------------------------

      else   ! bgc_flux_type = 'constant'

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (ice_growth(ij) > c0) PVt(ij) = -PVc
            do nn = 1,nbtrcr
            if (bgc_tracer_type(nn) >= c0 ) then
              if (ice_growth(ij) >= c0 .AND. &
                 cinit_v(ij,nn) < ocean_bio(i,j,nn)/phi_sk) then
                 congel_alg(ij,nn) = (ocean_bio(i,j,nn)*sk_l/phi_sk &
                                        - cinit(ij,nn))/dt
              elseif (ice_growth(ij) < c0) then
                 f_meltn(ij,nn) = min(c1, meltb(i,j)/sk_l)*cinit(ij,nn)/dt
             endif
           endif
           enddo !nn
         enddo   !ij

      endif  ! bgc_flux_type

      !-----------------------------------------------------------------------
      ! begin building biogeochemistry terms
      !-----------------------------------------------------------------------

      call algal_dyn (nx_block,        ny_block,        &
                      icells,          dt,              &
                      indxi,           indxj,           &
                      fswthru,         react,           & 
                      cinit_v,         nbtrcr,          &
                      grow_alg_skl,    n_algae,         &
                      iTin,                             &
                      upNOn,           upNHn,           &
                      Zoo_skl,                          &
                      Nerror,          conserve_N)

      !-----------------------------------------------------------------------
      ! compute new tracer concencentrations
      !-----------------------------------------------------------------------
  
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
     do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
        
          do nn = 1,nbtrcr
      !-----------------------------------------------------------------------
      ! if PVt(ij) > 0, ie melt, then ocean_bio term drops out (MJ2006)
      ! Combine boundary fluxes
      !-----------------------------------------------------------------------
           
         PVflag(nn) = SIGN(PVflag(nn),PVt(ij))
         cinit_tmp = max(c0, cinit_v(ij,nn) + react(ij,nn))
         flux_bio_temp(nn) = (PVflag(nn)*PVt(ij)*cinit_tmp &
                           -  PVflag(nn)*min(c0,PVt(ij))*ocean_bio(i,j,nn)) &
                           + f_meltn(ij,nn)*cling(nn) - congel_alg(ij,nn)

         if (cinit_tmp*sk_l < flux_bio_temp(nn)*dt) then
             flux_bio_temp(nn) = cinit_tmp*sk_l/dt*(c1-puny)
         endif

         cinit(ij,nn) = cinit_tmp*sk_l - flux_bio_temp(nn)*dt
         flux_bio(i,j,nn) = flux_bio(i,j,nn) + flux_bio_temp(nn)*phi_sk  

         ! Uncomment to update ocean concentration
         ! Currently not coupled with ocean biogeochemistry
         !ocean_bio(i,j,nn) = ocean_bio(i,j,nn) &
         !                  + flux_bio(i,j,nn)/hmix(i,j)*aicen(i,j)
         if (.not. conserve_N(i,j)) then
              write(nu_diag,*) 'N not conserved in skl_bgc, Nerror(i,j):',Nerror(i,j)
              write(nu_diag,*) 'sk_bgc < 0 after algal fluxes, ij,nn,cinit,flux_bio',&
                               ij,nn,cinit(ij,nn),flux_bio(i,j,nn)
              write(nu_diag,*) 'cinit_tmp,flux_bio_temp,f_meltn,congel_alg,PVt,PVflag: '
              write(nu_diag,*) cinit_tmp,flux_bio_temp(nn),f_meltn(ij,nn), &
                               congel_alg(ij,nn),PVt(ij),PVflag(nn)
              write(nu_diag,*) 'congel, meltb: ',congel(i,j),meltb(i,j)
              call abort_ice ('ice_bgc.F90: BGC error')
         elseif (cinit(ij,nn) < c0) then
              write(nu_diag,*) 'sk_bgc < 0 after algal fluxes, ij,nn,cinit,flux_bio',&
                               ij,nn,cinit(ij,nn),flux_bio(i,j,nn)
              write(nu_diag,*) 'cinit_tmp,flux_bio_temp,f_meltn,congel_alg,PVt,PVflag: '
              write(nu_diag,*) cinit_tmp,flux_bio_temp(nn),f_meltn(ij,nn), &
                               congel_alg(ij,nn),PVt(ij),PVflag(nn)
              write(nu_diag,*) 'congel, meltb: ',congel(i,j),meltb(i,j)
              call abort_ice ('ice_bgc.F90: BGC error3')
         endif
         
      !-----------------------------------------------------------------------
      ! reload tracer array:  Bulk tracer concentration (mmol or mg per m^3)
      !-----------------------------------------------------------------------
         trcrn(i,j,bio_index(nn))     = cinit(ij,nn) * phi_sk/sk_l
        
       enddo  !nbtrcr  
      !---------------------------------------------------------------------------
      ! compute the total chlorophyll for the delta-Eddington shortwave calculation
      !---------------------------------------------------------------------------
      trcrn_sw(i,j,:) = c0      
      if (dEdd_algae) then
       do nn = 1,n_algae
         trcrn_sw(i,j, nbtrcr_sw) =trcrn_sw(i,j,nbtrcr_sw) + &
                    F_abs_chl(nn)*R_chl2N(nn)*trcrn(i,j,nt_bgc_N(nn))*sk_l/hin(ij)* &
                    real(nilyr,kind=dbl_kind)
       enddo 
      endif
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
      subroutine z_biogeochemistry (nx_block,     ny_block,  &
                                    icells,       n_cat,     &
                                    dt,           indxii,    &
                                    indxjj,       nbtrcr,    & 
                                    n_algae,      n_zaero,   &   
                                    fcells,       findxi,    &
                                    findxj,       nfcells,   &
                                    nfindxi,      nfindxj,   &  
                                    findxij,      nfindxij,  & 
                                    aicen,        vicen,     & 
                                    hice_old,     ocean_bio, & 
                                    flux_bio,     bphin,     &
                                    iphin,        trcrn,     &  
                                    iDin,         sss,       &
                                    fswthrul,     grow_alg,  &
                                    upNOn,        upNHn,     &
                                    dh_top,       dh_bot,    &
                                    dh_top_chl,   dh_bot_chl,&
                                    zfswin,       TLAT,      &
                                    TLON,         hbri,      & 
                                    hbri_old,     darcy_V,   &
                                    darcy_V_chl,  bgrid,     &
                                    igrid,        icgrid,    &
                                    bphi_min,     zbgc_snow, &
                                    dhice,        zbgc_atm,  &
                                    iTin,         trcrn_sw,  &
                                    max_nsw,      swgrid,    &
                                    Zoo,          meltb)

      use ice_calendar,    only: istep1, time
      use ice_therm_shared, only: ktherm
      use ice_state, only: aice, nt_fbri, nbtrcr_sw, nt_zbgc_frac 
      use ice_constants, only: c0, c1, c2,rad_to_deg, p5, puny
      use ice_shortwave, only: hi_ssl

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells, n_cat,      & ! number of true cells with aicen > 0
         nbtrcr, n_algae,    & ! number of bgc tracers, number of autotrophs
         n_zaero,            & ! number of aerosols
         fcells, nfcells,    & ! first ice icells and not first ice icells 
         max_nsw             ! number of tracers active in dEdd  shortwave
                              
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxii, indxjj,   & ! compressed indices for icells with aicen > puny
         findxi, findxj,   & ! compressed indices for first ice icells with aicen > puny
         nfindxi,nfindxj,  & ! compressed indices for not first ice  icells with aicen > puny
         findxij,nfindxij    !

      real (kind=dbl_kind), intent(in) :: &
         dt             ! time step 

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bgrid          ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid          ! biology vertical interface points

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1), intent(in) :: &
         iTin           ! salinity vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         icgrid         ! CICE interface coordinate   

      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         swgrid         ! dEdd shortwave ice grid  

      real (kind=dbl_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         hbri        , & ! brine height  (m)
         dhice       , & ! change due to sublimation/condensation (m)
         bphi_min        ! surface porosity

      real (kind=dbl_kind), dimension (nx_block*ny_block), &
         intent(inout) :: &
         hbri_old        ! brine height  (m)

      real (kind=dbl_kind), dimension (nx_block*ny_block,nbtrcr), &
         intent(inout) :: &
         zbgc_snow  , &    ! tracer input from snow (mmol/m^3*m)
         zbgc_atm          ! tracer input from  atm (mmol/m^3 *m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice  (m)
         sss         , & ! ocean salinity (ppt)
         TLAT        , &
         TLON        , &
         hice_old    , & ! ice height (m)
         meltb           ! bottom melt in dt (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         darcy_V     , & ! darcy velocity
         darcy_V_chl , & ! darcy velocity for algae
         dh_bot      , & ! change in brine bottom (m)
         dh_bot_chl  , & ! change in brine bottom (m) felt by algae
         dh_top      , & ! change in brine top (m)
         dh_top_chl      ! change in brine top (m) felt by algae
           
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_nsw), &
         intent(inout) :: &
         trcrn_sw         ! tot chlorophyll (mg/m^3) and aerosols on 
                          ! the shortwave grid

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr+1), &
         intent(in) :: &
         fswthrul         ! visible short wave radiation on icgrid (W/m^2)  

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr + 1), &
         intent(out) :: & 
         zfswin           ! visible Short wave flux on igrid (W/m^2)  
       
      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr + 1), &
         intent(inout) :: & 
         Zoo           ! N losses to the system from reaction terms (ie. zooplankton/bacteria) 
                       ! (mmol/m^3)  

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr), intent(inout) :: &
         flux_bio      ! total ocean tracer flux (mmol/m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr), intent(in) :: &   
         !change to  inout when updating ocean fields
         ocean_bio       ! ocean concentrations (mmol/m^3) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1,n_algae), intent(inout) :: &
         upNOn       , & ! algal nitrate uptake rate  (mmol/m^3/s)
         upNHn           ! algal ammonium uptake rate (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1,n_algae), intent(inout) :: &
         grow_alg           ! algal growth rate          (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), intent(in) :: &
         bphin            ! Porosity on the bgrid

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1), intent(in) :: &
         iphin            ! Porosity on the igrid   

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(in) :: &
         iDin             ! Diffusivity/h on the igrid (1/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn           ! bulk tracer concentration (mmol/m^3)

      !-----------------------------------------------------------------------------
      ! algae absorption coefficient for 0.5 m thick layer
      ! Grenfell (1991): SA = specific absorption coefficient= 0.004 m^2/mg Chla
      ! generalizing kalg_bio(k) = SA*\sum R_chl2N(m)*trcrn(i,j,nt_bgc_N(m)+k-1)
      ! output kalg on the icgrid
      !-----------------------------------------------------------------------------
      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, jk      , & ! horizontal index, combines i and j loops
         k, m, mm, nn    ! vertical biology layer index 

      real (kind=dbl_kind), dimension(icells) :: &
         hin         , & ! ice thickness (m)        
         hin_old     , & ! ice thickness before current melt/growth (m)
         ice_conc        ! algal concentration in ice above hin > hinS

      real (kind=dbl_kind), dimension (icells,nblyr+2) :: &
         bphin_N         ! porosity for tracer model has minimum 
                         ! bphin_N >= bphimin

      real (kind=dbl_kind), dimension (icells,nblyr+1) :: &
         iphin_N         ! tracer porosity on the igrid

      real (kind=dbl_kind), dimension(nblyr+1):: &
         sbdiagz      , & ! sub-diagonal matrix elements
         diagz        , & ! diagonal matrix elements
         spdiagz      , & ! super-diagonal matrix elements
         rhsz         , & ! rhs of tri-diagonal matrix equation
         Ml_diag      , & ! lumped mass matrix
         D_spdiag, D_sbdiag ! artificial diffusion matrix

      real (kind=dbl_kind), dimension(icells,nblyr+1,nbtrcr):: &
         react            ! biological sources and sinks for equation matrix

      real (kind=dbl_kind), dimension(icells, nblyr+1,nbtrcr):: &
         in_init_cons     , & ! Initial bulk concentration*h (mmol/m^2)
         biomat_cons      , & ! Matrix output of (mmol/m^2)
         biomat_brine         ! brine concentration (mmol/m^3)

      real (kind=dbl_kind), dimension(nblyr+1):: &
         biomat_low           ! Low order solution

      real (kind=dbl_kind), dimension(icells,nbtrcr):: &
         C_top,            &  ! bulk top tracer source: h phi C(meltwater) (mmol/m^2)
         C_bot,            &  ! bulk bottom tracer source: h phi C(ocean) (mmol/m^2)
         Source_top,       &  ! For cons: (+) top tracer source into ice (mmol/m^2/s)
         Source_bot,       &  ! For cons: (+) bottom tracer source into ice (mmol/m^2/s)
         Sink_bot,         &  ! For cons: (+ or -) remaining bottom flux into ice(mmol/m^2/s)
         Sink_top,         &  ! For cons: (+ or -) remaining bottom flux into ice(mmol/m^2/s)
         ocean_b,          &  ! ocean_bio
         sum_react,        &
         rtau_ret,         &  ! retention frequency (s^-1)
         rtau_rel             ! release frequency   (s^-1)

      real (kind=dbl_kind):: &
         sum_old,     &  !
         sum_new,     &  !
         sum_tot,     &  !
         zspace          !1/nblyr

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: &
         trtmp0      , & ! temporary, remapped tracers
         trtmp           ! temporary, remapped tracers

      real (kind=dbl_kind), dimension(icells):: &
         darcyV,     &     !
         dhtop,      &     !
         dhbot             !
        
      real (kind=dbl_kind),dimension(nx_block,ny_block,nblyr+1) :: & 
         Nerror     ! Change in N after reactions

      real (kind=dbl_kind), dimension(icells, nbtrcr) :: & 
         atm_add_cons         !zbgc_snow+zbgc_atm (mmol/m^3*m)

      logical (kind=log_kind):: & 
         good_numerics

      logical (kind=log_kind), dimension(nx_block,ny_block,nblyr+1) :: &
         conserve_N

      real (kind=dbl_kind), dimension(nblyr+1):: &  ! temporary variables for
         Diff               , & ! diffusivity 
         initcons           , & ! initial concentration
         biocons            , & !  new concentration
         dmobile            , & !
         initcons_mobile    , & !
         initcons_stationary
 
      real (kind=dbl_kind), dimension (icells, nilyr+1):: &
         icegrid         ! correct for large ice surface layers

      real (kind=dbl_kind):: &
         top_conc     ! 1% (min_bgc) of surface concentration 
                      ! when hin > hbri:  just used in sw calculation

      real (kind=dbl_kind):: &
         bio_tmp      ! temporary variable

      real (kind=dbl_kind), dimension(nbtrcr) :: & 
         mobile       ! c1 if mobile, c0 otherwise

      integer (kind=int_kind) :: &
         tcells, ntcells         ! thin ice cells, not thin ice cells

      integer (kind=int_kind), dimension(icells) :: &
         tindxi,  tindxj,  tindxij, & ! compressed indices for tcells and ntcells 
         ntindxi, ntindxj, ntindxij   ! compressed indices for tcells and ntcells 

      ! local parameters
         
      real (kind=dbl_kind), parameter :: &
         accuracy = 1.0e-14_dbl_kind

      integer, parameter :: &
         nt_zfswin = 1    ! for interpolation of short wave to bgrid

  !-------------------------------------
  ! Initialize 
  !----------------------------------- 

        zspace = c1/real(nblyr,kind=dbl_kind)
        in_init_cons(:,:,:) = c0
        atm_add_cons(:,:) = c0
        sum_react(:,:) = c0
        dhtop(:) = c0
        dhbot(:) = c0
        darcyV(:) = c0
        C_top(:,:) = c0
        mobile(:) = c0
        conserve_N(:,:,:) = .true.
          
        do m = 1, nbtrcr
          do k  = 1, nblyr+1
          if (fcells > 0) then
          do jk = 1, fcells   
              i = findxi(jk)
              j = findxj(jk)
              ij = findxij(jk)  
              bphin_N(ij,nblyr+2) =c1
              bphin_N(ij,k) = bphin(i,j,k)
              iphin_N(ij,k) = iphin(ij,k)
              bphin_N(ij,1) = bphi_min(ij) 
     
              trcrn(i,j,bio_index(m) + k-1) = ocean_bio(i,j,m)*zbgc_init_frac(m)
              in_init_cons(ij,k,m) = trcrn(i,j,bio_index(m) + k-1)*hbri_old(ij)

              if (trcrn(i,j,bio_index(m) + k-1) < c0  ) then
                    write(nu_diag,*)'zbgc initialization error first ice'
                    write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,m:'&
                                    ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                     TLON(i,j)*rad_to_deg,istep1,m
                    write(nu_diag,*)'hbri(ij),hbri_old(ij)' 
                    write(nu_diag,*)hbri(ij),hbri_old(ij)
                    write(nu_diag,*) 'trcrn(i,j,bio_index(m) + k-1)'
                    write(nu_diag,*) trcrn(i,j,bio_index(m) + k-1)
                    call abort_ice ('ice_bgc.F90: first ice zbgc error')
              endif 
         enddo    !fcells
         endif   !fcells > 0
         if (nfcells > 0) then
         do jk = 1, nfcells   
              i = nfindxi(jk)
              j = nfindxj(jk)  
              ij = nfindxij(jk)
              bphin_N(ij,nblyr+2) =c1 
              bphin_N(ij,k) = bphin(i,j,k)
              iphin_N(ij,k) = iphin(ij,k)
              bphin_N(ij,1) = bphi_min(ij)
              in_init_cons(ij,k,m) = trcrn(i,j,bio_index(m) + k-1)* hbri_old(ij)     

              if ( trcrn(i,j,bio_index(m) + k-1) < c0  ) then
                    write(nu_diag,*)'zbgc initialization error, not first ice'
                    write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,m:'&
                                    ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                     TLON(i,j)*rad_to_deg,istep1,m
                    write(nu_diag,*)'hbri(ij),hbri_old(ij)' 
                    write(nu_diag,*)hbri(ij),hbri_old(ij)
                    write(nu_diag,*) 'trcrn(i,j,bio_index(m) + k-1)'
                    write(nu_diag,*) trcrn(i,j,bio_index(m) + k-1)
                    call abort_ice ('ice_bgc.F90: start of z_biogeochemistry')
              endif 
         enddo         !jk
         endif         ! nfcells > 0
         enddo         !k
       enddo           !m
      !-----------------------------------------------------------------
      !     boundary conditions
      !-----------------------------------------------------------------
       ice_conc(:) = c0

       do m = 1,nbtrcr 
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells
          i = indxii(ij)
          j = indxjj(ij)  
          hin(ij) = vicen(i,j)/aicen(i,j)
          hin_old(ij) = hice_old(i,j) 

      !-----------------------------------------------------------------
      !   time constants for mobile/stationary phase changes
      !-----------------------------------------------------------------
          if (hin_old(ij)  > hin(ij)) then  !melting
            rtau_rel(ij,m) = c1/tau_rel(m)
            rtau_ret(ij,m) = c0
          else                              !not melting
            rtau_ret(ij,m) = c1/tau_ret(m)
            rtau_rel(ij,m) = c0
          endif  
         
          ocean_b(ij,m) = ocean_bio(i,j,m)
          dhtop(ij)     = dh_top(i,j)
          dhbot(ij)     = dh_bot(i,j)
          darcyV(ij)    = darcy_V(i,j)
          C_top(ij,m)   = in_init_cons(ij,1,m)*trcrn(i,j,nt_zbgc_frac+m-1)!mobile fraction

          if (dhtop(ij)+darcyV(ij)/bphin_N(ij,1)*dt < -puny) then !snow/top ice melt
                 C_top(ij,m) = (zbgc_snow(ij,m)+zbgc_atm(ij,m))/abs(dhtop(ij) +&
                            darcyV(ij)/bphin_N(ij,1)*dt + puny)*hbri_old(ij)       
          elseif (dhtop(ij)+darcyV(ij)/bphin_N(ij,1)*dt >= -puny .and. &
                        abs((zbgc_snow(ij,m)+zbgc_atm(ij,m))*dt) >  puny) then
                 atm_add_cons(ij,m) =  abs(zbgc_snow(ij,m) + zbgc_atm(ij,m))
          else   !only positive fluxes 
                 flux_bio(i,j,m) = flux_bio(i,j,m) +  max(c0,(zbgc_atm(ij,m) + zbgc_snow(ij,m))/dt)
          endif
          C_bot(ij,m) = ocean_bio(i,j,m)*hbri_old(ij)*iphin_N(ij,nblyr+1)            
       enddo             !ij 
       enddo             !m
   
      !-----------------------------------------------------------------
      ! Interpolate shortwave flux, fswthrul (defined at top to bottom with nilyr+1 
      !  evenly spaced  with spacing = (1/nilyr) to grid variable zfswin:
      !-----------------------------------------------------------------
      
       trtmp(:,:,:) = c0 
       tcells = 0     !thin ice cells
       ntcells = 0    !not thin ice cells 
       tindxi(:) = 0
       tindxj(:) = 0
       tindxij(:) = 0   
       ntindxi(:) = 0
       ntindxj(:) = 0
       ntindxij(:) = 0   
    
       do k = 1, nilyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells    
           i = indxii(ij)
           j = indxjj(ij)               
           ! contains cice values (fswthrul(i,j,1) is surface value)
           ! and fwsthrul(i,j,nilyr+1) is output
           trtmp0(i,j,nt_zfswin+k-1) = fswthrul(i,j,k) 
        enddo   !icells
        enddo   !k

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells    
           i = indxii(ij)
           j = indxjj(ij)      
       
           call remap_zbgc  (ntrcr,  nilyr+1,                    &
                             nt_zfswin,                          &
                             trtmp0(i,j,1:ntrcr), trtmp(i,j,:),  &
                             0,                nblyr+1,          &
                             hin(ij),          hbri(ij),         &
                             icgrid(1:nilyr+1),                  &
                             igrid(1:nblyr+1), ice_conc(ij)) 

           if (hbri_old(ij) > thinS .and. hbri(ij) > thinS) then 
              ntcells = ntcells + 1
              ntindxi(ntcells) = i
              ntindxj(ntcells) = j
              ntindxij(ntcells) = ij   
           else
              tcells = tcells + 1 
              tindxi(tcells) = i
              tindxj(tcells) = j
              tindxij(tcells) = ij  
           endif 
       enddo  
       
       do k = 1,nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells  
           i = indxii(ij)
           j = indxjj(ij) 
           zfswin(i,j,k) = trtmp(i,j,nt_zfswin+k-1)
       enddo                    !ij
       enddo
  
      !-----------------------------------------------------------------
      ! Initialze Biology  
      !----------------------------------------------------------------- 
      do mm = 1, nbtrcr

         mobile(mm) = c0
         if (bgc_tracer_type(mm) .GE. c0) mobile(mm) = c1

      do k = 1, nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells   
           i = indxii(ij)
           j = indxjj(ij) 
           biomat_cons(ij,k,mm) = in_init_cons(ij,k,mm)
       enddo !ij
      enddo  !k
      enddo  !mm

      !-----------------------------------------------------------------
      ! Compute FCT
      !----------------------------------------------------------------- 
      do mm = 1, nbtrcr 

         if (ntcells > 0) then
         do jk = 1, ntcells  
            i = ntindxi(jk)
            j = ntindxj(jk)
            ij = ntindxij(jk)
            do k = 1,nblyr+1
              initcons_mobile(k) = in_init_cons(ij,k,mm)*trcrn(i,j,nt_zbgc_frac+mm-1)
              initcons_stationary(k) = mobile(mm)*(in_init_cons(ij,k,mm)-initcons_mobile(k))
              dmobile(k) = mobile(mm)*(initcons_mobile(k)*(exp(-dt*rtau_ret(ij, mm))-c1) + &
                                 initcons_stationary(k)*(c1-exp(-dt*rtau_rel(ij,mm))))
              initcons_mobile(k) = max(c0,initcons_mobile(k) + dmobile(k))
              initcons_stationary(k) = max(c0,initcons_stationary(k) - dmobile(k))

              Diff(k) = iDin(i,j,k) 
              initcons(k) = initcons_mobile(k)                          
              biocons(k) =  initcons_mobile(k)
            enddo

            call compute_FCT_matrix &
                                (initcons,sbdiagz, dt,               &
                                diagz, spdiagz, rhsz, bgrid,         & 
                                igrid, darcyV(ij),    dhtop(ij),     &
                                dhbot(ij),   iphin_N(ij,:),          &
                                Diff, hbri_old(ij),                  &
                                atm_add_cons(ij,mm), bphin_N(ij,:),  &
                                C_top(ij,mm), C_bot(ij,mm),          &
                                Source_bot(ij,mm), Source_top(ij,mm),&
                                Sink_bot(ij,mm),Sink_top(ij,mm),     &
                                D_sbdiag, D_spdiag, ML_diag)

            call tridiag_solverz &
                               (nblyr+1, sbdiagz,               &
                                diagz,   spdiagz,               &
                                rhsz,    biocons)

            call check_conservation_FCT &
                               (initcons,    &
                                biocons,     &
                                biomat_low,               &
                                Source_top(ij,mm),        &
                                Source_bot(ij,mm),        &
                                Sink_bot(ij,mm),          &
                                Sink_top(ij,mm),          &
                                dt, flux_bio(i,j,mm),     &
                                good_numerics)
                
            call compute_FCT_corr & 
                                (initcons,   &
                                 biocons, dt,&
                                 D_sbdiag, D_spdiag, ML_diag)  

            call regrid_stationary &
                                (initcons_stationary(:), hbri_old(ij),  &
                                 hbri(ij),            meltb(i,j),       &
                                 flux_bio(i,j,mm),    dt,               &
                                 nblyr,               igrid ) 

            biomat_cons(ij,:,mm) =  biocons(:) +  initcons_stationary(:)

            sum_old = (biomat_low(1) + biomat_low(nblyr+1))*zspace/c2
            sum_new = (biocons(1)+ biocons(nblyr+1))*zspace/c2
            sum_tot = (biomat_cons(ij,1,mm) + biomat_cons(ij,nblyr+1,mm))*zspace/c2
            do k = 2,nblyr
                sum_old = sum_old + biomat_low(k)*zspace
                sum_new = sum_new + biocons(k)*zspace
                sum_tot = sum_tot + biomat_cons(ij,k,mm)*zspace
            enddo
            trcrn(i,j,nt_zbgc_frac+mm-1) = zbgc_frac_init(mm)
            if (sum_tot > c0 .and. mobile(mm) > c0) trcrn(i,j,nt_zbgc_frac+mm-1) = sum_new/sum_tot

            if (abs(sum_new-sum_old) > accuracy*sum_old .or. &
                minval(biocons(:)) < c0  .or. minval(initcons_stationary(:)) < c0 &
                .or. .not. good_numerics) then
                write(nu_diag,*)'zbgc FCT tracer solution failed,nn', nn
                write(nu_diag,*)'sum_new,sum_old:',sum_new,sum_old
                write(nu_diag,*)'ij,mm,biocons(:):',ij,mm,biocons(:)
                write(nu_diag,*)'biomat_low:',biomat_low
                write(nu_diag,*)'Diff(:):',Diff(:)
                write(nu_diag,*)'dmobile(:):',dmobile(:)
                write(nu_diag,*)'mobile(mm):',mobile(mm)
                write(nu_diag,*)'initcons_stationary(:):',initcons_stationary(:)
                write (nu_diag, *) 'trcrn(i,j,nt_zbgc_frac+mm-1):',trcrn(i,j,nt_zbgc_frac+mm-1)
                write (nu_diag, *) 'in_init_cons(ij,:,mm):',in_init_cons(ij,:,mm)
                write (nu_diag, *) 'rtau_ret(ij, mm),rtau_rel(ij, mm)',rtau_ret(ij, mm),rtau_rel(ij, mm)
                write(nu_diag,*)'darcyV(ij),dhtop(ij),dhbot(ij)'
                write(nu_diag,*)darcyV(ij),dhtop(ij),dhbot(ij)
                         write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,mm:'&
                                         ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                         TLON(i,j)*rad_to_deg,istep1,mm
                call abort_ice ('z_biogeochemistry: FCT error')
            endif
        enddo  ! ntcells 
        endif  ! ntcells > 0
        if (tcells > 0) then
        do jk = 1, tcells  
            i = tindxi(jk)
            j = tindxj(jk)
            ij = tindxij(jk)

            Call thin_ice_flux(hbri(ij),hbri_old(ij),iphin_N(ij,:),biomat_cons(ij,:,mm), &
                            flux_bio(i,j,mm),zbgc_snow(ij,mm),zbgc_atm(ij,mm),igrid, dt)
        enddo !tcells   
        endif !tcells > 0
        do k = 1,nblyr+1 
        do ij = 1, icells    
           i = indxii(ij)
           j = indxjj(ij)
           biomat_brine(ij,k,mm) =  biomat_cons(ij,k,mm)/hbri(ij)/iphin_N(ij,k) 
        enddo
        enddo ! k
    enddo ! mm 

    react(:,:,:) = c0  
    grow_alg(:,:,:,:) = c0

    if (solve_zbgc) then

       do k = 1, nblyr+1   
         call algal_dyn (nx_block,      ny_block,      &
                         icells,        dt,            &
                         indxii,        indxjj,        &
                         zfswin(:,:,k), react(:,k,:),  & 
                         biomat_brine(:,k,:),          &
                         nbtrcr,                       &
                         grow_alg(:,:,k,:), n_algae,   &
                         iTin(:,k),                    &
                         upNOn(:,:,k,:),               &
                         upNHn(:,:,k,:), Zoo(:,:,k),   &
                         Nerror(:,:,k),                &
                         conserve_N(:,:,k))
       enddo       !k

     endif        !solve_zbgc

      !-----------------------------------------------------------------
      ! Update the tracer variable
      !-----------------------------------------------------------------
    
      do m = 1,nbtrcr
      do k = 1,nblyr+1                  !back to bulk quantity
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells    
            i = indxii(ij)
            j = indxjj(ij)
            bio_tmp = (biomat_brine(ij,k,m) + react(ij,k,m))*iphin_N(ij,k)*hbri(ij)
                       
            if (.not. conserve_N(i,j,k)) then  
                write (nu_diag, *) 'N in algal_dyn not conserved'
                write (nu_diag, *) 'Nerror(i,j,k):', Nerror(i,j,k)
                write (nu_diag, *) 'ij,k,m,hbri(ij),hbri_old(ij),bio_tmp,biomat_cons(ij,k,m),ocean_bio(i,j,m)'
                write (nu_diag, *) ij,k,m,hbri(ij),hbri_old(ij),bio_tmp,biomat_cons(ij,k,m),ocean_bio(i,j,m)
                write (nu_diag, *) 'react(ij,k,m),iphin_N(ij,k),biomat_brine(ij,k,m)'
                write (nu_diag, *) react(ij,k,m),iphin_N(ij,k),biomat_brine(ij,k,m)
                write (nu_diag, *) 'trcrn(i,j,nt_zbgc_frac+m-1):',trcrn(i,j,nt_zbgc_frac+m-1)
                write (nu_diag, *) 'in_init_cons(ij,k,m):',in_init_cons(ij,k,m)
                write (nu_diag, *) 'trcrn(i,j,bio_index(m) + k-1)'
                write (nu_diag, *) trcrn(i,j,bio_index(m) + k-1)
                write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,m:'&
                                 ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                 TLON(i,j)*rad_to_deg,istep1,m
                call abort_ice ('z_biogeochemistry: after integration ')                 
            elseif (abs(bio_tmp) < puny) then  
                     bio_tmp = c0
            elseif (bio_tmp > 1.0e6_dbl_kind) then
                write (nu_diag, *) 'very large bgc value'
                write (nu_diag, *) 'ij,k,m,hbri(ij),hbri_old(ij),bio_tmp,biomat_cons(ij,k,m),ocean_bio(i,j,m)'
                write (nu_diag, *) ij,k,m,hbri(ij),hbri_old(ij),bio_tmp,biomat_cons(ij,k,m),ocean_bio(i,j,m)
                write (nu_diag, *) 'react(ij,k,m),iphin_N(ij,k),biomat_brine(ij,k,m)'
                write (nu_diag, *) react(ij,k,m),iphin_N(ij,k),biomat_brine(ij,k,m)
                write (nu_diag, *) 'trcrn(i,j,nt_zbgc_frac+m-1):',trcrn(i,j,nt_zbgc_frac+m-1)
                write (nu_diag, *) 'in_init_cons(ij,k,m):',in_init_cons(ij,k,m)
                write (nu_diag, *) 'trcrn(i,j,bio_index(m) + k-1)'
                write (nu_diag, *) trcrn(i,j,bio_index(m) + k-1)
                write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,m:'&
                                 ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                 TLON(i,j)*rad_to_deg,istep1,m
                call abort_ice ('z_biogeochemistry: after integration ')
            elseif (bio_tmp < c0) then
                write (nu_diag, *) 'negative bgc'
                write (nu_diag, *) 'ij,k,m,nlt_bgc_Nit,hbri(ij),hbri_old(ij)'
                write (nu_diag, *) ij,k,m,nlt_bgc_Nit,hbri(ij),hbri_old(ij)
                write (nu_diag, *) 'bio_tmp,biomat_cons(ij,k,m),ocean_bio(i,j,m)'
                write (nu_diag, *) bio_tmp,biomat_cons(ij,k,m),ocean_bio(i,j,m)
                write (nu_diag, *) 'react(ij,k,m),iphin_N(ij,k),biomat_brine(ij,k,m)'
                write (nu_diag, *) react(ij,k,m),iphin_N(ij,k),biomat_brine(ij,k,m)
                write (nu_diag, *) 'trcrn(i,j,nt_zbgc_frac+m-1):',trcrn(i,j,nt_zbgc_frac+m-1)
                write (nu_diag, *) 'in_init_cons(ij,k,m):',in_init_cons(ij,k,m)
                write (nu_diag, *) 'rtau_ret(ij, m),rtau_ret(ij, m)',rtau_ret(ij, m),rtau_ret(ij, m)
                write (nu_diag, *) 'trcrn(i,j,bio_index(m) + k-1)'
                write (nu_diag, *) trcrn(i,j,bio_index(m) + k-1)
                write(nu_diag,*)'Category,i,j,ij,TLAT,TLON-- istep1,m:'&
                                 ,n_cat,i,j,ij,TLAT(i,j)*rad_to_deg,&
                                 TLON(i,j)*rad_to_deg,istep1,m
                call abort_ice ('z_biogeochemistry: after integration ')
           endif
           trcrn(i,j,bio_index(m)+k-1) = max(c0, bio_tmp/hbri(ij))
         enddo     !jk
      enddo        !k
      enddo        !m   
    

     if (dEdd_algae) &
         call compute_shortwave_trcr(nx_block,    ny_block,   &
                                    indxii,       indxjj,     &
                                    icells,       n_algae,    &
                                    trcrn(:,:,1:ntrcr),       &
                                    trcrn_sw(:,:,1:nbtrcr_sw),&
                                    swgrid,       hin,        &
                                    hbri,         ntrcr,      &
                                    nilyr,        nblyr,      &
                                    igrid,                    &
                                    nbtrcr_sw,    n_zaero )   
  
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
! authors: Scott Elliott, LANL
!          Nicole Jeffery, LANL

      subroutine algal_dyn (nx_block,     ny_block,     &
                            icells,       dt,           &
                            indxi,        indxj,        &
                            fswthru,      reactb,       & 
                            ltrcrn,       nbtrcr,       &
                            grow_alg,     n_algae,      &
                            T_bot,                      &
                            upNOn,        upNHn,        &
                            Zoo,                        &
                            Nerror,       conserve_N)      

      use ice_constants, only: p1, p5, c0, c1, secday, puny
      use ice_domain_size, only: n_zaero, n_doc, n_dic,  n_don, n_fed, n_fep

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of cells with aicen > puny
         nbtrcr , n_algae      ! number of layer tracers,
                                ! number of autotrophic types

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt         ! time step 

      real (kind=dbl_kind), dimension(nx_block*ny_block), intent(in) :: &
         T_bot      ! ice temperature (oC)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fswthru   ! average shortwave passing through current ice layer (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         Zoo , &     ! N losses from zooplankton/bacteria... (mmol/m^3)
         Nerror      ! Change in N after reactions (mmol/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,n_algae), intent(out) :: &
         grow_alg   !  algal growth rate   (mmol/m^3/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,n_algae), intent(out) :: &
         upNOn,  &  !  algal NO uptake rate   (mmol/m^3/s)
         upNHn      !  algal NH uptake rate   (mmol/m^3/s)

      real (kind=dbl_kind), dimension(icells,nbtrcr), &
         intent(inout) :: &
         reactb        ! biological reaction terms (mmol/m3)

      real (kind=dbl_kind), dimension(icells,nbtrcr), &
         intent(in) :: &
         ltrcrn     ! brine concentrations  in layer (mmol/m^3) 

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(inout) :: & 
         conserve_N

      !  local variables
      !------------------------------------------------------------------------------------
      ! Parameters for 3 possible autotrophs nt_bgc_N(1:3):  diatoms, flagellates, phaeocystis
      !                2 types of dissolved organic carbon nt_bgc_DOC(1:2): 
      !                        polysaccharids, lipids
      !                1 DON (proteins)
      !                1 particulate iron (nt_bgc_Fep) 
      ! Limiting macro/micro nutrients: nt_bgc_Nit -> nitrate, nt_bgc_NH -> ammonium, 
      !                        nt_bgc_Sil -> silicate, nt_bgc_Fed -> dissolved iron   
      ! --------------------------------------------------------------------------------------

      real (kind=dbl_kind),  parameter, dimension(max_algae) :: & 
         ! if dEdd_algae = F, then absorption is included in dEdd, and we don't include here
         chlabs     = (/ 0.03_dbl_kind, 0.01_dbl_kind, 0.05_dbl_kind/), & !these seem high
                                                                          !0.003 1/m/(mg/m^3)
         alpha2max_low  = (/ 0.67_dbl_kind, 0.67_dbl_kind, 0.67_dbl_kind/), & 
                         ! light limitation (1/(W/m^2))
                                                                       !low PAR value
         alpha2max_high  = (/ 0.25_dbl_kind, 0.25_dbl_kind, 0.25_dbl_kind/), & ! light limitation (1/(W/m^2))
                                                                       !high PAR value
         beta2max   = (/ 0.01_dbl_kind, 0.0025_dbl_kind, 0.01_dbl_kind/), & 
                              ! corresponding light inhibition (1/W/m^2)
         !Eppley 1972 curve mux_max = 0.851(1.066)^T , at T = 0
         mu_max     = (/ 0.851_dbl_kind, 0.851_dbl_kind, 0.851_dbl_kind/), & 
                      ! (1/day) maximum growth rate 'Jin2006'   
         grow_Tdep  = (/0.06_dbl_kind, 0.06_dbl_kind, 0.06_dbl_kind/),& ! (1/C)and its T dependence
         fr_graze   = (/c0, p1, p1/) ,  & ! A93 val for S, set to zero in Jin 06   Diatoms are not grazed
         mort_pre   = (/0.02_dbl_kind, 0.02_dbl_kind, 0.02_dbl_kind/),& ! (1/day) prefix to mortality
         mort_Tdep  = (/0.03_dbl_kind ,0.03_dbl_kind, 0.03_dbl_kind/) , & ! (1/C) T dependence of mortality
         k_exude    = (/c0, c0, c0/),    & ! algal carbon  exudation rate (1/d)
         K_Nit      = (/ c1, c1, c1/) ,  & ! nitrate half saturation (mmol/m^3) 
         K_Am       = (/ 0.3_dbl_kind, 0.3_dbl_kind, 0.3_dbl_kind/) ,  & ! ammonium half saturation (mmol/m^3) 
         K_Sil      = (/3.0_dbl_kind , c0, c0/), & ! silicon half saturation (mmol/m^3)
         K_Fe       = (/1.0_dbl_kind , 0.2_dbl_kind, 0.1_dbl_kind/) ! 0.2-1  (nM)
                     ! (nM) iron half saturation  or micromol/m^3
                     ! Timmermans et al 2004 for values 2e-4-1.1e-3 mmol/m^3 for Antarctic diatoms
            
      real (kind=dbl_kind), parameter, dimension(max_DON) :: & 
         f_don     = (/ 0.6_dbl_kind/), & !fraction of spilled grazing 
                     !that goes to each DON pool(proteins and amino acids)
         kn_bac    = (/0.03_dbl_kind/), &!(1/d) Bacterial degredation of DON
         f_don_Am  =  (/0.25_dbl_kind/)   !frac of remineralized DON to Am

      real (kind=dbl_kind),  parameter, dimension(max_DOC) :: &  
         f_doc     = (/ 0.4_dbl_kind, 0.4_dbl_kind, 0.2_dbl_kind /), &
                      ! fraction of mort_N that goes to each doc pool
         f_exude    = (/c1, c1, c1 /) , & ! fraction of exuded carbon to each DOC pool
         k_bac      = (/0.03_dbl_kind,0.03_dbl_kind,0.03_dbl_kind/)    
                      ! (1/d) Bacterial degredation of DOC (1/month)
      real (kind=dbl_kind), parameter :: & 
         T_max      = c0            , & 
         fsal       = c1            , & ! 0.2  Arrigo and Sullivan 1992 (only for T < -5oC or S_b > 60 ppt)
         op_dep_min = 0.1           , & ! 0.01 Light attenuates for optical depths exceeding
                                        ! op_dep_min (small effect unless chla > 100 mg/m^2)
         fr_graze_s = 0.5_dbl_kind  , & ! fraction of grazing spilled or slopped
         fr_graze_e = 0.5_dbl_kind  , & ! fraction of assimilation excreted 
         fr_mort2min= 0.5_dbl_kind  , & ! fractionation of mortality to Am
         k_nitrif   = c0            , & !(1/day) nitrification rate   0.015_dbl_kind 
         t_iron_conv= 65.0_dbl_kind , & ! desorption loss pFe to dFe (Parekh et al, 2004 use 61 days)
         max_loss   = 0.9_dbl_kind  , & ! restrict uptake to 90% of remaining value 
         max_dfe_doc1 = 0.2_dbl_kind    ! 0.1852_dbl_kind 
                                        ! (nM Fe/muM C) max ratio of dFe to saccharides allowed in the ice
                                        ! tuned for the Southern Ocean with 16.2 muM of saccharides at ISPOL
                                        ! and a maximum dFe of 3 nM in Southern Ocean waters. 

      real (kind=dbl_kind), parameter :: &
         fr_resp_s  = 0.75_dbl_kind, &  ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS   = 0.5_dbl_kind , &  ! and conversion given high yield
         t_sk_conv  = 3.0_dbl_kind , &  ! at a Stefels rate (d)
         t_sk_ox    = 10.0_dbl_kind     ! DMS in turn oxidizes slowly (d)

      integer (kind=int_kind) :: i, j, ij, k, n

      real (kind=dbl_kind), dimension(max_algae) :: &
         Nin        , &     ! algal nitrogen concentration on volume (mmol/m^3) 
         Cin        , &     ! algal carbon concentration on volume (mmol/m^3)
         chlin              ! algal chlorophyll concentration on volume (mg/m^3)

      real (kind=dbl_kind), dimension(max_doc) :: &
         Docin              ! dissolved organic carbon concentration on volume (mmolC/m^3) 

      real (kind=dbl_kind), dimension(max_dic) :: &
         Dicin              ! dissolved inorganic carbon concentration on volume (mmolC/m^3) 

      real (kind=dbl_kind), dimension(max_don) :: &  !proteins
         Donin              ! dissolved organic nitrogen concentration on volume (mmolN/m^3) 

      real (kind=dbl_kind), dimension(max_fe) :: &  !iron
         Fedin              ! dissolved iron concentration on volume (umol/m^3) 

      real (kind=dbl_kind), dimension(max_fe) :: &  !iron
         Fepin              ! algal nitrogen concentration on volume (umol/m^3) 

      real (kind=dbl_kind) :: &
         Nitin      , &     ! nitrate concentration on volume (mmol/m^3) 
         Amin       , &     ! ammonia/um concentration on volume (mmol/m^3) 
         Silin      , &     ! silicon concentration on volume (mmol/m^3) 
         DMSPpin    , &     ! DMSPp concentration on volume (mmol/m^3)
         DMSPdin    , &     ! DMSPd concentration on volume (mmol/m^3)
         DMSin      , &     ! DMS concentration on volume (mmol/m^3)
         PONin      , &     ! PON concentration on volume (mmol/m^3)
         op_dep     , &     ! bottom layer attenuation exponent (optical depth)
         Iavg_loc           ! bottom layer attenuated Fswthru (W/m^2)

      real (kind=dbl_kind), dimension(max_algae) :: &
         L_lim    , &  ! overall light limitation 
         Nit_lim  , &  ! overall nitrate limitation
         Am_lim   , &  ! overall ammonium limitation
         N_lim    , &  ! overall nitrogen species limitation
         Sil_lim  , &  ! overall silicon limitation
         Fe_lim   , &  ! overall iron limitation
         fr_Nit   , &  ! fraction of local ecological growth as nitrate
         fr_Am    , &  ! fraction of local ecological growth as ammonia
         growmax_N, &  ! maximum growth rate in N currency (mmol/m^3/s)
         grow_N   , &  ! true growth rate in N currency (mmol/m^3/s)
         potU_Nit , &  ! potential nitrate uptake (mmol/m^3/s)
         potU_Am  , &  ! potential ammonium uptake (mmol/m^3/s)
         U_Nit    , &  ! actual nitrate uptake (mmol/m^3/s)
         U_Am     , &  ! actual ammonium uptake (mmol/m^3/s)
         U_Sil    , &  ! actual silicon uptake (mmol/m^3/s)
         U_Fe     , &  ! actual iron uptake   (umol/m^3/s)
         U_Nit_f  , &  ! fraction of Nit uptake due to each algal species
         U_Am_f   , &  ! fraction of Am uptake due to each algal species
         U_Sil_f  , &  ! fraction of Sil uptake due to each algal species
         U_Fe_f        ! fraction of Fe uptake due to each algal species

      real (kind=dbl_kind) :: &
         dTemp        , &  ! sea ice temperature minus sst (oC) < 0
         U_Nit_tot    , &  ! actual nitrate uptake (mmol/m^3/s)
         U_Am_tot     , &  ! actual ammonium uptake (mmol/m^3/s)
         U_Sil_tot    , &  ! actual silicon uptake (mmol/m^3/s)
         U_Fe_tot     , &  ! actual iron uptake   (umol/m^3/s)
         nitrif       , &  ! nitrification (mmol/m^3/s)
         mort_N       , &  ! total algal mortality (mmol N/m^3)
         mort_C       , &  ! total algal mortality (mmol C/m^3)
         graze_N      , &  ! total algae grazed (mmol N/m^3)
         graze_C      , &  ! total algae grazed (mmol C/m^3)
         exude_C      , &  ! total carbon exuded by algae (mmol C/m^3)
         resp_N       , &  ! total N in respiration (mmol N/m^3)
         growth_N     , &  ! total algal growth (mmol N/m^3)
         fr_graze_p   , &  ! fraction of N grazed that becomes protein
                           !  (rest is assimilated) < (1-fr_graze_a)
                           !  and fr_graze_a*fr_graze_e becomes ammonia
         fr_mort_p         ! fraction of N mortality that becomes protein 
                           ! < (1-fr_mort2min)

      real (kind=dbl_kind), dimension(max_algae) :: &
         resp     , &  ! respiration (mmol/m^3/s)
         graze    , &  ! grazing (mmol/m^3/s)
         mort          ! sum of mortality and excretion (mmol/m^3/s)

!  source terms underscore s, removal underscore r

      real (kind=dbl_kind), dimension(max_algae) :: &
         N_s       , &  ! net algal nitrogen sources (mmol/m^3)
         N_r            ! net algal nitrogen removal (mmol/m^3)

      real (kind=dbl_kind), dimension(max_doc) :: &
        DOC_r      , &  ! net DOC removal (mmol/m^3)
        DOC_s           ! net DOC sources (mmol/m^3)

      real (kind=dbl_kind), dimension(max_don) :: &
        DON_r      , &  ! net DON removal (mmol/m^3)
        DON_s           ! net DON sources (mmol/m^3)

      real (kind=dbl_kind), dimension(max_fe) :: &
        Fed_r_l     , &  ! removal due to loss of binding saccharids (nM)
        Fed_r       , &  ! net Fed removal (nM)
        Fed_s       , &  ! net Fed sources (nM)
        Fep_r       , &  ! net Fep removal (nM)
        Fep_s       , &  ! net Fep sources (nM)
        rFep        , &  ! ratio of particulate Fe to tot Fep
        rFed             ! ratio of dissolved Fe to tot Fed

      real (kind=dbl_kind) :: &
         dN        , &  ! change in N (mmol/m^3)
         N_s_p     , &  ! algal nitrogen photosynthesis (mmol/m^3)
         N_r_g     , &  ! algal nitrogen losses to grazing (mmol/m^3)
         N_r_r     , &  ! algal nitrogen losses to respiration (mmol/m^3)
         N_r_mo    , &  ! algal nitrogen losses to mortality (mmol/m^3)
         Nit_s_n   , &  ! nitrate from nitrification (mmol/m^3)
         Nit_s_r   , &  ! nitrate from respiration (mmol/m^3)
         Nit_r_p   , &  ! nitrate uptake by algae (mmol/m^3)
         Nit_s     , &  ! net nitrate sources (mmol/m^3)
         Nit_r     , &  ! net nitrate removal (mmol/m^3)
         Am_s_e    , &  ! ammonium source from excretion (mmol/m^3)
         Am_s_r    , &  ! ammonium source from respiration (mmol/m^3)
         Am_s_mo   , &  ! ammonium source from mort/remin (mmol/m^3) 
         Am_r_p    , &  ! ammonium uptake by algae (mmol/m^3)
         Am_r_n    , &  ! ammonium removal to nitrification (mmol/m^3)
         Am_s      , &  ! net ammonium sources (mmol/m^3)
         Am_r      , &  ! net ammonium removal (mmol/m^3)
         Sil_s_r   , &  ! silicon from respiration (mmol/m^3)
         Sil_r_p   , &  ! silicon uptake by algae (mmol/m^3)
         Sil_s     , &  ! net silicon sources (mmol/m^3)
         Sil_r     , &  ! net silicon removal (mmol/m^3)
         Fe_r_p    , &  ! iron uptake by algae  (nM)
         DOC_r_c   , &  ! net doc removal from bacterial consumption (mmol/m^3)
         doc_s_m   , &  ! protein source due to algal mortality (mmol/m^3)
         doc_s_g        ! protein source due to grazing (mmol/m^3)         

      real (kind=dbl_kind) :: &
         DMSPd_s_r , &  ! skl dissolved DMSP from respiration (mmol/m^3)
         DMSPd_s_mo, &  ! skl dissolved DMSP from MBJ algal mortexc (mmol/m^3)
         DMSPd_r   , &  ! skl dissolved DMSP conversion (mmol/m^3) DMSPD_sk_r
         DMSPd_s   , &  ! net skl dissolved DMSP sources (mmol/m^3)
         DMS_s_c   , &  ! skl DMS source via conversion (mmol/m^3)
         DMS_r_o   , &  ! skl DMS losses due to oxidation (mmol/m^3)
         DMS_s     , &  ! net skl DMS sources (mmol/m^3)
         DMS_r     , &  ! net skl DMS removal (mmol/m^3)
         Fed_tot   , &  ! total dissolved iron from all sources (nM)
         Fed_tot_r , &  ! total dissolved iron losses (nM)
         Fed_tot_s , &  ! total dissolved iron sources (nM)
         Fep_tot   , &  ! total particulate iron from all sources (nM)
         Fep_tot_r , &  ! total particulate iron losses (nM)
         Fep_tot_s , &  ! total particulate iron sources (nM)
         Zoo_s_a   , &  ! N Losses due to zooplankton assimilation (mmol/m^3)
         Zoo_s_s   , &  ! N Losses due to grazing spillage (mmol/m^3)
         Zoo_s_m   , &  ! N Losses due to algal mortality (mmol/m^3)
         Zoo_s_b        ! N losses due to bacterial recycling of DON (mmol/m^3)
    
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
    do ij = 1, icells
       i = indxi(ij)
       j = indxj(ij)

      !-----------------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------------

       conserve_N(i,j) = .true.
       Nin(:)     = c0
       Cin(:)     = c0
       chlin(:)   = c0
       DOCin(:)   = c0
       DICin(:)   = c0
       DONin(:)   = c0
       Fedin(:)   = c0
       Fepin(:)   = c0
       Nitin      = c0
       Amin       = c0
       Silin      = c0
       DMSPpin    = c0
       DMSPdin    = c0
       DMSin      = c0
       PONin      = c0 
       U_Am_tot   = c0
       U_Nit_tot  = c0
       U_Sil_tot  = c0
       U_Fe_tot   = c0
       U_Am_f(:)  = c0
       U_Nit_f(:) = c0
       U_Sil_f(:) = c0
       U_Fe_f(:)  = c0
       DOC_s(:)   = c0
       DOC_r(:)   = c0
       DOC_r_c    = c0
       nitrif     = c0 
       mort_N     = c0
       mort_C     = c0
       graze_N    = c0
       graze_C    = c0
       exude_C    = c0
       resp_N     = c0
       growth_N   = c0
       Nit_r      = c0 
       Am_s       = c0
       Am_r       = c0 
       Sil_r      = c0
       Fed_r(:)   = c0
       Fed_s(:)   = c0
       Fep_r(:)   = c0
       Fep_s(:)   = c0
       DMSPd_s    = c0 
       dTemp      = min(T_bot(ij)-T_max,c0)
       Fed_tot    = c0
       Fed_tot_r  = c0
       Fed_tot_s  = c0
       rFed(:)    = c0
       Fep_tot    = c0
       Fep_tot_r  = c0
       Fep_tot_s  = c0
       rFep(:)    = c0
     
       Nitin     = ltrcrn(ij,nlt_bgc_Nit)
       op_dep = c0
       do k = 1, n_algae
          Nin(k)   = ltrcrn(ij,nlt_bgc_N(k))
          chlin(k) = R_chl2N(k)* Nin(k)  
          op_dep = op_dep + chlabs(k)*chlin(k)
       enddo
       if (tr_bgc_C)   then
        ! do k = 1, n_algae
        !     Cin(k)=  ltrcrn(ij,nlt_bgc_C(k))
        ! enddo
         do k = 1, n_doc
             DOCin(k)= ltrcrn(ij,nlt_bgc_DOC(k))
         enddo
         do k = 1, n_dic
             DICin(k)= ltrcrn(ij,nlt_bgc_DIC(k))
         enddo
       endif
       if (tr_bgc_Am)        Amin     = ltrcrn(ij,nlt_bgc_Am)
       if (tr_bgc_Sil)       Silin    = ltrcrn(ij,nlt_bgc_Sil)
       if (tr_bgc_DMS) then
        !       DMSPpin  = ltrcrn(ij,nlt_bgc_DMSPp)
             DMSPdin  = ltrcrn(ij,nlt_bgc_DMSPd)
             DMSin    = ltrcrn(ij,nlt_bgc_DMS) 
       endif
       if (tr_bgc_PON)       PONin    = ltrcrn(ij,nlt_bgc_PON) 
       if (tr_bgc_DON) then
         do k = 1, n_don
             DONin(k) = ltrcrn(ij,nlt_bgc_DON(k))
         enddo
       endif
       if (tr_bgc_Fe ) then
         do k = 1, n_fed 
             Fedin(k) = ltrcrn(ij,nlt_bgc_Fed(k))
         enddo
         do k = 1, n_fep 
             Fepin(k) = ltrcrn(ij,nlt_bgc_Fep(k))
         enddo
       endif

      !-----------------------------------------------------------------------
      ! Total iron from all pools
      !-----------------------------------------------------------------------

       do k = 1,n_fed
         Fed_tot = Fed_tot + Fedin(k)
       enddo
       do k = 1,n_fep
         Fep_tot = Fep_tot + Fepin(k)
       enddo
       if (Fed_tot > puny) then
       do k = 1,n_fed
         rFed(k) = Fedin(k)/Fed_tot
       enddo
       endif
       if (Fep_tot > puny) then
       do k = 1,n_fep
         rFep(k) = Fepin(k)/Fep_tot
       enddo
       endif

      !-----------------------------------------------------------------------
      ! Light limitation  (op_dep) defined above
      !-----------------------------------------------------------------------

       if (op_dep > op_dep_min .and. .not. dEdd_algae) then
         Iavg_loc = fswthru(i,j) * (c1 - exp(-op_dep)) / op_dep
       else
         Iavg_loc = fswthru(i,j)
       endif

       do k = 1, n_algae
          ! With light inhibition ! Maybe include light inhibition for diatoms but phaeocystis

           L_lim = (c1 - exp(-alpha2max_low(k)*Iavg_loc)) * exp(-beta2max(k)*Iavg_loc)      

          ! Without light inhibition
          !L_lim(k) = (c1 - exp(-alpha2max_low(k)*Iavg_loc)) 

      !-----------------------------------------------------------------------
      ! Nutrient limitation
      !-----------------------------------------------------------------------

          Nit_lim(k) = Nitin/(Nitin + K_Nit(k))
          Am_lim(k)  = c0
          N_lim(k) = Nit_lim(k)
          if (tr_bgc_Am) then
             Am_lim(k) = Amin/(Amin + K_Am(k))
             N_lim(k)  = min(c1, Nit_lim(k) + Am_lim(k))  
          endif
          Sil_lim(k) = c1
          if (tr_bgc_Sil .and. K_Sil(k) > c0) Sil_lim(k) = Silin/(Silin + K_Sil(k))

      !-----------------------------------------------------------------------
      ! Iron limitation
      !-----------------------------------------------------------------------

          Fe_lim(k) = c1         
          if (tr_bgc_Fe  .and. K_Fe (k) > c0) Fe_lim (k) = Fed_tot/(Fed_tot + K_Fe(k))
  
      !----------------------------------------------------------------------------
      ! Growth and uptake computed within the bottom layer 
      ! Note here per A93 discussions and MBJ model, salinity is a universal 
      ! restriction.  Comparison with available column nutrients inserted 
      ! but in tests had no effect.
      ! Primary production reverts to SE form, see MBJ below and be careful
      !----------------------------------------------------------------------------

          growmax_N(k) = mu_max(k) / secday * exp(grow_Tdep(k) * dTemp)* Nin(k) *fsal
          grow_N(k)    = min(L_lim(k), N_lim(k), Sil_lim(k), Fe_lim(k)) * growmax_N(k)
          potU_Nit(k)  = Nit_lim(k)* growmax_N(k)
          potU_Am(k)   = Am_lim(k)* growmax_N(k) 
          U_Am(k)      = min(grow_N(k), potU_Am(k))
          U_Nit(k)     = grow_N(k) - U_Am(k)
          U_Sil(k)     = R_Si2N(k) * grow_N(k)
          U_Fe (k)     = R_Fe2N(k) * grow_N(k)

          U_Am_tot     = U_Am_tot  + U_Am(k)
          U_Nit_tot    = U_Nit_tot + U_Nit(k)
          U_Sil_tot    = U_Sil_tot + U_Sil(k)
          U_Fe_tot     = U_Fe_tot  + U_Fe(k)
       enddo
       do k = 1, n_algae
          if (U_Am_tot > c0) U_Am_f(k) = U_Am(k)/U_Am_tot
          if (U_Nit_tot > c0) U_Nit_f(k) = U_Nit(k)/U_Nit_tot
          if (U_Sil_tot > c0) U_Sil_f(k) = U_Sil(k)/U_Sil_tot
          if (U_Fe_tot > c0) U_Fe_f(k) = U_Fe(k)/U_Fe_tot
       enddo

       if (tr_bgc_Sil) U_Sil_tot = min(U_Sil_tot, max_loss * Silin/dt)
       if (tr_bgc_Fe)  U_Fe_tot  = min(U_Fe_tot, max_loss * Fed_tot/dt)
       U_Nit_tot = min(U_Nit_tot, max_loss * Nitin/dt)  
       U_Am_tot  = min(U_Am_tot,  max_loss * Amin/dt)    

       do k = 1, n_algae
          U_Am(k)  = U_Am_f(k)*U_Am_tot
          U_Nit(k) = U_Nit_f(k)*U_Nit_tot
          U_Sil(k) = U_Sil_f(k)*U_Sil_tot
          U_Fe(k)  = U_Fe_f(k)*U_Fe_tot

          if (R_Si2N(k) > c0) then
             grow_N(k) = min(U_Sil(k)/R_Si2N(k),U_Nit(k) + U_Am(k), U_Fe(k)/R_Fe2N(k))
          else
             grow_N(k) = min(U_Nit(k) + U_Am(k),U_Fe(k)/R_Fe2N(k))
          endif

          fr_Am(k) = c0
          if (tr_bgc_Am) then
             fr_Am(k) = p5
             if (grow_N(k) > c0) fr_Am(k) = min(U_Am(k)/grow_N(k), c1)
          endif
          fr_Nit(k) = c1 - fr_Am(k)
          U_Nit(k)  = fr_Nit(k) * grow_N(k)
          U_Am(k)   = fr_Am(k)  * grow_N(k)
          U_Sil(k)  = R_Si2N(k) * grow_N(k)
          U_Fe (k)  = R_Fe2N(k) * grow_N(k)
    
      !-----------------------------------------------------------------------
      ! Define reaction terms
      !-----------------------------------------------------------------------

      ! Since the framework remains incomplete at this point,
      ! it is assumed as a starting expedient that 
      ! DMSP loss to melting results in 10% conversion to DMS
      ! which is then given a ten day removal constant.
      ! Grazing losses are channeled into rough spillage and assimilation
      ! then following ammonia there is some recycling.

      !--------------------------------------------------------------------
      ! Algal reaction term
      ! N_react = (grow_N*(c1 - fr_graze-fr_resp) - mort)*dt  
      !--------------------------------------------------------------------

          resp(k)   = fr_resp  * grow_N(k)  
          graze(k)  = fr_graze(k) * grow_N(k)
          mort(k)   = min(max_loss * Nin(k)/dt, mort_pre(k) * exp(mort_Tdep(k)*dTemp)  * Nin(k) / secday)
 
        ! history variables
          grow_alg(i,j,k) = grow_N(k)
          upNOn(i,j,k) = U_Nit(k)
          upNHn(i,j,k) = U_Am(k)

          N_s_p  = grow_N(k) * dt  
          N_r_g  = graze(k)  * dt 
          N_r_r  = resp(k)   * dt
          N_r_mo = mort(k)   * dt
          N_s(k)    = (c1- fr_resp - fr_graze(k)) * grow_N(k) *dt   !N_s_p
          N_r(k)    = mort(k) * dt                                  !N_r_g  + N_r_mo + N_r_r 

          graze_N   = graze_N + graze(k)
          graze_C   = graze_C + R_C2N(k)*graze(k)
          mort_N    = mort_N + mort(k)      
          mort_C    = mort_C + R_C2N(k)*mort(k)
          resp_N    = resp_N + resp(k)
          growth_N  = growth_N + grow_N(k)
 
      enddo ! n_algae
      !--------------------------------------------------------------------
      ! Ammonium source: algal grazing, respiration, and mortality
      !--------------------------------------------------------------------

          Am_s_e  = graze_N*(c1-fr_graze_s)*fr_graze_e*dt
          Am_s_mo = mort_N*fr_mort2min*dt
          Am_s_r  = resp_N*dt
          Am_s    = Am_s_r + Am_s_e + Am_s_mo

      !--------------------------------------------------------------------
      ! Nutrient net loss terms: algal uptake
      !--------------------------------------------------------------------
        
       do n = 1, n_algae
          Am_r_p  = U_Am(n)   * dt
          Am_r    = Am_r + Am_r_p 
          Nit_r_p = U_Nit(n)  * dt                
          Nit_r   = Nit_r + Nit_r_p 
          Sil_r_p = U_Sil(n) * dt
          Sil_r   = Sil_r + Sil_r_p 
          Fe_r_p  = U_Fe (n) * dt
          Fed_tot_r = Fed_tot_r + Fe_r_p  
          exude_C = exude_C + k_exude(n)* R_C2N(k)*Nin(k) / secday 
       enddo

      !--------------------------------------------------------------------
      ! nitrification
      !--------------------------------------------------------------------

       nitrif  = k_nitrif /secday * Amin 
       Am_r = Am_r +  nitrif*dt
       Nit_s_n = nitrif * dt  !source from NH4
       Nit_s   = Nit_s_n  

      !--------------------------------------------------------------------
      ! PON:  currently using PON to shadow nitrate
      !
      ! N Losses are counted in Zoo.  These arise from mortality not 
      ! remineralized (Zoo_s_m), assimilated grazing not excreted (Zoo_s_a), 
      !spilled N not going to DON (Zoo_s_s) and  bacterial recycling 
      ! of DON (Zoo_s_b). 
      !--------------------------------------------------------------------

       if (tr_bgc_Am) then
         Zoo_s_a = graze_N*(c1-fr_graze_e)*(c1-fr_graze_s) *dt
         Zoo_s_s = graze_N*fr_graze_s*dt
         Zoo_s_m = mort_N*dt  -  Am_s_mo
       else
         Zoo_s_a = graze_N*dt*(c1-fr_graze_s)
         Zoo_s_s = graze_N*fr_graze_s*dt
         Zoo_s_m = mort_N*dt 
       endif

         Zoo_s_b = c0

      !--------------------------------------------------------------------
      ! DON (n_don = 1)
      ! Proteins   
      !--------------------------------------------------------------------

       DON_r(:) = c0
       DON_s(:) = c0

       if (tr_bgc_DON) then
       do n = 1, n_don   
          DON_r(n) =  kn_bac(n)/secday * DONin(n) * dt
          DON_s(n) =  graze_N*f_don(n)*fr_graze_s * dt
          Zoo_s_s = Zoo_s_s - DON_s(n)
          Zoo_s_b = Zoo_s_b + DON_r(n)*(c1-f_don_Am(n))
        !  Am_s = Am_s + DON_r(n)*f_don_Am(n)
      enddo
      endif
     
       Zoo(i,j) = Zoo_s_a + Zoo_s_s + Zoo_s_m + Zoo_s_b

      !--------------------------------------------------------------------
      ! DOC
      ! polysaccharids, lipids
      !--------------------------------------------------------------------

       do n = 1, n_doc   
          
          DOC_r(n) =  k_bac(n)/secday * DOCin(n) * dt
          DOC_s(n) =  f_doc(n)*(fr_graze_s *graze_C + mort_C)*dt &
                      + f_exude(n)*exude_C
      enddo

      !--------------------------------------------------------------------
      ! Iron sources from remineralization  (follows ammonium)
      ! only Fed_s(1)  has remineralized sources
      !--------------------------------------------------------------------
      
      Fed_s(1) = Fed_s(1) + Am_s * R_Fe2N(1)    ! remineralization source

      !--------------------------------------------------------------------
      !  Conversion to dissolved Fe from Particulate requires DOC(1)
      !  Otherwise the only source of dFe is from remineralization
      !--------------------------------------------------------------------

      if (tr_bgc_C .and. DOCin(1) > c0) then
         
        if (Fed_tot/DOCin(1) > max_dfe_doc1) then             
          do n = 1,n_fed                                    ! low saccharid:dFe ratio leads to 
             Fed_r_l(n)  = Fedin(n)/t_iron_conv*dt/secday   ! loss of bioavailable Fe to particulate fraction
             Fep_tot_s   = Fep_tot_s + Fed_r_l(n)
             Fed_r(n)    = rFed(n) * Fed_tot_r + Fed_r_l(n) ! removal includes uptake and coagulation
          enddo  
          do n = 1,n_fep
             Fep_s(n) = rFep(n)* Fep_tot_s                  ! source from dissolved Fe 
          enddo
        elseif (Fed_tot/DOCin(1) < max_dfe_doc1) then  
          do n = 1,n_fep                                    ! high saccharid:dFe ratio leads to 
             Fep_r(n)  = Fepin(n)/t_iron_conv*dt/secday     ! gain of bioavailable Fe from particulate fraction
             Fed_tot_s = Fed_tot_s + Fep_r(n)
          enddo  
          do n = 1,n_fed
             Fed_s(n) = Fed_s(n) + rFed(n)* Fed_tot_s       ! source from particulate Fe
             Fed_r(n) = rFed(n)* Fed_tot_r                  ! algal uptake
          enddo    
       endif         
      endif
      do n = 1,n_fep
         Fep_s(n) = rFep(n)* Zoo(i,j) + Fep_s(n)           ! source from algal mortality/grazing
      enddo ! losses not direct to Fed 

      !--------------------------------------------------------------------
      ! Sulfur cycle begins here
      !--------------------------------------------------------------------
      ! Grazing losses are channeled into rough spillage and assimilation
      ! then onward and the MBJ mortality channel is included
      ! It is assumed as a starting expedient that 
      ! DMSP loss to melting gives partial conversion to DMS in product layer
      ! which then undergoes Stefels removal.

      !--------------------------------------------------------------------
      ! DMSPd  reaction term  (DMSPd conversion is outside the algal loop)
      ! DMSPd_react = R_S2N*((fr_graze_s+fr_excrt_2S*fr_graze_e*fr_graze_a)
      !                      *fr_graze*grow_N + fr_mort2min*mort)*dt
      !             - [\DMSPd]/t_sk_conv*dt
      !--------------------------------------------------------------------

       DMSPd_s_r = fr_resp_s  * R_S2N(k) * N_r_r   !respiration fraction to DMSPd
       DMSPd_s_mo= fr_mort2min * R_S2N(k)* N_r_mo  !mortality and extracellular excretion

       DMSPd_s = DMSPd_s + DMSPd_s_r + DMSPd_s_mo 
       DMSPd_r = (c1/t_sk_conv) * (c1/secday)  * (DMSPdin) * dt

      !--------------------------------------------------------------------
      ! DMS reaction term + DMSPd loss term 
      ! DMS_react = ([\DMSPd]*y_sk_DMS/t_sk_conv - c1/t_sk_ox *[\DMS])*dt
      !--------------------------------------------------------------------

       DMS_s_c = y_sk_DMS * DMSPd_r 
       DMS_r_o = DMSin * dt / (t_sk_ox * secday)
       DMS_s   = DMS_s_c
       DMS_r   = DMS_r_o

      !-----------------------------------------------------------------------
      ! Load reaction array
      !-----------------------------------------------------------------------

       dN = c0
       do k = 1,n_algae
              reactb(ij,nlt_bgc_N(k))  = N_s(k) - N_r(k)
              dN = dN + reactb(ij,nlt_bgc_N(k))
       enddo
       if (tr_bgc_C) then
        ! do k = 1,n_algae
        !      reactb(ij,nlt_bgc_C(k))  = R_C2N(k)*reactb(ij,nlt_bgc_N(k))
        ! enddo
         do k = 1,n_doc
              reactb(ij,nlt_bgc_DOC(k))= DOC_s(k) - DOC_r(k)  
         enddo
       endif
              reactb(ij,nlt_bgc_Nit)   = Nit_s   - Nit_r
              dN = dN + reactb(ij,nlt_bgc_Nit)
       if (tr_bgc_Am)  then
              reactb(ij,nlt_bgc_Am)    = Am_s    - Am_r
              dN = dN + reactb(ij,nlt_bgc_Am)
       endif
       if (tr_bgc_Sil) then
              reactb(ij,nlt_bgc_Sil)   = Sil_s   - Sil_r
       endif
       if (tr_bgc_DON) then
         do k = 1,n_don
              reactb(ij,nlt_bgc_DON(k))= DON_s(k) - DON_r(k)  
              dN = dN + reactb(ij,nlt_bgc_DON(k))
         enddo
       endif 
       if (tr_bgc_Fe ) then
        do k = 1,n_fed
              reactb(ij,nlt_bgc_Fed(k))= Fed_s (k) - Fed_r (k) 
        enddo
        do k = 1,n_fep
              reactb(ij,nlt_bgc_Fep(k))= Fep_s (k) - Fep_r (k) 
        enddo
       endif 
       if (tr_bgc_DMS) then
              reactb(ij,nlt_bgc_DMSPd) = DMSPd_s - DMSPd_r
              reactb(ij,nlt_bgc_DMS)   = DMS_s   - DMS_r
       endif
       Nerror(i,j) = dN + Zoo(i,j)
      ! if (abs(Nerror(i,j)) > max(reactb(ij,:))*1,0e-5) then
      !      conserve_N(i,j) = .false.
      !      write (nu_diag, *) 'Conservation error!'
      !      write (nu_diag, *) 'Nerror(i,j),dN, DONin(1),kn_bac(1),secday,dt,n_doc'
      !      write (nu_diag, *) Nerror(i,j),dN, DONin(1),kn_bac(1),secday,dt,n_doc
      !      write (nu_diag, *) 'reactb(ij,nlt_bgc_Nit),reactb(ij,nlt_bgc_N(1)),reactb(ij,nlt_bgc_N(2)'
      !      write (nu_diag, *) reactb(ij,nlt_bgc_Nit),reactb(ij,nlt_bgc_N(1)),reactb(ij,nlt_bgc_N(2))
      !      write (nu_diag, *) 'reactb(ij,nlt_bgc_Am),reactb(ij,nlt_bgc_DON(1)), DON_r(1),DON_s(1)'
      !      write (nu_diag, *) reactb(ij,nlt_bgc_Am),reactb(ij,nlt_bgc_DON(1)),DON_r(1),DON_s(1)
      !      write (nu_diag, *) 'Zoo(i,j):',Zoo(i,j)
      ! endif
          
      enddo   !ij

      end subroutine algal_dyn

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
!          Nicole Jeffery, LANL

      subroutine bgc_diags (dt)

      use ice_broadcast, only: broadcast_scalar
      use ice_constants, only: c0, mps_to_cmpdy, c100, p5, c1, secday
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, pbloc
      use ice_timers, only: timer_bgc, ice_timer_start, ice_timer_stop
      use ice_domain_size, only: ncat, nltrcr, n_algae, max_algae, n_zaero, max_aero, &
                   n_dic, max_dic, n_doc, max_doc, n_don, max_don, n_fed, n_fep, max_fe
      use ice_state, only:aice, vicen, vice, trcr, nt_bgc_N, nt_bgc_C, nt_bgc_chl, nt_bgc_Am, &
          nt_bgc_DMS, nt_bgc_DMSPd, nt_bgc_DMSPp, nt_bgc_Nit, nt_bgc_Sil, &
          nt_bgc_PON, nt_bgc_DON, nt_bgc_DIC, nt_bgc_DOC, nt_zaero, nt_bgc_Fed, &
          nt_bgc_Fep

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, nn, ii,jj, iblk,kk
      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         pNit_sk, pAm_sk, pSil_sk, &
         pDMSPp_sk, pDMSPd_sk, pDMS_sk, &
         pNit_ac, pAm_ac, pSil_ac, pDMSP_ac, pDMS_ac, &
         pflux_NO, pflux_Am,  &
         pflux_snow_NO, pflux_snow_Am,  &
         pflux_atm_NO, pflux_atm_Am,  pgrow_net

      real (kind=dbl_kind), dimension(npnt,max_algae) :: &
         pN_ac, pN_tot, pN_sk, pflux_N
      real (kind=dbl_kind), dimension(npnt,max_doc) :: &
         pDOC_ac, pDOC_sk
      real (kind=dbl_kind), dimension(npnt,max_don) :: &
         pDON_ac, pDON_sk
      real (kind=dbl_kind), dimension(npnt,max_fe ) :: &
         pFed_ac,  pFed_sk, pFep_ac, pFep_sk 
      real (kind=dbl_kind), dimension(npnt,max_aero) :: &
        pflux_zaero, pflux_snow_zaero, pflux_atm_zaero, &
        pflux_atm_zaero_s

      ! vertical  fields of category 1 at diagnostic points for bgc layer model
      real (kind=dbl_kind), dimension(npnt,2) :: &
         pNOs, pAms, pPONs
      real (kind=dbl_kind), dimension(npnt,2,max_algae) :: &
         pNs
      real (kind=dbl_kind), dimension(npnt,2,max_doc) :: &
         pDOCs
      real (kind=dbl_kind), dimension(npnt,2,max_don) :: &
         pDONs
      real (kind=dbl_kind), dimension(npnt,2,max_fe ) :: &
         pFeds, pFeps 
      real (kind=dbl_kind), dimension(npnt,2,max_aero) :: &
         pzaeros
      real (kind=dbl_kind), dimension(npnt,nblyr+1) :: &
         pNO, pAm, pPON, pzfswin, pZoo  
      real (kind=dbl_kind), dimension(npnt,nblyr+1,max_algae) :: &
         pN
      real (kind=dbl_kind), dimension(npnt,nblyr+1,max_aero) :: &
         pzaero
      real (kind=dbl_kind), dimension(npnt,nblyr+1,max_doc) :: &
         pDOC
      real (kind=dbl_kind), dimension(npnt,nblyr+1,max_don) :: &
         pDON
      real (kind=dbl_kind), dimension(npnt,nblyr+1,max_fe ) :: &
         pFed, pFep 
      real (kind=dbl_kind), dimension (nblyr+1) :: & 
         zspace

      zspace(:) = c1/real(nblyr,kind=dbl_kind)
      zspace(1) = zspace(1)*p5
      zspace(nblyr+1) = zspace(nblyr+1)*p5    

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
               pAm_ac(n)   = c0
               pSil_ac(n)  = c0
               pDMSP_ac(n) = c0
               pDMS_ac(n)  = c0
               pN_ac(n,:)  = c0
               pDOC_ac(n,:)= c0
               pDON_ac(n,:)= c0
               pFed_ac (n,:)= c0
               pFep_ac (n,:)= c0
               pNit_ac(n) = c0
               if (tr_bgc_N) then
                 do k = 1,n_algae
                    pN_ac(n,k)    = ocean_bio(i,j,nlt_bgc_N(k),iblk) 
                 enddo  !n_algae
               endif    !tr_bgc_N
               if (tr_bgc_C) then
                 do k = 1,n_doc
                    pDOC_ac(n,k)    = ocean_bio(i,j,nlt_bgc_DOC(k),iblk) 
                 enddo  !n_algae
               endif    !tr_bgc_N
               if (tr_bgc_DON) then
                 do k = 1,n_don
                    pDON_ac(n,k)    = ocean_bio(i,j,nlt_bgc_DON(k),iblk) 
                 enddo 
               endif
               if (tr_bgc_Fe ) then
                 do k = 1,n_fed 
                    pFed_ac (n,k)   = ocean_bio(i,j,nlt_bgc_Fed (k),iblk) 
                 enddo 
                 do k = 1,n_fep 
                    pFep_ac (n,k)   = ocean_bio(i,j,nlt_bgc_Fep (k),iblk) 
                 enddo 
               endif
               if (tr_bgc_Nit) &
               pNit_ac(n)  = ocean_bio(i,j,nlt_bgc_Nit,iblk)  ! nit(i,j,iblk)
               if (tr_bgc_Am) &
               pAm_ac(n)   = ocean_bio(i,j,nlt_bgc_Am,iblk)   ! amm(i,j,iblk)
               if (tr_bgc_Sil) &
               pSil_ac(n)  = ocean_bio(i,j,nlt_bgc_Sil,iblk)  ! sil(i,j,iblk)
               if (tr_bgc_DMS) then
               pDMSP_ac(n) = ocean_bio(i,j,nlt_bgc_DMSPp,iblk)! dmsp(i,j,iblk)
               pDMS_ac(n)  = ocean_bio(i,j,nlt_bgc_DMS,iblk)  ! dms(i,j,iblk)
               endif

               ! fluxes in mmol/m^2/d
               ! concentrations are bulk in mmol/m^3
               ! iron is in 10^-3 mmol/m^3  (nM)

               if (skl_bgc) then
                  pNit_sk(n)   = c0
                  pAm_sk(n)    = c0
                  pSil_sk(n)   = c0
                  pDMSPp_sk(n) = c0
                  pDMSPd_sk(n) = c0
                  pDMS_sk(n)   = c0
                  pN_sk(n,:)   = c0
                  pflux_N(n,:) = c0
                  pDOC_sk(n,:) = c0
                  pDON_sk(n,:) = c0
                  pFed_sk(n,:) = c0
                  pFep_sk(n,:) = c0
                  
                  do k = 1,n_algae            
                    pN_sk(n,k)       = trcr    (i,j,nt_bgc_N(k),   iblk)
                    pflux_N(n,k)     = flux_bio(i,j,nlt_bgc_N(k),  iblk)*mps_to_cmpdy/c100 
                  enddo
                  if (tr_bgc_C) then
                     do k = 1,n_doc
                        pDOC_sk(n,k)       = trcr    (i,j,nt_bgc_DOC(k),   iblk)
                     enddo
                  endif
                  if (tr_bgc_DON) then
                     do k = 1,n_don
                        pDON_sk(n,k)       = trcr    (i,j,nt_bgc_DON(k),   iblk)
                     enddo
                  endif
                  if (tr_bgc_Fe ) then
                     do k = 1,n_fed 
                        pFed_sk (n,k)       = trcr    (i,j,nt_bgc_Fed(k),   iblk)
                     enddo
                     do k = 1,n_fep 
                        pFep_sk (n,k)       = trcr    (i,j,nt_bgc_Fep(k),   iblk)
                     enddo
                  endif
                  if (tr_bgc_Nit) then
                     pNit_sk(n)  = trcr    (i,j, nt_bgc_Nit, iblk) 
                     pflux_NO(n) = flux_bio(i,j,nlt_bgc_Nit, iblk)*mps_to_cmpdy/c100 
                  endif
                  if (tr_bgc_Am) then
                     pAm_sk(n)   = trcr    (i,j, nt_bgc_Am,  iblk)
                     pflux_Am(n)= flux_bio(i,j,nlt_bgc_Am,   iblk)*mps_to_cmpdy/c100 
                  endif
                  if (tr_bgc_Sil) then
                     pSil_sk(n)  = trcr    (i,j, nt_bgc_Sil, iblk) 
                  endif
                  if (tr_bgc_DMS) then
                     pDMSPp_sk(n) = trcr   (i,j,nt_bgc_DMSPp,iblk)
                     pDMSPd_sk(n) = trcr   (i,j,nt_bgc_DMSPd,iblk)
                     pDMS_sk  (n) = trcr   (i,j,nt_bgc_DMS,  iblk)
                   endif

               elseif (z_tracers) then   ! zbgc
                 pflux_NO(n) = c0
                 pN_tot(n,:) = c0
                 pflux_Am(n) = c0
                 pflux_atm_Am(n) = c0
                 pflux_snow_Am(n) = c0
                 pflux_N(n,:) = c0
                 pflux_NO(n) = c0
                 pflux_atm_NO(n) = c0
                 pflux_snow_NO(n) = c0
                 pflux_zaero(n,:) = c0
                 pflux_atm_zaero_s(n,:) = c0
                 pflux_atm_zaero(n,:) = c0
                 pflux_snow_zaero(n,:) = c0
                 if (tr_bgc_Nit) then
                    pflux_NO(n)         =   flux_bio(i,j,nlt_bgc_Nit,iblk)*mps_to_cmpdy/c100 
                    pflux_atm_NO(n)     =   fbio_atmice(i,j,nlt_bgc_Nit,iblk)*mps_to_cmpdy/c100 
                    pflux_snow_NO(n)    =   fbio_snoice(i,j,nlt_bgc_Nit,iblk)*mps_to_cmpdy/c100
                 endif
                 if (tr_bgc_Am) then
                    pflux_Am(n)       = flux_bio(i,j,nlt_bgc_Am,iblk)*mps_to_cmpdy/c100 
                    pflux_atm_Am(n)   = fbio_atmice(i,j,nlt_bgc_Am,iblk)*mps_to_cmpdy/c100 
                    pflux_snow_Am(n)  = fbio_snoice(i,j,nlt_bgc_Am,iblk)*mps_to_cmpdy/c100
                 endif 
                 if (tr_bgc_N)  then
                   do k = 1,n_algae
                     pflux_N(n,k)      = flux_bio(i,j,nlt_bgc_N(k),iblk)*mps_to_cmpdy/c100 
                   enddo
                 endif
                 if (tr_zaero)  then
                   do k = 1,n_zaero
                     pflux_zaero(n,k)      = flux_bio(i,j,nlt_zaero(k),iblk)*mps_to_cmpdy/c100 
                     pflux_atm_zaero_s(n,k)= flux_bio_atm(i,j,nlt_zaero(k),iblk)*mps_to_cmpdy/c100 
                     pflux_atm_zaero(n,k)  = fbio_atmice(i,j,nlt_zaero(k),iblk)*mps_to_cmpdy/c100 
                     pflux_snow_zaero(n,k) = fbio_snoice(i,j,nlt_zaero(k),iblk)*mps_to_cmpdy/c100
                   enddo
                 endif

                 do k = 1, nblyr+1
                   pzfswin(n,k) = c0
                   pZoo(n,k) = c0
                   do nn = 1,ncat
                     pzfswin(n,k) = pzfswin(n,k) +  zfswin(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                     pZoo(n,k) = pZoo(n,k) + Zoo(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                   enddo !nn
                   if (vice(i,j,iblk) > c0) then
                     pzfswin(n,k) = pzfswin(n,k)/vice(i,j,iblk)
                     pZoo(n,k)    = pZoo(n,k)/vice(i,j,iblk)
                   endif !vice
                   pAm(n,k) = c0
                   pN(n,k,:) = c0
                   pDOC(n,k,:) = c0
                   pDON(n,k,:) = c0
                   pFed(n,k,:) = c0
                   pFep(n,k,:) = c0
                   pzaero(n,k,:) = c0
                   pPON(n,k) = c0
                   pNO(n,k) = c0
                   if (tr_bgc_Nit) pNO(n,k) =  trcr(i,j,nt_bgc_Nit+k-1,iblk)                   
                   if (tr_bgc_Am) pAm(n,k) =  trcr(i,j,nt_bgc_Am+k-1,iblk)     
                   if (tr_bgc_N) then
                     do nn = 1, n_algae
                        pN(n,k,nn)   =  trcr(i,j,nt_bgc_N(nn)+k-1,iblk)
                     enddo   
                   endif     
                   if (tr_bgc_C) then
                     do nn = 1, n_doc
                        pDOC(n,k,nn)   =  trcr(i,j,nt_bgc_DOC(nn)+k-1,iblk)
                     enddo   
                   endif   
                   if (tr_bgc_DON) then
                     do nn = 1, n_don
                        pDON(n,k,nn)   =  trcr(i,j,nt_bgc_DON(nn)+k-1,iblk)
                     enddo   
                   endif    
                   if (tr_bgc_Fe)  then
                     do nn = 1, n_fed
                        pFed(n,k,nn)   =  trcr(i,j,nt_bgc_Fed(nn)+k-1,iblk)
                     enddo   
                     do nn = 1, n_fep
                        pFep(n,k,nn)   =  trcr(i,j,nt_bgc_Fep(nn)+k-1,iblk)
                     enddo   
                   endif   
                   if (tr_zaero) then
                     do nn = 1, n_zaero
                        pzaero(n,k,nn)   =  trcr(i,j,nt_zaero(nn)+k-1,iblk)
                     enddo   
                   endif
                   if (tr_bgc_PON) pPON(n,k) =  trcr(i,j,nt_bgc_PON+k-1,iblk)
                 enddo  !k
                 if (tr_bgc_N) then
                   do nn = 1,n_algae
                      pN_tot(n,nn) = ice_bio_net(i,j,nlt_bgc_N(nn),iblk)
                   enddo
                   pgrow_net(n) =  grow_net(i,j,iblk)
                 endif !tr_bgc_N
                 do k = 1,2  !snow concentration
                   pAms(n,k) = c0
                   pNs(n,k,:) = c0
                   pDOCs(n,k,:) = c0
                   pDONs(n,k,:) = c0
                   pFeds (n,k,:)= c0
                   pFeps (n,k,:)= c0
                   pzaeros(n,k,:) = c0
                   pPONs(n,k) = c0
                   pNOs(n,k) = c0
                   if (tr_bgc_Nit) pNOs(n,k)  = trcr(i,j,nt_bgc_Nit+nblyr+k,iblk)  
                   if (tr_bgc_Am) pAms(n,k) = trcr(i,j,nt_bgc_Am+nblyr+k,iblk)
                   if (tr_bgc_N) then
                     do nn = 1, n_algae
                       pNs(n,k,nn) =  trcr(i,j,nt_bgc_N(nn)+nblyr+k,iblk)
                     enddo
                   endif
                   if (tr_bgc_C) then
                     do nn = 1, n_doc
                       pDOCs(n,k,nn) =  trcr(i,j,nt_bgc_DOC(nn)+nblyr+k,iblk)
                     enddo
                   endif
                   if (tr_bgc_DON) then
                     do nn = 1, n_don
                       pDONs(n,k,nn) =  trcr(i,j,nt_bgc_DON(nn)+nblyr+k,iblk)
                     enddo
                   endif
                   if (tr_bgc_Fe ) then
                     do nn = 1, n_fed 
                       pFeds(n,k,nn) =  trcr(i,j,nt_bgc_Fed(nn)+nblyr+k,iblk)
                     enddo
                     do nn = 1, n_fep 
                       pFeps(n,k,nn) =  trcr(i,j,nt_bgc_Fep(nn)+nblyr+k,iblk)
                     enddo
                   endif
                   if (tr_zaero) then
                     do nn = 1, n_zaero
                       pzaeros(n,k,nn) =  trcr(i,j,nt_zaero(nn)+nblyr+k,iblk)
                     enddo
                   endif
                   if (tr_bgc_PON)pPONs(n,k) =trcr(i,j,nt_bgc_PON+nblyr+k,iblk)
                 enddo   !k 
               endif
             endif                 ! my_task = pmloc
            
             call broadcast_scalar(pNit_ac  (n), pmloc(n))             
             call broadcast_scalar(pAm_ac   (n), pmloc(n))             
             call broadcast_scalar(pSil_ac  (n), pmloc(n))             
             call broadcast_scalar(pDMSP_ac (n), pmloc(n))             
             call broadcast_scalar(pDMS_ac  (n), pmloc(n))  
             call broadcast_scalar(pflux_NO (n), pmloc(n))             
             call broadcast_scalar(pflux_Am(n), pmloc(n))

            do k = 1,n_algae
             call broadcast_scalar(pN_ac    (n,k), pmloc(n)) 
             call broadcast_scalar(pflux_N  (n,k), pmloc(n))    
            enddo
            do k = 1,n_doc
             call broadcast_scalar(pDOC_ac    (n,k), pmloc(n))   
            enddo
            do k = 1,n_don
             call broadcast_scalar(pDON_ac    (n,k), pmloc(n))   
            enddo
            do k = 1,n_fed
             call broadcast_scalar(pFed_ac    (n,k), pmloc(n))   
            enddo
            do k = 1,n_fep
             call broadcast_scalar(pFep_ac    (n,k), pmloc(n))   
            enddo
            if (skl_bgc) then              ! skl_bgc
             do k = 1,n_algae
               call broadcast_scalar(pN_sk    (n,k), pmloc(n))            
             enddo
             do k = 1,n_doc
               call broadcast_scalar(pDOC_sk    (n,k), pmloc(n))            
             enddo
             do k = 1,n_don
               call broadcast_scalar(pDON_sk    (n,k), pmloc(n))            
             enddo
             do k = 1,n_fed
               call broadcast_scalar(pFed_sk    (n,k), pmloc(n))            
             enddo
             do k = 1,n_fep
               call broadcast_scalar(pFep_sk    (n,k), pmloc(n))            
             enddo
             call broadcast_scalar(pNit_sk  (n), pmloc(n))             
             call broadcast_scalar(pAm_sk   (n), pmloc(n))             
             call broadcast_scalar(pSil_sk  (n), pmloc(n))             
             call broadcast_scalar(pDMSPp_sk(n), pmloc(n))             
             call broadcast_scalar(pDMSPd_sk(n), pmloc(n))             
             call broadcast_scalar(pDMS_sk  (n), pmloc(n))        
            endif   !tr_bgc_sk

           if (z_tracers) then                   !  z_bgc
             do k = 1,n_algae 
                call broadcast_scalar(pN_tot  (n,k), pmloc(n))    
             enddo
             do k = 1,n_zaero 
                call broadcast_scalar(pflux_atm_zaero_s(n,k), pmloc(n)) 
                call broadcast_scalar(pflux_atm_zaero  (n,k), pmloc(n)) 
                call broadcast_scalar(pflux_snow_zaero (n,k), pmloc(n))
             enddo
             call broadcast_scalar(pflux_atm_NO (n), pmloc(n))            
             call broadcast_scalar(pflux_atm_Am(n), pmloc(n))    
             call broadcast_scalar(pflux_snow_NO (n), pmloc(n))             
             call broadcast_scalar(pflux_snow_Am(n), pmloc(n))
             call broadcast_scalar(pgrow_net (n), pmloc(n))
             do k = 1,nblyr+1
              call broadcast_scalar(pzfswin (n,k), pmloc(n))
              call broadcast_scalar(pZoo (n,k), pmloc(n))
              call broadcast_scalar(pNO (n,k), pmloc(n))
              call broadcast_scalar(pAm (n,k), pmloc(n))
              call broadcast_scalar(pPON (n,k), pmloc(n))
              do nn = 1,n_algae
                 call broadcast_scalar(pN (n,k,nn), pmloc(n))
              enddo
              do nn = 1,n_doc
                 call broadcast_scalar(pDOC (n,k,nn), pmloc(n))
              enddo
              do nn = 1,n_don
                 call broadcast_scalar(pDON (n,k,nn), pmloc(n))
              enddo
              do nn = 1,n_fed
                 call broadcast_scalar(pFed (n,k,nn), pmloc(n))
              enddo
              do nn = 1,n_fep
                 call broadcast_scalar(pFep (n,k,nn), pmloc(n))
              enddo
              do nn = 1,n_zaero
                 call broadcast_scalar(pzaero (n,k,nn), pmloc(n))
              enddo
             enddo !k
             do k = 1,2
              call broadcast_scalar(pNOs (n,k), pmloc(n))
              call broadcast_scalar(pAms (n,k), pmloc(n))
              call broadcast_scalar(pPONs(n,k), pmloc(n))
              do nn = 1,n_algae
                 call broadcast_scalar(pNs (n,k,nn), pmloc(n))   
              enddo
              do nn = 1,n_doc
                 call broadcast_scalar(pDOCs (n,k,nn), pmloc(n))   
              enddo
              do nn = 1,n_don
                 call broadcast_scalar(pDONs (n,k,nn), pmloc(n))   
              enddo
              do nn = 1,n_fed
                 call broadcast_scalar(pFeds (n,k,nn), pmloc(n))   
              enddo
              do nn = 1,n_fep
                 call broadcast_scalar(pFeps (n,k,nn), pmloc(n))   
              enddo
              do nn = 1,n_zaero
                 call broadcast_scalar(pzaeros (n,k,nn), pmloc(n))   
              enddo
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
       write(nu_diag,803) 'zfswin(1) PAR  ','zfswin(2) PAR '
       write(nu_diag,*) '---------------------------------------------------'
       write(nu_diag,802) ((pzfswin(n,k),n=1,2), k = 1,nblyr+1)              
       write(nu_diag,*) '      '          
       write(nu_diag,803) 'Losses: Zoo(1)(mmol/m^3)  ','Zoo(2)'
       write(nu_diag,803) '        Brine Conc.       ',' Brine Conc'
       write(nu_diag,*) '---------------------------------------------------'
       write(nu_diag,802) ((pZoo(n,k),n=1,2), k = 1,nblyr+1)              
       write(nu_diag,*) '      '          
       if (tr_bgc_Nit) then
         write(nu_diag,*) '    nitrate conc. (mmol/m^3) or flux (mmol/m^2/d)'
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,900) 'Ocean conc       = ',pNit_ac(1),pNit_ac(2)
         write(nu_diag,900) 'ice-ocean flux   = ',pflux_NO(1),pflux_NO(2)
         if (skl_bgc) then
           write(nu_diag,900) 'Bulk ice conc.   = ',pNit_sk(1),pNit_sk(2)
         elseif (z_tracers) then
           write(nu_diag,900) 'atm-ice flux     = ',pflux_atm_NO(1),pflux_atm_NO(2)
           write(nu_diag,900) 'snow-ice flux    = ',pflux_snow_NO(1),pflux_snow_NO(2)
           write(nu_diag,*) '             snow + ice conc'
           write(nu_diag,803) '    nitrate(1)','   nitrate(2)'
           write(nu_diag,*) '---------------------------------------------------'
           write(nu_diag,802) ((pNOs(n,k),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pNO(n,k),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '    '        
         endif
       endif
       if (tr_bgc_PON .and. z_tracers) then
           write(nu_diag,*) '    PON snow + ice conc. (mmol/m^3)'
           write(nu_diag,803) '    PON(1)','    PON(2)'
           write(nu_diag,*) '---------------------------------------------------'
           write(nu_diag,802) ((pPONs(n,k),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pPON(n,k),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) ' '
       endif
       if (tr_bgc_Am) then
         write(nu_diag,*) '    ammonium conc. (mmol/m^3) or flux (mmol/m^2/d)'
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,900) 'Ocean conc       = ',pAm_ac(1),pAm_ac(2)
         write(nu_diag,900) 'ice-ocean flux   = ',pflux_Am(1),pflux_Am(2)
         if (skl_bgc) then
           write(nu_diag,900) 'Bulk ice conc.   = ',pAm_sk(1),pAm_sk(2)
         elseif (z_tracers) then
           write(nu_diag,900) 'atm-ice flux     = ',pflux_atm_Am(1),pflux_atm_Am(2)
           write(nu_diag,900) 'snow-ice flux    = ',pflux_snow_Am(1),pflux_snow_Am(2)
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  ammonium(1)','  ammonium (2)'
           write(nu_diag,*) '---------------------------------------------------'
           write(nu_diag,802) ((pAms(n,k),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pAm(n,k),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '       '     
         endif
       endif
       if (tr_bgc_N) then
         write(nu_diag,900) 'tot algal growth (1/d) = ',pgrow_net(1),pgrow_net(2)
       do kk = 1,n_algae
         write(nu_diag,*) '  algal conc. (mmol N/m^3) or flux (mmol N/m^2/d)'
         write(nu_diag,1020) '  type:', kk
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,900) 'Ocean conc           = ',pN_ac(1,kk),pN_ac(2,kk)
         write(nu_diag,900) 'ice-ocean flux       = ',pflux_N(1,kk),pflux_N(2,kk)
         if (skl_bgc) then
           write(nu_diag,900) 'Bulk ice conc.   = ',pN_sk(1,kk),pN_sk(2,kk)
         elseif (z_tracers) then
           write(nu_diag,900) 'Tot ice (mmolN/m^2) = ',pN_tot(1,kk),pN_tot(2,kk)
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  algal N(1)','  algal N(2) '
           write(nu_diag,*) '---------------------------------------------------'
           write(nu_diag,802) ((pNs(n,k,kk),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pN(n,k,kk),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '         '   
         endif
       enddo
       endif
       if (tr_bgc_C) then
       do kk = 1,1 !n_doc
         write(nu_diag,*) '  DOC conc. (mmol C/m^3)'
         write(nu_diag,1020) '  type:', kk
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,900)  'Ocean conc       = ',(pDOC_ac(n,kk),n=1,2)
         if (skl_bgc) then
           write(nu_diag,900)'Bulk ice conc.   = ',pDOC_sk(1,kk),pDOC_sk(2,kk)
         elseif (z_tracers) then
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  DOC(1)','  DOC(2) '
           write(nu_diag,*) '---------------------------------------------------'
           write(nu_diag,802) ((pDOCs(n,k,kk),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pDOC(n,k,kk),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '      '      
         endif
       enddo
       endif
       if (tr_bgc_DON) then
       do kk = 1,n_don
         write(nu_diag,*) '  DON conc. (mmol N/m^3)'
         write(nu_diag,1020) '  type:', kk
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,900)  'Ocean conc       = ',(pDON_ac(n,kk),n=1,2)
         if (skl_bgc) then
           write(nu_diag,900)'Bulk ice conc.   = ',pDON_sk(1,kk),pDON_sk(2,kk)
         elseif (z_tracers) then
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  DON(1)','  DON(2) '
           write(nu_diag,*) '---------------------------------------------------'
           write(nu_diag,802) ((pDONs(n,k,kk),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pDON(n,k,kk),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '      '      
         endif
       enddo
       endif
       if (tr_bgc_Fe ) then
       do kk = 1,n_fed
         write(nu_diag,*) ' dFe  conc. (nM)'
         write(nu_diag,1020) '  type:', kk
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,900)  'Ocean conc       = ',(pFed_ac (n,kk),n=1,2)
         if (skl_bgc) then
           write(nu_diag,900)'Bulk ice conc.   = ',pFed_sk (1,kk),pFed_sk (2,kk)
         elseif (z_tracers) then
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  Fed (1)','  Fed (2) '
           write(nu_diag,*) '---------------------------------------------------'
           write(nu_diag,802) ((pFeds (n,k,kk),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pFed (n,k,kk),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '      '      
         endif
       enddo
       do kk = 1,n_fep
         write(nu_diag,*) ' pFe  conc. (nM)'
         write(nu_diag,1020) '  type:', kk
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,900)  'Ocean conc       = ',(pFep_ac (n,kk),n=1,2)
         if (skl_bgc) then
           write(nu_diag,900)'Bulk ice conc.   = ',pFep_sk (1,kk),pFep_sk (2,kk)
         elseif (z_tracers) then
           write(nu_diag,*) '             snow + ice conc.'
           write(nu_diag,803) '  Fep (1)','  Fep (2) '
           write(nu_diag,*) '---------------------------------------------------'
           write(nu_diag,802) ((pFeps (n,k,kk),n=1,2), k = 1,2)              
           write(nu_diag,802) ((pFep (n,k,kk),n=1,2), k = 1,nblyr+1)              
           write(nu_diag,*) '      '      
         endif
       enddo
       endif
       if (tr_bgc_DMS) then
         write(nu_diag,*) '    DMS (mmol/m^3)      '
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,900) 'Ocean DMSP    = ',pDMSP_ac(1),pDMSP_ac(2)
         write(nu_diag,900) 'Ocean DMS     = ',pDMS_ac(1),pDMS_ac(2)
         if (skl_bgc) then
          write(nu_diag,900) 'Ice DMSPp    = ',pDMSPp_sk(1),pDMSPp_sk(2)
          write(nu_diag,900) 'Ice DMSPd    = ',pDMSPd_sk(1),pDMSPd_sk(2)
          write(nu_diag,900) 'Ice DMS      = ',pDMS_sk(1),pDMS_sk(2)    
         endif
       endif
       if (tr_zaero .and. z_tracers) then
       do kk = 1,n_zaero
         write(nu_diag,*) '  aerosol conc. (kg/m^3) or flux (kg/m^2/d)'
         write(nu_diag,1020) '  type: ',kk
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,900) 'Atm source flux     = ',pflux_atm_zaero_s(1,kk),pflux_atm_zaero_s(2,kk)
         write(nu_diag,900) 'ice-ocean flux      = ',pflux_zaero(1,kk),pflux_zaero(2,kk)
         write(nu_diag,900) 'direct atm-ice flux = ',pflux_atm_zaero(1,kk),pflux_atm_zaero(2,kk)
         write(nu_diag,900) 'snow-ice flux       = ',pflux_snow_zaero(1,kk),pflux_snow_zaero(2,kk)
         write(nu_diag,*) '             snow + ice conc.'
         write(nu_diag,803) ' aerosol(1)','    aerosol(2) '
         write(nu_diag,*) '---------------------------------------------------'
         write(nu_diag,802) ((pzaeros(n,k,kk),n=1,2), k = 1,2)              
         write(nu_diag,802) ((pzaero(n,k,kk),n=1,2), k = 1,nblyr+1)   
         write(nu_diag,*) '            '
      enddo
      endif

      endif                   ! print_points
      endif                   ! my_task = master_task 

  802 format (f24.17,2x,f24.17)
  803 format (a25,2x,a25)
  900 format (a25,2x,f24.17,2x,f24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
 1020 format (a30,2x,i6)    ! integer

      end subroutine bgc_diags

!=======================================================================
!
! Find ice-ocean flux when ice is thin and internal dynamics/reactions are
!             assumed to be zero
!
! authors     Nicole Jeffery, LANL

      subroutine thin_ice_flux (hin, hin_old, phin, Cin, flux_o_tot, &
                                zbgc_snow, zbgc_atm, igrid,dt) 

      use ice_constants, only: c1, p5

      real (kind=dbl_kind), dimension(nblyr+1), intent(in) :: &
         phin

      real (kind=dbl_kind), dimension(nblyr+1), intent(inout) :: &
         Cin              ! initial concentration*hin_old*phin

      real (kind=dbl_kind), intent(in) :: &
         hin_old   , &     ! brine thickness (m) 
         hin       , &     ! new brine thickness (m)
         dt        , &     ! time step
         zbgc_atm  , &     ! atmospheric flux (mmol/m^3*h)
         zbgc_snow         ! snow flux (mmol/m^3*h)

      real (kind=dbl_kind), intent(inout) :: &
         flux_o_tot        ! tracer flux, gravity+molecular drainage flux ,
                           ! and boundary flux to ocean (mmol/m^2/s)  
                           ! positive into the ocean  

      real (kind=dbl_kind), dimension (nblyr + 1), intent(in) :: &
         igrid            ! biology nondimensional grid interface points 

     ! local variables

     integer (kind=int_kind) :: &
         k         ! vertical biology layer index
   
     real (kind=dbl_kind) :: &
         sum_bio, zspace  
   
         zspace = c1/real(nblyr,kind=dbl_kind)
         sum_bio = (Cin(1)+Cin(nblyr+1))/hin_old*zspace*p5
         do k = 2, nblyr
          sum_bio = sum_bio + Cin(k)/hin_old*zspace
          Cin(k) = Cin(k) *hin/hin_old
         enddo
         Cin(1) = Cin(1) * hin/hin_old
         Cin(nblyr+1) = Cin(nblyr+1)* hin/hin_old
         flux_o_tot = flux_o_tot -(hin-hin_old)*sum_bio/dt + zbgc_atm/dt + zbgc_snow/dt

     end subroutine thin_ice_flux

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
! Dumps all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_bgc()

      use ice_domain_size, only: ncat, n_algae, n_doc, n_dic, &
          n_don, n_zaero, n_fed, n_fep
      use ice_state, only: trcrn, nt_bgc_Am, &
          nt_bgc_DMS, nt_bgc_DMSPd, nt_bgc_C, nt_bgc_chl, &
          nt_bgc_DMSPp, nt_bgc_Nit, nt_bgc_Sil, &
          nt_bgc_PON, nt_bgc_DON, nt_bgc_DOC, nt_bgc_DIC, &
          nt_bgc_N, nt_zaero, nt_bgc_Fed, nt_bgc_Fep, &
          nt_zbgc_frac, nbtrcr
      use ice_restart,only: write_restart_field
      use ice_constants, only: field_loc_center, field_type_scalar

      ! local variables

      integer (kind=int_kind) :: &
         k, mm  ! vertical index or n_algae

      logical (kind=log_kind) :: diag

      character (len=3) :: nchar, ncharb

      diag = .true.

      !-----------------------------------------------------------------
      ! Skeletal layer BGC
      !-----------------------------------------------------------------
      if (skl_bgc) then
         do k = 1, n_algae
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_N(k),:,:), &
                                  'ruf8','bgc_N'//trim(nchar),ncat,diag)
            if (tr_bgc_chl) &
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_chl(k),:,:), &
                                  'ruf8','bgc_chl'//trim(nchar),ncat,diag)
         enddo
        if (tr_bgc_C)  then
          ! do k = 1, n_algae
          !  write(nchar,'(i3.3)') k
          !  call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_C(k),:,:), &
          !                        'ruf8','bgc_C'//trim(nchar),ncat,diag)
          ! enddo
           do k = 1, n_doc
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DOC(k),:,:), &
                                  'ruf8','bgc_DOC'//trim(nchar),ncat,diag)
           enddo
           do k = 1, n_dic
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DIC(k),:,:), &
                                  'ruf8','bgc_DIC'//trim(nchar),ncat,diag)
           enddo
         endif
         if (tr_bgc_Nit) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Nit,:,:), &
                                  'ruf8','bgc_Nit',ncat,diag)
         if (tr_bgc_Am) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Am,:,:), &
                                  'ruf8','bgc_Am',ncat,diag)
         if (tr_bgc_Sil) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Sil,:,:), &
                                  'ruf8','bgc_Sil',ncat,diag)
         if (tr_bgc_DMS) then
           call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPp,:,:), &
                                  'ruf8','bgc_DMSPp',ncat,diag)
           call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPd,:,:), &
                                  'ruf8','bgc_DMSPd',ncat,diag)
           call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMS,:,:), &
                                  'ruf8','bgc_DMS',ncat,diag)
         endif
         if (tr_bgc_PON) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_PON,:,:), &
                                  'ruf8','bgc_PON',ncat,diag)
      
        if (tr_bgc_DON)  then
           do k = 1, n_don
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DON(k),:,:), &
                                  'ruf8','bgc_DON'//trim(nchar),ncat,diag)
           enddo
         endif
        if (tr_bgc_Fe )  then
           do k = 1, n_fed 
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Fed (k),:,:), &
                                  'ruf8','bgc_Fed'//trim(nchar),ncat,diag)
           enddo
           do k = 1, n_fep 
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Fep (k),:,:), &
                                  'ruf8','bgc_Fep'//trim(nchar),ncat,diag)
           enddo
         endif
      !-----------------------------------------------------------------
      ! Z layer BGC
      !-----------------------------------------------------------------
      else 
         if (tr_bgc_Nit) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Nit+k-1,:,:),'ruf8', &
                                 'bgc_Nit'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_N) then
         do mm = 1,n_algae
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0, &
                                  trcrn(:,:,nt_bgc_N(mm)+k-1,:,:),'ruf8', &
                                 'bgc_N'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         if (tr_bgc_chl) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_chl(mm)+k-1,:,:),'ruf8', &
                                 'bgc_chl'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         endif
         enddo   !n_algae
         endif   ! tr_bgc_N
         if (tr_bgc_C) then
        ! do mm = 1,n_algae
        ! write(ncharb, '(i3.3)') mm
        ! do k = 1,nblyr+3
        !    write(nchar,'(i3.3)') k
        !    call write_restart_field(nu_dump_bgc,0,  &
        !                          trcrn(:,:,nt_bgc_C(mm)+k-1,:,:),'ruf8', &
        !                         'bgc_C'//trim(ncharb)//trim(nchar),ncat,diag)
        ! enddo
        ! enddo
         do mm = 1,n_doc
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_DOC(mm)+k-1,:,:),'ruf8', &
                                 'bgc_DOC'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         do mm = 1,n_dic
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_DIC(mm)+k-1,:,:),'ruf8', &
                                 'bgc_DIC'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif  !tr_bgc_C
         if (tr_bgc_Am) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Am+k-1,:,:),'ruf8', &
                                 'bgc_Am'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_Sil) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Sil+k-1,:,:),'ruf8', &
                                 'bgc_Sil'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_DMS) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPp+k-1,:,:),'ruf8', &
                                 'bgc_DMSPp'//trim(nchar),ncat,diag)
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPd+k-1,:,:),'ruf8', &
                                 'bgc_DMSPd'//trim(nchar),ncat,diag)
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMS+k-1,:,:),'ruf8', &
                                 'bgc_DMS'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_PON) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_PON+k-1,:,:),'ruf8', &
                                 'bgc_PON'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_DON) then
         do mm = 1,n_don
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_DON(mm)+k-1,:,:),'ruf8', &
                                 'bgc_DON'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif
         if (tr_bgc_Fe ) then
         do mm = 1,n_fed
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_Fed(mm)+k-1,:,:),'ruf8', &
                                 'bgc_Fed'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         do mm = 1,n_fep
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_Fep(mm)+k-1,:,:),'ruf8', &
                                 'bgc_Fep'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif
         if (tr_zaero) then
         do mm = 1,n_zaero
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_zaero(mm)+k-1,:,:),'ruf8', &
                                 'zaero'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif
         do mm = 1,nbtrcr
          write(nchar,'(i3.3)') mm
          call write_restart_field(nu_dump_bgc,0,   &
                                trcrn(:,:,nt_zbgc_frac+mm-1,:,:),'ruf8', &
                                'zbgc_frac'//trim(nchar),ncat,diag)
         enddo
      endif
      !-----------------------------------------------------------------
      ! Ocean BGC
      !-----------------------------------------------------------------
      if (tr_bgc_N) then
      do k = 1,n_algae
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,algalN(:,:,k,:),'ruf8','algalN'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_C) then
      do k = 1,n_doc
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,doc(:,:,k,:),'ruf8','doc'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_dic
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,dic(:,:,k,:),'ruf8','dic'//trim(nchar),1,diag)
      enddo  !k
      endif      
      if (tr_bgc_Nit) &
      call write_restart_field(nu_dump_bgc,0,nit,   'ruf8','nit',   1,diag)
      if (tr_bgc_Am) &
      call write_restart_field(nu_dump_bgc,0,amm,   'ruf8','amm',   1,diag)
      if (tr_bgc_Sil) &
      call write_restart_field(nu_dump_bgc,0,sil,   'ruf8','sil',   1,diag)
      if (tr_bgc_DMS) then
        call write_restart_field(nu_dump_bgc,0,dmsp,  'ruf8','dmsp',  1,diag)
        call write_restart_field(nu_dump_bgc,0,dms,   'ruf8','dms',   1,diag)
      endif
      if (tr_bgc_DON) then
      do k = 1,n_don
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,don(:,:,k,:),'ruf8','don'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_Fe ) then
      do k = 1,n_fed
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,fed(:,:,k,:),'ruf8','fed'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_fep
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,fep(:,:,k,:),'ruf8','fep'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_zaero) then
      do k = 1,n_zaero
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,zaeros(:,:,k,:),'ruf8','zaeros'//trim(nchar),1,diag)
      enddo  !k
      endif

      end subroutine write_restart_bgc

!=======================================================================
!
! Reads all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_bgc()

      use ice_constants, only: field_loc_center, field_type_scalar
      use ice_domain_size, only: ncat, n_algae, n_doc, n_dic,&
          n_don, n_zaero, n_fed, n_fep
      use ice_state, only: trcrn,nt_bgc_Am, &
          nt_bgc_DMS, nt_bgc_DMSPd, nt_bgc_C, &
          nt_bgc_DMSPp, nt_bgc_Nit, nt_bgc_Sil, &
          nt_bgc_PON, nt_bgc_DON, nt_bgc_DOC, nt_bgc_DIC, &
          nt_bgc_N, nt_zaero, nt_bgc_Fed, nt_bgc_Fep,  nt_bgc_chl, &
          nt_zbgc_frac, nbtrcr
      use ice_restart,only: read_restart_field

      ! local variables

      integer (kind=int_kind) :: &
          k, mm    ! counting indices
      logical (kind=log_kind) :: diag

      character (len=3) :: nchar, ncharb

      diag = .true.

      !-----------------------------------------------------------------
      ! Skeletal Layer BGC
      !-----------------------------------------------------------------
      if (skl_bgc) then
       if (my_task == master_task) write(nu_diag,*) 'skl bgc restart'

       do k = 1, n_algae
          write(nchar,'(i3.3)') k
          call read_restart_field(nu_restart_bgc,0, &
                 trcrn(:,:,nt_bgc_N(k),:,:), &
                 'ruf8','bgc_N'//trim(nchar),ncat,diag)
          if (tr_bgc_chl) &
          call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_chl(k),:,:), &
                 'ruf8','bgc_chl'//trim(nchar),ncat,diag)
       enddo   !k
       if (tr_bgc_C) then
         ! do k = 1, n_algae
         !     write(nchar,'(i3.3)') k
         !     call read_restart_field(nu_restart_bgc,0,  &
         !        trcrn(:,:,nt_bgc_C(k),:,:), &
         !        'ruf8','bgc_C'//trim(nchar),ncat,diag)
         ! enddo
          do k = 1, n_doc
              write(nchar,'(i3.3)') k
              call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_DOC(k),:,:), &
                 'ruf8','bgc_DOC'//trim(nchar),ncat,diag)
          enddo
          do k = 1, n_dic
              write(nchar,'(i3.3)') k
              call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_DIC(k),:,:), &
                 'ruf8','bgc_DIC'//trim(nchar),ncat,diag)
          enddo
       endif
       if (tr_bgc_Nit) &
       call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Nit,:,:), &
           'ruf8','bgc_Nit',ncat,diag)
       if (tr_bgc_Am) &
       call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Am,:,:), &
           'ruf8','bgc_Am',ncat,diag)
       if (tr_bgc_Sil) &
       call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Sil,:,:), &
           'ruf8','bgc_Sil',ncat,diag)
       if(tr_bgc_DMS) then
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPp,:,:), &
           'ruf8','bgc_DMSPp',ncat,diag)
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPd,:,:), &
           'ruf8','bgc_DMSPd',ncat,diag)
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMS,:,:), &
           'ruf8','bgc_DMS',ncat,diag)
       endif
       if (tr_bgc_PON) &
       call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_PON,:,:), &
           'ruf8','bgc_PON',ncat,diag)
       if (tr_bgc_DON) then
          do k = 1, n_don
              write(nchar,'(i3.3)') k
              call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_DON(k),:,:), &
                 'ruf8','bgc_DON'//trim(nchar),ncat,diag)
          enddo
       endif
       if (tr_bgc_Fe) then
          do k = 1, n_fed 
              write(nchar,'(i3.3)') k
              call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_Fed (k),:,:), &
                 'ruf8','bgc_Fed'//trim(nchar),ncat,diag)
          enddo
          do k = 1, n_fep 
              write(nchar,'(i3.3)') k
              call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_Fep (k),:,:), &
                 'ruf8','bgc_Fep'//trim(nchar),ncat,diag)
          enddo
       endif

      else
      !-----------------------------------------------------------------
      ! Z Layer BGC
      !-----------------------------------------------------------------
      if (tr_bgc_Nit) then
      if (my_task == master_task) write(nu_diag,*) 'z bgc restart: min/max Nitrate'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Nit+k-1,:,:),'ruf8', &
              'bgc_Nit'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif   !Nit
      if (tr_bgc_N) then
      do mm = 1,n_algae
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) ' min/max Algal N'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0, &
               trcrn(:,:,nt_bgc_N(mm)+k-1,:,:),'ruf8', &
              'bgc_N'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      if (tr_bgc_chl) then
      if (my_task == master_task) write(nu_diag,*) ' min/max Algal chla'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_chl(mm)+k-1,:,:),'ruf8', &
              'bgc_chl'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif   ! tr_bgc_chl
      enddo   !n_algae
      endif   ! tr_bgc_N
      if (tr_bgc_C) then
     ! do mm = 1,n_algae
     ! write(ncharb,'(i3.3)') mm
     ! if (my_task == master_task) write(nu_diag,*) ' min/max Algal C'
     ! do k=1,nblyr+3
     !    write(nchar,'(i3.3)') k
     !    call read_restart_field(nu_restart_bgc,0,  &
     !          trcrn(:,:,nt_bgc_C(mm)+k-1,:,:),'ruf8', &
     !         'bgc_C'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
     ! enddo
     ! enddo  !mm
      do mm = 1,n_doc
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) ' min/max DOC'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_DOC(mm)+k-1,:,:),'ruf8', &
              'bgc_DOC'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      do mm = 1,n_dic
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) ' min/max DIC'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_DIC(mm)+k-1,:,:),'ruf8', &
              'bgc_DIC'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif  ! tr_bgc_C
      if (tr_bgc_Am) then
      if (my_task == master_task) write(nu_diag,*) ' min/max ammonium'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Am+k-1,:,:),'ruf8', &
              'bgc_Am'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_Sil) then
      if (my_task == master_task) write(nu_diag,*) ' min/max silicate'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Sil+k-1,:,:),'ruf8', &
              'bgc_Sil'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_DMS) then
      do k=1,nblyr+3
      if (my_task == master_task) write(nu_diag,*) ' min/max DMSPp'
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPp+k-1,:,:),'ruf8', &
              'bgc_DMSPp'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      if (my_task == master_task) write(nu_diag,*) ' min/max DMSPd'
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPd+k-1,:,:),'ruf8', &
              'bgc_DMSPd'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      if (my_task == master_task) write(nu_diag,*) ' min/max DMS'
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMS+k-1,:,:),'ruf8', &
              'bgc_DMS'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_PON) then
      if (my_task == master_task) write(nu_diag,*) ' min/max PON'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_PON+k-1,:,:),'ruf8', &
              'bgc_PON'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_DON) then
      do mm = 1,n_don
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) ' min/max DON'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_DON(mm)+k-1,:,:),'ruf8', &
              'bgc_DON'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif
      if (tr_bgc_Fe) then
      do mm = 1,n_fed
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) ' min/max dFe '
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_Fed (mm)+k-1,:,:),'ruf8', &
              'bgc_Fed'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      do mm = 1,n_fep
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) ' min/max pFe '
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_Fep (mm)+k-1,:,:),'ruf8', &
              'bgc_Fep'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif
      if (tr_zaero) then
      do mm = 1,n_zaero
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) ' min/max z aerosols'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_zaero(mm)+k-1,:,:),'ruf8', &
              'zaero'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif
      do mm = 1,nbtrcr
          write(nchar,'(i3.3)') mm
          call read_restart_field(nu_restart_bgc,0, &
               trcrn(:,:,nt_zbgc_frac+mm-1,:,:),'ruf8', &
               'zbgc_frac'//trim(nchar),ncat,diag)
      enddo
      endif

      !-----------------------------------------------------------------
      ! Ocean BGC
      !-----------------------------------------------------------------

      if (my_task == master_task) write(nu_diag,*) 'mixed layer ocean bgc restart'
      if (tr_bgc_N) then
      do k = 1,n_algae
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,algalN(:,:,k,:),'ruf8','algalN'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_C) then
      do k = 1,n_doc
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,doc(:,:,k,:),'ruf8','doc'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_dic
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,dic(:,:,k,:),'ruf8','dic'//trim(nchar),1,diag)
      enddo  !k
      endif  !tr_bgc_C

      if (tr_bgc_Nit) &
      call read_restart_field(nu_restart_bgc,0,nit   ,'ruf8','nit'   ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_Am) &
      call read_restart_field(nu_restart_bgc,0,amm   ,'ruf8','amm'   ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_Sil) &
      call read_restart_field(nu_restart_bgc,0,sil   ,'ruf8','sil'   ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_DMS) then
        call read_restart_field(nu_restart_bgc,0,dmsp  ,'ruf8','dmsp'  ,1,diag)
        call read_restart_field(nu_restart_bgc,0,dms   ,'ruf8','dms'   ,1,diag)
      endif
      if (tr_bgc_DON) then
      do k = 1,n_don
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,don(:,:,k,:),'ruf8','don'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_Fe ) then
      do k = 1,n_fed
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,fed(:,:,k,:),'ruf8','fed'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_fep
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,fep(:,:,k,:),'ruf8','fep'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_zaero) then
      do k = 1,n_zaero
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,zaeros(:,:,k,:),'ruf8','zaeros'//trim(nchar),1,diag)
      enddo  !k
      endif

      end subroutine read_restart_bgc

!=======================================================================
!
! author: Nicole Jeffery, LANL

      subroutine get_atm_bgc (dt)

      use ice_blocks, only: nx_block, ny_block, block, get_block
      use ice_constants, only: p01, c0
      use ice_domain, only: nblocks, distrb_info, blocks_ice, nblocks
      use ice_state, only: nbtrcr
    !  use ice_zbgc_shared, only: flux_bio_atm, atm_bio, ocean_bio 

          
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step  

      !  local variables

      integer (kind=int_kind) :: &
         i, j, nn       , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         iblk               ! block index 

      type (block) :: &
         this_block      ! block information for current block

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
             
      do nn = 1,nbtrcr
           do j = jlo, jhi
           do i = ilo, ihi 
                flux_bio_atm(i,j,nn,iblk) = atm_bio_all(i,j,bio_index_o(nn),iblk)/dt*grid_o_t
                atm_bio(i,j,nn,iblk)      = atm_bio_all(i,j,bio_index_o(nn),iblk)
           enddo
           enddo
      enddo

      enddo ! iblk
      !$OMP END PARALLEL DO

      end subroutine get_atm_bgc

!=======================================================================

!  Increase aerosol in  snow surface due to deposition
!  and vertical cycling : after update_aerosol

      subroutine update_snow_bgc (dt, meltt, melts,    &
                                meltb,    congel,      &
                                snoice,   nbtrcr,      &
                                fsnow,    ntrcr,       &
                                trcrn,    bio_index,   &
                                aice_old, zbgc_snow,   &
                                vice_old, vsno_old,    &
                                vicen, vsnon, aicen,   &
                                flux_bio_atm, zbgc_atm)

      use ice_domain_size, only: nblyr, nslyr
      use ice_shortwave, only: hi_ssl, hs_ssl
      use ice_constants, only: c0, rhos, rhoi, hs_min, puny, &
                         c2, c1
      use ice_zbgc_shared, only: kscavz

      integer (kind=int_kind), intent(in) :: &
         nbtrcr,             & ! number of distinct snow tracers
         ntrcr                 ! number of tracers

      integer (kind=int_kind), dimension (nbtrcr), intent(in) :: &
         bio_index       

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), intent(in) :: &
         meltt,    & ! thermodynamic melt/growth rates
         melts,    &
         meltb,    &
         congel,   &
         snoice,   &
         fsnow,    &
         vicen,    & ! ice volume (m)
         vsnon,    & ! snow volume (m)
         aicen,    & ! ice area fraction
         aice_old, & ! values prior to thermodynamic changes
         vice_old, &
         vsno_old 

      real (kind=dbl_kind),dimension(nbtrcr), intent(inout) :: &
         zbgc_snow, & ! aerosol contribution from snow to ice
         zbgc_atm     ! and atm to ice concentration * volume (mmol/m^3*m)

      real (kind=dbl_kind), dimension(nbtrcr), &
         intent(in) :: &
         flux_bio_atm   ! aerosol deposition rate (mmol/m^2 s)

      real (kind=dbl_kind), dimension(ntrcr), &
         intent(inout) :: &
         trcrn       ! ice/snow tracer array

      !  local variables

      integer (kind=int_kind) :: i, j, ij, k, n

      real (kind=dbl_kind) :: &
         dzssl,  dzssl_new,      & ! snow ssl thickness
         dzint,  dzint_new,      & ! snow interior thickness
         hs,                     & ! snow thickness (m)
         dhs_evap,               & ! snow thickness change due to evap
         dhs_melts,              & ! ... due to surface melt
         dhs_snoice,             & ! ... due to snow-ice formation
         hslyr,                  & ! snow layer thickness (m)
         hslyr_old,              & ! old snow layer thickness (m)
         hs_old,                 & ! old snow thickness (m)
         dznew,                  & ! change in the snow sl (m)
         sloss1, sloss2,         & ! aerosol mass loss (kg/m^2)
         ar                        ! 1/aicen(i,j)

      real (kind=dbl_kind), dimension(nbtrcr) :: &
       !  kscav, kscavsi   , & ! scavenging by melt water
         aerotot, aerotot0, & ! for conservation check (mmol/m^3)
         aero_cons          ! for conservation check (mmol/m^2)


      real (kind=dbl_kind), dimension(nbtrcr,2) :: &
         aerosno,  & ! kg/m^2
         aerosno0  ! for diagnostic prints
     
     ! use in snow only: Flanner et al, 2007 for aerosols
     ! data kscav   / .03_dbl_kind, .20_dbl_kind,&
     !      .02_dbl_kind,.02_dbl_kind,.01_dbl_kind,.01_dbl_kind /
     ! data kscavsi / .03_dbl_kind, .20_dbl_kind,&
     !      .02_dbl_kind,.02_dbl_kind,.01_dbl_kind,.01_dbl_kind /

    !-------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------
         aerosno (:,:) = c0
         aerosno0(:,:) = c0
         aero_cons(:) = c0

         hs_old    = vsno_old/aice_old
         hslyr_old = hs_old/real(nslyr,kind=dbl_kind)

         dzssl  = min(hslyr_old/c2, hs_ssl)
         dzint  = hs_old - dzssl

         if (aicen > c0) then
            ar = c1/aicen
            hs = vsnon*ar
            dhs_melts  = -melts*ar
            dhs_snoice = snoice*ar*rhoi/rhos
         else ! ice disappeared during time step
            hs = vsnon/aice_old
            dhs_melts  = -melts/aice_old
            dhs_snoice = snoice/aice_old*rhoi/rhos
         endif

         dhs_evap = hs - (hs_old + dhs_melts - dhs_snoice &
                                 + fsnow/rhos*dt)

         ! trcrn() has units kg/m^3
         do k=1,nbtrcr
            aerosno (k,:) = &
               trcrn(bio_index(k)+ nblyr+1:bio_index(k)+ nblyr+2)*vsno_old
            aerosno0(k,:) = aerosno(k,:)
            aerotot0(k) = aerosno(k,2) + aerosno(k,1) 
         enddo

    !-------------------------------------------------------------------
    ! evaporation
    !-------------------------------------------------------------------
         dzint  = dzint  + min(dzssl  + dhs_evap, c0)
         dzssl  = max(dzssl  + dhs_evap, c0)
    !-------------------------------------------------------------------
    ! surface snow melt
    !-------------------------------------------------------------------
         if (-dhs_melts > puny) then
            do k = 1, nbtrcr
               sloss1 = c0
               sloss2 = c0
               if (dzssl > puny)  &
                  sloss1 = kscavz(k)*aerosno(k,1)  &
                                   *min(-dhs_melts,dzssl)/dzssl
               aerosno(k,1) = aerosno(k,1) - sloss1
               if (dzint > puny)  &
                   sloss2 = kscavz(k)*aerosno(k,2) &
                                    *max(-dhs_melts-dzssl,c0)/dzint
               aerosno(k,2) = aerosno(k,2) - sloss2
               zbgc_snow(k) = zbgc_snow(k) + (sloss1+sloss2)
            enddo  ! 

            ! update snow thickness
            dzint=dzint+min(dzssl+dhs_melts, c0)
            dzssl=max(dzssl+dhs_melts, c0)

            if ( dzssl <= puny ) then ! ssl melts away
               aerosno(:,2) = aerosno(:,1) + aerosno(:,2)
               aerosno(:,1) = c0
               dzssl = max(dzssl, c0)
            endif
            if (dzint <= puny ) then  ! all snow melts away
               zbgc_snow(:) = zbgc_snow(:) &
                                + aerosno(:,1) + aerosno(:,2)
               aerosno(:,:) = c0
               dzint = max(dzint, c0)
            endif
         endif

    !-------------------------------------------------------------------
    ! snowfall
    !-------------------------------------------------------------------
         if (fsnow > c0) dzssl = dzssl + fsnow/rhos*dt

    !-------------------------------------------------------------------
    ! snow-ice formation
    !-------------------------------------------------------------------
         if (dhs_snoice > puny) then
            do k = 1, nbtrcr
               sloss1 = c0
               sloss2 = c0
               if (dzint > puny)  &
                  sloss2 = min(dhs_snoice, dzint)  &
                           *aerosno(k,2)/dzint
               aerosno(k,2) = aerosno(k,2) - sloss2
               if (dzssl > puny)  &
                  sloss1 = max(dhs_snoice-dzint, c0)  &
                           *aerosno(k,1)/dzssl
               aerosno(k,1) = aerosno(k,1) - sloss1
               zbgc_snow(k)= zbgc_snow(k) &
                               + (sloss2+sloss1)
              !                  + kscavsi(k)*(sloss2+sloss1)
            enddo
            dzssl  = dzssl - max(dhs_snoice-dzint, c0)
            dzint  = max(dzint-dhs_snoice, c0)
         endif

    !-------------------------------------------------------------------
    ! aerosol deposition
    !-------------------------------------------------------------------
         if (aicen > c0) then
            hs = vsnon * ar
         else
            hs = c0
         endif
         if (hs >= hs_min)  then !should this really be hs_min or 0? 
                                  ! should use same hs_min value as in radiation
            do k=1,nbtrcr
               aerosno(k,1) = aerosno(k,1) &
                                + flux_bio_atm(k)*dt*aicen
            enddo
         else  
            do k=1,nbtrcr
               aerosno(k,1) = aerosno(k,1) &
                                + (hs/hs_min)*flux_bio_atm(k)*dt*aicen
               zbgc_atm(k) = zbgc_atm(k) &
                                + (c1-hs/hs_min)*flux_bio_atm(k)*dt*aicen
            enddo
         endif

    !-------------------------------------------------------------------
    ! redistribute aerosol within vertical layers
    !-------------------------------------------------------------------
         if (dzssl <= puny) then   ! nothing in SSL
            do k=1,nbtrcr
               aerosno(k,2) = aerosno(k,2) + aerosno(k,1)
               aerosno(k,1) = c0
            enddo
         endif
         if (dzint <= puny) then   ! nothing in Snow Int
            do k = 1, nbtrcr
               zbgc_snow(k) = zbgc_snow(k) + aerosno(k,2)
               aerosno(k,2) = c0
            enddo
         endif

         hslyr      = hs/real(nslyr,kind=dbl_kind)
         dzssl_new  = min(hslyr/c2, hs_ssl)
         dzint_new  = hs - dzssl_new

         if (hs > puny) then !should this really be hs_min or 0? 
            do k = 1, nbtrcr
               dznew = min(dzssl_new-dzssl, c0)
               sloss1 = c0
               if (dzssl > puny) &
                  sloss1 = dznew*aerosno(k,1)/dzssl ! not neccesarily a loss
                  dznew = max(dzssl_new-dzssl, c0)
               if (dzint > puny) &
                  sloss1 = sloss1 + aerosno(k,2)*dznew/dzint
               aerosno(k,1) = aerosno(k,1) + sloss1 
               aerosno(k,2) = aerosno(k,2) - sloss1
            enddo
         else
            zbgc_snow(:) = zbgc_snow(:)  &
                             + aerosno(:,1) + aerosno(:,2)
            aerosno(:,:) = c0
         endif

    !-------------------------------------------------------------------
    ! check conservation
    !-------------------------------------------------------------------
         do k = 1, nbtrcr
            aerotot(k) = aerosno(k,2) + aerosno(k,1) &
                       + zbgc_snow(k) + zbgc_atm(k)
            aero_cons(k) = aerotot(k)-aerotot0(k) &
                          - (    flux_bio_atm(k)*aicen )*dt 
            if (aero_cons(k)  > puny .or. zbgc_snow(k) + zbgc_atm(k) < c0) then             
               write(nu_diag,*) 'Conservation failure: aerosols in snow'
               write(nu_diag,*) 'aerosol tracer:  ',k
               write(nu_diag,*) 'aero_cons(k),puny:', aero_cons(k),puny
               write(nu_diag,*) 'aerotot,aerotot0 ',aerotot(k),aerotot0(k) 
               write(nu_diag,*) ' aerosno(k,2),aerosno(k,1) ', aerosno(k,2),aerosno(k,1)
               write(nu_diag,*) 'flux_bio_atm(k)*aicen*dt', &
               flux_bio_atm(k)*aicen*dt
               write(nu_diag,*) 'zbgc_snow(k)', &
               zbgc_snow(k)
               write(nu_diag,*) 'zbgc_atm(k)', &
               zbgc_atm(k)
            endif
         enddo

    !-------------------------------------------------------------------
    ! reload tracers
    !-------------------------------------------------------------------
         if (vsnon > puny) then
            aerosno(:,:) = aerosno(:,:)/vsnon
            zbgc_snow(:) = zbgc_snow(:)*hs/vsnon  !divide by aicen
            zbgc_atm(:)  = zbgc_atm(:)*hs/vsnon   !divide by aicen
         endif
         do k = 1,nbtrcr
         do n = 1,2
             trcrn(bio_index(k)+nblyr+n)=aerosno(k,n)
         enddo
         enddo
    !-------------------------------------------------------------------
    ! check for negative values
    !-------------------------------------------------------------------
         if (minval(aerosno(:,1)) < -puny  .or. &
            minval(aerosno(:,2)) < -puny) then

            write(nu_diag,*) 'Snow aerosol negative in update_snow_bgc'
            write(nu_diag,*) ' my_task = ',&
                              my_task, &
                              ' printing point n = ',n, &
                              ' i, j = ',i,j
            write(nu_diag,*) 'aicen= '      ,aicen   
            write(nu_diag,*) 'vicen= '      ,vicen   
            write(nu_diag,*) 'vsnon= '      ,vsnon   
            write(nu_diag,*) 'viceold= '    ,vice_old
            write(nu_diag,*) 'vsnoold= '    ,vsno_old
            write(nu_diag,*) 'melts= '      ,melts   
            write(nu_diag,*) 'meltt= '      ,meltt   
            write(nu_diag,*) 'meltb= '      ,meltb   
            write(nu_diag,*) 'congel= '     ,congel  
            write(nu_diag,*) 'snoice= '     ,snoice  
            write(nu_diag,*) 'aero evap sno?= '  ,dhs_evap
            write(nu_diag,*) 'fsnow= '      ,fsnow   
            do k = 1, nbtrcr
              write(nu_diag,*) 'NBTRCR value k = ', k 
              write(nu_diag,*) 'aero snowssl (k)= '    ,aerosno0(k,1)
              write(nu_diag,*) 'aero new snowssl (k)= ',aerosno(k,1)
              write(nu_diag,*) 'aero snowint (k)= '    ,aerosno0(k,2)
              write(nu_diag,*) 'aero new snowint(k)= ',aerosno(k,2)
              write(nu_diag,*) 'flux_bio_atm(k)= ' , flux_bio_atm(k)
              write(nu_diag,*) 'zbgc_snow(k)= '  ,zbgc_snow(k)
              write(nu_diag,*) 'zbgc_atm(k)= '  ,zbgc_atm(k)

              do n = 1,2
                trcrn(bio_index(k)+nblyr+n)=max(trcrn(bio_index(k)+nblyr+n), c0)
              enddo
            enddo
         endif

      end subroutine update_snow_bgc

!=======================================================================
!
! Compute matrix elements for the low order solution of FEM-FCT scheme
! Predictor step
!
! July, 2014 by N. Jeffery
!
      subroutine compute_FCT_matrix &
                                     (C_in, sbdiag, dt,           &
                                      diag, spdiag, rhs, bgrid,   &
                                      igrid, darcyV, dhtop, dhbot,&
                                      iphin_N, iDin, hbri_old,    &
                                      atm_add, bphin_N,           &
                                      C_top, C_bot, Qbot, Qtop,   &
                                      Sink_bot, Sink_top,         &
                                      D_sbdiag, D_spdiag, Ml)

       use ice_constants, only: c1, c0, p5, c2, puny


      real (kind=dbl_kind), dimension(nblyr+1), intent(in) :: &
         C_in            ! Initial (bulk) concentration*hbri_old (mmol/m^2)
                         ! conserved quantity on igrid

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step
     
      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         iDin            ! Diffusivity on the igrid (1/s)

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         iphin_N         ! Porosity with min condition on igrid

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bphin_N, &      ! Porosity with min condition on igrid
         bgrid 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid           ! biology nondimensional grid layer points 

      real (kind=dbl_kind), dimension (nblyr+1), &
         intent(out) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs         , & ! rhs of tri-diagonal matrix eqn.
         ML,           & ! lumped mass matrix
         D_spdiag, D_sbdiag     ! artificial diffusion matrix

      real (kind=dbl_kind), intent(in) :: &
         dhtop         , & ! Change in brine top (m)
         dhbot         , & ! Change in brine bottom (m)
         hbri_old      , & ! brine height (m)
         atm_add       , & ! atm-ice flux
         C_top         , & ! bulk surface source (mmol/m^2)
         C_bot         , & ! bulk bottom source (mmol/m^2)
         darcyV            ! Darcy velocity (m/s

      real (kind=dbl_kind), intent(inout) :: &   ! positive into ice
         Qtop         , & ! top  flux source (mmol/m^2/s)
         Qbot         , & ! bottom flux  source (mmol/m^2/s)
         Sink_bot     , &     ! rest of bottom flux (sink or source) (mmol/m^2/s)
         Sink_top          ! rest oftop flux (sink or source) (mmol/m^2/s)

      ! local variables

      real (kind=dbl_kind) :: &
         vel, vel2, dphi_dx, vel_tot, zspace, dphi

      integer (kind=int_kind) :: &
         k                ! vertical index

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         Q_top,  Q_bot,              & ! surface and bottom source
         K_diag, K_spdiag, K_sbdiag, & ! advection matrix
         S_diag, S_spdiag, S_sbdiag, & ! diffusion matrix
         D_diag, iDin_phi

      real (kind=dbl_kind), dimension (nblyr) :: &
         kvector1, kvectorn1

!---------------------------------------------------------------------
!  Diag (jj) solve for j = 1:nblyr+1
!  spdiag(j) == (j,j+1) solve for j = 1:nblyr otherwise 0
!  sbdiag(j) == (j,j-1) solve for j = 2:nblyr+1 otherwise 0
!---------------------------------------------------------------------
     kvector1(:) = c0
     kvector1(1) = c1 
     kvectorn1(:) = c1
     kvectorn1(1) = c0 

     zspace = c1/real(nblyr,kind=dbl_kind) 
     Qbot = c0
     Qtop = c0
     Sink_bot = c0
     Sink_top = c0

! compute the lumped mass matrix

      ML(:) = zspace
      ML(1) = zspace/c2
      ML(nblyr+1) = zspace/c2

! compute matrix K: K_diag , K_sbdiag, K_spdiag 
! compute matrix S: S_diag, S_sbdiag, S_spdiag

      K_diag(:) = c0
      D_diag(:) = c0
      D_spdiag(:) = c0
      D_sbdiag(:) = c0
      K_spdiag(:) = c0
      K_sbdiag(:) = c0
      S_diag(:) = c0
      S_spdiag(:) = c0
      S_sbdiag(:) = c0
      iDin_phi(:) = c0
      

      iDin_phi(1) = c0  !element 1
      iDin_phi(nblyr+1) = iDin(nblyr+1)/iphin_N(nblyr+1)  !outside ice at bottom
      iDin_phi(nblyr) = p5*(iDin(nblyr+1)/iphin_N(nblyr+1)+iDin(nblyr)/iphin_N(nblyr))

      vel = (bgrid(2)*dhbot - (bgrid(2)-c1)*dhtop)/dt+darcyV/bphin_N(2)
      K_diag(1) = p5*vel/hbri_old   
      dphi_dx = (iphin_N(nblyr+1) - iphin_N(nblyr))/(zspace)  
      vel = (bgrid(nblyr+1)*dhbot - (bgrid(nblyr+1)-c1)*dhtop)/dt  +darcyV/bphin_N(nblyr+1)  
      vel = vel/hbri_old   
      vel2 = (dhbot/hbri_old/dt +darcyV/hbri_old) 
      K_diag(nblyr+1) =   min(c0, vel2) - iDin_phi(nblyr+1)/(zspace+ grid_o/hbri_old) &
                   + p5*(-vel + iDin_phi(nblyr)/bphin_N(nblyr+1)*dphi_dx) 

      do k = 1, nblyr-1
         vel = (bgrid(k+1)*dhbot - (bgrid(k+1)-c1)*dhtop)/dt+darcyV/bphin_N(k+1)
         iDin_phi(k+1) = p5*(iDin(k+2)/iphin_N(k+2) + iDin(k+1)/iphin_N(k+1))
         dphi_dx =  (iphin_N(k+1) - iphin_N(k))/(zspace)
         K_spdiag(k)= p5*(vel/hbri_old - &
                         iDin_phi(k)/(bphin_N(k+1))*dphi_dx)

         vel = (bgrid(k+1)*dhbot - (bgrid(k+1)-c1)*dhtop)/dt  +darcyV/bphin_N(k+1)   
         dphi_dx = c0
         dphi_dx = kvectorn1(k)*(iphin_N(k+1) - iphin_N(k))/(zspace) 
         K_sbdiag(k+1)= -p5*(vel/hbri_old- &
                         iDin_phi(k)/bphin_N(k+1)*dphi_dx)
         K_diag(k) = K_diag(1)*kvector1(k) + (K_spdiag(k) + K_sbdiag(k))*kvectorn1(k)

         S_diag(k+1) =   -(iDin_phi(k)+ iDin_phi(k+1))/zspace
         S_spdiag(k)   = iDin_phi(k)/zspace 
         S_sbdiag(k+1) = iDin_phi(k)/zspace
      enddo

      !k = nblyr

      vel = (bgrid(nblyr+1)*dhbot - (bgrid(nblyr+1)-c1)*dhtop)/dt+darcyV/bphin_N(nblyr+1)
      dphi_dx =  (iphin_N(nblyr+1) - iphin_N(nblyr))/(zspace)
      K_spdiag(nblyr)= p5*(vel/hbri_old - &
                         iDin_phi(nblyr)/(bphin_N(nblyr+1))*dphi_dx)
      vel = (bgrid(nblyr+1)*dhbot - (bgrid(nblyr+1)-c1)*dhtop)/dt  +darcyV/bphin_N(nblyr+1)   
      dphi_dx = kvectorn1(nblyr)*(iphin_N(nblyr+1) - iphin_N(nblyr))/(zspace) 
      K_sbdiag(nblyr+1)= -p5*(vel/hbri_old- &
                         iDin_phi(nblyr)/bphin_N(nblyr+1)*dphi_dx)
      K_diag(nblyr) = K_spdiag(nblyr) + K_sbdiag(nblyr)
      S_diag(nblyr+1) = -iDin_phi(nblyr)/zspace 
      S_spdiag(nblyr)   = iDin_phi(nblyr)/zspace 
      S_sbdiag(nblyr+1) = iDin_phi(nblyr)/zspace
       
! compute matrix artificial D: D_spdiag, D_diag  (D_spdiag(k) = D_sbdiag(k+1))

      do k = 1,nblyr
         D_spdiag(k)    = max(-K_spdiag(k), c0, -K_sbdiag(k+1))
         D_sbdiag(k+1)  = D_spdiag(k)
      enddo
      do k = 1,nblyr+1
         D_diag(k) = D_diag(k) - D_spdiag(k) - D_sbdiag(k)
      enddo

! compute Q_top and Q_bot: top and bottom sources 

      vel2 = -(dhtop/hbri_old/dt +darcyV/bphin_N(1)/hbri_old)

      Q_top(:) = c0
      Q_top(1) = max(c0,vel2*C_top + atm_add/dt)
      Qtop = Q_top(1)

      vel = (dhbot/hbri_old/dt +darcyV/hbri_old)  ! going from iphin_N(nblyr+1) to c1 makes a difference

      Q_bot(:) = c0
      Q_bot(nblyr+1) = max(c0,vel*C_bot) + iDin_phi(nblyr+1)*C_bot&  
                      /(zspace + grid_o/hbri_old)
 
      Qbot = Q_bot(nblyr+1)
 
      Sink_bot = K_diag(nblyr+1) +  K_spdiag(nblyr)
      Sink_top = K_diag(1) + K_sbdiag(2)

!compute matrix elements  (1 to nblyr+1)

     spdiag = -dt * (D_spdiag + K_spdiag + S_spdiag)
     sbdiag = -dt * (D_sbdiag + K_sbdiag + S_sbdiag)
     diag = Ml - dt *  (D_diag + K_diag + S_diag)
     rhs = Ml * C_in + dt * Q_top + dt* Q_bot 
      
     end subroutine compute_FCT_matrix

!=======================================================================
!
! Compute matrices for final solution FCT for passive tracers
! Corrector step
!
! July, 2014 by N. Jeffery
!
      subroutine compute_FCT_corr &
                                     (C_in, C_low, dt,  &
                                      D_sbdiag, D_spdiag, ML)

       use ice_constants, only: c1, c0, c6, puny


      real (kind=dbl_kind), dimension(nblyr+1), intent(in) :: &
         C_in            ! Initial (bulk) concentration*hbri_old (mmol/m^2)
                         ! conserved quantity on igrid

      real (kind=dbl_kind), dimension(nblyr+1), intent(inout) :: &
         C_low           ! Low order solution (mmol/m^2) corrected

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), dimension (nblyr+1), &
         intent(in) :: &
         D_sbdiag      , & ! sub-diagonal artificial diffusion matrix elements
         ML            , & ! Lumped mass diagonal matrix elements
         D_spdiag          ! super-diagonal artificial diffusion matrix elements

      ! local variables

      real (kind=dbl_kind) :: &
        zspace

      integer (kind=int_kind) :: &
         k                ! vertical index

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         M_spdiag, M_sbdiag,         & ! mass matrix
         F_diag, F_spdiag, F_sbdiag, & ! anti-diffusive matrix
         Pplus, Pminus,              & !
         Qplus, Qminus,              & !
         Rplus, Rminus,              & !
         a_spdiag, a_sbdiag            ! weightings of F

!---------------------------------------------------------------------
!  Diag (jj) solve for j = 1:nblyr+1
!  spdiag(j) == (j,j+1) solve for j = 1:nblyr otherwise 0
!  sbdiag(j) == (j,j-1) solve for j = 2:nblyr+1 otherwise 0
!---------------------------------------------------------------------

     zspace = c1/real(nblyr,kind=dbl_kind) 

! compute the mass matrix

      M_spdiag(:) = zspace/c6
      M_spdiag(nblyr+1) = c0
      M_sbdiag(:) = zspace/c6
      M_sbdiag(1) = c0

! compute off matrix F:

      F_diag(:) = c0
      F_spdiag(:) = c0
      F_sbdiag(:) = c0

      do k = 1, nblyr 
           F_spdiag(k) = M_spdiag(k)*(C_low(k)-C_in(k) - C_low(k+1)+ C_in(k+1))/dt &
                       + D_spdiag(k)*(C_low(k)-C_low(k+1))
           F_sbdiag(k+1) =  M_sbdiag(k+1)*(C_low(k+1)-C_in(k+1) - C_low(k)+ C_in(k))/dt &
                       + D_sbdiag(k+1)*(C_low(k+1)-C_low(k))

           if (F_spdiag(k)*(C_low(k) - C_low(k+1)) > c0) F_spdiag(k) = c0
           if (F_sbdiag(k+1)*(C_low(k+1) - C_low(k)) > c0) F_sbdiag(k+1) = c0
     enddo

    if (maxval(abs(F_spdiag)) > c0) then

! compute the weighting factors: a_spdiag, a_sbdiag

      a_spdiag(:) = c0
      a_sbdiag(:) = c0

      Pplus(1)  = max(c0, F_spdiag(1))  
      Pminus(1) = min(c0, F_spdiag(1))
      Pplus(nblyr+1)  = max(c0, F_sbdiag(nblyr+1))
      Pminus(nblyr+1) = min(c0, F_sbdiag(nblyr+1))
      Qplus(1) = max(c0,C_low(2)-C_low(1))
      Qminus(1)= min(c0,C_low(2)-C_low(1))
      Qplus(nblyr+1) = max(c0,C_low(nblyr)-C_low(nblyr+1))
      Qminus(nblyr+1)= min(c0,C_low(nblyr)-C_low(nblyr+1))
      Rplus(1)  = min(c1, ML(1)*Qplus(1)/dt/(Pplus(1)+puny))
      Rminus(1) = min(c1, ML(1)*Qminus(1)/dt/(Pminus(1)-puny))
      Rplus(nblyr+1)  = min(c1, ML(nblyr+1)*Qplus(nblyr+1)/dt/(Pplus(nblyr+1)+puny))
      Rminus(nblyr+1) = min(c1, ML(nblyr+1)*Qminus(nblyr+1)/dt/(Pminus(nblyr+1)-puny))
      do k = 2,nblyr
         Pplus(k)  = max(c0,F_spdiag(k)) + max(c0,F_sbdiag(k))
         Pminus(k) = min(c0,F_spdiag(k)) + min(c0,F_sbdiag(k))
         Qplus(k)  = max(c0, max(C_low(k+1)-C_low(k),C_low(k-1)-C_low(k)))
         Qminus(k) = min(c0, min(C_low(k+1)-C_low(k),C_low(k-1)-C_low(k)))
         Rplus(k)  = min(c1, ML(k)*Qplus(k)/dt/(Pplus(k)+puny))
         Rminus(k) = min(c1, ML(k)*Qminus(k)/dt/(Pminus(k)-puny))
      enddo
     
      do k = 1, nblyr 
         a_spdiag(k) = min(Rminus(k),Rplus(k+1))
         if (F_spdiag(k) > c0) a_spdiag(k) = min(Rplus(k),Rminus(k+1))
         a_sbdiag(k+1) = min(Rminus(k+1),Rplus(k))
         if (F_sbdiag(k+1) > c0) a_sbdiag(k+1) = min(Rplus(k+1),Rminus(k))
      enddo

     ! if (minval(a_spdiag) < c0 .or. maxval(a_spdiag) > c1 .or. &
     !     minval(a_sbdiag) < c0 .or. maxval(a_sbdiag) > c1) then
     !     write(nu_diag,*) 'Error in compute_FCT_corr, alpha out of bounds'
     !     write(nu_diag,*) 'a_spdiag:', a_spdiag
     !     write(nu_diag,*) 'a_sbdiag:', a_sbdiag
     !     write(nu_diag,*) 'F_spdiag:', F_spdiag
     !     write(nu_diag,*) 'F_sbdiag:', F_sbdiag
     !     write(nu_diag,*) 'C_low:', C_low
     !     write(nu_diag,*) 'C_in:', C_in
     !  endif

!compute F_diag:

     F_diag(1) = a_spdiag(1)*F_spdiag(1)
     F_diag(nblyr+1) = a_sbdiag(nblyr+1)* F_sbdiag(nblyr+1) 
     C_low(1) = C_low(1) + dt*F_diag(1)/ML(1)
     C_low(nblyr+1) = C_low(nblyr+1) + dt*F_diag(nblyr+1)/ML(nblyr+1)

     do k = 2,nblyr
         F_diag(k) = a_spdiag(k)*F_spdiag(k) + a_sbdiag(k)*F_sbdiag(k)
         C_low(k) = C_low(k) + dt*F_diag(k)/ML(k)
     enddo
      
     endif  !F_spdiag is nonzero

     end subroutine compute_FCT_corr
  
!=======================================================================
!
! Tridiagonal matrix solver-- for salinity
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
      subroutine tridiag_solverz (nmat,     sbdiag,   &
                                 diag,      spdiag,   &
                                 rhs,       xout)

      integer (kind=int_kind), intent(in) :: &
         nmat            ! matrix dimension

      real (kind=dbl_kind), dimension (nmat), &
           intent(in) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      real (kind=dbl_kind), dimension (nmat), &
           intent(inout) :: &
         xout            ! solution vector

      ! local variables     

      integer (kind=int_kind) :: &
         k               ! row counter

      real (kind=dbl_kind) :: &
         wbeta           ! temporary matrix variable

      real (kind=dbl_kind), dimension(nmat):: &
         wgamma          ! temporary matrix variable

         wbeta = diag(1)
         xout(1) = rhs(1) / wbeta

      do k = 2, nmat
            wgamma(k) = spdiag(k-1) / wbeta
            wbeta = diag(k) - sbdiag(k)*wgamma(k)
            xout(k) = (rhs(k) - sbdiag(k)*xout(k-1)) &
                         / wbeta
      enddo                     ! k

      do k = nmat-1, 1, -1
            xout(k) = xout(k) - wgamma(k+1)*xout(k+1)
      enddo                     ! k

      end subroutine tridiag_solverz

!=======================================================================
!
! authors     Nicole Jeffery, LANL

      subroutine check_conservation_FCT &
                                     (C_init, C_new, C_low, S_top, &
                                      S_bot, L_bot, L_top, dt,     &
                                      fluxbio, good_numerics) 

      use ice_constants, only: p5, c1, c4, c0

      real (kind=dbl_kind), dimension(nblyr+1), intent(in) :: &
         C_init        , & ! initial bulk concentration * h_old (mmol/m^2)
         C_new             ! new bulk concentration * h_new (mmol/m^2)

      real (kind=dbl_kind), dimension(nblyr+1), intent(out) :: &
         C_low             ! define low order solution = C_new

      real (kind=dbl_kind),  intent(in) :: &
         S_top         , & ! surface flux into ice (mmol/m^2/s)
         S_bot         , & ! bottom flux into ice (mmol/m^2/s)
         L_bot         , & !remaining  bottom flux into ice (mmol/m^2/s)
         L_top         , & !remaining  top  flux into ice (mmol/m^2/s)
         dt         

      real (kind=dbl_kind), intent(inout) :: &
         fluxbio            ! (mmol/m^2/s)  positive down (into the ocean)

      logical (kind=log_kind), intent(inout) :: &   
         good_numerics    ! true if conservation satisfied within error

     ! local variables

    integer (kind=int_kind) :: &
         k

     real (kind=dbl_kind) :: &
         diff_dt     , &
         C_init_tot  , &
         C_new_tot   , &
         zspace      , &  !1/nblyr
         accuracy          ! centered difference is Order(zspace^2)
        
         zspace = c1/real(nblyr,kind=dbl_kind)
         good_numerics = .true.

     !-------------------------------------
     !  Ocean flux: positive into the ocean
     !-------------------------------------    
         C_init_tot = (C_init(1) + C_init(nblyr+1))*zspace*p5
         
         C_new_tot = (C_new(1) + C_new(nblyr+1))*zspace*p5
         C_low(1) = C_new(1)
         C_low(nblyr+1) = C_new(nblyr+1)

         do k = 2, nblyr
            C_init_tot = C_init_tot + C_init(k)*zspace
            C_new_tot = C_new_tot + C_new(k)*zspace
            C_low(k) = C_new(k)
         enddo

         accuracy = 1.0e-14_dbl_kind*max(c1, C_init_tot, C_new_tot)  
         fluxbio = fluxbio - S_bot - L_bot*C_new(nblyr+1)-L_top*C_new(1)
         diff_dt =C_new_tot - C_init_tot - (S_top+S_bot+ L_bot*C_new(nblyr+1)+L_top*C_new(1))*dt

         if (minval(C_low) < c0) then 
           write(nu_diag,*) 'Positivity of zbgc low order solution failed: C_low:',C_low
           good_numerics = .false.
         endif
           
         if (abs(diff_dt) > accuracy ) then
           good_numerics = .false.
           write(nu_diag,*) 'Conservation of zbgc low order solution failed: diff_dt:',&
                        diff_dt
           write(nu_diag,*) 'Total initial tracer', C_init_tot
           write(nu_diag,*) 'Total final1  tracer', C_new_tot
           write(nu_diag,*) 'bottom final tracer', C_new(nblyr+1)
           write(nu_diag,*) 'top final tracer', C_new(1)
           write(nu_diag,*) 'Near bottom final tracer', C_new(nblyr)
           write(nu_diag,*) 'Near top final tracer', C_new(2)
           write(nu_diag,*) 'Top flux*dt into ice:', S_top*dt
           write(nu_diag,*) 'Bottom flux*dt into ice:', S_bot*dt
           write(nu_diag,*) 'Remaining bot flux*dt into ice:', L_bot*C_new(nblyr+1)*dt
           write(nu_diag,*) 'Remaining top flux*dt into ice:', L_top*C_new(1)*dt
         endif
         
     end subroutine check_conservation_FCT

!=======================================================================
!
! authors     Nicole Jeffery, LANL

      subroutine compute_shortwave_trcr &
                                    (nx_block,    ny_block,  &
                                    indxii,       indxjj,    &
                                    icells,       n_algae,   &
                                    trcrn,        trcrn_sw,  &
                                    swgrid,       hin,       &
                                    hbri,         ntrcr,     &
                                    nilyr,        nblyr,     &
                                    igrid,                   &
                                    nbtrcr_sw,    n_zaero )   
      
      use ice_shortwave, only: hi_ssl
      use ice_constants, only: c0, c1, c2, p5
      use ice_domain_size, only: nslyr
      use ice_state, only: nt_bgc_N, nt_zaero

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells , n_zaero  , & ! number of cells with aicen > puny 
         nbtrcr_sw, n_algae, & ! nilyr+nslyr+2 for chlorophyll
         ntrcr

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxii, indxjj ! compressed indices for icells with aicen > puny

      integer (kind=int_kind), intent(in) :: &
         nblyr         , & ! number of bio layers
         nilyr             ! number of ice layers 

      real (kind=dbl_kind), dimension (nx_block, ny_block,ntrcr), &
         intent(in) ::       &
         trcrn          ! aerosol or chlorophyll

      real (kind=dbl_kind), dimension (nx_block, ny_block, nbtrcr_sw), &
         intent(inout) ::    &
         trcrn_sw       ! ice on shortwave grid tracers

      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         swgrid         ! 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid          ! CICE bio grid 
         
      real(kind=dbl_kind), dimension (icells), intent(in) :: &
         hin          , & ! CICE ice thickness
         hbri             ! brine height 

      !  local variables

      integer (kind=int_kind) :: k, n, ij, i, j

      real (kind=dbl_kind), dimension (icells, ntrcr+2) :: &
         trtmp0, &       ! temporary, remapped tracers
         trtmp

      real (kind=dbl_kind), dimension (icells,nilyr+1):: &
         icegrid         ! correct for large ice surface layers

      real (kind=dbl_kind):: &
         top_conc        ! 1% (min_bgc) of surface concentration 
                         ! when hin > hbri:  just used in sw calculation

      !-----------------------------------------------------------------
      ! Compute aerosols and algal chlorophyll on shortwave grid
      !-----------------------------------------------------------------
      trtmp0(:,:) = c0
      trtmp(:,:) = c0
    
      if (tr_bgc_N)  then
      do k = 1, nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
           do n = 1, n_algae
           do ij = 1, icells    
             i = indxii(ij)
             j = indxjj(ij)      
             trtmp0(ij,nt_bgc_N(1) + k-1) = trtmp0(ij,nt_bgc_N(1) + k-1) + &
                                R_chl2N(n)*F_abs_chl(n)*trcrn(i,j,nt_bgc_N(n)+k-1)
           enddo   !icells
           enddo   !n
      enddo   !k

!DIR$ CONCURRENT !Cray      
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells    
           i = indxii(ij)
           j = indxjj(ij)   
           icegrid(ij,:) = swgrid(:)
           if (swgrid(1)*hin(ij)*c2 > hi_ssl) then
                icegrid(ij,1) = hi_ssl/c2/hin(ij)
                icegrid(ij,2) = p5/real(nilyr,kind=dbl_kind)+ icegrid(ij,1)
           endif
   
           top_conc = trtmp0(ij,nt_bgc_N(1))*min_bgc
           call remap_zbgc  (ntrcr,  nilyr+1,            &
                             nt_bgc_N(1),                &
                             trtmp0(ij,1:ntrcr),         &
                             trtmp(ij,1:ntrcr+2),        &
                             1,                nblyr+1,  &
                             hin(ij),          hbri(ij), &
                             icegrid(ij,1:nilyr+1),      &
                             igrid(1:nblyr+1), top_conc) 

       enddo  
       
       do k = 1,nilyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          do ij = 1, icells  
             i = indxii(ij)
             j = indxjj(ij) 
             trcrn_sw(i,j,nlt_chl_sw+nslyr+k) = trtmp(ij,nt_bgc_N(1) + k-1)
          enddo    ! ij
       enddo       !k
       do n = 1,n_algae   ! snow contribution
       do ij = 1, icells  
            i = indxii(ij)
            j = indxjj(ij) 
            trcrn_sw(i,j,nlt_chl_sw)= trcrn_sw(i,j,nlt_chl_sw) + &
                     R_chl2N(n)*F_abs_chl(n)*trcrn(i,j,nt_bgc_N(n)+nblyr+1) 
                              !snow surface layer
            trcrn_sw(i,j,nlt_chl_sw+1:nlt_chl_sw+nslyr)= trcrn_sw(i,j,nlt_chl_sw+1:nlt_chl_sw+nslyr) + &
                     R_chl2N(n)*F_abs_chl(n)*trcrn(i,j,nt_bgc_N(n)+nblyr+2) 
                              !only 1 snow layer in zaero
       enddo       !ij
       enddo       !n
     endif         ! tr_bgc_N

     if (tr_zaero) then
    
     do n = 1, n_zaero
      trtmp0(:,:) = c0
      trtmp(:,:)  = c0
      do k = 1, nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
           do ij = 1, icells    
             i = indxii(ij)
             j = indxjj(ij)      
             trtmp0(ij,nt_zaero(n) + k-1) = trcrn(i,j,nt_zaero(n)+k-1)
           enddo   !icells
       enddo   !k

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells    
           i = indxii(ij)
           j = indxjj(ij) 
           icegrid(ij,:) = swgrid(:)
           if (swgrid(1)*hin(ij)*c2 > hi_ssl) then
                icegrid(ij,1) = hi_ssl/c2/hin(ij)
                icegrid(ij,2) = p5/real(nilyr,kind=dbl_kind)+ icegrid(ij,1)
           endif
           top_conc = trtmp0(ij,nt_zaero(n))*min_bgc
           call remap_zbgc  (ntrcr,  nilyr+1,            &
                             nt_zaero(n),                &
                             trtmp0(ij,1:ntrcr),        &
                             trtmp(ij,1:ntrcr+2),       &
                             1,                nblyr+1,  &
                             hin(ij),          hbri(ij), &
                             icegrid(ij,1:nilyr+1),         &
                             igrid(1:nblyr+1), top_conc) 

       enddo  
       
       do k = 1,nilyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          do ij = 1, icells  
             i = indxii(ij)
             j = indxjj(ij) 
             trcrn_sw(i,j,nlt_zaero_sw(n)+nslyr+k) = trtmp(ij,nt_zaero(n) + k-1)
          enddo    ! ij
       enddo       !k
       do ij = 1, icells  
            i = indxii(ij)
            j = indxjj(ij) 
            trcrn_sw(i,j,nlt_zaero_sw(n))= trcrn(i,j,nt_zaero(n)+nblyr+1) !snow ssl
            trcrn_sw(i,j,nlt_zaero_sw(n)+1:nlt_zaero_sw(n)+nslyr)= trcrn(i,j,nt_zaero(n)+nblyr+2)
       enddo       !ij
     enddo       !n
     endif        !tr_zaero
  

     end subroutine compute_shortwave_trcr

!=======================================================================

! read atmospheric aerosols
!
! authors: Elizabeth Hunke, LANL

      subroutine fzaero_data

      use ice_calendar, only: month, mday, istep, sec
      use ice_domain_size, only: max_blocks, n_zaero
      use ice_blocks, only: nx_block, ny_block
      use ice_zbgc_shared, only: flux_bio_atm
      use ice_forcing, only: interp_coeff_monthly, read_clim_data_nc, interpolate_data
      use ice_constants, only: c0, field_type_scalar, field_loc_center

#ifdef ncdf 
      ! local parameters

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), &
         save :: &
         aero_data    ! field values at 2 temporal data points

      character (char_len_long) :: & 
         aero_file,   &   ! netcdf filename
         fieldname        ! field name in netcdf file

      integer (kind=int_kind) :: & 
         ixm,ixp     , & ! record numbers for neighboring months
         maxrec      , & ! maximum record number
         recslot     , & ! spline slot for current record
         midmonth        ! middle day of month

      logical (kind=log_kind) :: readm

    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------

      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(month+maxrec-2,maxrec) + 1
      ixp  = mod(month,         maxrec) + 1
      if (mday >= midmonth) ixm = -99  ! other two points will be used
      if (mday <  midmonth) ixp = -99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values 
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

!      aero_file = trim(atm_data_dir)//'faero.nc'   
      ! Cam5 monthly total black carbon deposition on the gx1 grid"
      aero_file = '/usr/projects/climate/njeffery/DATA/CAM/Hailong_Wang/Cam5_bc_monthly_popgrid.nc'   

      fieldname='bcd'
      call read_clim_data_nc (readm, 0,  ixm, month, ixp, &
                              aero_file, fieldname, aero_data, &
                              field_loc_center, field_type_scalar)


      call interpolate_data (aero_data, flux_bio_atm(:,:,nlt_zaero(1),:)) ! kg/m^2/s

      where (flux_bio_atm(:,:,nlt_zaero(1),:) > 1.e20) flux_bio_atm(:,:,nlt_zaero(1),:) = c0

#endif

      end subroutine fzaero_data

!=======================================================================

! Initialize ocean iron from file
!
! authors: Nicole Jeffery, LANL

      subroutine init_bgc_data (fed1,fep1)

      use ice_domain_size, only: max_blocks
      use ice_blocks, only: nx_block, ny_block
      use ice_read_write, only: ice_open_nc, ice_read_nc, ice_close_nc
      use ice_constants, only: c0, p1

#ifdef ncdf
      use netcdf
#endif
           
      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(out) :: &
           fed1, &  ! first dissolved iron pool (nM)
           fep1    ! first particulate iron pool (nM)

      ! local parameters

      integer (kind=int_kind) :: &
         fid              , & ! file id for netCDF file 
         nbits

      logical (kind=log_kind) :: diag

      character (char_len_long) :: & 
         iron_file,   &   ! netcdf filename
         fieldname        ! field name in netcdf file

      nbits = 64              ! double precision data

    !-------------------------------------------------------------------
    ! Annual average data from Tagliabue, 2012 (top 50 m average
    ! poisson grid filled on gx1v6
    !-------------------------------------------------------------------
      fed1(:,:,:) = c0

      if (trim(fe_data_type) == 'clim') then
       	diag = .true.   ! write diagnostic information 
        iron_file = trim(bgc_data_dir)//'dFe_50m_annual_Tagliabue_gx1.nc'

        if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'Dissolved iron ocean concentrations from:'
            write (nu_diag,*) trim(iron_file)
            call ice_open_nc(iron_file,fid)
        endif

        fieldname='dFe'
        ! Currently only first fed  value is read
        call ice_read_nc(fid,1,fieldname,fed1,diag) 
        where ( fed1(:,:,:) > 1.e20) fed1(:,:,:) = p1  

        if (my_task == master_task) call ice_close_nc(fid)  

       	diag = .true.   ! write diagnostic information 
        iron_file = trim(bgc_data_dir)//'pFe_bathy_gx1.nc'

        if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'Particulate iron ocean concentrations from:'
            write (nu_diag,*) trim(iron_file)
            call ice_open_nc(iron_file,fid)
        endif

        fieldname='pFe'
        ! Currently only first fep value is read
        call ice_read_nc(fid,1,fieldname,fep1,diag) 
        where ( fep1(:,:,:) > 1.e20) fep1(:,:,:) = p1  

        if (my_task == master_task) call ice_close_nc(fid)  
    
      endif
    
      end subroutine init_bgc_data

!=======================================================================
!
! authors     Nicole Jeffery, LANL

      subroutine regrid_stationary &
                                    (C_stationary,           &
                                    hbri_old,                &
                                    hbri,                    &
                                    meltb,                   &
                                    flux_bio,                &
                                    dt,                      &
                                    nblyr,        igrid)
      
      use ice_constants, only: c0, c1, p5

      integer (kind=int_kind), intent(in) :: &
         nblyr             ! number of bio layers

      real (kind=dbl_kind), intent(inout) ::  &
         flux_bio          ! ocean tracer flux (mmol/m^2/s) positive into ocean
 
      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) ::  &     
         C_stationary      ! stationary bulk concentration*h (mmol/m^2)

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid          ! CICE bio grid 
         
      real(kind=dbl_kind),  intent(in) :: &
         dt           , & ! time step
         meltb        , & ! bottom melt
         hbri_old     , & ! previous timestep brine height
         hbri             ! brine height 

      !  local variables

      integer (kind=int_kind) :: k, n, nt

      real (kind=dbl_kind), dimension (nblyr+3) :: &
         trtmp0,   &    ! temporary, remapped tracers
         trtmp

      real (kind=dbl_kind):: &
         htemp,    &    ! ice thickness after melt (m)
         top_conc, &    ! 1% (min_bgc) of surface concentration 
         zspace,   &    ! bio grid spacing
         sum_old,  &    ! total tracer before melt loss
         sum_new        ! total tracer after melt

      zspace = c1/(real(nblyr,kind=dbl_kind))
      trtmp0(:) = c0
      trtmp(:) = c0
      nt = 1
      do k = 1, nblyr+1
         trtmp0(nblyr+2-k) = C_stationary(k)/hbri_old  ! reverse order
      enddo   !k
      htemp = min(c0,hbri - meltb)
      if (meltb > c0 .and. htemp > c0) then
      !-----------------------------------------------------------------
      ! Regrid C_stationary to remove bottom melt
      !-----------------------------------------------------------------
           top_conc =c0
           call remap_zbgc  (nblyr+1,  nblyr+1,          &
                             nt,                         &
                             trtmp0(1:nblyr+1),          &
                             trtmp(1:nblyr+3),           &
                             0,                nblyr+1,  &
                             hbri_old,         htemp,    &
                             igrid(1:nblyr+1),      &
                             igrid(1:nblyr+1), top_conc) 
    
           sum_old = (trtmp(1)+trtmp(nblyr+1))*htemp*p5*zspace
           sum_new = (C_stationary(1) + C_stationary(nblyr+1))*p5*zspace
           trtmp0(nblyr+1) = trtmp(1)
           trtmp0(1) = trtmp(nblyr+1)
           do k = 2,nblyr
              sum_old = sum_old + C_stationary(k)*zspace
              trtmp0(nblyr+2-k) = trtmp(nt + k-1)
              sum_new = sum_new + trtmp0(nblyr+2-k)*htemp*p5*zspace
           enddo       !k
           flux_bio = flux_bio + (sum_old - sum_new)/dt
       elseif (hbri > hbri_old) then
      !-----------------------------------------------------------------
      ! Regrid C_stationary to migrate if there's bottom growth
      !-----------------------------------------------------------------
           top_conc =c0
           call remap_zbgc  (nblyr+1,  nblyr+1,          &
                             nt,                         &
                             trtmp0(1:nblyr+1),          &
                             trtmp(1:nblyr+3),           &
                             0,                nblyr+1,  &
                             hbri_old,         hbri,     &
                             igrid(1:nblyr+1),           &
                             igrid(1:nblyr+1), top_conc) 

           do k = 1,nblyr+1
              trtmp0(nblyr+2-k) = trtmp(nt+k-1)
           enddo       !k
           flux_bio = flux_bio  ! no change
       endif 
      !----------------------------------------------------------------- 
      ! change in the brine height (liquid volume) rescales the tracer
      !-----------------------------------------------------------------
       do k = 1, nblyr+1
         C_stationary(k) = trtmp0(k)*hbri
       enddo   !k

     end subroutine regrid_stationary

!=======================================================================

      end module ice_algae

!=======================================================================
