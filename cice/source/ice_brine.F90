!=======================================================================
!
! Computes ice microstructural information for use in biogeochemistry
!
! authors: Nicole Jeffery, LANL
!
      module ice_brine

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: nilyr, nblyr, nblyr_hist, max_blocks, ncat
      use ice_fileunits, only: nu_diag, nu_rst_pointer, nu_dump_hbrine, &
          nu_restart_hbrine, flush_fileunit
      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_state, only: ntrcr, nt_qice, nt_sice
      use ice_zbgc_shared, only: cgrid, bgrid, igrid, exp_h, k_o, rhosi, &
           his_min, thinS, min_salin, igrid, remap_layers_bgc, &
           phi_snow, restart_hbrine, first_ice

      implicit none

      private
      public :: preflushing_changes, compute_microS_mushy, &
                update_hbrine, init_hbrine, write_restart_hbrine, &
                hbrine_diags
 
      real (kind=dbl_kind), parameter :: &   
         maxhinS = 1.25_dbl_kind  , & ! brine overflows if hinS > maxhinS*hin
         viscos  = 2.1e-6_dbl_kind, & ! kinematic viscosity (m^2/s) 
         ! Brine salinity as a cubic function of temperature
         a1      = -21.4_dbl_kind , & ! (psu/C)  
         a2      = -0.886_dbl_kind, & ! (psu/C^2)
         a3      = -0.012_dbl_kind, & ! (psu/C^3)
         ! Brine density as a quadratic of brine salinity
         b1      = 1000.0_dbl_kind, & ! (kg/m^3)  
         b2      = 0.8_dbl_kind       ! (kg/m^3/ppt)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat,max_blocks) :: &
         first_ice_real     ! .true. = c1, .false. = c0

!=======================================================================

      contains

!=======================================================================

!  Initialize brine height tracer

      subroutine init_hbrine 

      use ice_state, only: nt_fbri, trcrn

      integer (kind=int_kind) :: &
           k               ! vertical index

      real (kind=dbl_kind) :: & 
         zspace            ! grid spacing for CICE vertical grid

      !-----------------------------------------------------------------
      ! Calculate bio gridn: 0 to 1 corresponds to ice top to bottom 
      !-----------------------------------------------------------------

!echmod:  why nblyr_hist?  is this for history variables?

      bgrid(:)          = c0 ! bgc grid points         
      bgrid(nblyr_hist) = c1 ! bottom value
      igrid(:)          = c0 ! bgc interface grid points   
      igrid(1)          = c0 ! ice top
      igrid(nblyr+1)    = c1 ! ice bottom
      
      zspace = c1/(real(nblyr,kind=dbl_kind)) 
      do k = 2, nblyr+1
          bgrid(k) = zspace*(real(k,kind=dbl_kind) - c1p5)
      enddo
      
      do k = 2, nblyr
        igrid(k) = p5*(bgrid(k+1)+bgrid(k))
      enddo

      !-----------------------------------------------------------------
      ! Calculate CICE cgrid for interpolation ice top (0) to bottom (1) 
      !-----------------------------------------------------------------
       
      cgrid(1) = c0                           ! CICE vertical grid top point
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing
    
      do k = 2, nilyr+1
        cgrid(k) = zspace * (real(k,kind=dbl_kind) - c1p5) 
      enddo 

      !-----------------------------------------------------------------
      ! initialize restart variables
      !-----------------------------------------------------------------

      if (restart_hbrine) then
          call read_restart_hbrine
      else
          first_ice(:,:,:,:) = .true.            
          trcrn(:,:,nt_fbri,:,:) = c1
      endif

      end subroutine init_hbrine

!=======================================================================

! Computes the top and bottom brine boundary changes for flushing
! works for zsalinity and tr_salinity
!
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice with 
! dynamic salinity or the height ratio == hinS/vicen*aicen, where hinS is the 
! height of the brine surface relative to the bottom of the ice.  This volume fraction
! may be > 1 in which case there is brine above the ice surface (meltponds). 

      subroutine preflushing_changes (nx_block, ny_block,             &
                                      icells,   n_cat,                &
                                      indxii,   indxjj,               &
                                      aicen,    vicen,    vsnon,      &
                                      meltb,    meltt,    congel,     &
                                      snoice,   hice_old, fbri,       & 
                                      dh_top,   dh_bot,               &
                                      dhi_top,  dhi_bot,              &
                                      hinS_o,   hin,hsn,  firstice)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of cells with aicen > 0
         n_cat                 ! category
                              
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxii, indxjj ! compressed indices for icells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon       , & ! volume per unit area of snow         (m)
         meltb       , & ! bottom ice melt                      (m)
         meltt       , & ! top ice melt                         (m)
         congel      , & ! bottom ice growth                    (m)
         snoice          ! top ice growth from flooding         (m)
 
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: &
         fbri            ! trcrn(i,j,nt_fbri)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         dh_top       , & ! brine change in top for diagnostics (m)
         dh_bot       , & ! brine change in bottom for diagnostics (m)
         dhi_top      , & ! ice change in top for diagnostics (m)
         dhi_bot      , & ! ice change in bottom for diagnostics (m)
         hin          , & ! ice thickness (m) 
         hsn          , & ! snow thickness (m) 
         hice_old     , & ! old ice thickness (m)
         hinS_o           ! old brine height (m)

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(in) :: &
         firstice         ! if true, initialized values should be used     

      ! local variables

      integer (kind=int_kind) :: &
         i, j         , & ! horizontal indices
         ij               ! horizontal index, combines i and j loops
  
      real (kind=dbl_kind) :: &
         hin_old          ! ice thickness before current melt/growth (m)

      real (kind=dbl_kind):: &
         dhice            ! Change in hin due to subl/cond  (m)

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------

      dh_top    (:,:) = c0
      dh_bot    (:,:) = c0
      dhi_top   (:,:) = c0
      dhi_bot   (:,:) = c0
      hin       (:,:) = c0
      hsn       (:,:) = c0
      hinS_o    (:,:) = c0       

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
      
         i = indxii(ij)
         j = indxjj(ij)
!echmod would hin ever be <0? 
         hin(i,j) = max(c0, vicen(i,j) / aicen(i,j))
         if (hin(i,j) <= c0  .OR. fbri(i,j) <= c0) then
            write(nu_diag, *) 'preflushing: hin <= 0 or fbri <= c0:',i,j
            write(nu_diag, *) 'vicen, aicen', vicen(i,j), aicen(i,j)
            write(nu_diag, *) 'fbri, hice_old', fbri(i,j), hice_old(i,j)
            call abort_ice ('ice_brine error')
         endif
         hsn(i,j) = c0
!echmod would hsn ever be <0? 
         if (hin(i,j) > c0) hsn(i,j) = max(c0, vsnon(i,j) / aicen(i,j))
         hin_old = max(c0, hin(i,j) + meltb(i,j) + meltt(i,j) &
                                    - congel(i,j) - snoice(i,j))
         dhice = hin_old - hice_old(i,j)   ! change due to subl/cond
         dhi_top(i,j) = meltt(i,j) - dhice - snoice(i,j)
         dhi_bot(i,j) = congel(i,j) - meltb(i,j)   
         dh_top(i,j) = dhi_top(i,j)
         dh_bot(i,j) = dhi_bot(i,j)

         if ((hice_old(i,j) < puny) .OR. (hin_old < puny) &
                                    .OR. firstice(i,j)) then
            hin_old         = hin(i,j) 
            dh_top    (i,j) = c0
            dh_bot    (i,j) = c0
            dhi_bot   (i,j) = c0
            dhi_top   (i,j) = c0
            fbri(i,j)       = c1 
         endif

         hinS_o(i,j) = fbri(i,j)* hice_old(i,j)

      enddo  ! ij

      end subroutine preflushing_changes

!=======================================================================

! Computes ice microstructural properties for updating hbrine
!
! NOTE: This subroutine uses thermosaline_vertical output to compute
! average ice permeability and the surface ice porosity

      subroutine compute_microS_mushy (nx_block, ny_block,               &
                                       icells,   n_cat,                  &
                                       indxii,   indxjj,                 &
                                       trcrn,    hice_old,   hinS_o,     &
                                       sss,      sst,        bTin,       &
                                       bphin,    kperm,      zphi_min,   &
                                       bSin,     brine_sal,  brine_rho,  &
                                       iphin,    ibrine_rho, ibrine_sal, &
                                       sice_rho)

      use ice_therm_mushy, only: temperature_mush, liquid_fraction, permeability

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells , &            ! number of cells with aicen > 0
         n_cat                 ! ice category
                              
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxii, indxjj ! compressed indices for icells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         hice_old    , & ! previous timestep ice height (m)
         sss         , & ! ocean salinity (ppt)
         sst             ! ocean temperature (C)
       
      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(in) :: &
         trcrn           

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(out) :: & 
         zphi_min    , & ! surface porosity
         kperm           ! average ice permeability (m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         hinS_o      , & ! previous timestep brine height (m)
         sice_rho        ! average ice density

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), &
         intent(inout)  :: &
         iphin       , & ! porosity on the igrid 
         ibrine_rho  , & ! brine rho on interface  
         ibrine_sal      ! brine sal on interface   

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), &
         intent(inout)  :: &
         bSin        , & ! bulk salinity (ppt) on bgrid
         bTin        , & ! temperature on bgrid
         bphin       , & ! porosity on bgrid
         brine_sal   , & ! equilibrium brine salinity (ppt) 
         brine_rho       ! internal brine density (kg/m^3) 

      ! local variables

      real (kind=dbl_kind), dimension (nx_block*ny_block,nilyr) :: &
         cSin        , & ! bulk salinity (ppt)
         cqin            ! enthalpy ()

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2) :: &
         zTin        , & ! Temperature of ice layers on bgrid (C) 
         zSin        , & ! Salinity of ice layers on bgrid (C) 
         zqin            ! enthalpy on the bgrid ()

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k               ! vertical biology layer index 
      
      real (kind=dbl_kind), dimension(icells) :: &
         surface_S   , & ! salinity of ice above hin > hinS 
         hinc_old    , & ! mean ice thickness before current melt/growth (m)
         hinSc_old       ! mean brine thickness before current melt/growth (m)
      
!      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: &
      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2) :: &
         trtmp_s     , & ! temporary, remapped tracers   
         trtmp_q         ! temporary, remapped tracers   
     
      !-----------------------------------------------------------------
      ! Define ice salinity and temperature on bgrid
      !-----------------------------------------------------------------

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxii(ij)
            j = indxjj(ij)            
            cSin(ij,k) = trcrn(i,j,nt_sice+k-1)
            cqin(ij,k) = trcrn(i,j,nt_qice+k-1)
         enddo
      enddo
        
      trtmp_s(:,:,:) = c0
      trtmp_q(:,:,:) = c0
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells   ! map Sin and qin profiles to bgc grid 
         i = indxii(ij)
         j = indxjj(ij)
         hinS_o   (i,j) = min(hinS_o(i,j), maxhinS*hice_old(i,j))
         hinc_old (ij)  = hice_old(i,j)
         hinSc_old(ij)  = hinS_o  (i,j)

         call remap_layers_bgc (ntrcr,            nilyr,          &
                                nt_sice,                          &
                                trcrn(i,j,:),     trtmp_s(i,j,:), &
                                0,                nblyr+1,        &
                                hinc_old(ij),     hinc_old(ij),   &
                                cgrid(2:nilyr+1),                 &
                                bgrid(1:nblyr+1), surface_S(ij))
     
         call remap_layers_bgc (ntrcr,            nilyr,          &
                                nt_qice,                          &
                                trcrn(i,j,:),     trtmp_q(i,j,:), &
                                0,                nblyr+1,        &
                                hinc_old(ij),     hinc_old(ij),   &
                                cgrid(2:nilyr+1),                 &
                                bgrid(1:nblyr+1), surface_S(ij))
      enddo

      do k = 1,nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells   
            i = indxii(ij)
            j = indxjj(ij)

            zqin (ij,k)  = min(c0,        trtmp_q(i,j,k))

            zSin (ij,k)  = max(min_salin, trtmp_s(i,j,k))
            bSin (i,j,k) = zSin(ij,k)
            bSin (i,j,nblyr+2) = sss(i,j) 

            zTin (ij,k)  = temperature_mush(zqin(ij,k), zSin(ij,k))
            bTin (i,j,k) = zTin(ij,k)
            bTin (i,j,nblyr+2) = sst(i,j)

            bphin(i,j,k) = liquid_fraction (zTin(ij,k), zSin(ij,k))
            bphin(i,j,nblyr+2) = c1
         enddo ! ij                          
      enddo    ! k

      !-----------------------------------------------------------------
      ! Define ice multiphase structure
      !-----------------------------------------------------------------
     
      call prepare_hbrine (icells,     indxii, indxjj, &
                           bSin,       bTin,           &
                           brine_sal,  brine_rho,      &
                           ibrine_sal, ibrine_rho,     &
                           sice_rho,   bphin,  iphin,  &
                           kperm,      zphi_min,       &
                           igrid,      sss)
       
      end subroutine compute_microS_mushy

!=======================================================================

      subroutine prepare_hbrine (icells,     indxi,  indxj, &
                                 Sin,        zTin,           &
                                 brine_sal,  brine_rho,      &
                                 ibrine_sal, ibrine_rho,     &
                                 sice_rho,   zphin,  iphin,  &
                                 kperm,      zphi_min,       &
                                 igrid,      sss)

      use ice_therm_shared, only: calculate_Tin_from_qin

      integer (kind=int_kind), intent(in) :: &
         icells         ! number of cells with aicen > 0
                              
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj ! compressed indices for icells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), &
         intent(in) :: &
         Sin        , & ! salinity of ice layers on bio grid (ppt)
         zTin           ! temperature of ice layers on bio grid for history (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), &
         intent(inout) :: &
         zphin      , & ! porosity of layers
         brine_sal  , & ! equilibrium brine salinity (ppt)  
         brine_rho      ! internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), &
         intent(inout) :: &
         ibrine_rho , & ! brine density on interface (kg/m^3)
         ibrine_sal , & ! brine salinity on interface (ppt)
         iphin          ! porosity on interface

      real (kind=dbl_kind), dimension (nblyr+1), intent(in):: &
         igrid          ! biology grid interface points

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         sss            ! sea surface salinity (ppt)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         sice_rho   , & ! avg sea ice density (kg/m^3)
         kperm      , & ! harmonic average permeability (m^2)
         zphi_min       ! minimum porosity

      ! local variables

      real (kind=dbl_kind), dimension(icells, nblyr) :: &
          kin           !  permeability (m^2)
    
      real (kind=dbl_kind), dimension(icells) :: &
          k_min, ktemp

      real (kind=dbl_kind) :: &
          igrp, igrm, rigr
     
      integer (kind=int_kind) :: &
           k,m, i, j, ij   ! tracer indice

      !-----------------------------------------------------------------
      ! calculate equilibrium brine density and gradients 
      !-----------------------------------------------------------------

      sice_rho(:,:) = c0
      kperm   (:,:) = c0
      ktemp   (:)   = c0
      
      do k = 1, nblyr+2
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1,icells
            i = indxi(ij)
            j = indxj(ij) 

            brine_sal(i,j,k) = a1*zTin(i,j,k) + a2*zTin(i,j,k)**2 + a3*zTin(i,j,k)**3
            brine_rho(i,j,k) = b1 + b2*brine_sal(i,j,k)
            zphin    (i,j,k) = min(c1, max(puny, &
                               Sin(i,j,k)*rhosi/(brine_sal(i,j,k)*brine_rho(i,j,k)))) 
  
            if (k == nblyr+2) then 
               brine_sal(i,j,nblyr+2) = sss(i,j)  !max(brine_sal(i,j,nblyr+2),Sin(i,j,nblyr+2))
               brine_rho(i,j,nblyr+2) =  rhow  !b1 + b2*brine_sal(i,j,nblyr+2) 
               zphin(i,j,nblyr+2) = c1 !
               ibrine_sal(i,j,1) =brine_sal(i,j,2)! brine_sal(i,j,1)             !
               ibrine_sal(i,j,nblyr+1) =brine_sal(i,j,nblyr+2)! brine_sal(i,j,nblyr+2) !
               ibrine_rho(i,j,1) =brine_rho(i,j,2) !brine_rho(i,j,1)
               ibrine_rho(i,j,nblyr+1) = brine_rho(i,j,nblyr+2) !brine_rho(i,j,nblyr+2) !
               iphin(i,j,1) = zphin(i,j,2)
               iphin(i,j,nblyr+1) = zphin(i,j,nblyr+1)
                
            endif

            if (k > 1 .AND. k < nblyr+2) then

               kin(ij, k-1) = k_o*zphin(i,j,k)**exp_h 
               sice_rho(i,j) = sice_rho(i,j) + (rhoi*(c1-zphin(i,j,k)) + &
               brine_rho(i,j,k)*zphin(i,j,k))*(igrid(k)-igrid(k-1))
                  
            elseif (k == nblyr+2) then

               k_min(ij) = MINVAL(kin(ij,:))
               zphi_min(i,j) = zphin(i,j,2) 
                         
            endif
         enddo ! ij
      enddo      ! k         
    
      do k = 2, nblyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1,icells
            i = indxi(ij)
            j = indxj(ij) 
 
            if (k_min(ij) >  c0) then
               ktemp(ij) = ktemp(ij) + c1/kin(ij,k-1)
               if (k == nblyr) then
                  ktemp(ij)  = ktemp(ij) + c1/kin(ij,nblyr) 
                  kperm(i,j) = real(nblyr,kind=dbl_kind)/ktemp(ij)
               endif
            else
               kperm(i,j) = k_min(ij)
            endif

            igrp = igrid(k+1)-igrid(k  )
            igrm = igrid(k  )-igrid(k-1)
            rigr = c1/(igrid(k+1)-igrid(k-1))

            ibrine_sal(i,j,k) = (brine_sal(i,j,k+1)*igrp &
                              +  brine_sal(i,j,k  )*igrm) * rigr
            ibrine_rho(i,j,k) = (brine_rho(i,j,k+1)*igrp &
                              +  brine_rho(i,j,k  )*igrm) * rigr
            iphin(i,j,k) = min(c1, max(puny, &
                                    (zphin(i,j,k+1)*igrp &
                                  +  zphin(i,j,k  )*igrm) * rigr))
         enddo ! ij
      enddo    ! k         

      end subroutine prepare_hbrine

!=======================================================================

! Changes include brine height increases from ice and snow surface melt, 
! congelation growth, and upward pressure driven flow from snow loading.
!  
! Decreases arise from downward flushing and bottom melt.  
!
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice 
! with dynamic salinity or the height ratio == hinS/vicen*aicen, where 
! hinS is the height of the brine surface relative to the bottom of the 
! ice.  This volume fraction may be > 1 in which case there is brine 
! above the ice surface (ponds).

      subroutine update_hbrine (meltb,      meltt,       &
                                melts,      dt,          &
                                hin,        hsn,         &
                                hin_old,    firstice,    &
                                hinS,       hinS_old,    &
                                fbri,                    &
                                dhS_top,    dhS_bottom,  &
                                kperm,      zphi_min,    &
                                darcy_V)

      use ice_zbgc_shared, only: flood_frac

      real (kind=dbl_kind), intent(in) :: &
         dt  !timestep, tuning parameter
        
      real (kind=dbl_kind), intent(in):: &
         meltb,     & ! bottom melt over dt (m)
         meltt,     & ! true top melt over dt (m)
         melts,     & ! true snow melt over dt (m)
         hin,       & ! ice thickness (m)
         hsn,       & ! snow thickness (m)
         hin_old,   & ! past timestep ice thickness (m)
         hinS_old , & ! previous timestep hinS
         kperm        ! avg ice permeability

      real (kind=dbl_kind), intent(inout):: &
         darcy_V    , & ! Darcy velocity: m/s
         dhS_top    , & ! change in top hinS before darcy flow
         dhS_bottom , & ! change in bottom hinS initially before darcy flow
         hinS       , & ! thickness of brine (m) 
         fbri       , & ! brine height ratio tracer (hinS/hin) 
         zphi_min       ! surface porosity

      logical (kind=log_kind), intent(in) :: & 
         firstice

      ! local variables

      real (kind=dbl_kind) :: &  
         hinS_min  , & ! thinS or hin 
         dhinS_hin , & ! hinS-hin
         hbrn      , & ! brine height  (m) hinS-h_o
         dhinS     , & ! change in brine surface
         h_o       , & ! new ocean surface from ice bottom (m)
         darcy_coeff,& ! magnitude of the Darcy velocity/hbrn (1/s)
         hbrn_new      ! hbrn after flushing

      real (kind=dbl_kind), parameter :: &
         dh_min = 0.001_dbl_kind, &  !echmod WHAT IS THIS?
!echmod USE NAMELIST PARAMETERS rfracmin, rfracmax 
         run_off = c0   !fraction of melt that runs off directly to the ocean

         hbrn     = c0
         darcy_V  = c0
         hbrn_new = c0
         h_o = rhosi/rhow*hin + rhos/rhow*hsn 
       
         if (hinS_old > thinS .AND. hin_old > thinS) then

            dhS_top = -max(c0, min(hin_old-hinS_old, meltt)) * rhoi/rhow
            dhS_top = dhS_top - max(c0, melts) * rhos/rhow
            dhS_top = (c1 - run_off) * dhS_top
            dhinS   = dhS_bottom - dhS_top  
            hinS    = max(his_min, hinS_old + dhinS)
            hbrn    = hinS - h_o
            darcy_coeff = max(c0, kperm*gravit/(viscos*hinS_old))

            if (hbrn > c0 .AND. hinS > thinS ) then   
               hbrn_new = hbrn*exp(-darcy_coeff/zphi_min*dt)
               hinS = max(thinS, h_o + hbrn_new)    
            elseif (hbrn < c0) then
               if (hinS >= hin) zphi_min = phi_snow
               hbrn_new = hbrn*exp(-darcy_coeff/zphi_min*dt)
               hinS = max(his_min, h_o + hbrn_new)
            endif

            hbrn_new = hinS - h_o
            darcy_V  = -SIGN((hbrn-hbrn_new)/dt*zphi_min, hbrn)
            dhS_top  = dhS_top + SIGN((hbrn-hbrn_new), hbrn)

         else    ! very thin brine height 
            hinS_min  = min(thinS, hin)
            hinS = max(hinS_min, hinS_old+dhS_bottom-dhS_top)
            dhinS_hin = hinS - h_o
            if (abs(dhinS_hin) > dh_min) &
               hinS = max(his_min, h_o + SIGN(dh_min,dhinS_hin))
         endif 

         fbri = hinS/hin

      end subroutine update_hbrine

!=======================================================================

      subroutine read_restart_hbrine(filename_spec)

! Reads all values needed for hbrine
! author Elizabeth C. Hunke, LANL

      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init, &
                              istep0
      use ice_domain, only: nblocks
      use ice_flux, only: sss  
      use ice_state
      use ice_exit, only: abort_ice
      use ice_restart, only: lenstr, restart_dir, restart_file, pointer_file, runtype
      use ice_read_write, only: ice_open, ice_read

      character(len=char_len_long), intent(in), optional :: filename_spec

      ! local variables

      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday  , &   ! year, month, day
          nsize, nr0 , nt              

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
         if (n == 0) call abort_ice('hbrine restart: filename discrepancy')
         string1 = trim(filename0(1:n-1))
         string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
         write(filename,'(a,a,a,a)') &
            string1(1:lenstr(string1)), &
            restart_file(1:lenstr(restart_file)),'.hb', &
            string2(1:lenstr(string2))
      endif ! master_task

      call ice_open(nu_restart_hbrine,filename,0)

      if (my_task == master_task) then
        read(nu_restart_hbrine) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
        write(nu_diag,*) 'hbrine Restart read at istep=',istep0,time,time_forc
      endif

      diag = .true.
       
      do n = 1, ncat
          call ice_read(nu_restart_hbrine,0,trcrn(:,:,nt_fbri,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
          call ice_read(nu_restart_hbrine,0,first_ice_real(:,:,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
      enddo
      do iblk = 1, nblocks
         do n = 1,ncat
            do j = 1, ny_block
            do i = 1, nx_block
               if (first_ice_real(i,j,n,iblk) >= c1) then
                  first_ice (i,j,n,iblk) = .true.
               else
                  first_ice (i,j,n,iblk) = .false.
               endif
            enddo
            enddo
         enddo
      enddo

      if (my_task == master_task) close(nu_restart_hbrine)

      end subroutine read_restart_hbrine

!=======================================================================

      subroutine write_restart_hbrine(filename_spec)

! Dumps all values needed for a hbrine restart
! author Elizabeth C. Hunke, LANL

      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_domain, only: nblocks
      use ice_state
      use ice_flux, only: sss  
      use ice_restart, only: lenstr, restart_dir, restart_file, pointer_file, runtype
      use ice_read_write, only: ice_open, ice_write

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
              restart_file(1:lenstr(restart_file)),'.hb.', &
              iyear,'-',month,'-',mday,'-',sec
      endif

      ! begin writing restart data
      call ice_open(nu_dump_hbrine,filename,0)

      if (my_task == master_task) then
        write(nu_dump_hbrine) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
        write(nu_diag,*) 'hbrine Restart written ',istep1,time,time_forc
      endif

      diag = .true.

      do iblk = 1, nblocks
       do n = 1,ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (first_ice(i,j,n,iblk)) then
               first_ice_real(i,j,n,iblk) = c1
            else
               first_ice_real(i,j,n,iblk) = c0
            endif
         enddo
         enddo
       enddo
      enddo

      do n = 1, ncat
         call ice_write(nu_dump_hbrine,0,trcrn(:,:,nt_fbri,n,:),'ruf8',diag)
         call ice_write(nu_dump_hbrine,0,first_ice_real(:,:,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_dump_hbrine)

      end subroutine write_restart_hbrine

!=======================================================================
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW
!          Nicole Jeffery, LANL

      subroutine hbrine_diags (dt)
              
      use ice_broadcast, only: broadcast_scalar
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, pbloc, &
                                plat, plon
      use ice_domain_size, only: ncat, nltrcr
      use ice_state, only: aice, aicen, vicen, vice, trcr, nt_fbri, &
                          trcrn, hbrine, nt_sice
      use ice_zbgc_shared, only: darcy_V

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, iblk

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         phinS, phinS1, pdarcy_V, pfbri

      real (kind=dbl_kind), dimension(npnt,nilyr) :: &
         pSin

      !-----------------------------------------------------------------
      ! Dynamic brine height
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
               phinS1(n) = c0             
               phinS(n) = c0            
               pfbri(n) = trcrn(i,j,nt_fbri,1,iblk) 
               pdarcy_V(n) = darcy_V(i,j,1,iblk)
               if (aice(i,j,iblk) > c0) &
                       phinS (n) = trcr(i,j,nt_fbri,iblk) &
                                 * vice(i,j,iblk)/aice(i,j,iblk)
               if (aicen(i,j,1,iblk)> c0) &
                       phinS1(n) = trcrn(i,j,nt_fbri,1,iblk) &
                                 * vicen(i,j,1,iblk)/aicen(i,j,1,iblk)
               do k = 1,nilyr
                       pSin(n,k) = trcr(i,j,nt_sice+k-1,iblk)
               enddo
            endif                 ! my_task = pmloc
           
            do k = 1,nilyr
            call broadcast_scalar(pSin(n,k), pmloc(n))   
            enddo
            call broadcast_scalar(pfbri(n), pmloc(n))  
            call broadcast_scalar(phinS1(n), pmloc(n))  
            call broadcast_scalar(phinS(n), pmloc(n)) 
            call broadcast_scalar(pdarcy_V(n), pmloc(n))
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
        write(nu_diag,*) '-------- hbrine -------'
        write(nu_diag,900) 'hbrine, (m)            = ',phinS(1),phinS(2)
        write(nu_diag,900) 'fbri, cat1 (m)         = ',pfbri(1),pfbri(2)
        write(nu_diag,900) 'hbrine cat1, (m)       = ',phinS1(1),phinS1(2)  
        write(nu_diag,900) 'darcy_V cat1, (m/s)    = ',pdarcy_V(1),pdarcy_V(2)
        do k = 1, nilyr
        write(nu_diag,900) 'salinity profile (ppt) = ',pSin(1,k),pSin(2,k)
        enddo
      endif                   ! print_points
      endif                   ! my_task = master_task 

  900 format (a25,2x,f24.17,2x,f24.17)

      end subroutine hbrine_diags

!=======================================================================

      end module ice_brine

!=======================================================================
