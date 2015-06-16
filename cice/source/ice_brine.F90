!  SVN:$Id$
!=======================================================================
!
! Computes ice microstructural information for use in biogeochemistry
!
! authors: Nicole Jeffery, LANL
!
      module ice_brine

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: nilyr, nblyr, max_blocks, ncat
      use ice_fileunits, only: nu_diag, nu_rst_pointer, nu_dump_hbrine, &
          nu_restart_hbrine, flush_fileunit
      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_state, only: ntrcr, nt_qice, nt_sice, nt_bgc_S 
      use ice_zbgc_shared, only: cgrid, bgrid, igrid, exp_h, k_o, rhosi, &
           hbr_min, thinS, min_salin, remap_zbgc, &
           phi_snow, restart_hbrine, first_ice, icgrid, swgrid, Ra_c, bphimin 

      implicit none

      private
      public :: preflushing_changes, compute_microS_mushy, &
                update_hbrine, init_hbrine, write_restart_hbrine, &
                hbrine_diags, compute_microS, calculate_drho 
 
      real (kind=dbl_kind), parameter :: &   
         maxhbr  = 1.25_dbl_kind  , & ! brine overflows if hbr > maxhbr*hin
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

      bgrid(:)       = c0 ! bgc grid points         
      bgrid(nblyr+2) = c1 ! bottom value
      igrid(:)       = c0 ! bgc interface grid points   
      igrid(1)       = c0 ! ice top
      igrid(nblyr+1) = c1 ! ice bottom
      
      zspace = c1/max(c1,(real(nblyr,kind=dbl_kind)))
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
      ! Calculate CICE icgrid for ishortwave interpolation top(0) , bottom (1)
      !-----------------------------------------------------------------
       
      icgrid(1) = c0                        
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing
    
      do k = 2, nilyr+1
         icgrid(k) = zspace * (real(k,kind=dbl_kind)-c1)
      enddo 

      !-----------------------------------------------------------------
      ! Calculate CICE swgrid for dEdd ice: top of ice (0) , bottom of ice (1)
      !-----------------------------------------------------------------
      swgrid(1) = c1/60.0_dbl_kind         
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing
      swgrid(2) = zspace/c2 + swgrid(1)
      do k = 3, nilyr+1
         swgrid(k) = zspace * (real(k,kind=dbl_kind)-c1p5)
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
! dynamic salinity or the height ratio == hbr/vicen*aicen, where hbr is the 
! height of the brine surface relative to the bottom of the ice.  This volume fraction
! may be > 1 in which case there is brine above the ice surface (meltponds). 

      subroutine preflushing_changes (nx_block, ny_block,             &
                                      icells,   n_cat,                &
                                      indxii,   indxjj,               &
                                      aicen,    vicen,    vsnon,      &
                                      meltb,    meltt,    congel,     &
                                      snoice,   hice_old, dhice,      &
                                      fbri,     dhbr_top, dhbr_bot,   &
                                      hbr_old,  hin,hsn,  firstice)
 
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
  
      real (kind=dbl_kind), dimension(nx_block*ny_block), intent(inout) :: &
         hin          , & ! ice thickness (m) 
         hsn          , & ! snow thickness (m) 
         hbr_old      , & ! old brine height (m)
         dhice            ! change due to sublimation (<0)/condensation (>0) (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         fbri         , & ! trcrn(i,j,nt_fbri)
         dhbr_top     , & ! brine change in top for diagnostics (m)
         dhbr_bot     , & ! brine change in bottom for diagnostics (m)
         hice_old         ! old ice thickness (m)

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(in) :: &
         firstice         ! if true, initialized values should be used     

      ! local variables

      integer (kind=int_kind) :: &
         i, j         , & ! horizontal indices
         ij               ! horizontal index, combines i and j loops
  
      real (kind=dbl_kind) :: &
         hin_old          ! ice thickness before current melt/growth (m)

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------

      dhbr_top(:,:) = c0
      dhbr_bot(:,:) = c0
      dhice(:) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
      
         i = indxii(ij)
         j = indxjj(ij)
         if (fbri(i,j) <= c0) then
            write(nu_diag, *) 'preflushing: fbri <= c0:',i,j
            write(nu_diag, *) 'vicen, aicen', vicen(i,j), aicen(i,j)
            write(nu_diag, *) 'fbri, hice_old', fbri(i,j), hice_old(i,j)
            call abort_ice ('ice_brine error')
         endif

         hin(ij) = vicen(i,j) / aicen(i,j)
         hsn(ij) = vsnon(i,j) / aicen(i,j)
         hin_old = max(c0, hin(ij) + meltb (i,j) + meltt (i,j) &
                                   - congel(i,j) - snoice(i,j))
         dhice(ij) = hin_old - hice_old(i,j)   ! change due to subl/cond
         dhbr_top(i,j) = meltt (i,j)- snoice(i,j) - dhice(ij) 
         dhbr_bot(i,j) = congel(i,j) - meltb(i,j)

         if ((hice_old(i,j) < puny) .OR. (hin_old < puny) &
                                    .OR. firstice(i,j)) then
            hice_old  (i,j) = hin(ij)
            dhbr_top  (i,j) = c0
            dhbr_bot  (i,j) = c0
            dhice      (ij) = c0
            fbri      (i,j) = c1 
         endif

         hbr_old(ij) = fbri(i,j) * hice_old(i,j)

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
                                       trcrn,    hice_old,   hbr_old,    &
                                       sss,      sst,        bTin,       &
                                       iTin,                             &
                                       bphin,    kperm,      bphi_min,   &
                                       bSin,     brine_sal,  brine_rho,  &
                                       iphin,    ibrine_rho, ibrine_sal, &
                                       sice_rho, iDin)  

      use ice_therm_mushy, only: temperature_mush, liquid_fraction, permeability
      use ice_zbgc_shared, only: Dm, l_sk, viscos_dynamic  

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells , &            ! number of cells with aicen > 0
         n_cat                 ! ice category
                              
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxii, indxjj  ! compressed indices for icells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         hice_old    , & ! previous timestep ice height (m)
         sss         , & ! ocean salinity (ppt)
         sst             ! ocean temperature (C)
       
      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(in) :: &
         trcrn           

      real (kind=dbl_kind), dimension(nx_block*ny_block), intent(out) :: & 
         kperm       , & ! average ice permeability (m^2)
         bphi_min        ! surface porosity

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(inout) :: &
         hbr_old           ! previous timestep brine height (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), &
         intent(inout)  :: &
         iDin           ! tracer diffusivity/h^2 (1/s) includes gravity drainage/molecular

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1), &
         intent(inout)  :: &
         iphin       , & ! porosity on the igrid 
         ibrine_rho  , & ! brine rho on interface  
         ibrine_sal  , & ! brine sal on interface   
         iTin            ! Temperature on the igrid (oC)

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2), &
         intent(inout)  :: &
         bSin        , &    ! bulk salinity (ppt) on bgrid
         brine_sal   , & ! equilibrium brine salinity (ppt) 
         brine_rho       ! internal brine density (kg/m^3) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), intent(inout) :: &
         bTin        , & ! Temperature on bgrid
         bphin           ! porosity on bgrid

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         sice_rho        ! average ice density  

      ! local variables

      real (kind=dbl_kind), dimension (nx_block*ny_block,nilyr) :: &
         cSin        , & ! bulk salinity (ppt)
         cqin            ! enthalpy ()

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2) :: &
         zTin        , & ! Temperature of ice layers on bgrid (C) 
         zSin        , & ! Salinity of ice layers on bgrid (C) 
         bqin            ! enthalpy on the bgrid ()

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1) :: &
         ikin            ! permeability (m^2) 

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k               ! vertical biology layer index 
      
      real (kind=dbl_kind), dimension(icells) :: &
         surface_S   , & ! salinity of ice above hin > hbr 
         hinc_old    , & ! mean ice thickness before current melt/growth (m)
         hbrc_old        ! mean brine thickness before current melt/growth (m)
      
      real (kind=dbl_kind), dimension (icells,ntrcr+2) :: & ! nblyr+2)
         trtmp_s     , & ! temporary, remapped tracers   
         trtmp_q         ! temporary, remapped tracers   
      
      real (kind=dbl_kind), dimension(icells,nblyr+1) :: &   
         drho            ! brine density difference (kg/m^3)
     
    real(kind=dbl_kind), parameter :: &
         Smin = p01
    
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
        
      trtmp_s(:,:) = c0
      trtmp_q(:,:) = c0
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells   ! map Sin and qin (cice) profiles to bgc grid 
         i = indxii(ij)
         j = indxjj(ij)
         surface_S(ij) = min_salin
         hbr_old  (ij) = min(hbr_old(ij), maxhbr*hice_old(i,j))
         hinc_old (ij) = hice_old(i,j)
         hbrc_old (ij) = hbr_old (ij)

         call remap_zbgc(ntrcr,            nilyr,          &
                         nt_sice,                          &
                         trcrn(i,j,:),     trtmp_s(ij,:),  &
                         0,                nblyr,        &
                         hinc_old(ij),     hinc_old(ij),   &
                         cgrid(2:nilyr+1),                 &
                         bgrid(2:nblyr+1), surface_S(ij))
     
         call remap_zbgc(ntrcr,            nilyr,          &
                         nt_qice,                          &
                         trcrn(i,j,:),     trtmp_q(ij,:),  &
                         0,                nblyr,        &
                         hinc_old(ij),     hinc_old(ij),   &
                         cgrid(2:nilyr+1),                 &
                         bgrid(2:nblyr+1), surface_S(ij))
      enddo

      do k = 1, nblyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells   
            i = indxii(ij)
            j = indxjj(ij)
            bqin (ij, k+1) = min(c0,        trtmp_q(ij,nt_qice+k-1))
            bSin (ij, k+1) = max(Smin,      trtmp_s(ij,nt_sice+k-1))
            bTin (i,j,k+1) = temperature_mush(bqin(ij, k+1), bSin(ij,k+1))
            bphin(i,j,k+1) = liquid_fraction (bTin(i,j,k+1), bSin(ij,k+1))
         enddo ! ij                          
      enddo    ! k

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells   
            i = indxii(ij)
            j = indxjj(ij)
            bSin (ij, 1) = bSin(ij,2)
            bTin (i,j,1) = bTin(i,j,2)
            bphin(i,j,1) = bphin(i,j,2)
            bphin(i,j,nblyr+2) = c1
            bSin (ij, nblyr+2) = sss(i,j) 
            bTin (i,j,nblyr+2) = sst(i,j)
            bphin(i,j,nblyr+2) = c1
         enddo ! ij                          

      !-----------------------------------------------------------------
      ! Define ice multiphase structure
      !-----------------------------------------------------------------
     
      call prepare_hbrine (icells,     indxii, indxjj, &
                           bSin,       bTin,  iTin,    &
                           brine_sal,  brine_rho,      &
                           ibrine_sal, ibrine_rho,     &
                           sice_rho,                   &
                           bphin,      iphin,          &
                           kperm,      bphi_min,       &
                           igrid,      sss)
 
      call calculate_drho(icells,nx_block,ny_block,indxii,indxjj,igrid,bgrid,&
                                   brine_rho,ibrine_rho,drho)
      iDin(:,:,:) = c0

      do k= 2, nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxii(ij)
            j = indxjj(ij)  

            ikin(ij,k) = k_o*iphin(ij,k)**exp_h 
            iDin(i,j,k) =  iphin(ij,k)*Dm/hbr_old(ij)**2  
            if (hbr_old(ij) .GE. Ra_c) &
               iDin(i,j,k) =iDin(i,j,k) + l_sk*ikin(ij,k)*gravit/viscos_dynamic* &  
                           drho(ij,k)/hbr_old(ij)**2  
         enddo  !ij
      enddo    !k

 end subroutine compute_microS_mushy

!=======================================================================

      subroutine prepare_hbrine (icells,     indxi,  indxj,  &
                                 bSin,       bTin,    iTin,  &
                                 brine_sal,  brine_rho,      &
                                 ibrine_sal, ibrine_rho,     &
                                 sice_rho,   bphin,  iphin,  & 
                                 kperm,      bphi_min,       &
                                 igrid,      sss)

 
       use ice_therm_shared, only: calculate_Tin_from_qin

      integer (kind=int_kind), intent(in) :: &
         icells         ! number of cells with aicen > 0
                              
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj   ! compressed indices for icells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2), &
         intent(in) :: &
         bSin           ! salinity of ice layers on bio grid (ppt)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), &
         intent(in) :: &
         bTin           ! temperature of ice layers on bio grid for history (C)

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2), &
         intent(inout) :: &
         brine_sal  , & ! equilibrium brine salinity (ppt)  
         brine_rho      ! internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1), &
         intent(inout) :: &
         ibrine_rho , & ! brine density on interface (kg/m^3)
         ibrine_sal , & ! brine salinity on interface (ppt)
         iphin      , & ! porosity on interface
         iTin           ! Temperature on interface   

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), &
         intent(inout) :: &
         bphin          ! porosity of layers

      real (kind=dbl_kind), dimension (nblyr+1), intent(in):: &
         igrid          ! biology grid interface points

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         sss            ! sea surface salinity (ppt)

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(out) :: &
         kperm      , & ! harmonic average permeability (m^2)
         bphi_min       ! minimum porosity

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         sice_rho       ! avg sea ice density 

      ! local variables

      real (kind=dbl_kind), dimension(icells, nblyr+1) :: &
          kin       !  permeability  
    
      real (kind=dbl_kind), dimension(icells):: &
          k_min, ktemp
     
      real (kind=dbl_kind) :: &
          igrp, igrm, rigr  ! grid finite differences

      integer (kind=int_kind) :: &
           k, i, j, ij   ! tracer indice

     !-----------------------------------------------------------------
     !  calculate equilibrium brine density and gradients 
     !-----------------------------------------------------------------
    
      sice_rho(:,:) = c0
      
      do k = 1, nblyr+1
       
           if (k == 1) then 
              igrm = 0
           else
              igrm = igrid(k  ) - igrid(k-1)
           endif
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1,icells
            i = indxi(ij)
            j = indxj(ij) 

            brine_sal(ij, k)   = a1*bTin(i,j,k)    &
                               + a2*bTin(i,j,k)**2 &
                               + a3*bTin(i,j,k)**3
            brine_rho(ij, k)   = b1 + b2*brine_sal(ij,k)
            bphin    (i,j,k)   = min(c1, max(puny, bSin(ij,k)*rhosi &
                               /(brine_sal(ij,k)*brine_rho(ij,k)))) 
            kin      (ij, k)   = k_o*bphin(i,j,k)**exp_h 
            sice_rho(i,j)      = sice_rho(i,j) + (rhoi*(c1-bphin(i,j,k)) + &
                                 brine_rho(ij,k)*bphin(i,j,k))*igrm
         enddo ! ij
      enddo    ! k         

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij) 

            brine_sal (ij, nblyr+2) = sss       (i,j)
            brine_rho (ij, nblyr+2) = rhow
            bphin     (i,j,nblyr+2) = c1
            ibrine_sal(ij, 1)       = brine_sal (ij, 2)
            ibrine_sal(ij, nblyr+1) = brine_sal (ij, nblyr+2)
            ibrine_rho(ij, 1)       = brine_rho (ij, 2)
            ibrine_rho(ij, nblyr+1) = brine_rho (ij, nblyr+2)
            iTin      (ij, 1)       = bTin(i,j,2)
            iTin      (ij, nblyr+1) = bTin(i,j,nblyr+1)
            iphin     (ij, 1)       = bphin     (i,j,2)
            iphin     (ij, nblyr+1) = bphin     (i,j,nblyr+1)
            k_min     (ij)          = MINVAL(kin(ij, 2:nblyr+1))
            kperm     (ij)          = c0  ! initialize
            ktemp     (ij)          = c0
            bphi_min  (ij)          = max(bphin(i,j,1),bSin(ij,2)*rhosi/bphin(i,j,2)/ &  
                                  (brine_sal(ij,1)*brine_rho(ij,1))*phi_snow)
         enddo ! ij

      do k = 2, nblyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij) 
 
            if (k_min(ij) > c0) then
               ktemp(ij) = ktemp(ij) + c1/kin(ij,k)
               kperm(ij) = k_min(ij)
            endif

            igrp = igrid(k+1) - igrid(k  )
            igrm = igrid(k  ) - igrid(k-1)
            rigr = c1 / (igrid(k+1)-igrid(k-1))

            ibrine_sal(ij,k) = (brine_sal(ij,k+1)*igrp &
                             +  brine_sal(ij,k  )*igrm) * rigr
            ibrine_rho(ij,k) = (brine_rho(ij,k+1)*igrp &
                             +  brine_rho(ij,k  )*igrm) * rigr
            iphin     (ij,k) = min(c1, max(puny, &
                                  (bphin(i,j,k+1)*igrp &
                                +  bphin(i,j,k  )*igrm) * rigr))
         enddo ! ij
      enddo    ! k         

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij) 
 
            if (k_min(ij) > c0) then
                  ktemp(ij) = ktemp(ij) + c1/kin(ij,nblyr+1) 
                  kperm(ij) = real(nblyr,kind=dbl_kind)/ktemp(ij)
            endif
         enddo ! ij

      end subroutine prepare_hbrine

!=======================================================================

! Changes include brine height increases from ice and snow surface melt, 
! congelation growth, and upward pressure driven flow from snow loading.
!  
! Decreases arise from downward flushing and bottom melt.  
!
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice 
! with dynamic salinity or the height ratio == hbr/vicen*aicen, where 
! hbr is the height of the brine surface relative to the bottom of the 
! ice.  This volume fraction may be > 1 in which case there is brine 
! above the ice surface (ponds).

      subroutine update_hbrine (meltb,      meltt,       &
                                melts,      dt,          &
                                hin,        hsn,         &
                                hin_old,    hbr,         &
                                hbr_old,                 &
                                fbri,       snoice,      &
                                dhS_top,    dhS_bottom,  &
                                dh_top_chl, dh_bot_chl,  &
                                kperm,      bphi_min,    &
                                darcy_V, darcy_V_chl, bphin)


      real (kind=dbl_kind), intent(in) :: &
         dt             ! timestep
           
      real (kind=dbl_kind), intent(in):: &
         meltb,       & ! bottom melt over dt (m)
         meltt,       & ! true top melt over dt (m)
         melts,       & ! true snow melt over dt (m)
         hin,         & ! ice thickness (m)
         hsn,         & ! snow thickness (m)
         hin_old,     & ! past timestep ice thickness (m)
         hbr_old,     & ! previous timestep hbr
         kperm,       & ! avg ice permeability 
         bphin,       & ! upper brine porosity  
         snoice,      & ! snoice change (m)
         dhS_bottom     ! change in bottom hbr initially before darcy flow

      real (kind=dbl_kind), intent(inout):: &
         darcy_V    , & ! Darcy velocity: m/s
         darcy_V_chl, & ! Darcy velocity: m/s for bgc 
         dhS_top    , & ! change in top hbr before darcy flow
         dh_bot_chl , & ! change in bottom for algae  
         dh_top_chl , & ! change in bottom for algae  
         hbr        , & ! thickness of brine (m) 
         fbri       , & ! brine height ratio tracer (hbr/hin) 
         bphi_min       ! surface porosity   

      ! local variables

      real (kind=dbl_kind) :: &  
         hbr_min    , & ! thinS or hin 
         dhbr_hin   , & ! hbr-hin
         hbrocn     , & ! brine height above sea level (m) hbr-h_ocn
         dhbr       , & ! change in brine surface
         h_ocn      , & ! new ocean surface from ice bottom (m)
         darcy_coeff, & ! magnitude of the Darcy velocity/hbrocn (1/s)
         hbrocn_new     ! hbrocn after flushing

      real (kind=dbl_kind), parameter :: &
         dh_min = p001  ! brine remains within dh_min of sea level
                        ! when ice thickness is less than thinS
 
         hbrocn      = c0
         darcy_V     = c0
         darcy_V_chl = c0  
         hbrocn_new  = c0
         h_ocn = rhosi/rhow*hin + rhos/rhow*hsn 
         
         if (hbr_old > thinS .AND. hin_old > thinS .AND. hin > thinS ) then
            hbr_min = thinS
            dhS_top = -snoice -max(c0, min(hin_old-hbr_old, meltt)) * rhoi/rhow 
            dhS_top = dhS_top - max(c0, melts) * rhos/rhow
            dh_top_chl = dhS_top
            dhbr    = dhS_bottom - dhS_top  
            hbr     = max(puny, hbr_old+dhbr) !max(hbr_min, hbr_old + dhbr) 
            hbrocn  = hbr - h_ocn
            darcy_coeff = max(c0, kperm*gravit/(viscos*hbr_old))

            if (hbrocn > c0 .AND. hbr > thinS ) then 
               bphi_min   = bphin  
               hbrocn_new = hbrocn*exp(-darcy_coeff/bphi_min*dt)
               hbr = max(hbr_min, h_ocn + hbrocn_new)    
            elseif (hbrocn < c0 .AND. hbr > thinS) then
               if (hbr < hin .OR. hsn <= c0 .OR. phi_snow > c1) bphi_min = bphin
               hbrocn_new = hbrocn*exp(-darcy_coeff/bphi_min*dt)
               hbr = max(hbr_min, h_ocn + hbrocn_new)
            endif

            hbrocn_new = hbr - h_ocn
            darcy_V    = -SIGN((hbrocn-hbrocn_new)/dt*bphi_min, hbrocn)
            darcy_V_chl= darcy_V 
            dhS_top    = dhS_top - darcy_V*dt/bphi_min 
            dh_top_chl = dh_top_chl - darcy_V_chl*dt/bphi_min
            dh_bot_chl = dhS_bottom 
  
         else    ! very thin brine height 
            hbr_min  = min(thinS, hin)
            hbr = max(hbr_min, hbr_old+dhS_bottom-dhS_top)
            dhbr_hin = hbr - h_ocn
            if (abs(dhbr_hin) > dh_min) &
               hbr = max(hbr_min, h_ocn + SIGN(dh_min,dhbr_hin)) 
            dhS_top = hbr_old - hbr + dhS_bottom
            dh_top_chl = dhS_top
            dh_bot_chl = dhS_bottom
         endif 
        
         fbri = hbr/hin 

      end subroutine update_hbrine

!=======================================================================

      subroutine read_restart_hbrine()

! Reads all values needed for hbrine
! author Elizabeth C. Hunke, LANL

      use ice_communicate, only: my_task, master_task
      use ice_domain, only: nblocks
      use ice_fileunits, only: nu_diag, nu_restart_hbrine
      use ice_state, only: trcrn, nt_fbri
      use ice_restart,only: read_restart_field

      ! local variables

      integer (kind=int_kind) :: &
         i, j, n, iblk          ! counting indices

      logical (kind=log_kind) :: &
         diag

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) 'brine restart'

      call read_restart_field(nu_restart_hbrine,0,trcrn(:,:,nt_fbri,:,:),'ruf8', &
                              'fbrn',ncat,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_hbrine,0,first_ice_real(:,:,:,:),'ruf8', &
                              'first_ice',ncat,diag,field_loc_center,field_type_scalar)
       
      do iblk = 1, nblocks
         do n = 1,ncat
            do j = 1, ny_block
            do i = 1, nx_block
               if (first_ice_real(i,j,n,iblk) >= p5) then
                  first_ice (i,j,n,iblk) = .true.
               else
                  first_ice (i,j,n,iblk) = .false.
               endif
            enddo
            enddo
         enddo
      enddo

      end subroutine read_restart_hbrine

!=======================================================================

      subroutine write_restart_hbrine()

! Dumps all values needed for a hbrine restart
! author Elizabeth C. Hunke, LANL

      use ice_domain, only: nblocks
      use ice_state, only: trcrn, nt_fbri
      use ice_restart,only: write_restart_field

      ! local variables

      integer (kind=int_kind) :: &
         i, j, n, iblk

      logical (kind=log_kind) :: diag

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

      call write_restart_field(nu_dump_hbrine,0,trcrn(:,:,nt_fbri,:,:),'ruf8', &
                               'fbrn',ncat,diag)
      call write_restart_field(nu_dump_hbrine,0,first_ice_real(:,:,:,:),'ruf8', &
                               'first_ice',ncat,diag)

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
      use ice_state, only: aice, aicen, vicen, vice, trcr, nt_sice , nt_fbri, &
                          trcrn
      use ice_zbgc_shared, only: darcy_V
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, iblk

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         phinS, phinS1, pdarcy_V, pfbri

      real (kind=dbl_kind), dimension(npnt,nilyr) :: &
         pSin, pSin1

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
                       phinS(n) = trcr(i,j,nt_fbri,iblk)*vice(i,j,iblk)/aice(i,j,iblk)
               if (aicen(i,j,1,iblk)> c0) &
                       phinS1(n) = trcrn(i,j,nt_fbri,1,iblk)*vicen(i,j,1,iblk)/&
                                                aicen(i,j,1,iblk)                    
               do k = 1,nilyr
                  pSin1(n,k) = trcrn(i,j,nt_sice+k-1,1,iblk)
                  pSin(n,k)  = trcr(i,j,nt_sice+k-1,iblk)
               enddo
            endif                 ! my_task = pmloc
           
            do k = 1,nilyr
                 call broadcast_scalar(pSin(n,k), pmloc(n))   
                 call broadcast_scalar(pSin1(n,k), pmloc(n))   
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
          write(nu_diag,*) '------ hbrine ------'
          write(nu_diag,900) 'hbrine, (m)        = ',phinS(1),phinS(2)
          write(nu_diag,900) 'fbri, cat1 (m)     = ',pfbri(1),pfbri(2)
          write(nu_diag,900) 'hbrine cat1, (m)   = ',phinS1(1),phinS1(2)  
          write(nu_diag,900) 'darcy_V cat1, (m/s)= ',pdarcy_V(1),pdarcy_V(2)             
          write(nu_diag,*) '                         '
          write(nu_diag,*) '------ Thermosaline Salinity ------'
          write(nu_diag,803) 'Sice1(1) cat1 S (ppt)','Sice1(2) cat1 S'
          write(nu_diag,*) '---------------------------------------------------'
          write(nu_diag,802) ((pSin1(n,k),n=1,2), k = 1,nilyr)              
          write(nu_diag,*) '                         '
          write(nu_diag,803) 'Sice(1) bulk S (ppt) ','Sice(2) bulk S'
          write(nu_diag,*) '---------------------------------------------------'
          write(nu_diag,802) ((pSin(n,k),n=1,2), k = 1,nilyr)              
          write(nu_diag,*) '                         '
      endif                   ! print_points
      endif                   ! my_task = master_task 

  802 format (f24.17,2x,f24.17)
  803 format (a25,2x,a25)
  900 format (a25,2x,f24.17,2x,f24.17)

      end subroutine hbrine_diags

!=======================================================================
!
! Computes ice microstructural properties for zbgc and zsalinity 
!
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice with 
! dynamic salinity or the height ratio == hbr/vicen*aicen, where hbr is the 
! height of the brine surface relative to the bottom of the ice.  
! This volume fraction
! may be > 1 in which case there is brine above the ice surface (meltponds). 
! 
      subroutine compute_microS    (nx_block,  ny_block,             &
                                   icells,     n_cat,                &
                                   indxii,     indxjj,               &
                                   trcrn,      hice_old,             &
                                   hbr_old,    sss,      sst,        &
                                   bTin,       iTin,     bphin,      &
                                   kperm,  bphi_min,                 &
                                   Rayleigh_criteria,    firstice,   &
                                   bSin,                 brine_sal,  &
                                   brine_rho,  iphin,    ibrine_rho, &
                                   ibrine_sal, sice_rho, sloss,      &
                                   salinz)
 
      use ice_therm_shared, only: solve_zsal, calculate_Tin_from_qin
      use ice_calendar, only: istep1, time
      use ice_state, only: nt_fbri, nt_Tsfc

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells , &            ! number of true cells with aicen > 0
         n_cat                 ! ice category
                              
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxii, indxjj ! compressed indices for icells with aicen > puny


      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         hice_old    , & ! previous timestep ice height (m)
         sss         , & ! ocean salinity (ppt)
         sst             ! ocean temperature (oC)
 
      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn           

      real (kind=dbl_kind), dimension (nx_block*ny_block), &
         intent(inout) :: &
         hbr_old      ! old brine height

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), &
         intent(inout) :: &
         bTin          ! Temperature of ice layers on bio grid for history file (^oC) 

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2), &
         intent(out) :: &
         bSin

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1), &
         intent(out) :: &
         iTin          ! Temperature on the interface grid

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(in) :: &
         salinz           ! initial salinity profile for new ice (on cice grid)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), intent(inout) :: &
         bphin            ! Porosity of layers

      real (kind=dbl_kind), dimension(nx_block*ny_block), &
         intent(out) :: & 
         bphi_min    , & ! surface porosity
         kperm        , & ! average ice permeability (m^2)
         sloss           ! (g/m^2) salt from brine runoff for hbr > maxhbr*hin

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         sice_rho           ! average ice density

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         Rayleigh_criteria   ! .true. if ice exceeded a minimum thickness hin >= Ra_c 
        
      logical (kind=log_kind), dimension (nx_block,ny_block), & 
         intent(in) :: &
         firstice            ! .true. if ice is newly formed

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2), intent(inout)  :: &
         brine_sal       ,& ! equilibrium brine salinity (ppt) 
         brine_rho          ! Internal brine density (kg/m^3) 

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1), intent(inout)  :: &
         iphin        , & !  porosity on the igrid 
         ibrine_rho    ,& ! brine rho on interface  
         ibrine_sal       !brine sal on interface   

      ! local variables
 
      integer (kind=int_kind) :: &
         qcells, pcells, tcells    ! number of cells with first ice (qcells)

      integer (kind=int_kind), dimension(icells) :: &
         qndxi, qndxj , & ! compressed indices 
         pndxi, pndxj , & ! compressed indices 
         tndxi, tndxj,  &
         qndxij, pndxij, tndxij

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k           , & ! vertical biology layer index 
         pij, qij, tij   ! compressed indices
      
      real (kind=dbl_kind), dimension(icells) :: &
         surface_S   , & ! salinity of ice above hin > hbr 
         hinc_old    , & ! ice thickness (cell quantity) before current melt/growth (m)
         hbrc_old   , & ! brine thickness(cell quantity) before current melt/growth (m)
         h_o             ! freeboard height (m)

      logical (kind=log_kind), dimension (icells):: &
         Rayleigh        ! .true. if ice exceeded a minimum thickness hin >= Ra_c 

     real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: &
         trtmp0       , & ! temporary, remapped tracers  
         trtmp            ! temporary, remapped tracers   
     
      real (kind=dbl_kind) :: &
         Tmlts          ! melting temperature

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

       qcells = 0
       pcells = 0
       tcells = 0 
       sloss(:) = c0  
       bTin(:,:,:) = c0
       bSin(:,:) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells
          i = indxii(ij)
          j = indxjj(ij)
          hinc_old(ij) = hice_old(i,j)
          !--------------------------------------
          ! Rayleigh condition for salinity and bgc:
          !
          !  Implemented as a minimum thickness criteria
          !  for category 1 ice only.
          !  When hin >= Ra_c (m),  pressure flow
          !   is allowed. 
          ! 
          !  Turn off by putting Ra_c = 0 
          !   in ice_in namelist
          !--------------------------------------         
          Rayleigh(ij) = .true.
          if (n_cat == 1 .AND. hbr_old(ij) < Ra_c) then
             Rayleigh(ij) = Rayleigh_criteria(i,j) ! only category 1 ice can be false 
          endif
                     
          if (firstice(i,j)) then
             qcells = qcells + 1
             qndxi(qcells) = i
             qndxj(qcells) = j
             qndxij(qcells) = ij 
          elseif (hbr_old(ij) > maxhbr*hice_old(i,j)) then
             pcells = pcells + 1
             pndxi(pcells) = i
             pndxj(pcells) = j
             pndxij(pcells) = ij 
          else
             tcells = tcells + 1
             tndxi(tcells) = i
             tndxj(tcells) = j
             tndxij(tcells) = ij 
          endif         
                
       enddo  !ij

      !----------------------------------------------------------------
      ! Define ice salinity on Sin
      !---------------------------------------------------------------     

      trtmp(:,:,:) = c0    
      surface_S(:) = min_salin   
      if (solve_zsal) then  
         
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, qcells   !first_ice:  Map initial profile cice to bgc grid 
         i = qndxi(ij)
         j = qndxj(ij) 
         qij =qndxij(ij)
         do k = 1, nilyr
              trcrn(i,j,nt_sice+k-1) = salinz(i,j,k)
         enddo

         call remap_zbgc  (ntrcr,            nilyr,         &
                           nt_sice,                         &
                           trcrn(i,j,:),     trtmp(i,j,:),  &
                           0,                nblyr,         &
                           hinc_old(qij),    hinc_old(qij), &
                           cgrid(2:nilyr+1),                &
                           bgrid(2:nblyr+1), surface_S(qij))     
      enddo
        
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu   
      do ij = 1, pcells   !brine overflow: Salinity flux stored in sloss
         i = pndxi(ij)    ! map bgc to bgc grid
         j = pndxj(ij) 
         pij =pndxij(ij)
         call remap_zbgc (ntrcr,            nblyr,          &
                          nt_bgc_S,                        &
                          trcrn(i,j,:),     trtmp(i,j,:),  &
                          0,                nblyr,         &
                          hbr_old(pij),                   &
                          maxhbr*hinc_old(pij),            &
                          bgrid(2:nblyr+1),                &
                          bgrid(2:nblyr+1), surface_S(pij))
      
      enddo
 
      do k = 1, nblyr    
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu     
         do ij = 1, tcells    
            i = tndxi(ij)
            j = tndxj(ij) 
            tij = tndxij(ij)
            bSin(tij,k+1) = max(min_salin,trcrn(i,j,nt_bgc_S+k-1))  
            bSin(tij,1) = bSin(tij,2)
            bSin(tij,nblyr+2) =  sss(i,j)
            if (trcrn(i,j,nt_bgc_S+k-1) < min_salin-puny) then
                write(nu_diag, *) 'Bad value compute_microS, tcells'
                write(nu_diag, *) 'trcrn(i,j,nt_bgc_S+k-1),i,j,k',trcrn(i,j,nt_bgc_S+k-1),i,j,k
                write(nu_diag, *) 'hbr_old(tij),i,j,tij',hbr_old(tij),i,j,tij
                write(nu_diag, *) 'hinc_old(tij),i,j,tij',hinc_old(tij),i,j,tij
                write(nu_diag, *) 'trcrn(i,j,nt_fbri),ntrcr,nt_fbri:',trcrn(i,j,nt_fbri),ntrcr,nt_fbri
                write(nu_diag, *) 'istep1, time, n_cat',istep1,time,n_cat
                call abort_ice ('ice: ice_brine.F90 error')
            endif
         enddo ! ij

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, qcells   !first_ice(i,j)
            i = qndxi(ij)
            j = qndxj(ij) 
            qij =qndxij(ij)
            trcrn(i,j,nt_bgc_S+k-1) = max(min_salin,trtmp(i,j,nt_sice+k-1)) 
            bSin(qij,k+1) = max(min_salin,trcrn(i,j,nt_bgc_S+k-1))   
            bSin(qij,1) = bSin(qij,2) 
            bSin(qij,nblyr+2) =  sss(i,j) 
            if (trcrn(i,j,nt_bgc_S+k-1) < min_salin-puny) then
                write(nu_diag, *) 'Bad value compute_microS, qcells'
                write(nu_diag, *) 'trcrn(i,j,nt_bgc_S+k-1),i,j,k',trcrn(i,j,nt_bgc_S+k-1),i,j,k
                write(nu_diag, *) 'hbr_old(qij),i,j,qij',hbr_old(qij),i,j,qij
                write(nu_diag, *) 'hinc_old(qij),i,j,qij',hinc_old(qij),i,j,qij
                write(nu_diag, *) 'istep1, time, n_cat',istep1,time,n_cat
                call abort_ice ('ice:  ice_brine.F90 error')
            endif    
         enddo ! ij

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                            
         do ij = 1, pcells  
            i = pndxi(ij)
            j = pndxj(ij) 
            pij =pndxij(ij)
            bSin(pij,k+1) = max(min_salin,trtmp(i,j,nt_bgc_S+k-1)) 
            sloss(pij) = sloss(pij) + rhosi*(hbr_old(pij)*trcrn(i,j,nt_bgc_S+k-1) - &
                         maxhbr*hice_old(i,j)*bSin(pij,k+1))*(igrid(k+1)-igrid(k))
            trcrn(i,j,nt_bgc_S+k-1) = bSin(pij,k+1)                               
            bSin(pij,1) = bSin(pij,2)
            if (k == nblyr) then
               bSin(pij,nblyr+2) =  sss(i,j)
               hbr_old(pij) = maxhbr*hinc_old(pij)
            endif
            if (trcrn(i,j,nt_bgc_S+k-1) < min_salin-puny) then
               write(nu_diag, *) 'Bad value compute_microS, pcells'
               write(nu_diag, *) 'trcrn(i,j,nt_bgc_S+k-1),i,j,k',trcrn(i,j,nt_bgc_S+k-1),i,j,k
               write(nu_diag, *) 'hbr_old(pij),i,j,pij',hbr_old(pij),i,j,pij
               write(nu_diag, *) 'istep1, time, n_cat',istep1,time,n_cat
            endif
         enddo !ij
      enddo ! k
      
    else   ! not solve_zsal  

      do k = 1, nblyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells    
            i = indxii(ij)
            j = indxjj(ij) 
            bSin(ij,k+1) = trcrn(i,j,nt_bgc_S+k-1) !read from bio grid
            bSin(ij,1) =  bSin(ij,2) 
            bSin(ij,nblyr+2) =  sss(i,j) 
         enddo  !ij
      enddo     !k

    endif   !solve_zsal

    !-----------------------------------------------------------------
    ! sea ice temperature for bio grid
    !----------------------------------------------------------------- 
      
    do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells    
          i = indxii(ij)
          j = indxjj(ij)
          Tmlts = -trcrn(i,j,nt_sice+k-1)*depressT
          trtmp0(i,j,nt_qice+k-1) = calculate_Tin_from_qin(trcrn(i,j,nt_qice+k-1),Tmlts)
       enddo  !ij
    enddo   ! k
         
    trtmp(:,:,:) = c0                
      
    !CICE to Bio: remap temperatures

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
    do ij = 1, icells    
          i = indxii(ij)
          j = indxjj(ij)
          call remap_zbgc (ntrcr,                       &
                           nilyr,               nt_qice,     &
                           trtmp0(i,j,1:ntrcr), trtmp(i,j,:),&
                           0,                   nblyr,       &
                           hinc_old(ij),        hbr_old(ij), &
                           cgrid(2:nilyr+1),                 & 
                           bgrid(2:nblyr+1),    surface_S(ij))
    enddo   !ij

    do k = 1, nblyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
       do ij = 1, icells
          i = indxii(ij)
          j = indxjj(ij)
          Tmlts = -bSin(ij,k+1) * depressT
          bTin(i,j,k+1) = min(Tmlts,trtmp(i,j,nt_qice+k-1))
          Tmlts =  -min_salin* depressT 
          bTin(i,j,1) = min(Tmlts,(bTin(i,j,2) + trcrn(i,j,nt_Tsfc))*p5) 
          Tmlts = -bSin(ij,nblyr+2)* depressT  
          bTin(i,j,nblyr+2) = sst(i,j)
       enddo   !ij
    enddo !k
   
    !------------------------------------------------------------------
    ! Define ice multiphase structure
    !----------------------------------------------------------------
     
     call prepare_hbrine (icells,     indxii, indxjj,  &
                           bSin,       bTin, iTin,     &
                           brine_sal,  brine_rho,      &
                           ibrine_sal, ibrine_rho,     &
                           sice_rho,   bphin,   iphin, &
                           kperm,      bphi_min,       &
                           igrid,      sss)
       
    end subroutine compute_microS

!==========================================================================================
!
! Find density difference about interface grid points
! for gravity drainage parameterization
!  

      subroutine calculate_drho &
                                  (icells,nx_block,ny_block,indxi,indxj,igrid,bgrid,&
                                   brine_rho,ibrine_rho,drho)

       integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells       

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj  ! compressed indices 

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid            ! biology nondimensional grid layer points 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid         ! biology grid interface points 

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2), intent(in) :: &
         brine_rho           ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr + 1), intent(in) :: &
         ibrine_rho           ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (icells, nblyr+1), intent(out) :: & 
         drho                 ! brine difference about grid point (kg/m^3)

! local variables

    integer (kind=int_kind) :: &
         k, i,j, ij,m, mm ! indices

    integer (kind=int_kind), dimension (icells) :: &
         mstop, mstart

     real (kind=dbl_kind), dimension (icells, nblyr+2) :: &  !on the zbgc vertical grid
         diff              ! difference array between grid points and 
                           ! the density averaging boundary 

     real (kind=dbl_kind), dimension (icells, nblyr + 1) :: &  !on the zbgc vertical grid
         rho_a         ,&  ! average brine density  above grid point (kg/m^3)
         rho_2a            ! average brine density  above and below grid points (kg/m^3)

     real (kind=dbl_kind), dimension (icells, nblyr + 1) :: &  !on the zbgc vertical grid
         rho_b         ,&  ! brine density  above grid point (kg/m^3)
         rho_2b            ! brine density  above and below grid points (kg/m^3)

       rho_a(:,:) = c0
       rho_2a(:,:) = c0
       rho_b(:,:) = c0
       rho_2b(:,:) = c0
       drho(:,:) = c0  !surface is snow or atmosphere 

       do k = 1, nblyr+1   !igrid values
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij) 

         !----------------------------------------------
         ! h_avg(k,ij) = igrid(k)                        
         ! Calculate rho_a(k), ie  average rho above igrid(k) 
         ! first part is good
         !----------------------------------------------

            if (k == 2) then
               rho_a(ij,2) = (brine_rho(ij,2)*bgrid(2) + (ibrine_rho(ij,2) + brine_rho(ij,2))*&
                         p5*(igrid(2)-bgrid(2)) )/igrid(2)
               rho_b(ij,2) = brine_rho(ij,2)

            elseif (k > 2 .AND. k < nblyr+1) then 
               rho_a(ij,k) = (rho_a(ij,k-1)*igrid(k-1) +  (ibrine_rho(ij,k-1) + brine_rho(ij,k))*&
                         p5*(bgrid(k)-igrid(k-1)) + (ibrine_rho(ij,k) + brine_rho(ij,k))*p5*&
                         (igrid(k)-bgrid(k)))/igrid(k)
               rho_b(ij,k) = brine_rho(ij,k)
            else
               rho_a(ij,nblyr+1) = (rho_a(ij,nblyr)*igrid(nblyr) + (ibrine_rho(ij,nblyr) + &
                        brine_rho(ij,nblyr+1))*p5*(bgrid(nblyr+1)-igrid(nblyr)) + &
                        brine_rho(ij,nblyr+1)*(igrid(nblyr+1)-bgrid(nblyr+1)))/igrid(nblyr+1)
               rho_a(ij,1) = brine_rho(ij,2)   !for k == 1 use grid point value
               rho_b(ij,nblyr+1) = brine_rho(ij,nblyr+1)
               rho_b(ij,1) =  brine_rho(ij,2)
            endif
           enddo  !ij
       enddo     !k

       !----------------------------------------------
       ! Calculate average above and below k rho_2a
       !----------------------------------------------

       do k = 1, nblyr+1   !igrid values
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij) 
            if (k == 1) then
               rho_2a(ij,1) = (rho_a(ij,1)*bgrid(2) + p5*(brine_rho(ij,2) + ibrine_rho(ij,2))*&
                      (igrid(2)-bgrid(2)))/igrid(2)
               rho_2b(ij,1) = brine_rho(ij,2)
            else
               mstop(ij) = 2*(k-1) + 1
               if (mstop(ij) < nblyr+1) then
                  rho_2a(ij,k) = rho_a(ij,mstop(ij))
                  mstart(ij) = 2;
                  mstop(ij) = 1;
               else
                  mstart(ij) = nblyr+2
                  mstop(ij) = nblyr+3
               endif                     

               do mm = mstart(ij),mstop(ij)
                  rho_2a(ij,k) =(rho_a(ij,nblyr+1) + rhow*(c2*igrid(k)-c1))*p5/igrid(k)
               enddo
               rho_2b(ij,k) = brine_rho(ij,k+1)
            endif
           drho(ij,k) = max(rho_b(ij,k) - rho_2b(ij,k),max(c0,c2*(rho_a(ij,k)-rho_2a(ij,k)), &
              c2*(brine_rho(ij,k)-brine_rho(ij,k+1))/real(nblyr,kind=dbl_kind)))
         enddo
       enddo

     end subroutine calculate_drho

!=======================================================================

      end module ice_brine

!=======================================================================
