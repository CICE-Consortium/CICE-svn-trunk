!=======================================================================
!
!BOP
!
! !MODULE: ice_brine - brine height motions and microstructure
!
! !DESCRIPTION:
!
! Computes ice microstructural information for use in biogeochemistry and salinity
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors: Nicole Jeffery, LANL
!
!
!
!
! !INTERFACE:
!
      module ice_brine
!
! !USES:
!   
      use ice_kinds_mod
      use ice_constants
      use ice_domain_size, only: nilyr, nblyr, nblyr_hist, max_blocks, ncat
      use ice_fileunits, only: nu_diag, nu_rst_pointer, nu_dump_hbrine, &
          nu_restart_hbrine, flush_fileunit
      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_state, only: ntrcr, nt_qice, nt_sice
      use ice_zbgc_shared, only: cgrid, bgrid, igrid, exp_h, k_o, rhosi, his_min,&
           thinS, min_salin, igrid, remap_layers_bgc_plus_xy, &
           phi_snow, restart_hbrine, first_ice
!
!EOP
!
      implicit none

      private
      public :: preflushing_changes, compute_microS_mushy, &
                update_hbrine, init_hbrine, write_restart_hbrine, &
                hbrine_diags, calculate_drho
 
      real (kind=dbl_kind), parameter :: &   
         maxhinS = 1.25_dbl_kind, & ! brine overflow condition if hinS > maxhinS*hin
         viscos     = 2.1e-6_dbl_kind , & ! kinematic viscosity (m^2/s) 
         a1 = -21.4_dbl_kind, &    !(psu/C)  Brine salinity as a cubic function of T 
         a2 = -0.886_dbl_kind, &   !(psu/C^2)
         a3 = -0.012_dbl_kind, &   !(psu/C^3)
                                   ! cox and Weeks
         b1 = 1000.0_dbl_kind, &   !rhofresh (kg/m^3)  Brine density as a quadratic of brine salinity
         b2 = 0.8_dbl_kind         !(kg/m^3/ppt)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat,max_blocks) :: &
         first_ice_real     ! .true. = c1, .false. = c0

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_hbrine
!
! !DESCRIPTION:
!
!  Initialize brine height tracer
!  
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_hbrine 
!
! !USES:
!
   use ice_state, only: nt_fbri, trcrn
!
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!
!EOP
!
      integer (kind=int_kind) :: &
           k               ! vertical index

      real (kind=dbl_kind) :: & 
         zspace                 !grid spacing for CICE vertical grid

      !-----------------------------------------------------------------
      ! Calculate bio gridn: ice top to ice bottom corresponds to 0 to 1
      !-----------------------------------------------------------------

      bgrid(:) = c0     !bgc grid points         
      bgrid(nblyr_hist) = c1 ! bottom value
      igrid(:) = c0     !bgc interface grid points   
      igrid(1) = c0                 ! ice top
      igrid(nblyr+1) = c1           ! ice bottom
      
      zspace = c1/(real(nblyr,kind=dbl_kind)) 
      do k = 2, nblyr+1
          bgrid(k) = zspace*(real(k,kind=dbl_kind)-c1p5)
      enddo
      
      do k = 2, nblyr
        igrid(k) = p5*(bgrid(k+1)+bgrid(k))
      enddo

    !-----------------------------------------------------------------
    ! Calculate cgrid of CICE for interpolation ice top (0) to ice bottom (1) 
    !-----------------------------------------------------------------
       
      cgrid(1) = c0                           !CICE vertical grid top point
      zspace = c1/(real(nilyr,kind=dbl_kind)) !CICE  grid spacing
    
      do k = 2, nilyr+1
        cgrid(k) = zspace * (real(k,kind=dbl_kind)-c1p5) 
      enddo 

      if (restart_hbrine) then
          call read_restart_hbrine
      else
          first_ice(:,:,:,:) = .true.            
          trcrn(:,:,nt_fbri,:,:) = c1
      endif

     1060    format (a30,2x,2D13.2)   ! dbl precision

      end subroutine init_hbrine

!=======================================================================
!BOP
!
! !ROUTINE: preflushing_changes
!
! !DESCRIPTION:
!
! Computes the top and bottom brine boundary changes for flushing
! works for zsalinity and tr_salinity
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
      subroutine preflushing_changes  (nx_block, ny_block,         &
                                   icells, n_cat, indxii,indxjj,   &    
                                   aicen,  vicen,                  &
                                   vsnon,  meltb, meltt, congel,&
                                   snoice, hice_old,        &   
                                   fbri,  dh_top, dh_bot, dh_bot_chl, &
                                   dhi_top, dhi_bot,               &
                                   hinS_o,hin,hsn,                 &
                                   firstice)                          
! !USES:
!   
!
! !INPUT/OUTPUT PARAMETERS:                                
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of true cells with aicen > 0
         n_cat                 ! category
                              
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxii, indxjj ! compressed indices for icells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon       , & ! volume per unit area of snow         (m)
         meltb       , & ! bottom ice melt                      (m)
         meltt       , & ! top ice melt                         (m)
         congel      , & ! bottom ice growth                    (m)
         snoice          ! top ice growth from flooding         (m)
 
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         fbri           !trcrn(i,j,nt_fbri)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         dh_bot_chl     ! for chlorophyll (algae may not flush with meltwater)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         dh_top       , &   !  brine change in top and bottom for diagnostics (m)
         dh_bot       , &    !
         dhi_top      , &  !  ice change in top and bottom for diagnostics (m)
         dhi_bot      , &   !
         hin          , & ! ice thickness (m) 
         hsn          , & ! snow thickness (m) 
         hice_old     , &   !  old ice height
         hinS_o           !  old brine heights

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         firstice      ! .true. initialized values should be used     
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          ! horizontal index, combines i and j loops
  
      real (kind=dbl_kind) :: &
         hin_old      ! ice thickness before current melt/growth (m)

      real (kind=dbl_kind):: &
         dhice          ! Change in hin due to subl/cond  (m)

      !--------------------                           
      ! initialize
      !-------------------

       dh_top(:,:) = c0
       dh_bot(:,:) = c0
       dh_bot_chl(:,:) = c0
       dhi_top(:,:) = c0
       dhi_bot(:,:) = c0
           hin(:,:) = c0
           hsn(:,:) = c0
       hinS_o(:,:) = c0       

       do ij = 1, icells
      
           i = indxii(ij)
           j = indxjj(ij)
 
           hin(i,j) = max(c0,vicen(i,j) / (aicen(i,j)))
           if (hin(i,j) <= c0  .OR. fbri(i,j) <= c0) then
               write(nu_diag, *) 'preflushing: hin <= 0 or fbri(i,j) <= c0 at  i, j:',  i,j
               write(nu_diag, *) 'vicen(i,j), aicen(i,j)'
               write(nu_diag, *)  vicen(i,j), aicen(i,j)
               write(nu_diag, *) 'fbri(i,j),hice_old(i,j)'
               write(nu_diag, *) fbri(i,j), hice_old(i,j)
               call abort_ice ('ice:ice_brine error')
           endif
           hsn(i,j) = c0
           if (hin(i,j) > c0) hsn(i,j) = max(c0,vsnon(i,j) / (aicen(i,j)))
           hin_old = max(c0,hin(i,j) + meltb(i,j) + meltt(i,j) - &
                          congel(i,j)- snoice(i,j))
           dhice = hin_old - hice_old(i,j)   !change due to subl/cond
           dhi_top(i,j) = meltt(i,j)-dhice -snoice(i,j)
           dhi_bot(i,j) = congel(i,j) - meltb(i,j)   
           dh_top(i,j) = dhi_top(i,j)
           dh_bot(i,j) = dhi_bot(i,j)
           dh_bot_chl(i,j) = congel(i,j) - meltb(i,j)

           if ((hice_old(i,j) < puny) .OR. (hin_old < puny) .OR. firstice(i,j)) then

             hin_old = hin(i,j) 
             dh_top(i,j) = c0
             dh_bot(i,j) = c0
             dhi_bot(i,j) = c0
             dhi_top(i,j) = c0
             dh_bot_chl(i,j) = c0
             fbri(i,j) = c1 

           endif

           hinS_o(i,j) = fbri(i,j)* hice_old(i,j)

         enddo  !ij


 end subroutine preflushing_changes

!=======================================================================
!BOP
!
! !ROUTINE: compute_microS_mushy
!
! !DESCRIPTION:
!
! Computes ice microstructural properties for updating hbrine
!
! NOTE: This subroutine uses thermosaline_vertical output to compute
! average ice permeability and the surface ice porosity
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine compute_microS_mushy (nx_block, ny_block,        &
                                   icells,   n_cat,               &
                                   indxii,    indxjj,             &
                                   trcrn, hice_old, hinS_o,       &
                                   sss, sst,                      &
                                   bTin,  bphin,  kmin, zphi_min, &
                                   bSin, brine_sal, brine_rho,    &
                                   iphin, ibrine_rho, ibrine_sal, &
                                   sice_rho)
!                            
! !USES:
!
      use ice_therm_mushy, only: temperature_mush, liquid_fraction, permeability
!
! !INPUT/OUTPUT PARAMETERS:                                
!
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
       
      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), intent(inout) :: &
         bSin        , & ! Bulk salinity (ppt) on bgrid
         bTin        , & ! Temperature on bgrid
         bphin           ! porosity on bgrid

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(in) :: &
         trcrn           

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(out) :: & 
         zphi_min    , & ! surface porosity
         kmin            ! average ice permeability (m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         hinS_o      , & ! previous timestep brine height (m)
         sice_rho        ! average ice density

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), intent(inout)  :: &
         brine_sal   , & ! equilibrium brine salinity (ppt) 
         brine_rho       ! Internal brine density (kg/m^3) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(inout)  :: &
         iphin       , & ! porosity on the igrid 
         ibrine_rho  , & ! brine rho on interface  
         ibrine_sal      ! brine sal on interface   
!
!EOP
!    
      real (kind=dbl_kind), dimension (nx_block*ny_block,nilyr) :: &
         cSin        , & ! Bulk salinity (ppt)
         cqin            ! enthalpy 

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2) :: &
         zTin        , & ! Temperature of ice layers on bgrid (^oC) 
         zSin        , & ! Salinity of ice layers on bgrid (^oC) 
         zqin        , & ! enthalpy on the bgrid
         zphin           ! porosity on the bgrid

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k               ! vertical biology layer index 
      
      real (kind=dbl_kind), dimension(icells) :: &
         surface_S   , & ! salinity of ice above hin > hinS 
         hinc_old    , & ! ice thickness (cell quantity) before current melt/growth (m)
         hinSc_old       ! brine thickness(cell quantity) before current melt/growth (m)
      
      real (kind=dbl_kind), dimension(icells,nblyr+1) :: &
         drho            ! brine density difference (kg/m^3)

     real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: &
         trtmp        , & ! temporary, remapped tracers   
         trtmp_q          ! temporary, remapped tracers   
     
      !----------------------------------------------------------------
      ! Define ice salinity and temperature on bgrid
      !---------------------------------------------------------------     

       do k = 1, nilyr
         do ij = 1, icells
           i = indxii(ij)
           j = indxjj(ij)            
           cSin(ij,k) = trcrn(i,j,nt_sice+k-1)
           cqin(ij,k) = trcrn(i,j,nt_qice+k-1)
         enddo
        enddo
        
       trtmp(:,:,:) = c0      
       trtmp_q(:,:,:) = c0       
       do ij = 1, icells   !Map Sin and qin profiles to bgc grid 
           i = indxii(ij)
           j = indxjj(ij)  
           if (hinS_o(i,j) > maxhinS*hice_old(i,j)) hinS_o(i,j) = maxhinS*hice_old(i,j)
           hinc_old(ij) = hice_old(i,j)
           hinSc_old(ij) = hinS_o(i,j)

          call remap_layers_bgc_plus_xy (ntrcr, nilyr,  nt_sice,  &
                             trcrn(i,j,:),    trtmp(i,j,:),      &
                             0,        nblyr+1,                    &
                             hinc_old(ij), hinc_old(ij),         &
                             cgrid(2:nilyr+1),                   &
                             bgrid(1:nblyr+1), surface_S(ij))
     
          call remap_layers_bgc_plus_xy (ntrcr, nilyr,  nt_qice,  &
                             trcrn(i,j,:),    trtmp_q(i,j,:),      &
                             0,        nblyr+1,                    &
                             hinc_old(ij), hinc_old(ij),         &
                             cgrid(2:nilyr+1),                   &
                             bgrid(1:nblyr+1), surface_S(ij))
     
       enddo
       do k = 1,nblyr+1
        do ij = 1, icells   
           i = indxii(ij)
           j = indxjj(ij)

           zSin(ij,k) = max(min_salin,trtmp(i,j,nt_sice+k-1))
           zqin(ij,k) = min(c0,trtmp_q(i,j,nt_qice+k-1))
           bSin(i,j,k) = zSin(ij,k)
           bSin(i,j,nblyr+2) =  sss(i,j) 
           zTin(ij,k) = temperature_mush(zqin(ij,k),zSin(ij,k))
           bTin(i,j,k) = zTin(ij,k)
           bTin(i,j,nblyr+2) = sst(i,j)
           zphin(ij,k) = liquid_fraction(zTin(ij,k),zSin(ij,k))
           bphin(i,j,k) = zphin(ij,k)
           bphin(i,j,nblyr+2) = c1
         enddo ! ij                          
      enddo ! k

      !------------------------------------------------------------------
      ! Define ice multiphase structure
      !----------------------------------------------------------------
     
      call prepare_hbrine(icells, indxii, indxjj, &
                            bSin, bTin,&
                            brine_sal, brine_rho, &
                            ibrine_sal, ibrine_rho, &
                            sice_rho, bphin,iphin, kmin, zphi_min, &
                            igrid,sss)
       
      call calculate_drho(icells,nx_block,ny_block,indxii,indxjj,igrid,bgrid,&
                                   brine_rho,ibrine_rho,drho)

 end subroutine compute_microS_mushy

!==========================================================================================
!BOP
!
! !ROUTINE: calculate_drho --  find density difference about interface grid points
!                              for gravity drainage parameterization
!                              *except for k == 1 which is the first grid point value
!                          
!
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
! authors     Nicole Jeffery, LANL
!
!
!
! !INTERFACE:
!
      subroutine calculate_drho &
                                  (icells,nx_block,ny_block,indxi,indxj,igrid,bgrid,&
                                   brine_rho,ibrine_rho,drho)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
       integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells       

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj  ! compressed indices 

      real (kind=dbl_kind), dimension (nblyr_hist), intent(in) :: &
         bgrid            ! biology nondimensional grid layer points 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid         ! biology grid interface points 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(in) :: &
         brine_rho           ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr + 1), intent(in) :: &
         ibrine_rho           ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (icells, nblyr+1), intent(out) :: & 
         drho                 ! brine difference about grid point (kg/m^3)
!
!EOP
!
    integer (kind=int_kind) :: &
         k, i,j, ij,m, mm ! indices

    integer (kind=int_kind), dimension (icells) :: &
         mstop, mstart

     real (kind=dbl_kind), dimension (icells, nblyr_hist) :: &  !on the zbgc vertical grid
         diff              ! difference array between grid points and 
                           ! the density averaging boundary 

     real (kind=dbl_kind), dimension (icells, nblyr + 1) :: &  !on the zbgc vertical grid
         rho_a         ,&  ! average brine density  above grid point (kg/m^3)
         rho_2a            ! average brine density  above and below grid points (kg/m^3)

       rho_a(:,:) = c0
       rho_2a(:,:) = c0
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
               rho_a(ij,2) = (brine_rho(i,j,2)*bgrid(2) + (ibrine_rho(i,j,2) + brine_rho(i,j,2))*&
                         p5*(igrid(2)-bgrid(2)) )/igrid(2)

            elseif (k > 2 .AND. k < nblyr+1) then 
               rho_a(ij,k) = (rho_a(ij,k-1)*igrid(k-1) +  (ibrine_rho(i,j,k-1) + brine_rho(i,j,k))*&
                         p5*(bgrid(k)-igrid(k-1)) + (ibrine_rho(i,j,k) + brine_rho(i,j,k))*p5*&
                         (igrid(k)-bgrid(k)))/igrid(k)
            else
               rho_a(ij,nblyr+1) = (rho_a(ij,nblyr)*igrid(nblyr) + (ibrine_rho(i,j,nblyr) + &
                        brine_rho(i,j,nblyr+1))*p5*(bgrid(nblyr+1)-igrid(nblyr)) + &
                        brine_rho(i,j,nblyr+1)*(igrid(nblyr+1)-bgrid(nblyr+1)))/igrid(nblyr+1)
               rho_a(ij,1) = brine_rho(i,j,2)   !for k == 1 use grid point value
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
              rho_2a(ij,1) = (rho_a(ij,1)*bgrid(2) + p5*(brine_rho(i,j,2) + ibrine_rho(i,j,2))*&
                      (igrid(2)-bgrid(2)))/igrid(2)
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
           endif
           drho(ij,k) = max(c0,c2*(rho_a(ij,k)-rho_2a(ij,k)), &
              c2*(brine_rho(i,j,k)-brine_rho(i,j,k+1))/real(nblyr,kind=dbl_kind))
          enddo
        enddo

     end subroutine calculate_drho

!=======================================================================
!BOP
!
! !ROUTINE: update_hbrine   Find changes in brine height 
!
! !DESCRIPTION:
!
! Changes include brine height increases from ice and snow surface melt, 
! congelation growth, and upward pressure driven flow from snow loading.
!  
! Decreases arise from downward flushing and bottom melt.  
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
     subroutine update_hbrine (icells, nx_block,ny_block,indxi, indxj, &
                               meltb, meltt, melts, &
                               dt, hin,hsn,hin_old, firstice, &
                               hinS, hinS_old, fbri, dhS_top,&
                               dhS_bottom,dh_top_chl, dh_bot_chl,kmin,zphi_min, &
                               darcy_V,darcy_V_chl, flood_val,melt_frac) 

! !USES:
     use ice_zbgc_shared, only: flood_frac
!
! !INPUT/OUTPUT PARAMETERS:
!
     integer (kind=int_kind), intent(in) :: &
        nx_block,ny_block, icells

     integer (kind=int_kind), dimension(nx_block*ny_block), intent(in) :: &
        indxi,indxj 

     real (kind=dbl_kind), intent(in) :: &
         dt  !timestep, tuning parameter
        
     real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &
         meltb,     & ! bottom melt over dt (m)
         meltt,     & ! true top melt over dt (m)
         melts,     & ! true snow melt over dt (m)
         hin,       & ! ice thickness (m)
         hsn,       & ! snow thickness (m)
         hin_old,   & ! past timestep ice thickness (m)
         hinS_old , & ! previous timestep hinS
         kmin         ! avg ice permeability

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout):: &
         darcy_V    , & ! Darcy velocity: m/s
         darcy_V_chl, & ! Darcy velocity: m/s for bgc
         dhS_top    , & ! change in top hinS before darcy flow
         dhS_bottom , & ! change in bottom hinS initially before darcy flow
         dh_bot_chl , & ! change in bottom for algae
         dh_top_chl , & ! change in bottom for algae
         hinS       , & ! thickness of brine (m) 
         fbri       , & ! brine height ratio tracer (hinS/hin) 
         melt_frac  , & ! fraction of top accumulation due to meltwater
         zphi_min     ! surface porosity

      logical (kind=log_kind), dimension(nx_block,ny_block), intent(in) :: & 
         firstice

      logical (kind=log_kind), dimension(nx_block,ny_block), intent(inout) :: & 
         flood_val       ! .true. if ocean flooded the surface
!
!EOP
! 
      real (kind=dbl_kind) :: &  
         hinS_min  , & ! thinS or hin 
         dhinS_hin , & ! hinS-hin
         hbrn      , & ! brine height  (m) hinS-h_o
         dhinS     , & ! change in brine surface
         h_o       , & ! new ocean surface from ice bottom (m)
         darcy_coeff,& ! magnitude of the Darcy velocity/hbrn (1/s)
         hbrn_new   ,& ! hbrn after flushing
         dhb_temp   ,& !
         dht_temp

      integer (kind=int_kind) :: &
         ij, i, j     

      real (kind=dbl_kind), parameter :: &
         dh_min = 0.001_dbl_kind, &  !
         run_off= c0   !fraction of melt that runs off directly to the ocean

    !move these two parameters to input file

    flood_frac = c0 !0.05_dbl_kind   ! 1% percent of water floods the surface
    flood_val(:,:) = .false.
    melt_frac(:,:) = c1

    do ij = 1,icells
       i = indxi(ij)
       j = indxj(ij)
         
       dhS_top(i,j) = dhS_top(i,j)     
       dhS_bottom(i,j) = dhS_bottom(i,j)
       hbrn = c0
       darcy_V(i,j) = c0
       darcy_V_chl(i,j) = c0
       hbrn_new = c0
       
       if (hinS_old(i,j) > thinS .AND. hin_old(i,j) > thinS) then

         dhS_top(i,j) = c0            
         if (meltt(i,j) > c0 .AND. hin_old(i,j) > hinS_old(i,j)) then  
          dhS_top(i,j) = -min(hin_old(i,j) -hinS_old(i,j), meltt(i,j))*rhoi/rhow
         endif 
         dhS_top(i,j) = dhS_top(i,j) - max(c0,melts(i,j))*rhos/rhow
         dhS_top(i,j) = (c1-run_off)*dhS_top(i,j)
         dhinS = dhS_bottom(i,j) - dhS_top(i,j)  
         h_o = rhosi/rhow*hin(i,j) + rhos/rhow*hsn(i,j) 
         hinS(i,j)  = max(his_min,hinS_old(i,j) + dhinS)
         hbrn = hinS(i,j) - h_o

         if (hbrn > c0 .AND. hinS(i,j) > thinS ) then   
              darcy_coeff = max(c0,kmin(i,j)*gravit/viscos/hinS_old(i,j))    
              hbrn_new = hbrn*exp(-darcy_coeff/zphi_min(i,j)*dt)
              hinS(i,j) = max(thinS,h_o + hbrn_new)    
              hbrn_new = hinS(i,j)-h_o   
              darcy_V(i,j) = - SIGN((hbrn-hbrn_new)/dt*zphi_min(i,j),hbrn)
              dh_top_chl(i,j) = dhS_top(i,j)  - darcy_V_chl(i,j)*dt/zphi_min(i,j)
              dhS_top(i,j) = dhS_top(i,j) + SIGN((hbrn-hbrn_new),hbrn) 
              dht_temp = dhS_top(i,j)
         elseif (hbrn < c0) then
              darcy_coeff = max(c0,kmin(i,j)*gravit/viscos/hinS_old(i,j))  
              if (hinS(i,j) .GE. hin(i,j)) zphi_min(i,j)  = phi_snow ! c1-rhos/rhoi
              hbrn_new = hbrn*exp(-darcy_coeff/zphi_min(i,j)*dt)
              hinS(i,j) = max(his_min,h_o + hbrn_new)
              hbrn_new = hinS(i,j) - h_o
              darcy_V(i,j) = -SIGN((hbrn-hbrn_new)/dt*zphi_min(i,j),hbrn)
              darcy_V_chl(i,j) = c0 !darcy_V(i,j)
              dh_top_chl(i,j) = dhS_top(i,j) - darcy_V_chl(i,j)*dt/zphi_min(i,j)
              dhS_top(i,j) = dhS_top(i,j) + SIGN((hbrn-hbrn_new),hbrn)
              dht_temp = dhS_top(i,j)
              dhS_top(i,j) = dhS_top(i,j) + flood_frac*hbrn_new
              ! flood_val(i,j) = .true.
              if (dht_temp .LE. c0 .AND. dhS_top(i,j) < c0 .and. flood_val(i,j)) then
                   melt_frac(i,j) = abs(dht_temp/dhS_top(i,j))
              elseif (dht_temp > c0 .OR. dhS_top(i,j) > c0) then
                 write(nu_diag, *) 'ice_brine error: dht_temp, dhS_top(i,j)',dht_temp, dhS_top(i,j)
                 write(nu_diag, *) 'hinS(i,j),hbrn_new,hbrn,kmin(i,j)',&
                 hinS(i,j),hbrn_new,hbrn,kmin(i,j)
              endif
          endif
    
        else    ! very thin brine height 
            h_o = rhosi/rhow*hin(i,j) + rhos/rhow*hsn(i,j) 
            hinS_min = min(thinS,hin(i,j))  !max(his_min,min(h_o,maxhinS*hin(i,j)))
            hinS(i,j) =  max(hinS_min,hinS_old(i,j)+dhS_bottom(i,j)-dhS_top(i,j))
            dhinS_hin = hinS(i,j)-h_o
            if (abs(dhinS_hin) > dh_min) hinS(i,j) = max(his_min,h_o + SIGN(dh_min,dhinS_hin))
        endif 
        
        fbri(i,j) = hinS(i,j)/hin(i,j)

        if (dh_bot_chl(i,j) > c0) then
          dh_bot_chl(i,j) = p1*dhS_bottom(i,j)
        else
         dh_bot_chl(i,j) = dhS_bottom(i,j)
        endif
    enddo

      end subroutine update_hbrine

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_hbrine - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_hbrine(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for hbrine
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
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
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
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
      enddo

      do n = 1, ncat
          call ice_read(nu_restart_hbrine,0,first_ice_real(:,:,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
      enddo
      do iblk = 1, nblocks
          do n = 1,ncat
           do j = 1, ny_block
           do i = 1, nx_block
             if (first_ice_real(i,j,n,iblk) .GE. c1) then
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
!
!BOP
!
! !IROUTINE: write_restart_hbrine - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_hbrine(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for a hbrine restart
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_domain, only: nblocks
      use ice_state
      use ice_flux, only: sss  
      use ice_restart, only: lenstr, restart_dir, restart_file, pointer_file, runtype
      use ice_read_write, only: ice_open, ice_write
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
      !--------------------------
      !hbrine and first ice
      !--------------------------

      do n = 1, ncat
         call ice_write(nu_dump_hbrine,0,trcrn(:,:,nt_fbri,n,:),'ruf8',diag)
      enddo

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

      do n = 1,ncat
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
              
!        
! !USES:
!
      use ice_broadcast, only: broadcast_scalar
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, pbloc, &
                                plat, plon
      use ice_domain_size, only: ncat, nltrcr
      use ice_state, only: aice, aicen, vicen, vice, trcr, nt_sice , nt_fbri, &
                          trcrn, hbrine
      use ice_zbgc_shared, only: darcy_V
!
! !INPUT PARAMETERS:
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
                       phinS(n) = trcr(i,j,nt_fbri,iblk)*vice(i,j,iblk)/aice(i,j,iblk)
               if (aicen(i,j,1,iblk)> c0)&
                       phinS1(n) = trcrn(i,j,nt_fbri,1,iblk)*vicen(i,j,1,iblk)/&
                                                aicen(i,j,1,iblk)
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
          write(nu_diag,*) '------ hbrine ------'
          write(nu_diag,900) 'hbrine, (m)        = ',phinS(1),phinS(2)
          write(nu_diag,900) 'fbri, cat1 (m)     = ',pfbri(1),pfbri(2)
          write(nu_diag,900) 'hbrine cat1, (m)   = ',phinS1(1),phinS1(2)  
          write(nu_diag,900) 'darcy_V cat1, (m/s)= ',pdarcy_V(1),pdarcy_V(2)             
          write(nu_diag,*) '                         '
          write(nu_diag,*) '------ Thermosaline Salinity ------'
          write(nu_diag,803) 'Sice(1) bulk S (ppt) ','Sice(2) bulk S'
          write(nu_diag,*) '---------------------------------------------------'
          write(nu_diag,802) ((pSin(n,k),n=1,2), k = 1,nilyr)              
          write(nu_diag,*) '                         '
      endif                   ! print_points
      endif                   ! my_task = master_task 

  802 format (f24.17,2x,f24.17)
  803 format (a25,2x,a25)
  900 format (a25,2x,f24.17,2x,f24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
  903 format (a25,5x,i4,1x,i4,1x,i4,1x,i4,7x,i4,1x,i4,1x,i4,1x,i4)

      end subroutine hbrine_diags

!=======================================================================
!BOP
!
! !ROUTINE: prepare_hbrine -- calculate quantities to use in solve hbrn
!
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
! authors     Nicole Jeffery, LANL
!
!
! June 2011 by N. Jeffery
!
! !INTERFACE:
!
      subroutine prepare_hbrine &
                                     (icells, indxi, indxj, Sin, zTin, brine_sal, &
                                      brine_rho, ibrine_sal, ibrine_rho, &
                                       sice_rho, zphin,iphin,kmin,zphi_min,&
                                       igrid,sss)
!
! !USES:

       use ice_therm_shared, only: calculate_Tin_from_qin

! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
        icells

      integer (kind=int_kind), dimension(nx_block*ny_block), intent(in) :: &
        indxi,indxj 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(in) :: &
         Sin           , &! Salinity of ice layers on bio grid 
         zTin              ! Temperature of ice layers on bio grid for history file (^oC) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(inout) :: &
         zphin          ,& ! Porosity of layers
         brine_sal     ,& ! equilibrium brine salinity (ppt)  
         brine_rho        ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(inout) :: &
         ibrine_rho    ,& ! brine rho on interface  
         ibrine_sal   , &      !brine sal on interface   
         iphin            ! porosity on interface

   real (kind=dbl_kind), dimension (nblyr+1), intent(in):: &
         igrid         ! biology grid interface points

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         sss        ! sea surface salinity

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         sice_rho    , &    ! avg sea ice density
         kmin        , &    ! harmonic average permeability (m^2)
         zphi_min           !minimum porosity
!
!EOP
! 
      real (kind=dbl_kind), dimension(icells, nblyr) :: &
          kin       !  permeability  
    
      real (kind=dbl_kind), dimension(icells):: &
          k_min, ktemp
     
      integer (kind=int_kind) :: &
           k,m, i, j, ij   ! tracer indice

     ! integer (kind=int_kind), dimension(icells) :: &
     !     loc

     !-----------------------------------------------------------------
     !  calculate equilibrium brine density and gradients 
     !-----------------------------------------------------------------

     ! zphi_avg(:) = c0
    
      sice_rho(:,:) = c0
      kmin(:,:) = c0
      ktemp(:) = c0
      
      do k = 1, nblyr+2
           do ij = 1,icells
             i = indxi(ij)
             j = indxj(ij) 

               brine_sal(i,j,k) = a1*zTin(i,j,k) + a2*zTin(i,j,k)**2 + a3*zTin(i,j,k)**3
               brine_rho(i,j,k) = b1 + b2*brine_sal(i,j,k)  !  + b3*brine_sal(i,j,k)**2
               zphin(i,j,k) = min(c1,max(puny, Sin(i,j,k)*rhosi/(brine_sal(i,j,k)*brine_rho(i,j,k)))) 
  
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
      !            zphi_avg(ij) = zphi_avg(ij) + zphin(i,j,k)*(igrid(k)-igrid(k-1))
                  sice_rho(i,j) = sice_rho(i,j) + (rhoi*(c1-zphin(i,j,k)) + &
                              brine_rho(i,j,k)*zphin(i,j,k))*(igrid(k)-igrid(k-1))
                  
               elseif (k .GE. nblyr+2) then

                  k_min(ij) = MINVAL(kin(ij,:))
                  zphi_min(i,j) = zphin(i,j,2) 
                         
               endif
           enddo  !ij   
      enddo                !k  
    
      do k = 2, nblyr
           do ij = 1,icells
             i = indxi(ij)
             j = indxj(ij) 
 
             if (k_min(ij) >  c0) then
               ktemp(ij) = ktemp(ij) + c1/kin(ij,k-1)
               if (k == nblyr) then
                 ktemp(ij) = ktemp(ij) + c1/kin(ij,nblyr) 
                 kmin(i,j) = real(nblyr,kind=dbl_kind)/ktemp(ij)
               endif
              else
               kmin(i,j) = k_min(ij)
             endif
             ibrine_sal(i,j,k) = (brine_sal(i,j,k+1)*(igrid(k+1)-igrid(k)) + &
                                 brine_sal(i,j,k)*(igrid(k)-igrid(k-1)))/(igrid(k+1)-igrid(k-1))
             ibrine_rho(i,j,k) = (brine_rho(i,j,k+1)*(igrid(k+1)-igrid(k)) + &
                                 brine_rho(i,j,k)*(igrid(k)-igrid(k-1)))/(igrid(k+1)-igrid(k-1))
             iphin(i,j,k) = min(c1, max(puny, (zphin(i,j,k+1)*(igrid(k+1)-igrid(k)) + &
                                 zphin(i,j,k)*(igrid(k)-igrid(k-1)))/(igrid(k+1)-igrid(k-1))))
           enddo !ij
       enddo                !k         

      end subroutine prepare_hbrine

      end module ice_brine

!=======================================================================
