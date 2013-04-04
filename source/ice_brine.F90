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
      use ice_domain_size, only: nilyr, nblyr, nblyr_hist
      use ice_fileunits, only: nu_diag
      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_state, only: ntrcr, nt_bgc_S, nt_qice, nt_sice
      use ice_zbgc_public, only: cgrid, bgrid, exp_h, k_o, rhosi, &
           thin, min_salin, Ra_c, igrid, remap_layers_bgc_plus_xy, remap_layers_bgc_plus
!
!EOP
!
      implicit none

      private
      public :: preflushing_changes, compute_micros, update_hbrine
 
      real (kind=dbl_kind), parameter :: &   
         maxhinS = 1.90_dbl_kind, & ! brine overflow condition if hinS > maxhinS*hin
         viscos     = 2.1e-6_dbl_kind , & ! kinematic viscosity (m^2/s) More like 2 for cold brine of 90 PSU
         a1 = -21.4_dbl_kind, &    !(psu/C)  Brine salinity as a cubic function of T 
         a2 = -0.886_dbl_kind, &   !(psu/C^2)
         a3 = -0.012_dbl_kind, &   !(psu/C^3)
                                   ! cox and Weeks
         b1 = 1000.0_dbl_kind, &   !rhofresh (kg/m^3)  Brine density as a quadratic of brine salinity
         b2 = 0.8_dbl_kind         !(kg/m^3/ppt)

!=======================================================================

      contains

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
                                   vsnon,  meltb, melt_top, congel,&
                                   snoice, dsnow, hice_old,        &   
                                   fbri,  dh_top, dh_bot, dh_bot_chl, &
                                   dhi_top, dhi_bot,               &
                                   hinS_o,hin,hsn,                 &
                                   first_ice)


                            
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
         hinS_o       , & !  old brine heights
         melt_top     , & ! top ice melt     (m)           
         dsnow         ! change in snow top (m)


      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         first_ice      ! for first category ice only .true. initialized values should be used

      
!
!EOP
!

     
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          ! horizontal index, combines i and j loops
  

      real (kind=dbl_kind) :: &
         hin_old     , & ! ice thickness before current melt/growth (m)
         hsn_old      ! snow thickness before current melt/growth (m)

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
           if (hin(i,j) .LE. c0 ) then
               write(nu_diag, *) 'preflushing: hin <= 0 at  i, j:',  i, j
               write(nu_diag, *) 'vicen(i,j), aicen(i,j)'
               write(nu_diag, *) 'fbri(i,j),hice_old(i,j)'
               write(nu_diag, *) fbri(i,j), hice_old(i,j)
               call abort_ice ('ice:ice_brine error')
           endif
           hsn(i,j) = c0
           if (hin(i,j) > c0) hsn(i,j) = max(c0,vsnon(i,j) / (aicen(i,j)))
          

           hin_old = max(c0,hin(i,j) + meltb(i,j) + melt_top(i,j) - &
                          congel(i,j)  - snoice(i,j))            
           hsn_old = max(c0,hsn(i,j) - dsnow(i,j))

           dsnow(i,j) = hsn(i,j) - hsn_old       
           dhice = hin_old - hice_old(i,j)   !change due to subl/cond
           dh_top(i,j) = melt_top(i,j)   
           dh_bot(i,j) = congel(i,j) - meltb(i,j)   
           dh_bot_chl(i,j) = congel(i,j) - meltb(i,j)

           if ((hice_old(i,j) < puny) .OR. (hin_old < puny) .OR. first_ice(i,j)) then

             hin_old = hin(i,j)  
             dh_top(i,j) = c0
             dh_bot(i,j) = c0
             dhi_bot(i,j) = c0
             dhi_top(i,j) = c0
             dh_bot_chl(i,j) = c0
             dsnow(i,j) = c0
             melt_top(i,j) = c0  
             fbri(i,j) = c1 

           else

             dhi_bot(i,j) = congel(i,j) - meltb(i,j)
             dhi_top(i,j) = melt_top(i,j) - snoice(i,j) 
             dhi_top(i,j) = dhi_top(i,j) - dhice  
                     !(dhice > 0 is cond (ice gain), dhice < 0 is sub (ice loss)

           endif

           hinS_o(i,j) = fbri(i,j)* hin_old
           hice_old(i,j) = hin_old   

         enddo  !ij


 end subroutine preflushing_changes

!=======================================================================
!BOP
!
! !ROUTINE: compute_microS
!
! !DESCRIPTION:
!
! Computes ice microstructural properties for zbgc and zsalinity 
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
      subroutine compute_microS    (nx_block, ny_block,           &
                                   icells,   n_cat,               &
                                   indxii,    indxjj,             &
                                  ! tcells, tndxi, tndxj, tndxij, &
                                  ! pcells, pndxi, pndxj, pndxij, &
                                  ! qcells, qndxi, qndxj, qndxij, &
                                   trcrn, hice_old,               &
                                   hinS_o,  sss, sst,             &
                                   zTin,  zphin,  kmin, zphi_min,    &
                                   Rayleigh_criteria, first_ice,  &
                                   Sin, brine_sal, brine_rho,     &
                                   iphin, ibrine_rho, ibrine_sal, &
                                   sice_rho, sloss)


                            
! !USES:
!
      use ice_therm_shared, only: solve_Sin, calculate_Tin_from_qin
      use ice_calendar, only: istep1, time
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
 

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn           

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         hinS_o      ! old brine height

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), &
         intent(inout) :: &
         zTin          ! Temperature of ice layers on bio grid for history file (^oC) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), &
         intent(out) :: &
         Sin

       
      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(inout) :: &
         zphin            ! Porosity of layers

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(out) :: & 
         zphi_min     , & ! surface porosity
         kmin        , & ! average ice permeability (m^2)
         sloss           ! (g/m^2) salt from brine runoff for hinS > maxhinS*hin

     
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         sice_rho           ! average ice density


      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         Rayleigh_criteria   ! .true. if ice exceeded a minimum thickness hin >= Ra_c 
        
      logical (kind=log_kind), dimension (nx_block,ny_block), & 
         intent(in) :: &
         first_ice            ! .true. if ice is newly formed (Category 1 ice only)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), intent(inout)  :: &
         brine_sal       ,& ! equilibrium brine salinity (ppt) 
         brine_rho          ! Internal brine density (kg/m^3) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(inout)  :: &
         iphin        , & !  porosity on the igrid 
         ibrine_rho    ,& ! brine rho on interface  
         ibrine_sal       !brine sal on interface   
!
!EOP
!

     
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
         surface_S   , & ! salinity of ice above hin > hinS 
         hinc_old    , & ! ice thickness (cell quantity) before current melt/growth (m)
         hinSc_old       ! brine thickness(cell quantity) before current melt/growth (m)

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
       sloss(:,:) = c0  
       zTin(:,:,:) = c0
       Sin(:,:,:) = c0

        do ij = 1, icells
      
           i = indxii(ij)
           j = indxjj(ij)
 
           hinc_old(ij) = hice_old(i,j)
           hinSc_old(ij) = hinS_o(i,j)

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
           if (n_cat == 1 .AND. hinSc_old(ij) < Ra_c) then
              Rayleigh(ij) = Rayleigh_criteria(i,j) ! only category 1 ice can be false 
           endif
           
           
            if (first_ice(i,j)) then
               qcells = qcells + 1
               qndxi(qcells) = i
               qndxj(qcells) = j
               qndxij(qcells) = ij 
            elseif (hinS_o(i,j) > maxhinS*hice_old(i,j)) then
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

     if (solve_Sin) then  
         
        do ij = 1, qcells   !first_ice:  Map initial profile to bgc grid 
          i = qndxi(ij)
          j = qndxj(ij) 
          qij =qndxij(ij)


          call remap_layers_bgc_plus_xy (ntrcr, nilyr,  nt_sice,  &
                             trcrn(i,j,:),    trtmp(i,j,:),      &
                             0,        nblyr,                    &
                             hinc_old(qij), hinc_old(qij),         &
                             cgrid(2:nilyr+1),                   &
                             bgrid(2:nblyr+1), surface_S(qij))
     
        enddo
        
        do ij = 1, pcells   !brine overflow: Salinity flux stored in sloss
         i = pndxi(ij)
         j = pndxj(ij) 
         pij =pndxij(ij)

         call remap_layers_bgc_plus_xy ( ntrcr, nblyr, nt_bgc_S, &
                             trcrn(i,j,:),    trtmp(i,j,:),      &
                             0,        nblyr,                    &
                             hinSc_old(pij), maxhinS*hinc_old(pij),&
                             bgrid(2:nblyr+1),         &
                             bgrid(2:nblyr+1), surface_S(pij))
      
        enddo
 
       do k = 1, nblyr 
        
         do ij = 1, tcells    
            i = tndxi(ij)
            j = tndxj(ij) 
            tij = tndxij(ij)

            Sin(i,j,k+1) = max(min_salin,trcrn(i,j,nt_bgc_S+k-1))   
            Sin(i,j,1) = Sin(i,j,2)
            Sin(i,j,nblyr+2) =  sss(i,j)

            if (trcrn(i,j,nt_bgc_S+k-1) < min_salin-puny) then
                write(nu_diag, *) 'Bad value compute_microS, tcells'
                write(nu_diag, *) 'trcrn(i,j,nt_bgc_S+k-1),i,j,k',trcrn(i,j,nt_bgc_S+k-1),i,j,k
                write(nu_diag, *) 'hinS_o(i,j), hinSc_old(tij),i,j,tij',hinS_o(i,j), hinSc_old(tij),i,j,tij
                write(nu_diag, *) 'istep1, time, n_cat',istep1,time,n_cat
               ! write(nu_diag, *) 'TLAT(i,j),TLON(i,j)'             
               ! write(nu_diag, *) TLAT(i,j)*rad_to_deg,TLON(i,j)*rad_to_deg
                call abort_ice ('ice: ice_brine.F90 error')
            endif
         enddo ! ij

         do ij = 1, qcells   !first_ice(i,j)
            i = qndxi(ij)
            j = qndxj(ij) 
            qij =qndxij(ij)

            trcrn(i,j,nt_bgc_S+k-1) = max(min_salin,trtmp(i,j,nt_sice+k-1))
            Sin(i,j,k+1) = max(min_salin,trcrn(i,j,nt_bgc_S+k-1))   
            Sin(i,j,1) = Sin(i,j,2) 
            Sin(i,j,nblyr+2) =  sss(i,j) 
  
            if (trcrn(i,j,nt_bgc_S+k-1) < min_salin-puny) then
                write(nu_diag, *) 'Bad value compute_microS, qcells'
                write(nu_diag, *) 'trcrn(i,j,nt_bgc_S+k-1),i,j,k',trcrn(i,j,nt_bgc_S+k-1),i,j,k
                write(nu_diag, *) 'hinS_o(i,j), hinSc_old(qij),i,j,qij',hinS_o(i,j), hinSc_old(qij),i,j,qij
                write(nu_diag, *) 'istep1, time, n_cat',istep1,time,n_cat
                !write(nu_diag, *) 'TLAT(i,j),TLON(i,j)'             
                !write(nu_diag, *) TLAT(i,j)*rad_to_deg,TLON(i,j)*rad_to_deg
                call abort_ice ('ice:  ice_brine.F90 error')
            endif
             
         enddo ! ij

                             !elseif (hinS_o(i,j) > maxhinS*hice_old(i,j)) then
         do ij = 1, pcells   !hinS_old(qij) > maxhinS*hin_old(qij)
             i = pndxi(ij)
             j = pndxj(ij) 
             pij =pndxij(ij)

             Sin(i,j,k+1) = max(min_salin,trtmp(i,j,nt_bgc_S+k-1))
             sloss(i,j) = sloss(i,j) + rhosi*(hinS_o(i,j)*trcrn(i,j,nt_bgc_S+k-1) - &
                          maxhinS*hice_old(i,j)*Sin(i,j,k+1))*(igrid(k+1)-igrid(k))
             trcrn(i,j,nt_bgc_S+k-1) = Sin(i,j,k+1)
                               
             Sin(i,j,1) = Sin(i,j,2)
             if (k == nblyr) then
               Sin(i,j,nblyr+2) =  sss(i,j)
               hinSc_old(pij) = maxhinS*hinc_old(pij)
               hinS_o(i,j) = hinSc_old(pij)
             endif

            if (trcrn(i,j,nt_bgc_S+k-1) < min_salin-puny) then
                write(nu_diag, *) 'Bad value compute_microS, pcells'
                write(nu_diag, *) 'trcrn(i,j,nt_bgc_S+k-1),i,j,k',trcrn(i,j,nt_bgc_S+k-1),i,j,k
                write(nu_diag, *) 'hinS_o(i,j), hinSc_old(pij),i,j,pij',hinS_o(i,j), hinSc_old(pij),i,j,pij
                write(nu_diag, *) 'istep1, time, n_cat',istep1,time,n_cat
                !write(nu_diag, *) 'TLAT(i,j),TLON(i,j)'             
                !write(nu_diag, *) TLAT(i,j)*rad_to_deg,TLON(i,j)*rad_to_deg
                call abort_ice ('ice:d diffuse_bio error')
            endif

          enddo !ij

      enddo ! k
      
    else   ! not solve_Sin  

    do k = 1, nblyr

         do ij = 1, icells    
            i = indxii(ij)
            j = indxjj(ij) 
       
             Sin(i,j,k+1) = trcrn(i,j,nt_bgc_S+k-1) !read from bio grid
             Sin(i,j,1) =  Sin(i,j,2) 
             Sin(i,j,nblyr+2) =  sss(i,j) 

          enddo  !ij
       enddo     !k

    endif   !solve_Sin

!echmod reset Sin for testing
!      do k = 1, nilyr
!       do ij = 1, icells    
!          i = indxii(ij)
!          j = indxjj(ij) 
!              Sin(i,j,k+1) = trcrn(i,j,nt_sice+k-1)
!       enddo
!      enddo
!echmod   

  
      !-----------------------------------------------------------------
      ! sea ice temperature  for bio grid
      !----------------------------------------------------------------- 

      
       do k = 1, nilyr
         do ij = 1, icells    
            i = indxii(ij)
            j = indxjj(ij) 

             Tmlts = -trcrn(i,j,nt_sice+k-1)*depressT
             trtmp0(i,j,nt_qice+k-1) = calculate_Tin_from_qin(trcrn(i,j,nt_qice+k-1),Tmlts)

          enddo  !ij
         enddo   ! k
         
        trtmp(:,:,:) = c0                
      
        !CICE to Bio: remap temperatures

        call remap_layers_bgc_plus (nx_block,ny_block,        &
                             indxii,   indxjj,           &
                             icells,                   &
                             ntrcr,                    &
                             nilyr,                    &
                             nt_qice,                  &
                             trtmp0,    trtmp,          &
                             0,        nblyr,          &
                             hinc_old, hinSc_old,         &
                             cgrid(2:nilyr+1),         &
                             bgrid(2:nblyr+1), surface_S)

        do k = 1, nblyr
         do ij = 1, icells
            i = indxii(ij)
            j = indxjj(ij)

            Tmlts = -Sin(i,j,k+1) * depressT
            zTin(i,j,k+1) = min(Tmlts,trtmp(i,j,nt_qice+k-1)) 
            !--------------------------
            !boundary points
            !-------------------------

            Tmlts =  -Sin(i,j,1)* depressT 
            zTin(i,j,1) = min(Tmlts,trtmp(i,j,nt_qice)) 
            Tmlts = -Sin(i,j,nblyr_hist)* depressT  
            zTin(i,j,nblyr+2) = sst(i,j)
           enddo   !ij
      enddo !k

    
      !------------------------------------------------------------------
      ! Define ice multiphase structure
      !----------------------------------------------------------------
     
        call prepare_hbrine(icells, indxii, indxjj, &
                            Sin, zTin,&
                            brine_sal, brine_rho, &
                            ibrine_sal, ibrine_rho, &
                            sice_rho, zphin,iphin, kmin, zphi_min, &
                            igrid,sss)
       

 end subroutine compute_microS

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
     subroutine update_hbrine &
                               (icells, nx_block,ny_block,indxi, indxj,&
                               melt_top, meltb, dsnow, &
                               dt, hin,hsn,hin_old, first_ice, &
                               hinS, hinS_old, fbri, dhS_top,&
                               dhS_bottom,kmin,zphi_min, &
                               darcy_V) 

! !USES:
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
         melt_top,  &  ! top melt over dt (m)
         meltb,     &  ! bottom melt over dt (m)
         dsnow,     & ! change in snow thickness (m)
         hin,       & ! ice thickness (m)
         hsn,       & ! snow thickness (m)
         hin_old,   & ! past timestep ice thickness (m)
         hinS_old   , & ! previous timestep hinS
         kmin,      & ! avg ice permeability
         zphi_min      ! surface porosity

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout):: &
         darcy_V    , & ! Darcy velocity: m/s
         dhS_top    , & ! change in top hinS before darcy flow
         dhS_bottom , & ! change in bottom hinS initially before darcy flow
         hinS       , & ! thickness of brine (m) 
         fbri           ! brine height ratio tracer (hinS/hin) 


      logical (kind=log_kind), dimension(nx_block,ny_block), intent(in) :: & 
         first_ice
!
!EOP
! 

      real (kind=dbl_kind) :: &  
         hbrn      , & ! brine height  (m) hinS-h_o
         snowm     , & ! snow melt  (dsnow < 0)
         snowf     , & ! snow accumulation (dsnow > 0)
         dhinS     , & ! change in brine surface
         h_o       , & ! new ocean surface from ice bottom (m)
         darcy_coeff,& ! magnitude of the Darcy velocity/hbrn (1/s)
         hbrn_new      ! nbrn after flushing

 
      integer (kind=int_kind) :: &
         ij, i, j     

     do ij = 1,icells
       i = indxi(ij)
       j = indxj(ij)
         
       dhS_top(i,j) = dhS_top(i,j)     
       dhS_bottom(i,j) = dhS_bottom(i,j)
      
       hbrn = c0
       darcy_V(i,j) = c0
       

       if (hinS_old(i,j) > thin ) then

         if (dsnow(i,j) > c0) then
           snowf = dsnow(i,j)
           snowm = c0
         else
           snowf = c0
           snowm = abs(dsnow(i,j))
         endif

      
         dhS_top(i,j) = c0            
         if (melt_top(i,j) > c0 .AND. hin_old(i,j) > hinS_old(i,j)) then  
          dhS_top(i,j) = -min(hin_old(i,j) -hinS_old(i,j), melt_top(i,j))*rhoi/rhosi 
         endif 
     

         dhS_top(i,j) = dhS_top(i,j) - snowm*rhos/rhosi         
         dhinS = dhS_bottom(i,j) - dhS_top(i,j)  
         h_o = rhosi/rhow*hin(i,j) + rhos/rhow*hsn(i,j) 

         hinS(i,j)  = hinS_old(i,j) + dhinS
         hbrn = hinS(i,j) - h_o

         if (abs(hbrn) > c0 ) then   

              darcy_coeff = max(c0,kmin(i,j)*gravit/viscos/hin(i,j))    
              hbrn_new = hbrn*exp(-darcy_coeff/zphi_min(i,j)*dt)
              hinS(i,j) = h_o + hbrn_new       
              darcy_V(i,j) =  -(hbrn-hbrn_new)/dt*zphi_min(i,j)
              dhS_bottom(i,j) = dhS_bottom(i,j)-(hbrn-hbrn_new)       
        
         endif

        else    ! very thin brine height 
            hinS(i,j) = hin(i,j)
        endif 


        fbri(i,j) = hinS(i,j)/hin(i,j)
 
       enddo   


 end subroutine update_hbrine

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
 
             if (k_min(ij) > puny) then
               ktemp(ij) = ktemp(ij) +  c1/kin(ij,k-1)
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
