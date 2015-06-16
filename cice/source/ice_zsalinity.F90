!=======================================================================
!
! Vertical salinity (trcrn(nt_bgc_S)) is solved on the bio grid (bgrid and igrid)
! with domain defined by the dynamic brine height (trcrn(nt_fbri) * vicen/aicen).
! The CICE Bitz and Lipscomb thermodynamics is solve in the cgrid with height
! vicen/aicen.
! Gravity drainage is parameterized as nonlinear advection
! Flushing is incorporated in the boundary changes and a darcy flow. 
! (see Jeffery et al., JGR, 2011).  
!
! authors: Nicole Jeffery, LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_zsalinity

      use ice_kinds_mod
      use ice_constants
      use ice_domain_size
      use ice_fileunits, only: nu_diag, nu_dump_S, nu_restart_S, &
           nu_rst_pointer, flush_fileunit
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_state
      use ice_zbgc_shared
      use ice_blocks, only: nx_block, ny_block
!
      implicit none

      private
      public :: init_zsalinity, S_diags, write_restart_S, solve_zsalinity, &
           column_sum_S, merge_S_fluxes

      real (kind=dbl_kind), parameter :: & 
         max_salin = 200.0_dbl_kind   , & !(ppt) maximum bulk salinity
         lapidus_g = 0.3_dbl_kind     , & ! constant for artificial viscosity/diffusion during growth
         lapidus_m = 0.007_dbl_kind !0.0035_dbl_kind  ! constant for artificial diffusion during melt

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         Rayleigh_real      ! .true. = c1, .false. = c0
    
!=======================================================================

      contains

!=======================================================================
!
!  Initialize vertical profile of biogeochemistry
!
      subroutine init_zsalinity (sss)

      use ice_domain, only: nblocks

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         sss      ! sea surface salinity (ppt)
      
      ! local variables

      integer (kind=int_kind) :: &
           i, j, iblk       , & ! horizontal indices
           k,m              , & ! vertical index
           n                    ! category index

      real (kind=dbl_kind) :: &  
         zspace                 !grid spacing for CICE vertical grid

      if (nblyr .LE. 7) then
          dts_b = 300.0_dbl_kind
      else
          dts_b = 50.0_dbl_kind 
      endif

     !-----------------------------------------------------------------------------   
     !     BGC Layer Model
     !-----------------------------------------------------------------------------   

      bphi(:,:,:,:,:) = c0   ! initial porosity for no ice 
      iDi(:,:,:,:,:)  = c0   !interface diffusivity
      bTiz(:,:,:,:,:) = c0   !initial bio grid ice temperature

      if (restart_S) then
          call read_restart_S
      else

      Rayleigh_criteria(:,:,:) = .false.    ! no ice initial condition 
      iki(:,:,:,:,:) = c0                   ! permeability

      do n = 1,ncat
           do iblk = 1, nblocks
              do j = 1, ny_block
                 do i = 1, nx_block 
                    do k = 1,nblyr
                       trcrn(i,j,nt_bgc_S+k-1,n,iblk) = sss(i,j,iblk)*salt_loss
                    enddo      !k
                 enddo        !i
              enddo        !j
           enddo           !iblk
       enddo           !n

      endif    !restart_S

      end subroutine init_zsalinity

!=======================================================================
!
! update vertical salinity 
! 
      subroutine solve_zsalinity(nx_block,             ny_block,     &
                                   icells,     n_cat,  dt,           &
                                   indxii,             indxjj,       & 
                                   trcrn_S,            trcrn_q,      &
                                   trcrn_Si,           aicen,        &
                                   vicen, &  
                                   bSin,               bTin,         &
                                   bphin,              iphin,        &
                                   ikin,               hbr_old,      &
                                   hbrin,              hin,          &
                                   hin_old,            iDin,         &
                                   darcy_V,            brine_sal,    &
                                   brine_rho,          ibrine_sal,   &
                                   ibrine_rho,                       &
                                   Rayleigh_criteria,                &
                                   first_ice,          sss,          &
                                   sst,                dh_top,       &
                                   dh_bot,                           &  
                                   TLAT,               TLON,         &
                                   l_stop,             istop,        &
                                   jstop,              fzsaln,       &
                                   fzsaln_g,           bphi_min) 
 

      use ice_therm_shared, only: calculate_Tin_from_qin, solve_zsal
      use ice_calendar, only: istep1, time

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells              ! number of true cells with aicen > 0
          
      integer (kind=int_kind), intent(in) :: &
         n_cat           ! category number 
                    
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxii, indxjj ! compressed indices for icells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt             ! time step

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         sss         , & ! ocean salinity (ppt)
         sst         , & ! ocean temperature (oC)
         hin_old     , & ! old ice thickness (m)
         dh_top      , & ! brine change in top and bottom for diagnostics (m)
         dh_bot      , &
         TLAT        , &
         TLON        , & 
         darcy_V    

      real (kind=dbl_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         hbr_old     , & ! old brine height  (m)
         hin         , & ! new ice thickness (m)
         hbrin       , & ! new brine height  (m)
         bphi_min   
 
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         fzsaln        , & ! total flux of salt out of ice over timestep(kg/m^2/s)
         fzsaln_g          ! gravity drainage flux of salt  over timestep(kg/m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block, nblyr+2), &
         intent(inout) :: &
         bTin          , & ! Ice Temperature ^oC (on bio grid)
         bphin             ! Ice porosity (on bio grid)

      real (kind=dbl_kind), dimension (nx_block*ny_block, nblyr+2), &
         intent(inout) :: &
         bSin          , & ! Ice salinity ppt (on bio  grid)
         brine_sal     , & ! brine salinity (ppt)
         brine_rho         ! brine density  (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block, nblyr), &
         intent(inout) :: &
         trcrn_S           ! salinity tracer ppt (on bio grid)

      real (kind=dbl_kind), dimension (nx_block,ny_block, nilyr), &
         intent(inout) :: &
         trcrn_q       , & ! enthalpy tracer 
         trcrn_Si          ! salinity on CICE grid

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         Rayleigh_criteria ! .true. if minimun ice thickness (Ra_c)  was reached 
      
      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         first_ice         ! for first category ice only .true. 
                           !initialized values should be used 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(out) :: &
         iDin          , & ! Diffusivity on the igrid   (1/s)
         ikin              ! permeability on the igrid 

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1), intent(inout) :: &
         iphin         , & ! porosity on the igrid 
         ibrine_rho    , & ! brine rho on interface  
         ibrine_sal        ! brine sal on interface   
         
      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop      ! indices of grid cell where code aborts
      
      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, zcells  , & !
         pcells,fcells,& ! horizontal index, combines i and j loops
         k, m, nint      ! vertical biology layer index 

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         zndxi, zndxj, zndxij, & ! compressed indices for icells with hbrin>thinS, hbr_old>thinS
         pndxi, pndxj, pndxij, & ! compressed indices for the rest of the icells
         fndxi, fndxj, fndxij    ! compressed indices for the rest of the icells
    
      real (kind=dbl_kind), dimension(icells) :: &
         surface_S        ! salinity of ice above hin > hbrin
      
      real (kind=dbl_kind), dimension(icells,2) :: &
         S_bot           

      real (kind=dbl_kind) :: &
         Tmlts, &             ! melting temperature
         dts                  ! local timestep (s)

      logical (kind=log_kind), dimension(icells) :: &
         Rayleigh
     
      real (kind=dbl_kind):: &
         Ttemp     ! initial temp profile on the CICE grid

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: &
         trtmp0       , &  ! temporary, remapped tracers     !need extra 
         trtmp            ! temporary, remapped tracers     !

! local parameters
          
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

       dts = dts_b
       nint = max(1,INT(dt/dts))  
       dts = dt/nint

       l_stop = .false.
       istop = 0
       jstop = 0
     
      !----------------------------------------------------------------
      ! Update boundary conditions
      !----------------------------------------------------------------
         
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxii(ij)
            j = indxjj(ij)    

            surface_S(ij) = min_salin

            Rayleigh(ij) = .true.
            if (n_cat == 1 .AND. hbr_old(ij) < Ra_c) then
               Rayleigh(ij) = Rayleigh_criteria(i,j) ! only category 1 ice can be false 
            endif

            if (dh_bot(i,j) + darcy_V(i,j)*dt > c0) then  
            
                 bSin(ij,nblyr+2) = sss(i,j)
                 bTin(i,j,nblyr+2) = sst(i,j)
                 brine_sal(ij,nblyr+2) = sss(i,j) 
                 brine_rho(ij,nblyr+2) = rhow
                 bphin(i,j,nblyr+2) = c1 
                 S_bot(ij,1) = c0
                 S_bot(ij,2) = c1  
   
      !-------------------------------
      ! bottom melt
      !--------------------------------

            else  
                 bSin(ij,nblyr+2) = bSin(ij,nblyr+1)  
                 Tmlts =  -bSin(ij,nblyr+2)* depressT 
                 bTin(i,j,nblyr+2) = bTin(i,j,nblyr+1)
                 bphin(i,j,nblyr+2)= iphin(ij,nblyr+1)
                 S_bot(ij,1) = c1
                 S_bot(ij,2) = c0
            endif

            if (abs(dh_top(i,j)) > puny .AND. abs(darcy_V(i,j)) > puny) then 
                 bSin(ij,1) = max(min_salin,-(brine_rho(ij,2)*brine_sal(ij,2)/rhosi * &
                       darcy_V(i,j)*dt - &
                       (dh_top(i,j) +darcy_V(i,j)*dt/bphi_min(ij))*min_salin)/dh_top(i,j))
                 brine_sal(ij,1) = brine_sal(ij,2)
                 brine_rho(ij,1) = brine_rho(ij,2)
                 bphin(i,j,1) = bphi_min(ij)
            else
                 bSin(ij,1) = min_salin
            endif
            
         enddo

      !----------------------------------------------------------------
      ! Solve for S using CICE T.  If solve_zsal = .true., then couple back
      ! to the thermodynamics
      !----------------------------------------------------------------

          call solve_S_dt (icells, nx_block,         ny_block,      &                    
                                      indxii,        indxjj,        &
                                      zcells, zndxi, zndxj, zndxij, &
                                      pcells, pndxi, pndxj, pndxij, &
                                      nint,   dts,   bSin,  bTin,   &
                                      aicen,  bphin, iphin, igrid,  &
                                      bgrid,  ikin,  hbr_old,       &
                                      hbrin,  hin,   hin_old,       &
                                      iDin,          darcy_V,       &
                                      brine_sal,     Rayleigh,      &
                                      first_ice,     sss,           &
                                      dt,     n_cat, dh_top,        &
                                      dh_bot,        brine_rho,     &
                                      ibrine_sal,    ibrine_rho,    &
                                      TLAT,          TLON,          &
                                      fzsaln,        fzsaln_g,      &
                                      istep1,        S_bot) 
  
     if (n_cat == 1)then
        do ij = 1, icells 
         i = indxii(ij)
         j = indxjj(ij)       
           Rayleigh_criteria(i,j) = Rayleigh(ij)
        enddo
     endif
         
     trtmp0(:,:,:) = c0
     trtmp (:,:,:) = c0
       
      do k = 1,nblyr                  !back to bulk quantity 
         do ij = 1, zcells 
            i = zndxi(ij)
            j = zndxj(ij)  
            m = zndxij(ij)  
            trcrn_S(i,j,k) =   bSin(m,k+1) 
            trtmp0(i,j,nt_sice+k-1) = trcrn_S(i,j,k)
         enddo
         do ij = 1,pcells
            i = pndxi(ij)
            j = pndxj(ij) 
            m = pndxij(ij) 
            trcrn_S(i,j,k) =   bSin(m,k+1) 
            trtmp0(i,j,nt_sice+k-1) = trcrn_S(i,j,k)
         enddo
       enddo           !ij
 
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
    do ij = 1, icells    
          i = indxii(ij)
          j = indxjj(ij)

          call remap_zbgc   (ntrcr, nilyr,   &
                             nt_sice,                 &
                             trtmp0(i,j,1:ntrcr),     &
                             trtmp(i,j,:),            &
                             1,         nblyr,        &
                             hin(ij),   hbrin(ij),    &
                             cgrid(2:nilyr+1),        &
                             bgrid(2:nblyr+1),        &
                             surface_S(ij))
    enddo

         do k = 1, nilyr

            do ij = 1, zcells
               i = zndxi(ij)
               j = zndxj(ij)      
               
               Tmlts = -trcrn_Si(i,j,k)*depressT
               Ttemp = min(-(min_salin+puny)*depressT, &
                       calculate_Tin_from_qin(trcrn_q(i,j,k),Tmlts)) 
               trcrn_Si(i,j,k) = min(-Ttemp/depressT, max(min_salin, &
                                 trtmp(i,j,nt_sice+k-1)))
               Tmlts = - trcrn_Si(i,j,k)*depressT 
               trcrn_q(i,j,k) = calculate_qin_from_Sin(Ttemp,Tmlts)             
          
            enddo !ij
            do ij = 1, pcells
               i = pndxi(ij)
               j = pndxj(ij)        
               
               Tmlts = -trcrn_Si(i,j,k)*depressT
               Ttemp = min(-(min_salin+puny)*depressT, &
                       calculate_Tin_from_qin(trcrn_q(i,j,k),Tmlts))
               trcrn_Si(i,j,k) = min(-Ttemp/depressT, max(min_salin, &
                                 trtmp(i,j,nt_sice+k-1)))
               Tmlts = -trcrn_Si(i,j,k)*depressT
               trcrn_q(i,j,k) = calculate_qin_from_Sin(Ttemp,Tmlts)             
          
            enddo !ij
         enddo !k

      end subroutine solve_zsalinity

!=======================================================================
!
!  solves salt continuity explicitly using 
!  Lax-Wendroff-type scheme (MacCormack)
!  (Mendez-Nunez and Carroll,  Monthly Weather Review, 1993)
!
! authors     Nicole Jeffery, LANL
!
      subroutine solve_S_dt(icells,   nx_block,      ny_block,      &                    
                                      indxi,         indxj,         &
                                      qcells, qndxi, qndxj, qndxij, &
                                      pcells, pndxi, pndxj, pndxij, &
                                      nint,   dts,   bSin,  bTin,   &
                                      aicen,  bphin, iphin, igrid,  &
                                      bgrid,  ikin,  hbri_old,      &
                                      hbrin,  hice,  hice_old,      &
                                      iDin,          darcy_V,       &
                                      brine_sal,     Rayleigh,      &
                                      first_ice,     sss,           &
                                      dt,     n_cat, dht,           &
                                      dhb,           brine_rho,     &
                                      ibrine_sal,    ibrine_rho,    &
                                      TLAT,          TLON,          &
                                      fzsaln,        fzsaln_g,      &
                                      istep1,        S_bot)    


      use ice_brine, only: calculate_drho

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block   , & ! block dimensions
         icells               , &
         n_cat,istep1         , &
         nint                     ! number of interations     
        
      integer (kind=int_kind), intent(out) :: &
         qcells, pcells
                
      integer (kind=int_kind), dimension(icells), intent(out) :: &
         qndxi, qndxj     , & ! compressed indices 
         qndxij           , &
         pndxi, pndxj     , & ! compressed indices 
         pndxij     

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj         ! compressed indices 

     real (kind=dbl_kind), intent(in) :: &
         dt            , & ! timestep (s)
         dts               ! local timestep (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         fzsaln        , & ! salt flux +ive to  ocean (kg/m^2/s)
         fzsaln_g          ! gravity drainage salt flux +ive to ocean (kg/m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         TLAT,TLON    , &
         aicen          

      logical (kind=log_kind), dimension (icells), &
         intent(inout) :: &
         Rayleigh   !if .true. convection is allowed.  if .false. not yet

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         first_ice

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2), intent(in) :: &
         brine_sal       , & ! Internal brine salinity (ppt)
         brine_rho           ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), intent(inout) :: &
         bphin             ! Porosity of layers

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1), intent(in) :: &
         ibrine_rho  , &    ! brine rho on interface 
         ibrine_sal         ! brine sal on interface 

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+1), intent(inout) :: &
         iphin          ! Porosity of layers on interface

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(out) :: &
         iDin        , &    ! Diffusivity on the igrid (1/s) with minimum bphi condition
         ikin               !permeability on interface

       real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         sss       , &   ! sea surface salinity
         dht       , &   ! change in the ice top  (positive for melting)
         dhb       , &   ! change in the ice bottom (positive for freezing)
         hice_old        ! old ice thickness (m)

       real (kind=dbl_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         hbri_old    , & ! brine thickness (m) 
         hbrin       , & ! new brine thickness (m)
         hice            ! ice thickness (m

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid           ! biology nondimensional grid layer points 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid           ! biology grid interface points 

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2), &
         intent(inout) :: &
         bSin            ! Bulk Salinity (ppt) contains previous timestep
                         ! and ocean ss

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+2), &
         intent(in) :: &
         bTin             ! Temperature of ice layers on bio grid for history file (^oC) 

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         darcy_V       ! Darcy velocity due to a pressure head (m/s) or melt      
  
      real (kind=dbl_kind), dimension (icells, 2), intent(in) :: &
         S_bot
      
      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k, m , mm   , & ! vertical biology layer index 
         ij_fault        ! location of conservation failure

      integer (kind=int_kind) :: &
         pij, qij            ! and h > thinS or h <= thinS

      real (kind=dbl_kind), dimension (icells,nblyr+1) :: &
         iDin_p          ! Diffusivity on the igrid (1/s)/bphi^3 

      real (kind=dbl_kind), dimension (icells,nblyr+2) :: &
         Din_p          ! Diffusivity on the igrid (1/s)/bphi^3  

      real (kind=dbl_kind), dimension (icells,nblyr+2) :: &
         Sintemp       ,& ! initial salinity
         pre_sin       ,&! estimate of  salinity of layers
         pre_sinb        ! estimate of  salinity of layers
 
      real (kind=dbl_kind), dimension (icells, nblyr+2) :: &
         bgrid_temp       ! biology nondimensional grid layer points 
                          ! with boundary values 
      real (kind=dbl_kind), dimension (icells) :: &
         dh       , &       ! (m) change in hbrine over dts
         dbgrid   , &       ! ratio of grid space to spacing across boundary 
                            ! i.e. 1/nilyr/(dbgrid(2)-dbgrid(1))
         S_surface          ! min_salin 

      real (kind=dbl_kind), dimension (nblyr+1) :: &  !on the zbgc vertical grid 2:nblyr+1
         dSbdx          ! gradient of brine rho on grid

      real (kind=dbl_kind), dimension (icells, nblyr+1) :: &  !on the zbgc vertical grid
         drho               ! brine difference rho_a-rho_b  (kg/m^3)

      real (kind=dbl_kind), dimension (icells, nblyr+2) :: &  !on the zbgc vertical grid 2:nblyr+1
         Q_s, C_s, &   ! Functions in continuity equation
         V_s, U_s, F_s

       real (kind=dbl_kind), dimension (icells, nblyr+1) :: &  !on the zbgc vertical igrid 1:nblyr+1
         Ci_s, &  !
         Ui_s, &  ! interface function
         Vi_s     ! for conservation check

      real (kind=dbl_kind), dimension (icells,nblyr) :: &
         vel, &               ! advective velocity times dt (m)
         lapidus_diff     , & ! lapidus term and 
         flux_corr

      real (kind=dbl_kind), dimension (icells,nblyr+1) :: &
         ivel

      real (kind=dbl_kind), dimension (icells):: &
         lapidus     ! artificial viscosity:  use lapidus_g for growth
   
      real (kind=dbl_kind), dimension (icells) :: &
         Ssum_old,Ssum_new, & ! depth integrated salt before and after timestep
         fluxcorr,          & ! flux correction to prevent S < min_salin
         Ssum_corr,        & ! numerical boundary flux correction
         fluxb, &  !bottom, top and total boundary flux (g/kg/m^2)
         fluxg, &  !bottom, top and total gravity drainage flux (g/kg/m^2)
         fluxm     !bottom, top and total molecular diffusion flux (g/kg/m^2)

      real (kind=dbl_kind) :: &
         sum_old,sum_new          , &  ! integrated salinity at t and t+dt
         dh_dt, dS_dt  
                        
      logical (kind=log_kind) :: &   
         write_flag       , &    ! set to true at each timestep        
         good_numerics    , &    ! output of check_conservation 
         stable,            &    ! if false, redo with smaller timestep
         test_conservation       ! test that salt change is balanced by fluxes 

       real  (kind=dbl_kind), dimension (nblyr):: &
          lapA    , &
          lapB    
   
     !--------------------------------------
     !  Initialize
     !--------------------------------------

       write_flag = .true.
       test_conservation = .false. 
       iDin_p(:,:) = c0   
       Din_p(:,:) = c0 
       lapA(:) = c1
       lapB(:) = c1
       lapA(nblyr) = c0
       lapB(1) = c0
       qcells = 0
       pcells = 0
       V_s(:,:) = c0
       U_s(:,:) = c0
       Q_s(:,:) = c0
       C_s(:,:) = c0
       Ci_s(:,:) = c0
       F_s(:,:) = c0
       Ui_s(:,:) = c0
       Vi_s(:,:) = c0
       ivel(:,:) = c0
       vel(:,:) = c0
       dh(:) = c0
       dbgrid(:) = c2
       S_surface(:) = min_salin

     !-----------------------------------------------------------------
     ! Find brine density gradient for gravity drainage parameterization
     !-----------------------------------------------------------------

         call calculate_drho(icells,nx_block,ny_block,indxi,indxj,igrid,bgrid,&
                                   brine_rho,ibrine_rho,drho)

     !-------------------------------------------------------------------------------
     ! Calculate bphi diffusivity on the grid points
     ! rhosi = 919-974 kg/m^2  set in bio_in
     ! rhow = 1026.0 density of sea water: uses kinematic viscosity (m^2/s) in Q18
     ! dynamic viscosity  divided by density = kinematic. 
     !------------------------------------------------------------------------------- 

        do k = 2, nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
           do ij = 1, icells
               iDin_p(ij,k) =k_o*gravit*l_skS/viscos_dynamic* &  
                     drho(ij,k)/(hbri_old(ij)**2)
 
               Din_p(ij,k) = (iDin_p(ij,k)*(igrid(k)-bgrid(k)) + iDin_p(ij,k-1)*&
                            (bgrid(k)-igrid(k-1)))/(igrid(k)-igrid(k-1))           
          enddo   !ij
         enddo                !k

     !---------------------------------------------------------------------------------------
     ! Critical Ra_c value is only for the onset of convection in thinS ice and not throughout
     !  therefore I need a flag to indicate the Ra_c was reached for a particular ice block
     ! Using a thickness minimum (Ra_c) for simplicity.
     !----------------------------------------------------------------------------------------    

        do ij = 1, icells
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          i = indxi(ij)
          j = indxj(ij) 
          bgrid_temp(ij,:) = bgrid(:)
          Din_p(ij,nblyr+2) =  iDin_p(ij,nblyr+1)
          if (.NOT. Rayleigh(ij) .AND. hbrin(ij) < Ra_c) then
                  Din_p(ij,:) =  c0  
                  iDin_p(ij,:) = c0  
          else
              Rayleigh(ij) = .true.
          endif
          if (hbri_old(ij) > thinS .AND. hbrin(ij) > thinS &
              .and.  hice_old(i,j) > thinS .AND. .NOT. first_ice(i,j)) then

               qcells = qcells + 1
               qndxi(qcells) = i
               qndxj(qcells) = j
               qndxij(qcells) = ij 
               bgrid_temp(ij,1) = c0
               bgrid_temp(ij,nblyr+2) = c1
               dbgrid(ij) = igrid(2)/(bgrid_temp(ij,2)-bgrid_temp(ij,1))

               !-----------------------------------
               ! surface boundary terms
               !-----------------------------------
                             
               lapidus(ij) = lapidus_g/real(nblyr,kind=dbl_kind)**2
               ivel(ij,1) = dht(i,j)/hbri_old(ij)
               U_s(ij,1) =  ivel(ij,1)/dt*dts 
               Ui_s(ij,1) = U_s(ij,1) 
               Ci_s(ij,1) = c0
               F_s(ij,1) =  brine_rho(ij,2)*brine_sal(ij,2)/rhosi*darcy_V(i,j)*dts/hbri_old(ij)/bSin(ij,1)

               !-----------------------------------
               ! bottom boundary terms
               !-----------------------------------

               ivel(ij,nblyr+1) =  dhb(i,j)/hbri_old(ij)           
               Ui_s(ij,nblyr+1) = ivel(ij,nblyr+1)/dt*dts  
               dSbdx(nblyr) = (ibrine_sal(ij,nblyr+1)*ibrine_rho(ij,nblyr+1) - &
                             ibrine_sal(ij,nblyr)*ibrine_rho(ij,nblyr))/&
                            (igrid(nblyr+1)-igrid(nblyr))	              
               C_s(ij,nblyr+1) = Dm/brine_sal(ij,nblyr+1)/brine_rho(ij,nblyr+1)*dts/hbri_old(ij)**2*&
	         (ibrine_sal(ij,nblyr+1)*ibrine_rho(ij,nblyr+1) - &
	         ibrine_sal(ij,nblyr)*ibrine_rho(ij,nblyr))/(igrid(nblyr+1)-igrid(nblyr))
               F_s(ij,nblyr+1) = darcy_V(i,j)*dts/hbri_old(ij)/bphin(i,j,nblyr+1)
               F_s(ij,nblyr+2) = darcy_V(i,j)*dts/hbri_old(ij)/bphin(i,j,nblyr+2)  
               vel(ij,nblyr) =(bgrid(nblyr+1)*(dhb(i,j)) -(bgrid(nblyr+1) - c1)* dht(i,j) )/hbri_old(ij)
               U_s(ij,nblyr+1) = vel(ij,nblyr)/dt*dts  
               V_s(ij,nblyr+1) = Din_p(ij,nblyr+1)/rhosi&
                       *(rhosi/brine_sal(ij,nblyr+1)/brine_rho(ij,nblyr+1))**exp_h&
                       *dts*dSbdx(nblyr) 
               dSbdx(nblyr+1) =  (brine_sal(ij,nblyr+2)*brine_rho(ij,nblyr+2) - &
                                brine_sal(ij,nblyr+1)*brine_rho(ij,nblyr+1))/&
                               (bgrid(nblyr+2)-bgrid(nblyr+1)+ grid_oS/hbri_old(ij) )  
               C_s(ij, nblyr+2) = Dm/brine_sal(ij,nblyr+2)/brine_rho(ij,nblyr+2)*dts/hbri_old(ij)**2*&
                              (brine_sal(ij,nblyr+2)*brine_rho(ij,nblyr+2) - & 
                              brine_sal(ij,nblyr+1)*brine_rho(ij,nblyr+1))/&
                              (bgrid(nblyr+2)-bgrid(nblyr+1)+ grid_oS/hbri_old(ij) ) 
               U_s(ij,nblyr+2) = ivel(ij,nblyr+1)/dt*dts 
               V_s(ij,nblyr+2) = Din_p(ij,nblyr+2)/rhosi &
                   *(bphin(i,j,nblyr+1)/bSin(ij,nblyr+2))**exp_h&
                   *dts*dSbdx(nblyr+1)
               Ci_s(ij,nblyr+1) = C_s(ij,nblyr+2)
               Vi_s(ij,nblyr+1) = V_s(ij,nblyr+2) 
               dh(ij) =(dhb(i,j)-dht(i,j))/dt*dts
          else                           !thinS ice solve or first_ice
               pcells = pcells + 1
               pndxi(pcells) =  i
               pndxj(pcells) = j
               pndxij(pcells) =  ij
           endif
      enddo

      if (qcells > 0) then
          do k = 2, nblyr  
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu  
           do ij = 1, qcells   ! full solve: 
             qij = qndxij(ij)
             i = qndxi(ij)
             j = qndxj(ij) 
              ivel(qij,k) =  (igrid(k)*dhb(i,j) - (igrid(k)-c1)*dht(i,j))/hbri_old(qij)
              Ui_s(qij,k) = ivel(qij,k)/dt*dts   
              Vi_s(qij,k) = iDin_p(qij,k)/rhosi &
                            *(rhosi/ibrine_rho(qij,k)/ibrine_sal(qij,k))**exp_h*dts*&
                            (brine_sal(qij,k+1)*brine_rho(qij,k+1)-&
                             brine_sal(qij,k)*brine_rho(qij,k))&
                            /(bgrid(k+1)-bgrid(k)) 
              dSbdx(k-1) = (ibrine_sal(qij,k)*ibrine_rho(qij,k) - &
                            ibrine_sal(qij,k-1)*ibrine_rho(qij,k-1))/(igrid(k)-igrid(k-1))
              F_s(qij,k) = darcy_V(i,j)*dts/hbri_old(qij)/bphin(i,j,k)
              C_s(qij,k) = Dm/brine_sal(qij,k)/brine_rho(qij,k)*dts/hbri_old(qij)**2*&
                         (ibrine_sal(qij,k)*ibrine_rho(qij,k) - &
                         ibrine_sal(qij,k-1)*ibrine_rho(qij,k-1))/(igrid(k)-igrid(k-1))
              Ci_s(qij,k) = Dm/ibrine_sal(qij,k)/ibrine_rho(qij,k)*dts/hbri_old(qij)**2*&
                          (brine_sal(qij,k+1)*brine_rho(qij,k+1)- &
                          brine_sal(qij,k)*brine_rho(qij,k))/(bgrid(k+1)-bgrid(k))
              vel(qij,k-1) = (bgrid(k)*(dhb(i,j)) - (bgrid(k) - c1)* dht(i,j))/hbri_old(qij)
              U_s(qij,k) = vel(qij,k-1)/dt*dts 
              V_s(qij,k) = Din_p(qij,k)/rhosi &
                       *(rhosi/brine_sal(qij,k)/brine_rho(qij,k))**exp_h&
                       *dts*dSbdx(k-1) 
              C_s(qij,2) = c0
              V_s(qij,2) = c0
          enddo  !ij
         enddo !k

      !-----------------------------------------------------------------
      ! Solve
      !--------------------------------------------------------  

       do m = 1, nint     
           lapidus_diff(:,:) = c0
           flux_corr(:,:) = c0
           Sintemp(:,:) = bSin(:,:)
           pre_sin(:,:) = bSin(:,:)  
           pre_sinb(:,:) = bSin(:,:)
           Ssum_old(:) = c0
           Ssum_new(:) = c0
           Ssum_corr(:) = c0
           fluxcorr(:) = c0
           fluxg(:) = c0
           fluxb(:) = c0
           fluxm(:) = c0
          do k = 2, nblyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
             do ij = 1, qcells
               qij = qndxij(ij) 
 !forward-difference 
               Ssum_old(qij) = Ssum_old(qij) + bSin(qij,k)*(igrid(k)-igrid(k-1))

               pre_sin(qij,k) =bSin(qij,k) + (Ui_s(qij,k)*(bSin(qij,k+1) - bSin(qij,k)) + &
                 V_s(qij,k+1)*bSin(qij,k+1)**3 - V_s(qij,k)*bSin(qij,k)**3 + &
                 (C_s(qij,k+1)+F_s(qij,k+1))*bSin(qij,k+1)-&
                 (C_s(qij,k)+F_s(qij,k))*bSin(qij,k))/(bgrid_temp(qij,k+1)-bgrid_temp(qij,k)) 

               pre_sin(qij,nblyr+1) = bSin(qij,nblyr+1) + (Ui_s(qij,nblyr+1)*(bSin(qij,nblyr+2) - &
                 bSin(qij,nblyr+1)) +  V_s(qij,nblyr+2)*bSin(qij,nblyr+2)**3 - &
                 V_s(qij,nblyr+1)*bSin(qij,nblyr+1)**3+ (C_s(qij,nblyr+2)+F_s(qij,nblyr+2))*&
                 bSin(qij,nblyr+2)-(C_s(qij,nblyr+1)+F_s(qij,nblyr+1))*bSin(qij,nblyr+1) )/&
                 (bgrid_temp(qij,nblyr+2)- bgrid_temp(qij,nblyr+1))
              
              enddo  !qcells
          enddo    !k

          do k = nblyr+1, 3, -1  !nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
              do ij = 1, qcells
                qij = qndxij(ij)
  !backward-difference 
                pre_sinb(qij,k) =p5*(bSin(qij,k) + pre_sin(qij,k) +  (Ui_s(qij,k-1)&
                  *(pre_sin(qij,k) - pre_sin(qij,k-1)) + &
                  V_s(qij,k)*pre_sin(qij,k)**3 - &
                  V_s(qij,k-1)*pre_sin(qij,k-1)**3 + (C_s(qij,k)+F_s(qij,k))*pre_sin(qij,k)-&
                  (C_s(qij,k-1)+F_s(qij,k-1))*pre_sin(qij,k-1))/(bgrid_temp(qij,k)-bgrid_temp(qij,k-1)) )

                 pre_sinb(qij,2) = p5*(bSin(qij,2) + pre_sin(qij,2) +  (Ui_s(qij,1)&
                  *(pre_sin(qij,2) -pre_sin(qij,1)) + &
                  V_s(qij,2)*pre_sin(qij,2)**3 - &
                  V_s(qij,1)*pre_sin(qij,1)**3 + (C_s(qij,2)+F_s(qij,2))*pre_sin(qij,2)-&
                  (C_s(qij,1)+F_s(qij,1))*pre_sin(qij,1) )/(bgrid_temp(qij,2)-bgrid_temp(qij,1)) )
              
                 Q_s(qij,k) = V_s(qij,k)*pre_sin(qij,k)**2 + U_s(qij,k) + C_s(qij,k) + F_s(qij,k) 
                 Q_s(qij,2) = V_s(qij,2)*pre_sin(qij,2)**2 + U_s(qij,2) + C_s(qij,2) + F_s(qij,2)
                 Q_s(qij,1) = V_s(qij,1)*pre_sin(qij,2)**2 + Ui_s(qij,1) + C_s(qij,1)+ F_s(qij,1)
                 Q_s(qij,nblyr+2) = V_s(qij,nblyr+2)*pre_sin(qij,nblyr+1)**2 + & 
                   Ui_s(qij,nblyr+1) + C_s(qij,nblyr+2) +  F_s(qij,nblyr+2)
             enddo  !qcells
          enddo   !k

        ! Add artificial viscosity   [Lapidus,1967] [Lohner et al, 1985]
        ! * more for melting ice
        !--------------------------------------------------------------------- 

         do k = 2, nblyr+1  
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
             do ij = 1, qcells
               i = qndxi(ij)
               j = qndxj(ij) 
               qij = qndxij(ij)
                  
               lapidus_diff(qij,k-1) =    lapidus(qij)/& ! lapidus(qij)/real(nblyr,kind=dbl_kind)**2/&
                  (igrid(k)-igrid(k-1))* &
                  ( lapA(k-1)*ABS(Q_s(qij,k+1)-Q_s(qij,k))*(abs(pre_sinb(qij,k+1))-abs(pre_sinb(qij,k)))/&
                  (bgrid_temp(qij,k+1)-bgrid_temp(qij,k) )**2 - &
                  lapB(k-1)*ABS(Q_s(qij,k)-Q_s(qij,k-1))*(abs(pre_sinb(qij,k))-abs(pre_sinb(qij,k-1)))/&
                  (bgrid_temp(qij,k)-bgrid_temp(qij,k-1))**2)
                          
               bSin(qij,k) = pre_sinb(qij,k) + lapidus_diff(qij,k-1)

               if (bSin(qij,k) < min_salin) then
                 flux_corr(qij,k-1) = min_salin - bSin(qij,k) !  flux into the ice
                 bSin(qij,k) = min_salin 
               elseif (bSin(qij,k) > -bTin(i,j,k)/depressT) then
                 flux_corr(qij,k-1) = bSin(qij,k)+bTin(i,j,k)/depressT !  flux into the ice
                 bSin(qij,k) = -bTin(i,j,k)/depressT
               elseif (bSin(qij,k) > max_salin) then
                 stable = .false.
                 write(nu_diag, *) 'Very Bad value in solve_dt_S-- istep1:',istep1
                 write(nu_diag, *) 'TLAT,TLON:',TLAT(i,j)*rad_to_deg,TLON(i,j)*rad_to_deg
                 write(nu_diag, *) 'bSin(qij,k),i,j,k',bSin(qij,k),i,j,k
                 write(nu_diag, *) 'V_s(qij,k),Ui_s(qij,k),U_s(qij,k),k',V_s(qij,k),&
                                   Ui_s(qij,k),U_s(qij,k),k
                 write(nu_diag, *) 'V_s(qij,k-1),Ui_s(qij,k-1),U_s(qij,k-1)',V_s(qij,k-1),&
                                   Ui_s(qij,k-1),U_s(qij,k-1)
                 write(nu_diag, *) 'V_s(qij,k+1),Ui_s(qij,k+1),U_s(qij,k+1)',V_s(qij,k+1),&
                                   Ui_s(qij,k+1),U_s(qij,k+1)
                 write(nu_diag, *) 'F_s(qij,k),F_s(qij,k+1),F_s(qij,k-1)', &
                                    F_s(qij,k),F_s(qij,k+1),F_s(qij,k-1)
                 write(nu_diag, *)'brine_rho(i,j,k),brine_sal(i,j,k)'
                 write(nu_diag,*)brine_rho(qij,k),brine_sal(qij,k)
                 write(nu_diag,*)'ivel(qij,k-1),ivel(qij,k),vel(qij,k-1),Din_p(qij,k),darcy_V(i,j)'
                 write(nu_diag,*)ivel(qij,k-1),ivel(qij,k),vel(qij,k-1),Din_p(qij,k),darcy_V(i,j)
                 write(nu_diag,*)'pre_sin(qij,k),pre_sin(qij,k-1),pre_sin(qij,k+1)'
                 write(nu_diag,*)pre_sin(qij,k),pre_sin(qij,k-1),pre_sin(qij,k+1)
                 write(nu_diag,*)'pre_sinb(qij,k),pre_sinb(qij,k-1),pre_sinb(qij,k+1)'
                 write(nu_diag,*)pre_sinb(qij,k),pre_sinb(qij,k-1),pre_sinb(qij,k+1)
                 write(nu_diag,*)'bTin(i,j,k),bTin(i,j,k-1)'
                 write(nu_diag,*)bTin(i,j,k),bTin(i,j,k-1)
                 write(nu_diag,*)'Q_s(qij,k),Q_s(qij,k+1),Q_s(qij,k-1)'
                 write(nu_diag,*)Q_s(qij,k),Q_s(qij,k+1),Q_s(qij,k-1)
                 write(nu_diag,*)'bgrid_temp(qij,k),bgrid_temp(qij,k+1),bgrid_temp(qij,k-1)'
                 write(nu_diag,*)bgrid_temp(qij,k),bgrid_temp(qij,k+1),bgrid_temp(qij,k-1)
                 write(nu_diag,*)'hbri_old(qij),hbrin(qij),dhb(i,j),dht(i,j),&
                                 dts, dt'
                 write(nu_diag,*)hbri_old(qij),hbrin(qij),dhb(i,j),dht(i,j),&
                                dts,dt
                 write(nu_diag,*)'darcy_V(i,j)'
                 write(nu_diag,*)darcy_V(i,j)
                 write(nu_diag,*)'bSin(qij,1),bSin(qij,nblyr+2)'
                 write(nu_diag,*)bSin(qij,1),bSin(qij,nblyr+2)
                 write(nu_diag,*)bSin(qij,1),bSin(qij,nblyr+2)
                 call abort_ice ('ice: Solve_S_dt error')                 
              endif
            
              if (k == nblyr+1) bSin(qij,nblyr+2) = S_bot(qij,1)*bSin(qij,nblyr+1) + &
                                                   S_bot(qij,2)*bSin(qij,nblyr+2) 

              Ssum_new(qij) = Ssum_new(qij) + bSin(qij,k)*(igrid(k)-igrid(k-1))
              fluxcorr(qij) = fluxcorr(qij) + (flux_corr(qij,k-1) + &
                               lapidus_diff(qij,k-1))*(igrid(k)-igrid(k-1))
        enddo  !qcells
       enddo   !k

       call calc_salt_fluxes (m,qcells,qndxi,qndxj,qndxij,icells,nx_block,ny_block,&
                             Ui_s,dh,dbgrid,hbri_old,Sintemp,pre_sin,fluxb,fluxg,fluxm,V_s,&
                             C_s,F_s,Ssum_corr,fzsaln_g,fzsaln,Ssum_old,fluxcorr,dts,S_surface)

       if (test_conservation) then
         good_numerics = .true.
         ij_fault = 0

         call check_conserve_salt(m,qcells,qndxij,icells,nx_block,ny_block, dt,dts,&
                                Ssum_old,Ssum_new,&
                                fluxcorr,dh,hbri_old,fluxb,fluxg,fluxm, Ssum_corr,&
                                good_numerics, ij_fault)

         if (.NOT. good_numerics) then
              i = qndxi(ij_fault)
              j = qndxj(ij_fault)
              qij = qndxij(ij_fault)
              write(nu_diag,*)',nint, Category,i,j,ij,TLAT,TLON-- istep1:'&
                              ,nint,n_cat,i,j,qij,TLAT(i,j)*rad_to_deg,&
                              TLON(i,j)*rad_to_deg,istep1
              write(nu_diag,*)'dhb(i,j),dht(i,j),Rayleigh(qij):'&
                              ,dhb(i,j),dht(i,j),Rayleigh(qij) 
              write(nu_diag,*) 'fzsaln(i,j),fzsaln_g(i,j):'
              write(nu_diag,*) fzsaln(i,j),fzsaln_g(i,j)
              write(nu_diag,*) 'Ssum_old(qij),Ssum_new(qij):'
              write(nu_diag,*) Ssum_old(qij),Ssum_new(qij)
              write(nu_diag,*) 'fluxb(qij):', fluxb(qij)
              write(nu_diag,*) 'fluxg(qij):',fluxg(qij)
              write(nu_diag,*) 'fluxm(qij):',fluxm(qij)
              write(nu_diag,*) 'bSin(qij,1),bSin(qij,nblyr+2)'
              write(nu_diag,*) bSin(qij,1),bSin(qij,nblyr+2)
              write(nu_diag,*) 'pre_sin(qij,1),pre_sin(qij,nblyr+2)'
              write(nu_diag,*) pre_sin(qij,1),pre_sin(qij,nblyr+2)
              do  k = 2,nblyr+1
                write(nu_diag, *) 'bSin(qij,k),i,j,k',bSin(qij,k),i,j,k
                write(nu_diag, *) 'V_s(qij,k),Ui_s(qij,k),U_s(qij,k),k',V_s(qij,k),&
                                  Ui_s(qij,k),U_s(qij,k),k
                write(nu_diag, *) 'V_s(qij,k-1),Ui_s(qij,k-1),U_s(qij,k-1)',V_s(qij,k-1),&
                                  Ui_s(qij,k-1),U_s(qij,k-1)
                write(nu_diag, *) 'V_s(qij,k+1),U_s(qij,k+1)',V_s(qij,k+1),&
                                  U_s(qij,k+1)
                write(nu_diag, *)'brine_rho(qij,k),brine_sal(qij,k)'
                write(nu_diag,*)brine_rho(qij,k),brine_sal(qij,k)
                write(nu_diag,*)'ivel(qij,k-1),ivel(qij,k),vel(qij,k-1),Din_p(qij,k),darcy_V(i,j)'
                write(nu_diag,*)ivel(qij,k-1),ivel(qij,k),vel(qij,k-1),Din_p(qij,k),darcy_V(i,j)
                write(nu_diag,*)'pre_sin(qij,k),pre_sin(qij,k-1),pre_sin(qij,k+1)'
                write(nu_diag,*)pre_sin(qij,k),pre_sin(qij,k-1),pre_sin(qij,k+1)
                write(nu_diag,*)'pre_sinb(qij,k),pre_sinb(qij,k-1),pre_sinb(qij,k+1)'
                write(nu_diag,*)pre_sinb(qij,k),pre_sinb(qij,k-1),pre_sinb(qij,k+1)
                write(nu_diag,*)'bTin(i,j,k),bTin(i,j,k-1)'
                write(nu_diag,*)bTin(i,j,k),bTin(i,j,k-1)
                write(nu_diag,*)'Q_s(qij,k),Q_s(qij,k+1),Q_s(qij,k-1)'
                write(nu_diag,*)Q_s(qij,k),Q_s(qij,k+1),Q_s(qij,k-1)
                write(nu_diag,*)'bgrid_temp(qij,1),bgrid_temp(qij,nblyr+2),bgrid_temp(qij,2)'
                write(nu_diag,*)bgrid_temp(qij,1),bgrid_temp(qij,nblyr+2),bgrid_temp(qij,2)
                write(nu_diag,*)'hbri_old(qij),hbrin(qij),dhb(i,j),dht(i,j),&
                                 dts, dt,dbgrid(qij)'
                write(nu_diag,*)hbri_old(qij),hbrin(qij),dhb(i,j),dht(i,j),&
                                dts,dt,dbgrid(qij)
              enddo
              call abort_ice ('ice: check_conserve_salt error')
          endif  !good_numerics
        endif  !test_conservation
       enddo !m
     endif   !qcells > 0 

     if (pcells > 0) then  !  add/melt ice only 

       do ij = 1, pcells
          i = pndxi(ij)
          j = pndxj(ij) 
          pij = pndxij(ij)

          sum_old = c0
          sum_new = c0
          dh_dt = hbrin(pij)-hbri_old(pij)
          dS_dt = c0
          if (dh_dt > c0) dS_dt = sss(i,j)*dh_dt/real(nblyr,kind=dbl_kind)
          do k = 2, nblyr+1 
            sum_old = sum_old + bSin(pij,k)*hbri_old(pij)*(igrid(k)-igrid(k-1))
            bSin(pij,k) = bSin(pij,k) + dS_dt
            sum_new = sum_new + bSin(pij,k)*hbrin(pij)*(igrid(k)-igrid(k-1))
          enddo  !k
          fzsaln(i,j) = fzsaln(i,j) -rhosi*(sum_new-sum_old)*p001/dt
          fzsaln_g(i,j) = c0
        enddo   !pcells
     endif    ! add/melt ice only

      !-----------------------------------------------------------------
      ! Move this to bgc calculation if using tr_salinity
      ! Calculate bphin, iphin, ikin, iDin and iDin_N 
      !-----------------------------------------------------------------

     iDin(:,:,:) = c0
     iphin(:,:)  = c1
     ikin(:,:,:) = c0    

       do k= 1, nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)  

              if (k < nblyr+1) bphin(i,j,k+1) = min(c1,max(puny, bSin(ij,k+1)*rhosi/&
                                             (brine_sal(ij,k+1)*brine_rho(ij,k+1)))) 
              if (k == 1) then     
                bphin(i,j,k) = min(c1,max(puny, bSin(ij,k)*rhosi/(brine_sal(ij,k)*brine_rho(ij,k))))  
                iphin(ij,k) = bphin(i,j,2)
              elseif (k == nblyr+1) then
                iphin(ij,nblyr+1) = bphin(i,j,nblyr+1)
              else
                iphin(ij,k) =min(c1, max(c0,(bphin(i,j,k+1) - bphin(i,j,k))/(bgrid(k+1) - &
                             bgrid(k))*(igrid(k)-bgrid(k)) + bphin(i,j,k)))
              endif 

              ikin(i,j,k) = k_o*iphin(ij,k)**exp_h 
           enddo  !ij
         enddo    !k

         do k = 2, nblyr+1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, qcells
               i = qndxi(ij)
               j = qndxj(ij) 
               qij = qndxij(ij)
               
               iDin(i,j,k) =  iphin(qij,k)*Dm/hbri_old(qij)**2  
               if (Rayleigh(qij) .AND. hbrin(qij) .GE. Ra_c) &
               iDin(i,j,k) =iDin(i,j,k) + l_sk*ikin(i,j,k)*gravit/viscos_dynamic* &  
                           drho(qij,k)/hbri_old(qij)**2 
            enddo  !ij

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, pcells
                i = pndxi(ij)
                j = pndxj(ij) 
                pij = pndxij(ij)
                iDin(i,j,k) = iphin(pij,k)*Dm/hbri_old(pij)**2

            enddo         !ij
        enddo       !k

770 format (I6,D16.6)        
781 format (I6,I6,I6)
790 format (I6,I6)
791 format (f24.17)
792 format (2D16.6)
793 format (3D16.6)
794 format (4D15.5)
800 format (F10.4)

     end subroutine solve_S_dt

!=======================================================================
!
! Calcuate salt fluxes
! 
      subroutine calc_salt_fluxes &
                                    (mint,qcells,qndxi,qndxj,qndxij,icells,&
                                     nx_block,ny_block,&
                                     Ui_s,dh,dbgrid,hbri_old,Sintemp,pre_sin,&
                                     fluxb,fluxg,fluxm,V_s,&
                                     C_s,F_s,Ssum_corr,fzsaln_g,fzsaln,Ssum_old,&
                                     fluxcorr,dts,S_surface)

      integer(kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         qcells, icells,mint      

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         qndxi, qndxj, qndxij ! compressed indices for icells with hin > thinS

      real (kind=dbl_kind), intent(in) :: &
         dts               !  halodynamic timesteps (s)

      real (kind=dbl_kind), dimension (icells,nblyr+2), intent(in) :: &
         Sintemp      , &  ! initial salinity
         pre_sin           ! estimate of  salinity of layers

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: &
         fzsaln        , & ! total salt flux  out of ice over timestep(kg/m^2/s)
         fzsaln_g          ! gravity drainage flux of salt  over timestep(kg/m^2/s)
                   
      real (kind=dbl_kind), dimension(icells), intent(in) :: &
         dh        , &  ! (m) change in hbrine over dts
         dbgrid    , &  ! ratio of grid space to spacing across boundary 
                        ! ie. 1/nilyr/(dbgrid(2)-dbgrid(1))
         S_surface      ! min_salin

      real (kind=dbl_kind), dimension(icells), intent(inout) :: &
         Ssum_old       , &  ! initial integrated salt content
         Ssum_corr      , &  ! boundary flux correction due to numerics
         fluxb          , &  ! total boundary salt flux into the ice (+ into the ice)
         fluxg          , &  ! total gravity drainage salt flux into the ice (+ into the ice)
         fluxm               ! total molecular diffusive salt flux into the ice (+ into the ice)

      real (kind=dbl_kind), dimension(icells), intent(in) :: &
         fluxcorr        ! flux correction to ensure S >= min_salin

      real (kind=dbl_kind), dimension (icells, nblyr+1), intent(in) :: & 
         Ui_s ! interface function

      real (kind=dbl_kind), dimension (icells, nblyr+2), intent(in) :: & 
         C_s, &   ! Functions in continuity equation
         F_s, &
         V_s

      real (kind=dbl_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         hbri_old          ! (m) initial brine height

     ! local variables

     integer (kind=int_kind) :: &
         i, j, ij, qij   ! vertical and horizontal indices

      real (kind=dbl_kind), dimension (icells) :: &
         Ssum_corr_flux,         & ! numerical boundary flux correction
         fluxb_b, fluxb_t, &  !bottom, top and total boundary flux (g/kg/m^2)
         fluxg_b, fluxg_t, &  !bottom, top and total gravity drainage flux (g/kg/m^2)
         fluxm_b, fluxm_t    !bottom, top and total molecular diffusion flux (g/kg/m^2)

      real (kind=dbl_kind) :: hin_next

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
           do ij = 1, qcells
               i = qndxi(ij)
               j = qndxj(ij) 
               qij = qndxij(ij)

               Ssum_old(qij) = Ssum_old(qij) + Sintemp(qij,nblyr+1)* &
                                 (igrid(nblyr+1)-igrid(nblyr))

     !----------------------------------------------------------------------------
     ! boundary fluxes (positive into the ice)
     !----------------------------------------------------------------------------

               fluxb_b(qij) = p5*Ui_s(qij,nblyr+1)*((c1-dh(qij)/hbri_old(qij))* &
                              Sintemp(qij,nblyr+2)*dbgrid(qij) + pre_sin(qij,nblyr+1)+&
                              (c1-dh(qij)/hbri_old(qij))*(c1-dbgrid(qij))*Sintemp(qij,nblyr+1)) + &
                              p5*((c1-dh(qij)/hbri_old(qij))*dbgrid(qij)*F_s(qij,nblyr+2)* &
                              Sintemp(qij,nblyr+2)  +  F_s(qij,nblyr+1)*(pre_sin(qij,nblyr+1) - &
                              (c1-dh(qij)/hbri_old(qij))*(dbgrid(qij)-c1)*Sintemp(qij,nblyr+1)))

               fluxb_t(qij) = -p5*Ui_s(qij,1)*(pre_sin(qij,1)*dbgrid(qij) +  &
                              (c1-dh(qij)/hbri_old(qij))*Sintemp(qij,2) -  &
                              (dbgrid(qij)-c1)*pre_sin(qij,2)) + &
                              -p5*(dbgrid(qij)*F_s(qij,1)*pre_sin(qij,1) + &
                              F_s(qij,2)*((c1-dh(qij)/hbri_old(qij))*Sintemp(qij,2) &
                              +(c1-dbgrid(qij))*pre_sin(qij,2)))

               fluxb(qij) = fluxb_b(qij) + fluxb_t(qij)

    !----------------------------------------------------------------------------------------
    ! gravity drainage fluxes (positive into the ice)
    !----------------------------------------------------------------------------------------

               fluxg_b(qij) =  p5*((c1-dh(qij)/hbri_old(qij))* dbgrid(qij)* &
                               V_s(qij,nblyr+2)*Sintemp(qij,nblyr+1)**3  +  &
                               V_s(qij,nblyr+1)*(pre_sin(qij,nblyr+1)**3 - &
                               (c1-dh(qij)/hbri_old(qij))*(dbgrid(qij) - c1)* &
                               Sintemp(qij,nblyr+1)**3))
                
               fluxg_t(qij) =  -p5*(dbgrid(qij)*V_s(qij,1)*pre_sin(qij,1)**3 + &
                               V_s(qij,2)*((c1-dh(qij)/hbri_old(qij))*Sintemp(qij,2)**3- &
                               (dbgrid(qij)-c1)*pre_sin(qij,2)**3))
                
               fluxg(qij) =  fluxg_b(qij) + fluxg_t(qij)

     !------------------------------------------------------------------------------------
     ! diffusion fluxes (positive into the ice)
     !----------------------------------------------------------------------------------

               fluxm_b(qij) = p5*((c1-dh(qij)/hbri_old(qij))*dbgrid(qij)*C_s(qij,nblyr+2)* &
                              Sintemp(qij,nblyr+2)  +  C_s(qij,nblyr+1)*(pre_sin(qij,nblyr+1) - &
                              (c1-dh(qij)/hbri_old(qij))*(dbgrid(qij)-c1)*Sintemp(qij,nblyr+1)))

	       fluxm_t(qij) =-p5*(dbgrid(qij)*C_s(qij,1)*pre_sin(qij,1) + &
                              C_s(qij,2)*((c1-dh(qij)/hbri_old(qij))*Sintemp(qij,2)+ &
                              (c1-dbgrid(qij))*pre_sin(qij,2)))
           
               fluxm(qij) =  fluxm_b(qij) + fluxm_t(qij)
              

               Ssum_corr(qij) = (-dh(qij)/hbri_old(qij) + p5*(dh(qij)/hbri_old(qij))**2)*&
                                Ssum_old(qij) 
               hin_next = hbri_old(qij) + real(mint,kind=dbl_kind)*dh(qij)
               Ssum_corr_flux(qij) =  dh(qij)*Ssum_old(qij)/hin_next +Ssum_corr(qij)
          
               fzsaln_g(i,j) = fzsaln_g(i,j) - hin_next*fluxg(qij)*rhosi*p001/dts  

               fzsaln(i,j) = fzsaln(i,j) - hin_next*(fluxb(qij) + &
                             fluxg(qij) + fluxm(qij)+ Ssum_corr_flux(qij) + &
                             fluxcorr(qij))*rhosi*p001/dts

       enddo

     end subroutine calc_salt_fluxes

!=======================================================================
! 
! Test salt conservation:   flux conservative form dSin/dt = -dF(x,Sin)/dx 
!  
      subroutine check_conserve_salt (mint,qcells, qndxij, icells,nx_block,&
                                      ny_block, dt,dts,&
                                      Ssum_old,Ssum_new,&
                                      fluxcorr,dh,hbri_old,fluxb,fluxg, &
                                      fluxm,Ssum_corr,& 
                                      good_numerics, ij_fault) 

      integer(kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         qcells, icells,mint      

      integer(kind=int_kind), intent(inout) :: &
         ij_fault          !ij value where conservation fails 

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         qndxij            ! compressed index for icells with hin > thinS

      real (kind=dbl_kind), intent(in) :: &
         dt, dts       ! thermodynamic and halodynamic timesteps (s)

      real (kind=dbl_kind), dimension(icells), intent(in) :: &
         dh       ! (m) change in hbrine over dts

      real (kind=dbl_kind), dimension(icells), intent(in) :: &
         Ssum_old       , &  ! initial integrated salt content
         Ssum_new       , &  ! final integrated salt content
         fluxcorr       , &  ! flux correction to ensure S >= min_salin
         Ssum_corr      , &  ! boundary flux correction due to numerics
         fluxb          , &  ! total boundary salt flux into the ice (+ into the ice)
         fluxg          , &  ! total gravity drainage salt flux into the ice (+ into the ice)
         fluxm              ! total gravity drainage salt flux into the ice (+ into the ice)

      logical (kind=log_kind), intent(inout) :: &   
         good_numerics          ! true if conservation satisfied within error

      real (kind=dbl_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         hbri_old          ! (m) initial brine height

     ! local variables

     integer (kind=int_kind) :: &
         ij, qij   ! vertical and horizontal indices

     real (kind=dbl_kind):: &
         diff2 , & !
         dsum  , & ! salt for conservation (g/kg/m^3)
         dsum_flux , & ! salt change in kg/m^2/s
         flux_tot  , &          ! fluxg+fluxb
         fzsaln_old, &  !
         order

     real (kind=dbl_kind), parameter :: &
         accuracy = 1.0e-14_dbl_kind ! g/kg/m^3 difference between boundary fluxes 

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
        do ij = 1,qcells
         qij = qndxij(ij)
         dsum = Ssum_new(qij) -  Ssum_old(qij) 
         flux_tot = fluxb(qij) + fluxg(qij) + fluxm(qij) + fluxcorr(qij) 
         dsum_flux =(Ssum_new(qij)*(hbri_old(qij) + (real(mint,kind=dbl_kind))*dh(qij)) - &
                     Ssum_old(qij)*(hbri_old(qij) + (real(mint,kind=dbl_kind)-c1)* &
                     dh(qij)) )*rhosi*p001/dts 
         order = abs(dh(qij)/hbri_old(qij))
         if (abs(dsum) > c0) then
           diff2 = abs(dsum- Ssum_corr(qij)- flux_tot)
           if (diff2 >  puny .AND. diff2 > order ) then 
              ij_fault = ij
              good_numerics = .false.
              write(nu_diag,*) 'Poor salt conservation: check_conserve_salt'
              write(nu_diag,*) 'qij, mint:', qij,mint
              write(nu_diag,*) 'Ssum_corr(qij)',Ssum_corr(qij)
              write(nu_diag,*) 'fluxb(qij),fluxg(qij),fluxm(qij),flux_tot,fluxcorr(qij):'
              write(nu_diag,*) fluxb(qij),fluxg(qij),fluxm(qij),flux_tot,fluxcorr(qij)
              write(nu_diag,*) 'fluxg(qij),',fluxg(qij)
              write(nu_diag,*) 'dsum,dsum_flux,',dsum,dsum_flux
              write(nu_diag,*) 'Ssum_new(qij),Ssum_old(qij),hbri_old(qij),dh(qij):'
              write(nu_diag,*) Ssum_new(qij),Ssum_old(qij),hbri_old(qij),dh(qij)
              write(nu_diag,*) 'diff2,order,puny',diff2,order,puny
           endif
         endif
        enddo

     end subroutine check_conserve_salt

!=======================================================================
!
! Aggregate flux information from all ice thickness categories
!
      subroutine merge_S_fluxes (nx_block, ny_block, &
                               icells,               &
                               indxi,      indxj,    &
                               aicen,                &
                               zsal_totn,  zsal_tot, & 
                               dsnown,     dsnow,    &
                               fzsal,      fzsaln,   &
                               fzsal_g,    fzsaln_g)

      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block, & ! block dimensions
          icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
          intent(in) :: &
          indxi, indxj    ! compressed indices for cells with aicen > puny

      ! single category fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &          
          aicen   , & ! concentration of ice
          dsnown  , & ! snow  growth                    (m)
          fzsaln  , & ! salt flux                       (kg/m**2/s)
          fzsaln_g    ! gravity drainage salt flux      (kg/m**2/s)

      real (kind=dbl_kind), dimension(nx_block*ny_block), intent(in):: &
          zsal_totn   ! tot salinity in category (psu*volume*rhosi)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &          
          zsal_tot, & ! tot salinity (psu*rhosi*total vol ice)
          dsnow   , & ! snow growth                     (m)
          fzsal   , & ! salt flux                       (kg/m**2/s)
          fzsal_g     ! gravity drainage salt flux      (kg/m**2/s)

      ! local variables

      integer (kind=int_kind) :: &
          ij, i, j, &   ! horizontal indices
          k             ! tracer indice

      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         zsal_tot(i,j)  = zsal_tot   (i,j) + zsal_totn  (ij)  !already *aicen

         ! ocean tot and gravity drainage salt fluxes
         fzsal   (i,j) = fzsal     (i,j) + fzsaln  (i,j)*aicen(i,j)
         fzsal_g (i,j) = fzsal_g   (i,j) + fzsaln_g(i,j)*aicen(i,j)

         ! ice/snow thickness
         dsnow    (i,j) = dsnow    (i,j) + dsnown  (i,j)*aicen(i,j)
      enddo                     ! ij

      end subroutine merge_S_fluxes

!==========================================================================
!
! For each grid cell, sum field over all ice layers.  "Net" refers to  the column
! integration while "avg"  is normalized by the ice thickness

      subroutine column_sum_S (nx_block, ny_block,       &
                             icells,   indxi,   indxj,   &
                             ntrcr, nblyr, vicen, trcrn, &
                             zsal_totn)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         nblyr, ntrcr              , & ! number of layers/tracers
         icells                 ! number of ice/ocean grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj          ! compressed i/j indices

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &          
         vicen       !volume of ice

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
         intent(in) :: &
         trcrn         ! input fields

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(out) :: &
          zsal_totn  ! avg salinity (psu*rhosi*vol of ice)

      ! local variables

      integer (kind=int_kind) :: &
           i, j, ij     , & ! horizontal indices
           k                ! layer index

     zsal_totn(:) = c0
      do k = 1, nblyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            zsal_totn(ij) = zsal_totn(ij) + rhosi*trcrn(i,j,nt_bgc_S+k-1) &
                         * vicen(i,j)*trcrn(i,j,nt_fbri)/real(nblyr,kind=dbl_kind)
                             
         enddo                  ! ij
      enddo                     ! n

      end subroutine column_sum_S

!=======================================================================
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW
!          Nicole Jeffery, LANL

      subroutine S_diags (dt)

      use ice_blocks, only: nx_block, ny_block
      use ice_broadcast, only: broadcast_scalar
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, pbloc, plat, plon
      use ice_domain_size, only: max_blocks
      use ice_grid, only: lmask_n, lmask_s, tarean, tareas, grid_type

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, nn, ii,jj, iblk

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         phinS, phinS1, phinS2,&
         phbrn,pdh_top1,pdh_bot1, pdh_top2,pdh_bot2, psice_rho, pfzsal, & 
         pfzsal_g, pdarcy_V1, pdarcy_V2 

      ! vertical  fields of category 1 at diagnostic points for bgc layer model
      real (kind=dbl_kind), dimension(npnt,nblyr+2) :: &
         pphin, pgrid, pphin1
      real (kind=dbl_kind), dimension(npnt,nblyr) :: &
        pSin, pSice, pSin1, pSin2

      real (kind=dbl_kind), dimension(npnt,nblyr+1) :: &
         pbTiz, piDin

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1, work2

      !-----------------------------------------------------------------
      ! salinity and microstructure  of the ice
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

               pfzsal(n) = fzsal(i,j,iblk)   
               pfzsal_g(n) = fzsal_g(i,j,iblk)            
               phinS(n) = c0            
               phinS1(n) = c0            
               phinS2(n) = c0
               phbrn(n) = c0
               psice_rho(n) = c0
               pdh_top1(n) = c0
               pdh_bot1(n) = c0
               pdh_top2(n) = c0
               pdh_bot2(n) = c0
               pdarcy_V1(n) = c0
               pdarcy_V2(n) = c0
               do nn = 1,ncat
                  psice_rho(n) = psice_rho(n) + sice_rho(i,j,nn,iblk)*aicen(i,j,nn,iblk)
               enddo
               if (aice(i,j,iblk) > c0) then
                 psice_rho(n) = psice_rho(n)/aice(i,j,iblk)
               endif
                if (tr_brine .and. aice(i,j,iblk) > c0) &
                   phinS(n) = trcr(i,j,nt_fbri,iblk)*vice(i,j,iblk)/aice(i,j,iblk)

                if (aicen(i,j,1,iblk)> c0) then
                  if (tr_brine) phinS1(n) = trcrn(i,j,nt_fbri,1,iblk)*vicen(i,j,1,iblk)/&
                                                aicen(i,j,1,iblk)
                  pdh_top1(n) = dhbr_top(i,j,1,iblk)
                  pdh_bot1(n) = dhbr_bot(i,j,1,iblk)
                  pdarcy_V1(n) = darcy_V(i,j,1,iblk)
                endif  
                if (aicen(i,j,2,iblk)> c0) then
                  if (tr_brine) phinS2(n) = trcrn(i,j,nt_fbri,2,iblk) *vicen(i,j,2,iblk)/&
                                                    aicen(i,j,2,iblk)
                  pdh_top2(n) = dhbr_top(i,j,2,iblk)
                  pdh_bot2(n) = dhbr_bot(i,j,2,iblk)
                  pdarcy_V2(n) = darcy_V(i,j,2,iblk)
                endif
                if (tr_brine .AND. aice(i,j,iblk) > c0) phbrn(n) = vice(i,j,iblk)/aice(i,j,iblk) - &
                       rhosi/rhow*vice(i,j,iblk)/aice(i,j,iblk) - & 
                       rhos/rhow*vsno(i,j,iblk)/aice(i,j,iblk)
                 do k = 1, nblyr+1
                   pbTiz(n,k) = c0
                   piDin(n,k) = c0
                  do nn = 1,ncat
                   pbTiz(n,k) = pbTiz(n,k) + bTiz(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                   piDin(n,k) = piDin(n,k) +  iDi(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                  enddo
                   if (vice(i,j,iblk) > c0)then
                    pbTiz(n,k) = pbTiz(n,k)/vice(i,j,iblk)
                    piDin(n,k) = piDin(n,k)/vice(i,j,iblk) 
                   endif
                 enddo                 !k
                 do k = 1, nblyr+2
                   pphin(n,k) = c0
                   pphin1(n,k) = c0
                   if (aicen(i,j,1,iblk) > c0) pphin1(n,k) = bphi(i,j,k,1,iblk)
                   do nn = 1,ncat
                     pphin(n,k) = pphin(n,k) + bphi(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                   enddo
                   if (vice(i,j,iblk) > c0) then
                    pphin(n,k) = pphin(n,k)/vice(i,j,iblk)
                   endif
                 enddo
                 do k = 1,nblyr
                       pSin(n,k) = c0
                       pSin1(n,k) = c0
                       pSin2(n,k) = c0                      
                       pSin(n,k)= trcr(i,j,nt_bgc_S+k-1,iblk)     
                       if (aicen(i,j,1,iblk) > c0) pSin1(n,k) = trcrn(i,j,nt_bgc_S+k-1,1,iblk)      
                       if (aicen(i,j,2,iblk)> c0) pSin2(n,k)  = trcrn(i,j,nt_bgc_S+k-1,2,iblk)
                 enddo 
                 do k = 1,nilyr
                   pSice(n,k) = trcr(i,j,nt_sice+k-1,iblk)   
                 enddo
            endif                 ! my_task = pmloc

           call broadcast_scalar(phinS(n), pmloc(n)) 
           call broadcast_scalar(phinS1(n), pmloc(n)) 
           call broadcast_scalar(phinS2(n), pmloc(n)) 
           call broadcast_scalar(phbrn(n), pmloc(n)) 
           call broadcast_scalar(pdh_top1(n), pmloc(n)) 
           call broadcast_scalar(pdh_bot1(n), pmloc(n)) 
           call broadcast_scalar(pdh_top2(n), pmloc(n)) 
           call broadcast_scalar(pdh_bot2(n), pmloc(n)) 
           call broadcast_scalar(psice_rho(n), pmloc(n))  
           call broadcast_scalar(pfzsal_g(n), pmloc(n))  
           call broadcast_scalar(pdarcy_V1(n), pmloc(n))  
           call broadcast_scalar(pdarcy_V2(n), pmloc(n)) 
           call broadcast_scalar(pfzsal(n), pmloc(n)) 

           do k = 1,nblyr+1
              call broadcast_scalar(pbTiz (n,k), pmloc(n))
              call broadcast_scalar(piDin (n,k), pmloc(n))
              call broadcast_scalar(pphin (n,k), pmloc(n))
              call broadcast_scalar(pphin1 (n,k), pmloc(n))
           enddo
           do k = 1,nblyr
              call broadcast_scalar(pSin (n,k), pmloc(n))
              call broadcast_scalar(pSin1 (n,k), pmloc(n))
              call broadcast_scalar(pSin2 (n,k), pmloc(n))
           enddo

           do k = 1,nilyr
             call broadcast_scalar(pSice (n,k), pmloc(n))
           enddo
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
        write(nu_diag,902) '      Brine height       '
        write(nu_diag,900) 'hbrin                   = ',phinS(1),phinS(2)
        write(nu_diag,900) 'hbrin cat 1             = ',phinS1(1),phinS1(2)
        write(nu_diag,900) 'hbrin cat 2             = ',phinS2(1),phinS2(2)
        write(nu_diag,900) 'Freeboard               = ',phbrn(1),phbrn(2)
        write(nu_diag,900) 'dhbrin cat 1 top        = ',pdh_top1(1),pdh_top1(2)
        write(nu_diag,900) 'dhbrin cat 1 bottom     = ',pdh_bot1(1),pdh_bot1(2)
        write(nu_diag,900) 'dhbrin cat 2 top        = ',pdh_top2(1),pdh_top2(2)
        write(nu_diag,900) 'dhbrin cat 2 bottom     = ',pdh_bot2(1),pdh_bot2(2)
        write(nu_diag,*) '                         '
        write(nu_diag,902) '     zSalinity         '
        write(nu_diag,900) 'Avg density (kg/m^3)   = ',psice_rho(1),psice_rho(2)
        write(nu_diag,900) 'Salt flux (kg/m^2/s)   = ',pfzsal(1),pfzsal(2)
        write(nu_diag,900) 'Grav. Drain. Salt flux = ',pfzsal_g(1),pfzsal_g(2)
        write(nu_diag,900) 'Darcy V cat 1 (m/s)    = ',pdarcy_V1(1),pdarcy_V1(2)
        write(nu_diag,900) 'Darcy V cat 2 (m/s)    = ',pdarcy_V2(1),pdarcy_V2(2)
        write(nu_diag,*) '                         '
        write(nu_diag,*) ' Top down bgc Layer Model'
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'bTiz(1) ice temp',' bTiz(2) ice temp  '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pbTiz(n,k),n = 1,2), k = 1,nblyr+1)
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'iDi(1) diffusivity  ','iDi(2) diffusivity  '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((piDin(n,k),n=1,2), k = 1,nblyr+1)
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'bphi(1) porosity   ','bphi(2) porosity  '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pphin(n,k),n=1,2), k = 1,nblyr+1)
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'phi1(1) porosity   ','phi1(2) porosity  '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pphin1(n,k),n=1,2), k = 1,nblyr+1)
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'Sin1(1) cat 1  ','Sin1(2) cat 1 '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pSin1(n,k),n=1,2), k = 1,nblyr)                         
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'Sin2(1) cat 2   ','Sin(2) cat 2'
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pSin2(n,k),n=1,2), k = 1,nblyr)                         
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'Sin(1) Avg S   ','Sin(2) Avg S '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pSin(n,k),n=1,2), k = 1,nblyr)            
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'Sice(1) Ice S   ','Sice(2) Ice S '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pSice(n,k),n=1,2), k = 1,nilyr)            
        write(nu_diag,*) '                         '
       endif                    ! print_points
      endif                     ! my_task = master_task

  802 format (f24.17,2x,f24.17)
  803 format (a25,2x,a25)
  900 format (a25,2x,f24.17,2x,f24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
  903 format (a25,5x,i4,1x,i4,1x,i4,1x,i4,7x,i4,1x,i4,1x,i4,1x,i4)

      end subroutine S_diags

!=======================================================================
!
! Dumps all values needed for a zsalinity  restart
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_S()

      use ice_domain_size, only: ncat
      use ice_domain, only: nblocks
      use ice_state, only: trcrn
      use ice_flux, only: sss  
      use ice_restart, only:  write_restart_field

      ! local variables

      integer (kind=int_kind) :: &
       i, j,  k, iblk  ! indices

      logical (kind=log_kind) :: diag

      character (len=3) :: nchar

      diag = .true.

      !------------------------------------------------------------------------------
      ! Salinity and extras
      !------------------------------------------------------------------------------

      do k = 1,nblyr
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_S,0,trcrn(:,:,nt_bgc_S+k-1,:,:),'ruf8', &
                   'zSalinity'//trim(nchar),ncat,diag)
      enddo
    
      call write_restart_field(nu_dump_S,0,sss,'ruf8','sss',1,diag)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (Rayleigh_criteria(i,j,iblk)) then
               Rayleigh_real(i,j,iblk) = c1
            elseif (.NOT. Rayleigh_criteria(i,j,iblk)) then
               Rayleigh_real(i,j,iblk) = c0
            endif
         enddo
         enddo
      enddo

      call write_restart_field(nu_dump_S,0,Rayleigh_real,'ruf8','Rayleigh',1,diag)

      end subroutine write_restart_S

!=======================================================================
!
! Reads all values needed for a zsalinity restart
!
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_S()

      use ice_domain_size, only: ncat
      use ice_domain, only: nblocks
      use ice_flux, only: sss  
      use ice_state, only: trcrn
      use ice_restart, only: read_restart_field
    !  use ice_constants, only: field_loc_center, field_type_scalar

      ! local variables

      integer (kind=int_kind) :: &
          i, j, k, iblk ! counting indices
      logical (kind=log_kind) :: diag
      character (len=3) :: nchar

      diag = .true.

      if (my_task == master_task) write(nu_diag,*)'zSalinity restart'
      do k = 1,nblyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_S,0,trcrn(:,:,nt_bgc_S+k-1,:,:),'ruf8', &
              'zSalinity'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
 
     if (my_task == master_task) write(nu_diag,*) 'sea surface salinity'
     call read_restart_field(nu_restart_S,0,sss,'ruf8','sss',1,diag)
     call read_restart_field(nu_restart_S,0,Rayleigh_real,'ruf8','Rayleigh',1,diag)
      
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (Rayleigh_real(i,j,iblk) .GE. c1) then
               Rayleigh_criteria (i,j,iblk) = .true.
            elseif (Rayleigh_real(i,j,iblk) < c1) then
               Rayleigh_criteria (i,j,iblk) = .false.
            endif
         enddo
         enddo
      enddo

      end subroutine read_restart_S

!=======================================================================

      end module ice_zsalinity

!=======================================================================
