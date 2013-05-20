!=======================================================================
!
!BOP
!
! !MODULE: ice_zsalinity - prognostic salinity using B&L thermodynamics
!
! !DESCRIPTION:
!
!  Halodynamics based on Jeffery and Hunke, 2013, the CICE sea ice model
!  with prognostic salinity: Arctic simulations, JGR,  submitted
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors: Nicole Jeffery, LANL
!          Elizabeth C. Hunke, LANL
!
! Vertical salinity (trcrn(nt_bgc_S)) is solved on the bio grid (bgrid and igrid)
! with domain defined by the dynamic brine height (trcrn(nt_fbri) * vicen/aicen).
! The CICE Bitz and Lipscomb thermodynamics is solve in the cgrid with height
! vicen/aicen.
! Gravity drainage is parameterized as nonlinear advection
! Flushing is incorporated in the boundary changes. 
! (see Jeffery et al., JGR, 2011 and 2012).  
!
!
! !INTERFACE:
!
      module ice_zsalinity
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_domain_size
      use ice_fileunits, only: nu_diag, nu_dump_S, nu_restart_S, &
           nu_rst_pointer, nu_dump_S, flush_fileunit
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_state
      use ice_zbgc_public
      use ice_blocks, only: nx_block, ny_block
!
!EOP
!
      implicit none

      private
      public :: init_zsalinity, first_ice, S_diags, write_restart_S, solve_zsalinity, &
           column_sum_S, merge_S_fluxes

      integer (kind=int_kind), parameter :: &
         restart_n  = 7       ! value  of nblyr in restart files

      real (kind=dbl_kind), parameter :: & 
         viscos_dynamic = 2.2_dbl_kind, & !1.8e-3_dbl_kind (pure water at 0^oC) (kg/m/s)
         max_salin = 80.0_dbl_kind  ! (ppt) maximum bulk salinity

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         Rayleigh_real , &  ! .true. = c1, .false. = c0
         first_ice_real     ! .true. = c1, .false. = c0
       
    real (kind=dbl_kind), parameter :: & 
         dts_b   = 120_dbl_kind, &!
         dts_s   = c10  !Not used 
        
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_zsalinity
!
! !DESCRIPTION:
!
!  Initialize vertical profile of biogeochemistry
!  
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_zsalinity (sss)
!
! !USES:
!
      use ice_domain, only: nblocks
!      use ice_zbgc_public, only: zTin, ocean_bio, cgrid, &
!                           bgrid, igrid, zphi, iDi, iki
     ! use ice_exit
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         sss      ! sea surface salinity (ppt)
!
!
!EOP
!
      integer (kind=int_kind) :: &
           i, j, iblk       , & ! horizontal indices
           ij               , & ! horizontal index, combines i and j loops
           l                , & ! month index
           k,m               , & ! vertical index
           n                , & ! category index
           nbits

      real (kind=dbl_kind), dimension(nilyr+2) :: &
         Temp      

      real (kind=dbl_kind) :: &
         xin                , & !transformed coordinate point       
         zspace                 !grid spacing for CICE vertical grid

      logical (kind=log_kind) :: &
         dbug             ! prints debugging output if true

       
      real (kind=dbl_kind), dimension(7), parameter :: & 
         Ps  =  (/0.8828e3_dbl_kind, -2.0570e3_dbl_kind, &  !polynomial fit
                  1.8255e3_dbl_kind, -0.7523e3_dbl_kind, &  !for initial
                  0.1456e3_dbl_kind, -0.0299e3_dbl_kind, &  !Salinity
                  0.0223_dbl_kind/)  

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

     !-----------------------------------------------------------------------------   
     !     BGC Layer Model
     !-----------------------------------------------------------------------------   

      dbug = .true.


      zphi(:,:,:,:,:) = c0    ! initial porosity for no ice 
      iki(:,:,:,:,:) = c0    ! permeability
      iDi(:,:,:,:,:) = c0    !interface diffusivity

      zTin(:,:,:,:,:) = c0   !initial bio grid ice temperature

      if (restart_S) then
          call read_restart_S
      else

       Rayleigh_criteria(:,:,:) = .false.    ! a no ice condition              
       first_ice(:,:,:) = .true.             ! if true, use initialization
      !-----------------------------------------------------------------
      ! Initial profile for S with hfrazilmin = 0.05 m and growing from no ice
      ! conditions:  Polynomial fit to high resolution run dt = 1s, nblyr = 40
      !              and hfrazilmin = 0.025 m and sss = 32.0 PSU 
      !------------------------------------------------------------------
      trcrn(:,:,nt_fbri,:,:) = c1
 
      if (tr_bgc_S)  then        ! take salinity from init_state: Error
        do n = 1,ncat
           do iblk = 1, nblocks
              do j = 1, ny_block
                 do i = 1, nx_block 
                    do k = 1,nblyr
                     ! if (k == 1) then
                     !   trcrn(i,j,nt_bgc_S+k-1,n,iblk) = salinz(i,j,1,iblk)
                     ! else
                       trcrn(i,j,nt_bgc_S+k-1,n,iblk) = sss(i,j,iblk)*salt_loss
                     ! endif
                    enddo      !k
                 enddo        !i
              enddo        !j
           enddo           !iblk
        enddo           !n
      endif

      endif    !restart_S


     1060    format (a30,2x,2D13.2)   ! dbl precision

      end subroutine init_zsalinity


!=======================================================================

!BOP
!
! !ROUTINE: solve_zsalinity
!
! !DESCRIPTION:
!
! update vertical salinity 
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine solve_zsalinity(nx_block, ny_block, &
                                   icells, n_cat, dt,&
                                   indxii, indxjj, & 
                                  ! tcells, tndxi, tndxj, tndxij, &
                                  ! pcells, pndxi, pndxj, pndxij, &
                                  ! qcells, qndxi, qndxj, qndxij, &
                                   trcrn_S, trcrn_q, trcrn_Si, &
                                   aicen, vicen, &  
                                   Sin, zTin,          &
                                   zphin, iphin,               &
                                   ikin, hinS_old, hinS, hin, &
                                   hin_old, iDin,  darcy_V,   &
                                   brine_sal, brine_rho,      &
                                   ibrine_sal,  ibrine_rho,   &
                                   Rayleigh_criteria,&
                                   first_ice, sss, sst,  &
                                   dh_top,dh_bot,  &  
                                   TLAT,TLON,        &
                                   l_stop, istop, jstop, fsicen, &
                                   fsicen_g, zphi_min,sloss) 

 
! !USES:
!
      use ice_therm_shared, only: solve_Sin
      use ice_therm_shared, only: calculate_Tin_from_qin
      use ice_calendar, only: istep1, time
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bgc1
!
! !INPUT/OUTPUT PARAMETERS:                                
!
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
         hin         , & ! new ice thickness (m)
         hinS_old    , & ! old brine height  (m)
         hinS        , & ! new brine height  (m)
         dh_top      , &  !  brine change in top and bottom for diagnostics (m)
         dh_bot      , &     !
         zphi_min     , &
         TLAT        , &
         sloss        , & !salinity flux to ocean from brine overflow  (g/m^2)
         TLON       
 

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         darcy_V      , &
         fsicen       , & ! total flux of salt out of ice over timestep(kg/m^2/s)
         fsicen_g         ! gravity drainage flux of salt  over timestep(kg/m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block, nblyr+2), &
         intent(inout) :: &
         Sin          , & ! Ice salinity ppt (on bio  grid)
         zTin          , & ! Ice Temperature ^oC (on bio grid)
         zphin         , & ! Ice porosity (on bio grid)
         brine_sal    , & ! brine salinity (ppt)
         brine_rho        ! brine density  (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block, nblyr), &
         intent(inout) :: &
         trcrn_S          ! salinity tracer ppt (on bio grid)

      real (kind=dbl_kind), dimension (nx_block,ny_block, nilyr), &
         intent(inout) :: &
         trcrn_q, &          ! enthalpy tracer 
         trcrn_Si           ! salinity on CICE grid

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         Rayleigh_criteria    ! .true. if minimun ice thickness (Ra_c)  was reached 
      
      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         first_ice    ! for first category ice only .true. initialized values should be used 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(out) :: &
         iDin          , & !  Diffusivity on the igrid   (1/s)
         ikin             !  permeability on the igrid 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(inout) :: &
         iphin         , & !  porosity on the igrid 
         ibrine_rho    , & ! brine rho on interface  
         ibrine_sal        !brine sal on interface   
       
    
      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop    ! indices of grid cell where code aborts

!
!EOP
!

     
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, zcells,pcells  , & ! horizontal index, combines i and j loops
         k, m, mm,nint     ! vertical biology layer index 
         
      !real (kind=dbl_kind), dimension (icells,nblyr) :: &       
      !   tracer_loss_z   ! loss term (tracer that remains in ice during flushing)   
                 !
      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         zndxi, zndxj, zndxij, & ! compressed indices for icells with hinS>thin, hinS_old>thin
         pndxi, pndxj, pndxij    ! compressed indices for the rest of the icells
    
      real (kind=dbl_kind), dimension(icells) :: &
         hinc         , & ! ice thickness (m)   
         hinSc        , & ! brine thickness (m)
         surface_S        ! salinity of ice above hin > hinS
  

      real (kind=dbl_kind) :: &
         Tmlts, &             ! melting temperature
         dts                  ! local timestep (s)

      logical (kind=log_kind), dimension(icells) :: & 
         Rayleigh      !

     
      real (kind=dbl_kind):: &
         Ttemp     ! initial temp profile on the CICE grid

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr+2) :: &
         trtmp0       , &  ! temporary, remapped tracers     !need extra 
         trtmp            ! temporary, remapped tracers     ! for nblyr_hist = nblyr+2


! local parameters
                           
      !-----------------------------------------------------------------------------
      !  Data from Cottier et al, 1999, JGR for fitting zbgc physical parameters
      !-----------------------------------------------------------------------------
    
      real (kind=log_kind), parameter :: &
          accuracy = 1.0e-5   !1.0e-5 kg/m^2/s difference between boundary fluxes 


      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      dts = dts_b
      nint = max(1,INT(dt/dts))  

       call ice_timer_start(timer_bgc1)  


       l_stop = .false.
       istop = 0
       jstop = 0
       fsicen_g(:,:) = c0 !only defined here
      ! tracer_loss_z(:,:) = c0

     
      !----------------------------------------------------------------
      ! Update boundary conditions
      !----------------------------------------------------------------
    
     
         do ij = 1, icells
            i = indxii(ij)
            j = indxjj(ij)    

            hinc(ij) = hin(i,j)
            hinSc(ij) = hinS(i,j)
            surface_S(ij) = min_salin
 
            Rayleigh(ij) = .true.
            if (n_cat == 1 .AND. hinS_old(i,j) < Ra_c) then
               Rayleigh(ij) = Rayleigh_criteria(i,j) ! only category 1 ice can be false 
            endif

            if ((dh_bot(i,j) > c0)) then
            
                 Sin(i,j,nblyr+2) = sss(i,j)
                 zTin(i,j,nblyr_hist) = sst(i,j)
                 brine_sal(i,j,nblyr_hist) = sss(i,j) 
                 brine_rho(i,j,nblyr_hist) = rhow 
                 zphin(i,j,nblyr_hist) = c1 
                      
         !-------------------------------
         ! !bottom melt
         !--------------------------------

            else  
                 Sin(i,j,nblyr+2) = Sin(i,j,nblyr+1)  
                 Tmlts =  -Sin(i,j,nblyr+2)* depressT 
                 zTin(i,j,nblyr_hist) =   zTin(i,j,nblyr+1)
            endif           
          
         enddo
      !----------------------------------------------------------------
      ! Solve for S using CICE T.  If solve_Sin = .true., then couple back
      ! to the thermodynamics
      !----------------------------------------------------------------


          call solve_S_dt (icells, nx_block, ny_block, &                    
                                      indxii,   indxjj,  &
                                      zcells, zndxi, zndxj, zndxij, &
                                      pcells, pndxi, pndxj, pndxij, &
                                      nint, dts, Sin, zTin,  aicen,  &
                                      zphin, iphin,      &
                                      igrid, bgrid,              &
                                      ikin, hinS_old, hinS, hin, hin_old,&
                                      iDin,darcy_V,      &
                                      brine_sal, &
                                      Rayleigh, &
                                      sss, dt, n_cat ,&
                                      dh_top,dh_bot, brine_rho,&
                                      ibrine_sal, ibrine_rho,TLAT,TLON,&
                                      fsicen,fsicen_g,istep1,  &
                                      zphi_min,sloss) 
       
      
     
     

  
     if (n_cat == 1)then

        do ij = 1, icells 
         i = indxii(ij)
         j = indxjj(ij)
          
           Rayleigh_criteria(i,j) = Rayleigh(ij)

        enddo

     endif

         trtmp0(:,:,:) = c0
         trtmp (:,:,:) = c0
       
        do k = 1,nblyr                  !back to bulk quantity  (do in solve_S_dt)

         do ij = 1, zcells 
            i = zndxi(ij)
            j = zndxj(ij)

               trcrn_S(i,j,k) =   Sin(i,j,k+1) 
               trtmp0(i,j,nt_sice+k-1) = trcrn_S(i,j,k)

            enddo
            do ij = 1,pcells
              i = pndxi(ij)
              j = pndxj(ij)

              trtmp0(i,j,nt_sice+k-1) = trcrn_S(i,j,k)

            enddo
         enddo           !ij
       
        
         call remap_layers_bgc_plus (nx_block,ny_block,        &
                             indxii,   indxjj,           &
                             icells,                   &
                             ntrcr,                    &
                             nilyr,                    &
                             nt_sice,                  &
                             trtmp0,    trtmp,          &
                             1,        nblyr,          &
                             hinc, hinSc,         &
                             cgrid(2:nilyr+1),         &
                             bgrid(2:nblyr+1), surface_S)
    
       

         do k = 1, nilyr

            do ij = 1, zcells
               i = zndxi(ij)
               j = zndxj(ij)        
               
               Tmlts = -trcrn_Si(i,j,k)*depressT
               Ttemp = calculate_Tin_from_qin(trcrn_q(i,j,k),Tmlts)
  
               trcrn_Si(i,j,k) = max(min_salin,trtmp(i,j,nt_sice+k-1))
               Tmlts = -trtmp(i,j,nt_sice+k-1)*depressT 

               
               Ttemp = min(Ttemp,Tmlts) 
               trcrn_q(i,j,k) = calculate_qin_from_Sin(Ttemp,Tmlts)
               
          
            enddo !ij
         enddo !k
        

 call ice_timer_stop(timer_bgc1) 
 
 
770 format (I6,D16.6)        
781 format (I6,I6,I6)
790 format (I6,I6)
791 format (f24.17)
792 format (2D16.6)
793 format (3D16.6)
794 format (4D15.5)
800 format (F10.4)

      end subroutine solve_zsalinity

!=======================================================================
!BOP
!
! !ROUTINE: solve_S_dt-- solve for bulk S  using T and  continuity
!                        but use own time-step
! !DESCRIPTION:
!
!  solves continuity explicitly using Lax-Wendroff-type scheme (MacCormack)
!  (Mendez-Nunez and Carroll,  Monthly Weather Review, 1993)
!
! !REVISION HISTORY:
!
! authors     Nicole Jeffery, LANL
!
!
!
! !INTERFACE:
!
      subroutine solve_S_dt &
                                     (icells, nx_block, ny_block,   &
                                      indxi,   indxj,               &
                                      qcells, qndxi, qndxj, qndxij, &   
                                      pcells, pndxi, pndxj, pndxij, &     
                                      nint, dts, Sin, zTin, aicen,  &
                                      zphin, iphin,                 &
                                      igrid, bgrid,                 &
                                      ikin, hin_old, hin, hice, hice_old, &             
                                      iDin,darcy_V,        &
                                      brine_sal,  Rayleigh, &
                                      sss, dt, n_cat, &
                                      dht, dhb, brine_rho, ibrine_sal, &
                                      ibrine_rho,TLAT,TLON,      &
                                      fsicen, fsicen_g,istep1, &
                                      zphi_min, sloss) !, tracer_loss_z)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block   , & ! block dimensions
         icells,      &
         n_cat,istep1, &
         nint                   ! number of interations     
        
      integer (kind=int_kind), intent(out) :: &
         qcells, pcells
                
      integer (kind=int_kind), dimension(icells), intent(out) :: &
         qndxi, qndxj ,& ! compressed indices 
         qndxij,       &
         pndxi, pndxj ,& ! compressed indices 
         pndxij

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj  ! compressed indices 

     real (kind=dbl_kind), intent(in) :: &
         dt,          &   ! timestep (s)
         dts              ! local timestep (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         fsicen,       & ! salt flux to (+ive) ocean (kg/m^2/s) or -ive if from the ocean
         fsicen_g        ! gravity drainage salt flux to (+ive) ocean (kg/m^2/s) or -ive if from the ocean

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         TLAT,TLON,    &
         aicen   ,     &
         zphi_min ,     &  ! surface porosity
         sloss           ! g/m^2   loss to ocean from brine runoff

     !real (kind=dbl_kind), dimension (icells,nblyr), intent(inout) :: &
     !    tracer_loss_z      !loss from each grid cell 
    

      logical (kind=log_kind), dimension (icells), &
         intent(inout) :: &
         Rayleigh    !if .true. convection is allowed.  if .false. not yet

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(in) :: &
         brine_sal       , & ! Internal brine salinity (ppt)
         brine_rho           ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(inout) :: &
         zphin             ! Porosity of layers

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(in) :: &
         ibrine_rho  , &    ! brine rho on interface 
         ibrine_sal         !brine sal on interface 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(inout) :: &
         iphin          ! Porosity of layers on interface

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1), intent(out) :: &
         iDin        , &    ! Diffusivity on the igrid (1/s) with minimum zphi condition
         ikin               !permeability on interface

       real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         sss       , &   ! sea surface salinity
         dht       , &   ! change in the ice top  (positive for melting)
         dhb       , &   ! change in the ice bottom (positive for freezing)
         hin_old   , &   ! brine thickness (m) 
         hin       , &   ! new brine thickness (m)
         hice      , &   ! ice thickness (m)
         hice_old        ! old ice thickness (m)


    
      real (kind=dbl_kind), dimension (nblyr_hist), intent(in) :: &
         bgrid            ! biology nondimensional grid layer points 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid         ! biology grid interface points 


      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), &
         intent(inout) :: &
         Sin             ! Bulk Salinity (ppt) contains previous timestep
                         ! and ocean ss

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), &
         intent(in) :: &
         zTin             ! Temperature of ice layers on bio grid for history file (^oC) 

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         darcy_V       ! Darcy velocity due to a pressure head (m/s) or melt
!
!EOP
!
     integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, ll      , & ! horizontal index, combines i and j loops
         k, m , mm, &   ! vertical biology layer index 
         qij_fault     ! location of conservation failure

      integer (kind=int_kind) :: &
         pij, qij            ! and h > thin or h <= thin


      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1) :: &
         iDin_p          ! Diffusivity on the igrid (1/s)/zphi^3 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist) :: &
         Din_p          ! Diffusivity on the igrid (1/s)/zphi^3  

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist) :: &
         Sintemp       ,& ! initial salinity
         pre_sin       ,&! estimate of  salinity of layers
         pre_sinb        ! estimate of  salinity of layers


      real (kind=dbl_kind), dimension (icells) :: &
         dh              ! (m) change in hbrine over dts

      real (kind=dbl_kind), dimension (nblyr+1) :: &  !on the zbgc vertical grid 2:nblyr+1
         dSbdx          ! gradient of brine rho on grid

      real (kind=dbl_kind), dimension (icells, nblyr+1) :: &  !on the zbgc vertical grid
         drho               ! brine difference rho_a-rho_b  (kg/m^3)

      real (kind=dbl_kind), dimension (icells, nblyr_hist) :: &  !on the zbgc vertical grid 2:nblyr+1
         Q_s, C_s, &   ! Functions in continuity equation
         V_s, U_s

       real (kind=dbl_kind), dimension (icells, nblyr+1) :: &  !on the zbgc vertical igrid 1:nblyr+1
         Ci_s, &  ! 
         Ui_s, &  ! interface function
         Vi_s     ! for conservation check

      real (kind=dbl_kind), dimension (icells,nblyr) :: &
         vel, &               ! advective velocity times dt (m)
         lapidus_diff               , & ! lapidus term and 
         flux_corr

      real (kind=dbl_kind), dimension (icells,nblyr+1) :: &
         ivel   

      real (kind=dbl_kind), dimension (icells):: &
         lapidus        ! artificial viscosity:  use lapidus_g for growth
 
   
      real (kind=dbl_kind), dimension (icells) :: &
         Ssum_old,Ssum_new, & ! depth integrated salt before and after timestep
         fluxcorr,          & ! flux correction to prevent S < min_salin
         Ssum_corr,        & ! numerical boundary flux correction
         fluxb, &  !bottom, top and total boundary flux (g/kg/m^2)
         fluxg, &  !bottom, top and total gravity drainage flux (g/kg/m^2)
         fluxm     !bottom, top and total molecular diffusion flux (g/kg/m^2)

      real (kind=dbl_kind) :: &
         sum_old,sum_new, hin_next  ! integrated salinity at t and t+dt
   
! local parameters

      logical (kind=log_kind) :: &   
         write_flag       , &    ! set to true at each timestep        
         good_numerics    , &    ! output of check_conservation 
         stable,            &    ! if false, redo with smaller timestep
         test_conservation       ! test that salt change is balanced by fluxes 

       real  (kind=dbl_kind), dimension (nblyr):: &
          lapA    , &
          lapB 
   
      real (kind=dbl_kind), dimension (nblyr+2) :: &
         vector_top, vector_bottom
     !--------------------------------------
     !  Initialize
     !--------------------------------------

     ! to remove  if-then statements    
      vector_top(:) = c0
      vector_bottom(:) = c0
      vector_top(2) = c1
      vector_bottom(nblyr+1) = c1

      write_flag = .true.
      test_conservation = .true.

      iDin_p(:,:,:) = c0   
      Din_p(:,:,:) = c0    
      lapA(:) = c1
      lapB(:) = c1
      lapA(nblyr) = c0
      lapB(1) = c0
      
       qcells = 0
       pcells = 0
       !tracer_loss(:) = c0   !should modify the bottom boundary flux
       !tracer_loss_z(:,:) = c0
       V_s(:,:) = c0
       U_s(:,:) = c0
       Q_s(:,:) = c0
       C_s(:,:) = c0
       Ci_s(:,:) = c0
       Ui_s(:,:) = c0
       Vi_s(:,:) = c0
       ivel(:,:) = c0
       vel(:,:) = c0
       dh(:) = c0
          
     !-----------------------------------------------------------------
     ! Find brine density gradient for gravity drainage parameterization
     !-----------------------------------------------------------------

         call calculate_drho(icells,nx_block,ny_block,indxi,indxj,igrid,bgrid,&
                                   brine_rho,ibrine_rho,drho)
     
     
     !-----------------------------------------------------------------
     ! Calculate zphi diffusivity on the grid points
     ! rhosi = 919-974 kg/m^2  set in bio_in
     ! rhow = 1026.0 density of sea water: uses kinematic viscosity (m^2/s)  in Q18
     ! dynamic viscosity  divided by density = kinematic. 
     !
     ! old representation of flushing
     !if (k == nblyr+1 .AND. Rayleigh(ij)) then
     !   k_min(ij) = MINLOC(kmin,DIM = 1,MASK = kmin .GE. c0)   
     !   darcy_V(i,j) = MINVAL(kmin) * gravit / viscos * hbrn(ij)/hin_old(ij)  !positive down
     ! endif
     !---------------------------------------------------------------- 

        
        do k = 2, nblyr+1
          
           do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij) 
          
               iDin_p(i,j,k) =k_o*gravit*l_skS/viscos_dynamic* &  
                     drho(ij,k)/(hin_old(i,j)**2)
 
               Din_p(i,j,k) = (iDin_p(i,j,k)*(igrid(k)-bgrid(k)) + iDin_p(i,j,k-1)*&
                            (bgrid(k)-igrid(k-1)))/(igrid(k)-igrid(k-1))           
             
          
          enddo   !ij
         enddo                !k


 
    !-------------------------------------------------
    ! Critical Ra_c value is only for the onset of convection in thin ice and not throughout
    !  therefore I need a flag to indicate the Ra_c was reached for a particular ice block
    ! Using a thickness minimum (Ra_c) for simplicity.
    !-------------------------------------------------    
        

        do ij = 1, icells
          i = indxi(ij)
          j = indxj(ij) 

          Din_p(i,j,nblyr+2) =  iDin_p(i,j,nblyr+1)

          if (.NOT. Rayleigh(ij) .AND. hin(i,j) < Ra_c) then
                    
                  Din_p(i,j,:) =  c0  
                  iDin_p(i,j,:) = c0  !
                  darcy_V(i,j) = c0

          else
              Rayleigh(ij) = .true.
          endif
 

          if (hin_old(i,j) > thin .AND.  hin(i,j) > thin) then
       
               qcells = qcells + 1
               qndxi(qcells) = i
               qndxj(qcells) = j
               qndxij(qcells) = ij 

               !-----------------------------------
               ! tracer_loss balanced by  the ocean flux
               !-----------------------------------

               ! if (darcy_V(i,j) < c0 .AND. hice(i,j) > hin(i,j)) then
               !    tracer_loss(ij) =  (min(-darcy_V(i,j)*dt/zphi_min(i,j), &
               !    (hice(i,j)-hin(i,j)))/hin_old(i,j)*& 
               !        min_salin)
                

               ! elseif  (darcy_V(i,j) > c0 .AND. hice(i,j) > hin(i,j)) then   
               !    tracer_loss(ij) =  -min(darcy_V(i,j)*dt/zphi_min(ij), &
               !    (hice(i,j)-hin(i,j)))/hin_old(i,j)*&
               !         min_salin
                
               ! endif

               !-----------------------------------
               ! surface boundary terms
               !-----------------------------------

               lapidus(ij) = lapidus_g
               if (dhb(i,j) < c0 ) lapidus(ij) = lapidus_m
               ivel(ij,1) =  dht(i,j)/hin_old(i,j) 
               U_s(ij,1) =  ivel(ij,1)/dt*dts 
               Ui_s(ij,1) = U_s(ij,1) 
             !  dSbdx(1) = c0 
               Ci_s(ij,1) = c0
               !-----------------------------------
               ! bottom boundary terms
               !-----------------------------------

               ivel(ij,nblyr+1) =  dhb(i,j)/hin_old(i,j)  
                                  !(igrid(nblyr+1)*dhb(i,j) - (igrid(nblyr+1)-c1)*dht(i,j))/hin_old(i,j)             
               Ui_s(ij,nblyr+1) = ivel(ij,nblyr+1)/dt*dts   !defined on interface 1 to n+1
               dSbdx(nblyr) = (ibrine_sal(i,j,nblyr+1) - ibrine_sal(i,j,nblyr))/(igrid(nblyr+1)-igrid(nblyr))
               C_s(ij,nblyr+1) = Dm/brine_sal(i,j,nblyr+1)*dts/hin_old(i,j)**2*&
                (ibrine_sal(i,j,nblyr+1) - ibrine_sal(i,j,nblyr))/(igrid(nblyr+1)-igrid(nblyr)) 
               vel(ij,nblyr) =(bgrid(nblyr+1)*(dhb(i,j)) - &
                    (bgrid(nblyr+1) - c1)* dht(i,j) )/hin_old(i,j)                   
               U_s(ij,nblyr+1) = vel(ij,nblyr)/dt*dts              
               V_s(ij,nblyr+1) = Din_p(i,j,nblyr+1)*brine_rho(i,j,nblyr+1)/rhosi&
                       *(rhosi/brine_sal(i,j,nblyr+1)/brine_rho(i,j,nblyr+1))**exp_h&
                       *dts*dSbdx(nblyr)   
               dSbdx(nblyr+1) =  (brine_sal(i,j,nblyr+2) - brine_sal(i,j,nblyr+1))/&
                               (bgrid(nblyr+2)-bgrid(nblyr+1)+ grid_oS/hin_old(i,j) )  
               C_s(ij, nblyr+2) = Dm/brine_sal(i,j,nblyr+2)*dts/hin_old(i,j)**2*&
                              (brine_sal(i,j,nblyr+2) - brine_sal(i,j,nblyr+1))/&
                               (bgrid(nblyr+2)-bgrid(nblyr+1)+ grid_oS/hin_old(i,j) )  
               U_s(ij,nblyr_hist) = ivel(ij,nblyr+1)/dt*dts 
               V_s(ij,nblyr_hist) = Din_p(i,j,nblyr+2)*brine_rho(i,j,nblyr+1)/rhosi &
                   *(rhosi/brine_sal(i,j,nblyr+1)/brine_rho(i,j,nblyr+1))**exp_h&
                   *dts*dSbdx(nblyr+1) 
          
               Ci_s(ij,nblyr+1) = C_s(ij,nblyr_hist)
               Vi_s(ij,nblyr+1) = V_s(ij,nblyr_hist)
 
               dh(ij) =(dhb(i,j)-dht(i,j))/dt*dts

              ivel(ij,2) =  (igrid(2)*dhb(i,j) - (igrid(2)-c1)*dht(i,j))/hin_old(i,j) 

              Ui_s(ij,2) = ivel(ij,2)/dt*dts   
          
              Vi_s(ij,2) = iDin_p(i,j,2)*ibrine_rho(i,j,2)/rhosi &
                            *(rhosi/ibrine_rho(i,j,2)/ibrine_sal(i,j,2))**exp_h*dts*&
                            (brine_sal(i,j,3)-brine_sal(i,j,2))&
                            /(bgrid(3)-bgrid(2))     
              dSbdx(1) = (ibrine_sal(i,j,2) - ibrine_sal(i,j,1))/(igrid(2)-igrid(1)+grid_oS/hin_old(i,j))

              C_s(ij,2) = Dm/brine_sal(i,j,2)*dts/hin_old(i,j)**2*&
                         (ibrine_sal(i,j,2) - ibrine_sal(i,j,1))/(igrid(2)-igrid(1)+grid_oS/hin_old(i,j)) 
              Ci_s(ij,2) = Dm/ibrine_sal(i,j,2)*dts/hin_old(i,j)**2*&
                           (brine_sal(i,j,3)-brine_sal(i,j,2))/(bgrid(3)-bgrid(2))
              vel(ij,1) = (bgrid(2)*(dhb(i,j)) - &
                    (bgrid(2) - c1)* dht(i,j) )/hin_old(i,j)
                  

              U_s(ij,2) = vel(ij,1)/dt*dts  
           
              V_s(ij,2) = Din_p(i,j,2)*brine_rho(i,j,2)/rhosi &
                       *(rhosi/brine_sal(i,j,2)/brine_rho(i,j,2))**exp_h&
                       *dts*dSbdx(1)   


          else                           !thin ice solve or first_ice
               pcells = pcells + 1
               pndxi(pcells) =  i
               pndxj(pcells) = j
               pndxij(pcells) =  ij
           endif
      enddo


      if (qcells > 0) then
           
          do k = 3, nblyr    
 
           do ij = 1, qcells   ! full solve: 
          
             qij = qndxij(ij)
             i = qndxi(ij)
             j = qndxj(ij) 

              ivel(qij,k) =  (igrid(k)*dhb(i,j) - (igrid(k)-c1)*dht(i,j))/hin_old(i,j) 

              Ui_s(qij,k) = ivel(qij,k)/dt*dts   
          
              Vi_s(qij,k) = iDin_p(i,j,k)*ibrine_rho(i,j,k)/rhosi &
                            *(rhosi/ibrine_rho(i,j,k)/ibrine_sal(i,j,k))**exp_h*dts*&
                            (brine_sal(i,j,k+1)-brine_sal(i,j,k))&
                            /(bgrid(k+1)-bgrid(k))     
              dSbdx(k-1) = (ibrine_sal(i,j,k) - ibrine_sal(i,j,k-1))/(igrid(k)-igrid(k-1))

              C_s(qij,k) = Dm/brine_sal(i,j,k)*dts/hin_old(i,j)**2*&
                         (ibrine_sal(i,j,k) - ibrine_sal(i,j,k-1))/(igrid(k)-igrid(k-1)) 
              Ci_s(qij,k) = Dm/ibrine_sal(i,j,k)*dts/hin_old(i,j)**2*&
                           (brine_sal(i,j,k+1)-brine_sal(i,j,k))/(bgrid(k+1)-bgrid(k))
              vel(qij,k-1) = (bgrid(k)*(dhb(i,j)) - &
                    (bgrid(k) - c1)* dht(i,j) )/hin_old(i,j)
                  

              U_s(qij,k) = vel(qij,k-1)/dt*dts  
           
              V_s(qij,k) = Din_p(i,j,k)*brine_rho(i,j,k)/rhosi &
                       *(rhosi/brine_sal(i,j,k)/brine_rho(i,j,k))**exp_h&
                       *dts*dSbdx(k-1)   

         

          enddo  !ij

         enddo !k
      
      
      !-----------------------------------------------------------------
      ! Solve
      !--------------------------------------------------------  
 
       do m = 1, nint
        
          lapidus_diff(:,:) = c0
          flux_corr(:,:) = c0
           Sintemp(:,:,:) = Sin(:,:,:)
           pre_sin(:,:,:) = Sin(:,:,:)  
           pre_sinb(:,:,:) = Sin(:,:,:)
           Ssum_old(:) = c0
           Ssum_new(:) = c0
           Ssum_corr(:) = c0
           fluxcorr(:) = c0
           fluxg(:) = c0
           fluxb(:) = c0
           fluxm(:) = c0
          

          do k = 2, nblyr
                                       
             do ij = 1, qcells
               i = qndxi(ij)
               j = qndxj(ij)    
               qij = qndxij(ij) 
            
       
 !forward-difference 

              Ssum_old(qij) = Ssum_old(qij) + Sin(i,j,k)*(igrid(k)-igrid(k-1))

              pre_sin(i,j,k) = Sin(i,j,k) + (Ui_s(qij,k)*(Sin(i,j,k+1) - Sin(i,j,k)) + &
                 V_s(qij,k+1)*Sin(i,j,k+1)**3 - V_s(qij,k)*Sin(i,j,k)**3 + &
                 C_s(qij,k+1)*Sin(i,j,k+1)-&
                 C_s(qij,k)*Sin(i,j,k))/(bgrid(k+1)-bgrid(k))

               pre_sin(i,j,nblyr+1) = Sin(i,j,nblyr+1) + (Ui_s(qij,nblyr+1)*(Sin(i,j,nblyr+2) - Sin(i,j,nblyr+1)) + &
                 V_s(qij,nblyr+2)*Sin(i,j,nblyr+1)**3 - V_s(qij,nblyr+1)*Sin(i,j,nblyr+1)**3+ C_s(qij,nblyr+2)*Sin(i,j,nblyr+2)-&
                 C_s(qij,nblyr+1)*Sin(i,j,nblyr+1) ) /(bgrid(nblyr+2)-bgrid(nblyr+1))             
                
              enddo  !qcells

          enddo    !k
      
          do k = nblyr+1, 3, -1  !nblyr+1

             do ij = 1, qcells
               i = qndxi(ij)
               j = qndxj(ij) 
               qij = qndxij(ij)
                      
  !backward-difference 
             

               pre_sinb(i,j,k) = p5*(Sin(i,j,k) + pre_sin(i,j,k) +  (Ui_s(qij,k-1)&
                  *(pre_sin(i,j,k) - pre_sin(i,j,k-1)) + &
                   V_s(qij,k)*pre_sin(i,j,k)**3 - &
                  V_s(qij,k-1)*pre_sin(i,j,k-1)**3 + C_s(qij,k)*pre_sin(i,j,k)-&
                  C_s(qij,k-1)*pre_sin(i,j,k-1))/(bgrid(k)-bgrid(k-1))) 

               pre_sinb(i,j,2) = p5*(Sin(i,j,2) + pre_sin(i,j,2) +  (Ui_s(qij,1)&
                  *(pre_sin(i,j,2) - min_salin) + &
                   V_s(qij,2)*pre_sin(i,j,2)**3 - &
                  V_s(qij,1)*pre_sin(i,j,1)**3 + C_s(qij,2)*pre_sin(i,j,2)-&
                   C_s(qij,1)*pre_sin(i,j,1))/(bgrid(2)-bgrid(1))) ! - &
              
                Q_s(qij,k) = V_s(qij,k)*pre_sin(i,j,k)**2 + U_s(qij,k) + C_s(qij,k)
        
                
 
                Q_s(qij,2) = V_s(qij,2)*pre_sin(i,j,2)**2 + U_s(qij,2) + C_s(qij,2)
                Q_s(qij,1) = V_s(qij,1)*pre_sin(i,j,2)**2 + Ui_s(qij,1) 
                Q_s(qij,nblyr_hist) = V_s(qij,nblyr_hist)*pre_sin(i,j,nblyr+1)**2 + & 
                        Ui_s(qij,nblyr+1) 
                 
            
             enddo  !qcells
          enddo   !k

        ! Add artificial viscosity   [Lapidus,1967] [Lohner et al, 1985]
        ! * more for melting ice
        !--------------------------------------------------------------------- 
    
         do k = 2, nblyr+1  
 
             do ij = 1, qcells
               i = qndxi(ij)
               j = qndxj(ij) 
               qij = qndxij(ij)
                  
               lapidus_diff(qij,k-1) =    lapidus(qij)/real(nblyr,kind=dbl_kind)**2/&
                  (igrid(k)-igrid(k-1))* &
                  ( lapA(k-1)*ABS(Q_s(qij,k+1)-Q_s(qij,k))*(abs(pre_sinb(i,j,k+1))-abs(pre_sinb(i,j,k)))/&
                  (bgrid(k+1)-bgrid(k) )**2 - &
                   lapB(k-1)*ABS(Q_s(qij,k)-Q_s(qij,k-1))*(abs(pre_sinb(i,j,k))-abs(pre_sinb(i,j,k-1)))/&
                  (bgrid(k)-bgrid(k-1))**2)
            
              
              Sin(i,j,k) = pre_sinb(i,j,k) + lapidus_diff(qij,k-1)

              if (Sin(i,j,k) < min_salin) then
                  flux_corr(qij,k-1) = min_salin - Sin(i,j,k) !  flux into the ice
                  Sin(i,j,k) = min_salin  
              elseif (Sin(i,j,k) > max_salin) then
                stable = .false.
                write(nu_diag, *) 'Very Bad value in solve_dt_S-- istep1:',istep1
                write(nu_diag, *) 'TLAT,TLON:',TLAT(i,j)*rad_to_deg,TLON(i,j)*rad_to_deg
                write(nu_diag, *) 'Sin(i,j,k),i,j,k',Sin(i,j,k),i,j,k
                write(nu_diag, *) 'V_s(qij,k),Ui_s(qij,k),U_s(qij,k),k',V_s(qij,k),&
                                  Ui_s(qij,k),U_s(qij,k),k
                write(nu_diag, *) 'V_s(qij,k-1),Ui_s(qij,k-1),U_s(qij,k-1)',V_s(qij,k-1),&
                                  Ui_s(qij,k-1),U_s(qij,k-1)
                write(nu_diag, *) 'V_s(qij,k+1),Ui_s(qij,k+1),U_s(qij,k+1)',V_s(qij,k+1),&
                                  Ui_s(qij,k+1),U_s(qij,k+1)
                write(nu_diag, *) 'C_s(qij,k),C_s(qij,k+1),C_s(qij,k-1)',C_s(qij,k),C_s(qij,k+1),C_s(qij,k-1)
                write(nu_diag, *)'brine_rho(i,j,k),brine_sal(i,j,k)'
                write(nu_diag,*)brine_rho(i,j,k),brine_sal(i,j,k)
                write(nu_diag,*)'ivel(qij,k-1),ivel(qij,k),vel(qij,k-1),Din_p(i,j,k),darcy_V(i,j)'
                write(nu_diag,*)ivel(qij,k-1),ivel(qij,k),vel(qij,k-1),Din_p(i,j,k),darcy_V(i,j)
                write(nu_diag,*)'pre_sin(i,j,k),pre_sin(i,j,k-1),pre_sin(i,j,k+1)'
                write(nu_diag,*)pre_sin(i,j,k),pre_sin(i,j,k-1),pre_sin(i,j,k+1)
                write(nu_diag,*)'pre_sinb(i,j,k),pre_sinb(i,j,k-1),pre_sinb(i,j,k+1)'
                write(nu_diag,*)pre_sinb(i,j,k),pre_sinb(i,j,k-1),pre_sinb(i,j,k+1)
                write(nu_diag,*)'zTin(i,j,k),zTin(i,j,k-1)'
                write(nu_diag,*)zTin(i,j,k),zTin(i,j,k-1)
                write(nu_diag,*)'Q_s(qij,k),Q_s(qij,k+1),Q_s(qij,k-1)'
                write(nu_diag,*)Q_s(qij,k),Q_s(qij,k+1),Q_s(qij,k-1)
                write(nu_diag,*)'bgrid(k),bgrid(k+1),bgrid(k-1)'
                write(nu_diag,*)bgrid(k),bgrid(k+1),bgrid(k-1)
                write(nu_diag,*)'hin_old(i,j),hin(i,j),dhb(i,j),dht(i,j),&
                                 dts, dt'
                write(nu_diag,*)hin_old(i,j),hin(i,j),dhb(i,j),dht(i,j),&
                                dts,dt
                write(nu_diag,*)'darcy_V(i,j)'
                write(nu_diag,*)darcy_V(i,j)
                write(nu_diag,*)'Sin(i,j,1),Sin(i,j,nblyr_hist)'
                write(nu_diag,*)Sin(i,j,1),Sin(i,j,nblyr_hist)
                write(nu_diag,*)Sin(i,j,1),Sin(i,j,nblyr_hist)

                call abort_ice ('ice: Solve_S_dt error')                 
              endif
            
              Ssum_new(qij) = Ssum_new(qij) + Sin(i,j,k)*(igrid(k)-igrid(k-1))
              fluxcorr(qij) = fluxcorr(qij) + flux_corr(qij,k-1)*(igrid(k)-igrid(k-1))
             
        enddo  !qcells
       enddo   !k

       call calc_salt_fluxes (m,qcells,qndxi,qndxj,qndxij,icells,nx_block,ny_block,&
                             Ui_s,dh,hin_old,Sintemp,pre_sin,fluxb,fluxg,fluxm,V_s,&
                             C_s,Ssum_corr,fsicen_g,fsicen,Ssum_old,fluxcorr,dts)
 

       if (test_conservation) then
         good_numerics = .true.
         qij_fault = 0

         call check_conserve_salt(m,qcells,qndxi,qndxj,qndxij,icells,nx_block,ny_block, dt,dts,&
                                Ssum_old,Ssum_new,&
                                fluxcorr,dh,hin_old,fluxb,fluxg,fluxm, Ssum_corr,&
                                good_numerics,qij_fault)

         if (.NOT. good_numerics) then
              i = qndxi(qij_fault)
              j = qndxj(qij_fault)
              write(nu_diag,*)',nint, Category,i,j,ij,TLAT,TLON-- istep1:'&
                              ,nint,n_cat,i,j,qij,TLAT(i,j)*rad_to_deg,&
                              TLON(i,j)*rad_to_deg,istep1
              write(nu_diag,*)'dhb(i,j),dht(i,j),Rayleigh(qij):'&
                              ,dhb(i,j),dht(i,j),Rayleigh(qij) 
              write(nu_diag,*) 'fsicen(i,j),fsicen_g(i,j):'
              write(nu_diag,*) fsicen(i,j),fsicen_g(i,j)
              write(nu_diag,*) 'Ssum_old(qij),Ssum_new(qij):'
              write(nu_diag,*) Ssum_old(qij),Ssum_new(qij)
              write(nu_diag,*) 'fluxb(qij):', fluxb(qij)
              write(nu_diag,*) 'fluxg(qij):',fluxg(qij)
              write(nu_diag,*) 'fluxm(qij):',fluxm(qij)
              write(nu_diag,*) 'Sin(i,j,1),Sin(i,j,nblyr+2)'
              write(nu_diag,*) Sin(i,j,1),Sin(i,j,nblyr+2)
              write(nu_diag,*) 'pre_sin(i,j,1),pre_sin(i,j,nblyr+2)'
              write(nu_diag,*) pre_sin(i,j,1),pre_sin(i,j,nblyr+2)
             do  k = 2,nblyr+1
                write(nu_diag, *) 'Sin(i,j,k),i,j,k',Sin(i,j,k),i,j,k
                write(nu_diag, *) 'V_s(qij,k),Ui_s(qij,k),U_s(qij,k),k',V_s(qij,k),&
                                  Ui_s(qij,k),U_s(qij,k),k
                write(nu_diag, *) 'V_s(qij,k-1),Ui_s(qij,k-1),U_s(qij,k-1)',V_s(qij,k-1),&
                                  Ui_s(qij,k-1),U_s(qij,k-1)
                write(nu_diag, *) 'V_s(qij,k+1),Ui_s(qij,k+1),U_s(qij,k+1)',V_s(qij,k+1),&
                                  Ui_s(qij,k+1),U_s(qij,k+1)
                write(nu_diag, *) 'C_s(qij,k),C_s(qij,k+1),C_s(qij,k-1)',C_s(qij,k),C_s(qij,k+1),C_s(qij,k-1)
                write(nu_diag, *)'brine_rho(i,j,k),brine_sal(i,j,k)'
                write(nu_diag,*)brine_rho(i,j,k),brine_sal(i,j,k)
                write(nu_diag,*)'ivel(qij,k-1),ivel(qij,k),vel(qij,k-1),Din_p(i,j,k),darcy_V(i,j)'
                write(nu_diag,*)ivel(qij,k-1),ivel(qij,k),vel(qij,k-1),Din_p(i,j,k),darcy_V(i,j)
                write(nu_diag,*)'pre_sin(i,j,k),pre_sin(i,j,k-1),pre_sin(i,j,k+1)'
                write(nu_diag,*)pre_sin(i,j,k),pre_sin(i,j,k-1),pre_sin(i,j,k+1)
                write(nu_diag,*)'pre_sinb(i,j,k),pre_sinb(i,j,k-1),pre_sinb(i,j,k+1)'
                write(nu_diag,*)pre_sinb(i,j,k),pre_sinb(i,j,k-1),pre_sinb(i,j,k+1)
                write(nu_diag,*)'zTin(i,j,k),zTin(i,j,k-1)'
                write(nu_diag,*)zTin(i,j,k),zTin(i,j,k-1)
                write(nu_diag,*)'Q_s(qij,k),Q_s(qij,k+1),Q_s(qij,k-1)'
                write(nu_diag,*)Q_s(qij,k),Q_s(qij,k+1),Q_s(qij,k-1)
                write(nu_diag,*)'bgrid(k),bgrid(k+1),bgrid(k-1)'
                write(nu_diag,*)bgrid(k),bgrid(k+1),bgrid(k-1)
                write(nu_diag,*)'hin_old(i,j),hin(i,j),dhb(i,j),dht(i,j),&
                                 dts, dt'
                write(nu_diag,*)hin_old(i,j),hin(i,j),dhb(i,j),dht(i,j),&
                                dts,dt
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
          do k = 2, nblyr+1 
            sum_old = sum_old + Sin(i,j,k)*hin_old(i,j)*(igrid(k)-igrid(k-1))
            Sin(i,j,k) = max(min_salin,Sin(i,j,k))
            sum_new = sum_new + Sin(i,j,k)*hin(i,j)*(igrid(k)-igrid(k-1))
          enddo  !k
          
          fsicen(i,j) = fsicen(i,j) -rhosi*(sum_new-sum_old)*p001/dt
          fsicen_g(i,j) = c0
        enddo   !pcells
     endif    ! add/melt ice only
   

      !-----------------------------------------------------------------
      ! Move this to bgc calculation if using tr_salinity
      ! Calculate zphin, iphin, ikin, iDin and iDin_N 
      !-----------------------------------------------------------------

     iDin(:,:,:) = c0
     iphin(:,:,:) = c1
     ikin(:,:,:) = c0    
             
       do k= 1, nblyr+1
          do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)  
            
              if (k < nblyr+1) zphin(i,j,k+1) = min(c1,max(puny, Sin(i,j,k+1)*rhosi/(brine_sal(i,j,k+1)*brine_rho(i,j,k+1)))) 
              if (k == 1) then     
                zphin(i,j,k) = min(c1,max(puny, Sin(i,j,k)*rhosi/(brine_sal(i,j,k)*brine_rho(i,j,k))))  
                iphin(i,j,k) = zphin(i,j,2)
              elseif (k == nblyr+1) then
                iphin(i,j,nblyr+1) = zphin(i,j,nblyr+1)
              else
                iphin(i,j,k) =min(c1, max(c0,(zphin(i,j,k+1) - zphin(i,j,k))/(bgrid(k+1) - &
                             bgrid(k))*(igrid(k)-bgrid(k)) + zphin(i,j,k)))
              endif 
              
              ikin(i,j,k) = k_o*iphin(i,j,k)**exp_h 
           
           enddo  !ij
         enddo    !k
       
 
         do k = 2, nblyr+1
 
            do ij = 1, qcells
               i = qndxi(ij)
               j = qndxj(ij) 
               qij = qndxij(ij)
               
               iDin(i,j,k) =  (iphin(i,j,k)*Dm)/hin_old(i,j)**2  
               if (Rayleigh(qij) .AND. hin(i,j) .GE. Ra_c) &
               iDin(i,j,k) =iDin(i,j,k) + l_sk*ikin(i,j,k)*gravit/viscos_dynamic* &  
                           drho(qij,k)/hin_old(i,j)**2   
            enddo  !ij
  
            do ij = 1, pcells
                i = pndxi(ij)
                j = pndxj(ij) 
                pij = pndxij(ij)
                iDin(i,j,k) = iphin(i,j,k)*Dm/hin_old(i,j)**2

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
           drho(ij,k) = max(c0,c2*(rho_a(ij,k)-rho_2a(ij,k)))
          enddo
        enddo


     end subroutine calculate_drho

!=======================================================================
!BOP
!
! !ROUTINE: calc_salt_fluxes
!                          
!
! !DESCRIPTION:
!
! calcuate fluxes
!                                                       
! !REVISION HISTORY:
!
! authors     Nicole Jeffery, LANL
!
!
!
! !INTERFACE:
!
      subroutine calc_salt_fluxes &
                                    (mint,qcells,qndxi,qndxj,qndxij,icells,nx_block,ny_block,&
                                    Ui_s,dh,hin_old,Sintemp,pre_sin,fluxb,fluxg,fluxm,V_s,&
                                    C_s,Ssum_corr,fsicen_g,fsicen,Ssum_old,fluxcorr,dts)
!
! !USES:
!
! 
!
! !INPUT/OUTPUT PARAMETERS:
!
  
      integer(kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         qcells, icells,mint      
 
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         qndxi, qndxj,qndxij ! compressed indices for icells with hin > thin

      real (kind=dbl_kind), intent(in) :: &
        dts       !  halodynamic timesteps (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr_hist), intent(in) :: &
         Sintemp       ,& ! initial salinity
         pre_sin          ! estimate of  salinity of layers

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: &
        fsicen, fsicen_g

      real (kind=dbl_kind), dimension(icells), intent(in) :: &
         dh      ! (m) change in hbrine over dts

      real (kind=dbl_kind), dimension(icells), intent(inout) :: &
         Ssum_old       , &  ! initial integrated salt content
         Ssum_corr      , &  ! boundary flux correction due to numerics
         fluxb, &  ! total boundary salt flux into the ice (+ into the ice)
         fluxg, &  ! total gravity drainage salt flux into the ice (+ into the ice)
         fluxm     ! total gravity drainage salt flux into the ice (+ into the ice)

      real (kind=dbl_kind), dimension(icells), intent(in) :: &
         fluxcorr        ! flux correction to ensure S >= min_salin

    real (kind=dbl_kind), dimension (icells, nblyr+1), intent(in) :: &  !on the zbgc vertical igrid 1:nblyr+1
         Ui_s ! interface function

      real (kind=dbl_kind), dimension (icells, nblyr_hist), intent(in) :: &  !on the zbgc vertical grid 2:nblyr+1
         C_s, &   ! Functions in continuity equation
         V_s

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         hin_old          ! (m) initial brine height
                           
!
!EOP
!
     integer (kind=int_kind) :: &
         i, j, ij, qij   ! vertical and horizontal indices

      real (kind=dbl_kind), dimension (icells) :: &
         Ssum_corr_flux,         & ! numerical boundary flux correction
         fluxb_b, fluxb_t, &  !bottom, top and total boundary flux (g/kg/m^2)
         fluxg_b, fluxg_t, &  !bottom, top and total gravity drainage flux (g/kg/m^2)
         fluxm_b, fluxm_t    !bottom, top and total molecular diffusion flux (g/kg/m^2)

      real (kind=dbl_kind) :: hin_next
              
           do ij = 1, qcells
               i = qndxi(ij)
               j = qndxj(ij) 
               qij = qndxij(ij)

                 Ssum_old(qij) = Ssum_old(qij) + Sintemp(i,j,nblyr+1)*(igrid(nblyr+1)-igrid(nblyr))
               !----------------------------------
               ! boundary fluxes (positive into the ice)
               !--------------------------------

                  fluxb_b(qij) = Ui_s(qij,nblyr+1)*((c1-dh(qij)/hin_old(i,j))* &
                                             Sintemp(i,j,nblyr+2) + p5*(pre_sin(i,j,nblyr+1)-&
                                             (c1-dh(qij)/hin_old(i,j))*Sintemp(i,j,nblyr+1)))
             
                  fluxb_t(qij) = -Ui_s(qij,1)*(min_salin +  &
                              p5*((c1-dh(qij)/hin_old(i,j))*Sintemp(i,j,2) - pre_sin(i,j,2)))

                  fluxb(qij) = fluxb_b(qij) + fluxb_t(qij)

               !----------------------------------
               ! gravity drainage fluxes (positive into the ice)
               !--------------------------------
               

                  fluxg_b(qij) =  ((c1-dh(qij)/hin_old(i,j))*  V_s(qij,nblyr+2)*Sintemp(i,j,nblyr+1)**3  +  &
                                  p5*V_s(qij,nblyr+1)*(pre_sin(i,j,nblyr+1)**3 - &
                                  (c1-dh(qij)/hin_old(i,j))*Sintemp(i,j,nblyr+1)**3))
                
                  fluxg_t(qij) =  -(V_s(qij,1)*pre_sin(i,j,1)**3 + &
                                  p5*V_s(qij,2)*((c1-dh(qij)/hin_old(i,j))*Sintemp(i,j,2)**3- &
                                      pre_sin(i,j,2)**3))

                
                  fluxg(qij) =  fluxg_b(qij) + fluxg_t(qij)
               !----------------------------------
               ! diffusion fluxes (positive into the ice)
               !--------------------------------
                  fluxm_b(qij) =  ((c1-dh(qij)/hin_old(i,j))* C_s(qij,nblyr+2)* &
                                 Sintemp(i,j,nblyr+2)  +  p5*C_s(qij,nblyr+1)*(pre_sin(i,j,nblyr+1) - &
                                  (c1-dh(qij)/hin_old(i,j))*Sintemp(i,j,nblyr+1)))
   
                  fluxm_t(qij) =-(C_s(qij,1)*pre_sin(i,j,1) + &
                                   p5*C_s(qij,2)*((c1-dh(qij)/hin_old(i,j))*Sintemp(i,j,2)-pre_sin(i,j,2)))
                
                  fluxm(qij) =  fluxm_b(qij) + fluxm_t(qij)

              Ssum_corr(qij) = (-dh(qij)/hin_old(i,j) + p5*(dh(qij)/hin_old(i,j))**2)*Ssum_old(qij) 
                               !add to boundary flux!!

              hin_next = hin_old(i,j) + real(mint,kind=dbl_kind)*dh(qij)
              Ssum_corr_flux(qij) =  dh(qij)*Ssum_old(qij)/hin_next +Ssum_corr(qij)
          
              fsicen_g(i,j) = fsicen_g(i,j) - hin_next*fluxg(qij)*rhosi*p001/dts  

              fsicen(i,j) = fsicen(i,j) - hin_next*(fluxb(qij) +  fluxg(qij) + fluxm(qij)+ &
                    Ssum_corr_flux(qij)+ fluxcorr(qij))*rhosi*p001/dts 
            
       enddo
        

     end subroutine calc_salt_fluxes


!=======================================================================
!BOP
!
! !ROUTINE: check_conserve_salt
!                          
!
! !DESCRIPTION:
!
!  Written in flux conservative form dSin/dt = -dF(x,Sin)/dx 
!                                                       
! !REVISION HISTORY:
!
! authors     Nicole Jeffery, LANL
!
!
!
! !INTERFACE:
!
      subroutine check_conserve_salt &
                                     (mint,qcells, qndxi,qndxj,qndxij,icells,nx_block,ny_block, dt,dts,&
                                      Ssum_old,Ssum_new,&
                                      fluxcorr,dh,hin_old,fluxb,fluxg,fluxm,Ssum_corr,& 
                                      good_numerics, qij_fault) 
!
! !USES:
!
! 
!
! !INPUT/OUTPUT PARAMETERS:
!
  
      integer(kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         qcells, icells,mint      
  
      integer(kind=int_kind), intent(inout) :: &
         qij_fault          !qndxij value where conservation fails 

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         qndxi, qndxj,qndxij ! compressed indices for icells with hin > thin

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

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         hin_old          ! (m) initial brine height
                           
!
!EOP
!
     integer (kind=int_kind) :: &
         i, j, ij, qij   ! vertical and horizontal indices

     real (kind=dbl_kind):: &
         diff2 , & !
         dsum  , & ! salt for conservation (g/kg/m^3)
         dsum_flux , & ! salt change in kg/m^2/s
         flux_tot  , &          ! fluxg+fluxb
         fsicen_old, &  !
         order
      
     real (kind=dbl_kind), parameter :: &
         accuracy = 1.0e-14_dbl_kind ! g/kg/m^3 difference between boundary fluxes 
              
     
        do ij = 1,qcells
         i =  qndxi(ij)
         j = qndxj(ij) 
         qij = qndxij(ij)

         dsum = Ssum_new(qij) -  Ssum_old(qij) 
 
         flux_tot = fluxb(qij) + fluxg(qij) + fluxm(qij) + fluxcorr(qij) 

         dsum_flux =(Ssum_new(qij)*(hin_old(i,j) + (real(mint,kind=dbl_kind))*dh(qij)) - &
                 Ssum_old(qij)*(hin_old(i,j) + (real(mint,kind=dbl_kind)-c1)*dh(qij)) )*rhosi*p001/dts 

         order = abs(dh(qij)/hin_old(i,j))

         if (abs(dsum) > c0) then
           diff2 = abs(dsum- Ssum_corr(qij)- flux_tot)
           if (diff2 >  puny .AND. diff2 > order ) then 
              qij_fault = qij
              good_numerics = .false.
              write(nu_diag,*) 'Poor salt conservation: check_conserve_salt'
              write(nu_diag,*) 'i,j,qij, mint:', i,j,qij,mint
              write(nu_diag,*) 'Ssum_corr(qij)',Ssum_corr(qij)
              write(nu_diag,*) 'fluxb(qij),fluxg(qij),fluxm(qij),flux_tot,fluxcorr(qij):'
              write(nu_diag,*) fluxb(qij),fluxg(qij),fluxm(qij),flux_tot,fluxcorr(qij)
              write(nu_diag,*) 'fluxg(qij),',fluxg(qij)
              write(nu_diag,*) 'dsum,dsum_flux,',dsum,dsum_flux
              write(nu_diag,*) 'Ssum_new(qij),Ssum_old(qij),hin_old(i,j),dh(qij):'
              write(nu_diag,*) Ssum_new(qij),Ssum_old(qij),hin_old(i,j),dh(qij)
              write(nu_diag,*) 'diff2,order,puny',diff2,order,puny
           endif
        
         endif

        enddo
        

     end subroutine check_conserve_salt

!=======================================================================
!
!BOP
!
! !IROUTINE: merge_fluxes - aggregate flux information over ITD
!
! !INTERFACE:
!
      subroutine merge_S_fluxes (nx_block, ny_block,   &
                               icells,               &
                               indxi,    indxj,      &
                               aicen,                &
                               S_totn, S_tot, &
                               hbrin, hbri,          &    
                               dsnown,               &
                               dsnow,                &
                               fsice,  fsicen,       &
                               fsice_g, fsicen_g)
!
! !DESCRIPTION:
!
! Aggregate flux information from all ice thickness categories
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
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
          indxi, indxj    ! compressed indices for cells with aicen > puny

      ! single category fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &          
          aicen   , & ! concentration of ice
          S_totn  , & ! tot salinity in category (psu*volume*rhosi)
          hbrin    , & ! brine height of cat n (m)
          dsnown  , & ! snow  growth                    (m)
          fsicen  , & ! salt flux                       (kg/m**2/s)
          fsicen_g    ! gravity drainage salt flux      (kg/m**2/s)
       

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &          
          S_tot  , & ! tot salinity (psu*rhosi*total vol ice)
          hbri    , & ! tot brine height over all categories (m)
          dsnow   , & ! snow growth                     (m)
          fsice   , & ! salt flux                       (kg/m**2/s)
          fsice_g     ! gravity drainage salt flux      (kg/m**2/s)
     

!EOP
!
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
         S_tot    (i,j)  = S_tot   (i,j) + S_totn  (i,j)  !already *aicen
         hbri     (i,j)  = hbri    (i,j) + hbrin(i,j)*aicen(i,j)  

         ! ocean tot and gravity drainage salt fluxes

         fsice   (i,j) = fsice     (i,j) + fsicen  (i,j)*aicen(i,j)
         fsice_g (i,j) = fsice_g   (i,j) + fsicen_g(i,j)*aicen(i,j)

         ! ice/snow thickness

         dsnow    (i,j) = dsnow    (i,j) + dsnown  (i,j)*aicen(i,j)
      enddo                     ! ij
      
      end subroutine merge_S_fluxes

!==========================================================================
!BOP
!
! !IROUTINE: column_sum_S - find total salinity
!
! !INTERFACE:
!
      subroutine column_sum_S (nx_block, ny_block,   &
                             icells,   indxi,   indxj, &
                             ntrcr, nblyr, vicen, trcrn, &
                             S_totn)
!
! !DESCRIPTION:
!
! For each grid cell, sum field over all ice layers.  "Net" refers to  the column
! integration while "avg"  is normalized by the ice thickness
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
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

     

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
          S_totn  ! avg salinity (psu*rhosi*vol of ice)
           
!
!EOP
!
      integer (kind=int_kind) :: &
           i, j, ij     , & ! horizontal indices
           k                ! layer index

     S_totn(:,:) = c0

   
      do k = 1, nblyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            S_totn(i,j) = S_totn(i,j) + rhosi*trcrn(i,j,nt_bgc_S+k-1) &
                         * vicen(i,j)*trcrn(i,j,nt_fbri)/real(nblyr,kind=dbl_kind)
                             
         enddo                  ! ij
      enddo                     ! n
      
      end subroutine column_sum_S


!=======================================================================
!BOP
!
! !IROUTINE: S_diags - writes max,min,global sums to standard out
!
! !INTERFACE:
!
      subroutine S_diags (dt)
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
      use ice_blocks, only: nx_block, ny_block
      use ice_broadcast, only: broadcast_scalar
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, pbloc, plat, plon
      use ice_domain_size, only: max_blocks
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
         phinS, phinS1, phinS2,&
         phbrn, pdh_top1,pdh_bot1, pdh_top2,pdh_bot2, psice_rho, pfsice, &
         pfsice_g, pdarcy_V1, pdarcy_V2, &
         pdhi_top1, pdhi_bot1,pdhi_top2, pdhi_bot2

      ! vertical  fields of category 1 at diagnostic points for bgc layer model
      real (kind=dbl_kind), dimension(npnt,nblyr+2) :: &
         pphin, pgrid, pphin1
      real (kind=dbl_kind), dimension(npnt,nblyr) :: &
        pSin, pSice, pSin1, pSin2

      real (kind=dbl_kind), dimension(npnt,nblyr+1) :: &
         pzTin, piDin

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

               pfsice(n) = fsice(i,j,iblk)   
               pfsice_g(n) = fsice_g(i,j,iblk)            
               phinS(n) = c0            
               phinS1(n) = c0            
               phinS2(n) = c0
               phbrn(n) = c0
               psice_rho(n) = c0
               pdh_top1(n) = c0
               pdh_bot1(n) = c0
               pdhi_top1(n) = c0
               pdhi_bot1(n) = c0
               pdh_top2(n) = c0
               pdh_bot2(n) = c0
               pdhi_top2(n) = c0
               pdhi_bot2(n) = c0
               pdarcy_V1(n) = c0
               pdarcy_V2(n) = c0
               do nn = 1,ncat
                  psice_rho(n) = psice_rho(n) + sice_rho(i,j,nn,iblk)*aicen(i,j,nn,iblk)
               enddo
               if (aice(i,j,iblk) > c0) then
                 psice_rho(n) = psice_rho(n)/aice(i,j,iblk)
               endif
              
               
                if (hbrine .and. aice(i,j,iblk) > c0) &
                   phinS(n) = trcr(i,j,nt_fbri,iblk)*vice(i,j,iblk)/aice(i,j,iblk)
                                     
                if (aicen(i,j,1,iblk)> c0) then
                  if (hbrine) phinS1(n) = trcrn(i,j,nt_fbri,1,iblk)*vicen(i,j,1,iblk)/&
                                                aicen(i,j,1,iblk)
                  pdh_top1(n) = dh_top(i,j,1,iblk)
                  pdhi_top1(n) = dhi_top(i,j,1,iblk)
                  pdh_bot1(n) = dh_bot(i,j,1,iblk)
                  pdhi_bot1(n) = dhi_bot(i,j,1,iblk)
                  pdarcy_V1(n) = darcy_V(i,j,1,iblk)
                endif  
                if (aicen(i,j,2,iblk)> c0) then
                  if (hbrine) phinS2(n) = trcrn(i,j,nt_fbri,2,iblk) *vicen(i,j,2,iblk)/&
                                                    aicen(i,j,2,iblk)
                  pdh_top2(n) = dh_top(i,j,2,iblk)
                  pdhi_top2(n) = dhi_top(i,j,2,iblk)
                  pdh_bot2(n) = dh_bot(i,j,2,iblk)
                  pdhi_bot2(n) = dhi_bot(i,j,2,iblk)
                  pdarcy_V2(n) = darcy_V(i,j,2,iblk)
                endif
                if (hbrine .AND. aice(i,j,iblk) > c0) phbrn(n) = phinS(n) - &
                       rhosi/rhow*vice(i,j,iblk)/aice(i,j,iblk) - & 
                       rhos/rhow*vsno(i,j,iblk)/aice(i,j,iblk)
               
              
                 do k = 1, nblyr+1
                   pzTin(n,k) = c0
                   piDin(n,k) = c0
                  do nn = 1,ncat
                   pzTin(n,k) = pzTin(n,k) + zTin(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                   piDin(n,k) = piDin(n,k) +  iDi(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
                  enddo
                   if (vice(i,j,iblk) > c0)then
                    pzTin(n,k) = pzTin(n,k)/vice(i,j,iblk)
                    piDin(n,k) = piDin(n,k)/vice(i,j,iblk) 
                   endif
                 enddo                 !k
                 do k = 1, nblyr+2
                   pphin(n,k) = c0
                   pphin1(n,k) = c0
                   if (aicen(i,j,1,iblk) > c0) pphin1(n,k) = zphi(i,j,k,1,iblk)
                   do nn = 1,ncat
                     pphin(n,k) = pphin(n,k) + zphi(i,j,k,nn,iblk)*vicen(i,j,nn,iblk)
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
                       if (aicen(i,j,2,iblk)> c0) pSin2(n,k) = trcrn(i,j,nt_bgc_S+k-1,2,iblk)
                   
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
           call broadcast_scalar(pdhi_top1(n), pmloc(n)) 
           call broadcast_scalar(pdhi_bot1(n), pmloc(n)) 
           call broadcast_scalar(pdhi_top2(n), pmloc(n)) 
           call broadcast_scalar(pdhi_bot2(n), pmloc(n))  
           call broadcast_scalar(psice_rho(n), pmloc(n))  
           call broadcast_scalar(pfsice_g(n), pmloc(n))  
           call broadcast_scalar(pdarcy_V1(n), pmloc(n))  
           call broadcast_scalar(pdarcy_V2(n), pmloc(n)) 
           call broadcast_scalar(pfsice(n), pmloc(n)) 

           do k = 1,nblyr+1
              call broadcast_scalar(pzTin (n,k), pmloc(n))
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
        write(nu_diag,902) '       Lat, Long         ',plat(1),plon(1), &
                                                       plat(2),plon(2)
        write(nu_diag,903) '  my_task, iblk, i, j     ', &
                              pmloc(1),pbloc(1),piloc(1),pjloc(1), &
                              pmloc(2),pbloc(2),piloc(2),pjloc(2)

      
       
        write(nu_diag,900) 'hinS                   = ',phinS(1),phinS(2)
        write(nu_diag,900) 'hinS cat 1             = ',phinS1(1),phinS1(2)
        write(nu_diag,900) 'hinS cat 2             = ',phinS2(1),phinS2(2)
        write(nu_diag,900) 'hbrn                   = ',phbrn(1),phbrn(2)
        write(nu_diag,900) 'dhinS cat 1 top        = ',pdh_top1(1),pdh_top1(2)
        write(nu_diag,900) 'dhin cat 1 top         = ',pdhi_top1(1),pdhi_top1(2)
        write(nu_diag,900) 'dhinS cat 1 bottom     = ',pdh_bot1(1),pdh_bot1(2)
        write(nu_diag,900) 'dhin  cat 1 bottom     = ',pdhi_bot1(1),pdhi_bot1(2)
        write(nu_diag,900) 'dhinS cat 2 top        = ',pdh_top2(1),pdh_top2(2)
        write(nu_diag,900) 'dhin  cat 2 top        = ',pdhi_top2(1),pdhi_top2(2)
        write(nu_diag,900) 'dhinS cat 2 bottom     = ',pdh_bot2(1),pdh_bot2(2)
        write(nu_diag,900) 'dhin  cat 2 bottom     = ',pdhi_bot2(1),pdhi_bot2(2)
        write(nu_diag,900) 'Avg sea ice density    = ',psice_rho(1),psice_rho(2)
        write(nu_diag,900) 'Total Salt flux        = ',pfsice(1),pfsice(2)
        write(nu_diag,900) 'Gravity Drainage Salt flux   = ',pfsice_g(1),pfsice_g(2)
        write(nu_diag,900) 'Darcy V cat 1 (m/s)    = ',pdarcy_V1(1),pdarcy_V1(2)
        write(nu_diag,900) 'Darcy V cat 2 (m/s)    = ',pdarcy_V2(1),pdarcy_V2(2)

      
        write(nu_diag,*) '                         '
        write(nu_diag,*) ' Top down bgc Layer Model'
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'zTin(1) ice temp',' zTin(2) ice temp  '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((pzTin(n,k),n = 1,2), k = 1,nblyr+1)
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'iDi(1) diffusivity  ','iDi(2) diffusivity  '
        write(nu_diag,*) '---------------------------------------------------'
        write(nu_diag,802) ((piDin(n,k),n=1,2), k = 1,nblyr+1)
        write(nu_diag,*) '                         '
        write(nu_diag,803) 'zphi(1) porosity   ','zphi(2) porosity  '
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
!BOP
!
! !IROUTINE: write_restart_S - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_S(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for a S restart
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
              restart_file(1:lenstr(restart_file)),'.S.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_S,filename,0)

      if (my_task == master_task) then
        write(nu_dump_S) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
        write(nu_diag,*) 'S Restart written ',istep1,time,time_forc
      endif

      diag = .true.
      !--------------------------
      !Salinity and extras
      !--------------------------

      do n = 1, ncat
         if (hbrine) call ice_write(nu_dump_S,0,trcrn(:,:,nt_fbri,n,:),'ruf8',diag)

         do k = 1,nblyr
         if (tr_bgc_S) call ice_write(nu_dump_S,0,trcrn(:,:,nt_bgc_S+k-1,n,:),'ruf8',diag)
         enddo

      enddo

      if (tr_bgc_S)  then
         call ice_write(nu_dump_S,0,sss,'ruf8',diag)


      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (Rayleigh_criteria(i,j,iblk)) then
               Rayleigh_real(i,j,iblk) = c1
            else
               Rayleigh_real(i,j,iblk) = c0
            endif
            if (first_ice(i,j,iblk)) then
               first_ice_real(i,j,iblk) = c1
            else
               first_ice_real(i,j,iblk) = c0
            endif
         enddo
         enddo
      enddo

      call ice_write(nu_dump_S,0,Rayleigh_real,'ruf8',diag)
      call ice_write(nu_dump_S,0,first_ice_real,'ruf8',diag)
      endif

      if (my_task == master_task) close(nu_dump_S)

      end subroutine write_restart_S

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_S - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_S(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for a Salinity restart
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
      use ice_read_write, only: ice_open, ice_read, ice_write
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday  , &   ! year, month, day
          nsize, nr0 , nt              !restart_n + nblyr, nr0 = 0 (interp restart_n to nblyr)

      character(len=char_len_long) :: &
         filename, filename0, string1, string2

      logical (kind=log_kind) :: &
         diag, hit_eof

     real (kind=dbl_kind), dimension (:,:,:,:), allocatable :: &
          read_S

     real (kind=dbl_kind), dimension (:,:,:), allocatable :: &
          interp_S
          

      if (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         read(nu_rst_pointer,'(a)') filename0
         filename = trim(filename0)
         close(nu_rst_pointer)

         ! reconstruct path/file
         n = index(filename0,trim(restart_file))
         if (n == 0) call abort_ice('S restart: filename discrepancy')
         string1 = trim(filename0(1:n-1))
         string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
         write(filename,'(a,a,a,a)') &
            string1(1:lenstr(string1)), &
            restart_file(1:lenstr(restart_file)),'.S', &
            string2(1:lenstr(string2))
      endif ! master_task

      call ice_open(nu_restart_S,filename,0)

      if (my_task == master_task) then
        read(nu_restart_S) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
        write(nu_diag,*) 'Salinity Restart read at istep=',istep0,time,time_forc
      endif

      diag = .true.

      if (tr_bgc_S .AND. restart_n /= nblyr) then
         nsize = restart_n + nblyr
         nr0 = 0
         nt = 1
         allocate (read_S(nx_block,ny_block,nsize,nblocks), &
                   interp_S(nx_block,ny_block,nsize))
      endif

      do n = 1, ncat
         if (hbrine ) then
              write(nu_diag,*) 'cat ',n, &
                               ' fbri'
              call ice_read(nu_restart_S,0,trcrn(:,:,nt_fbri,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
         endif
             
         if (tr_bgc_S ) then
           if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' S for each bgc layer'
           if (restart_n == nblyr) then

             do k = 1,nblyr
              call ice_read(nu_restart_S,0,trcrn(:,:,nt_bgc_S+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
             
             enddo

           else !restart_n .NE. nblyr
           
             do k = 1,restart_n
               call ice_read(nu_restart_S,0,read_S(:,:,k,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
             enddo
         
             do iblk = 1, nblocks
               call interp_restart(nx_block,ny_block, &
                                 nsize, restart_n, read_S(:,:,:,iblk), &
                                 interp_S, nr0, nt, nblyr)


               do k = 1,nblyr
                 do j = 1, ny_block
                 do i = 1, nx_block
                   trcrn(i,j,nt_bgc_S+k-1,n,iblk) = interp_S(i,j,k)
                 enddo !j
                 enddo !i
               enddo !k
         
             enddo  !iblk
                    
           endif  !restart_n

         endif  !tr_bgc_S 
        enddo
      

      if (my_task == master_task) &
         write(nu_diag,*) 'ocean  sss'
      
        if (tr_bgc_S )call ice_read(nu_restart_S,0,sss,'ruf8',diag, &
                       field_loc_center, field_type_scalar)

        call ice_write(nu_dump_S,0,Rayleigh_real,'ruf8',diag)
        call ice_write(nu_dump_S,0,first_ice_real,'ruf8',diag)
        do iblk = 1, nblocks
          do j = 1, ny_block
          do i = 1, nx_block
            if (Rayleigh_real(i,j,iblk) .GE. c1) then
               Rayleigh_criteria (i,j,iblk) = .true.
            else
               Rayleigh_criteria (i,j,iblk) = .false.
            endif
            if (first_ice_real(i,j,iblk) .GE. c1) then
               first_ice (i,j,iblk) = .true.
            else
               first_ice (i,j,iblk) = .false.
            endif
         enddo
         enddo
      enddo
     

      if (tr_bgc_S .AND. restart_n /= nblyr) then
         deallocate (read_S, interp_S)
      endif


      if (my_task == master_task) close(nu_restart_S)

      end subroutine read_restart_S

!=======================================================================
!BOP
!
! !IROUTINE: interp_restart
!
! !INTERFACE:
!
      subroutine interp_restart (nx_block,ny_block,        &
                                  ! indxi,   indxj,           &
                                  ! icells,                   &
                                   ntrcr,                    &
                                   nlyrn,                    &
                                   read_S,   dummy_S, &! trtmp,          &
                                   nr0,  it, &    
                                   nblyr)
!
! !DESCRIPTION:
!
! Interpolate bio restart file with  n_nd in donor and n_nr in receiver
!
! !REVISION HISTORY:
!
! author: njeffery
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
       !  icells            , & ! number of ocean/ice cells
         ntrcr             , & ! number of tracers in use
         it                , & ! tracer index in top layer
         nr0               , & ! receiver category
         nlyrn             , & ! number of ice layers
         nblyr                 ! number of biology layers

      !integer (kind=int_kind), dimension (icells), intent(in) :: &
      !   indxi             , & ! indices for i/j directions
      !   indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
         intent(in) ::       &
         read_S    !trcrn                 ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr), &
         intent(inout) ::    &
         dummy_S   !trtmp                 ! temporary, remapped ice tracers
!
!EOP
!
      integer (kind=int_kind) :: &
           i, j, ij      , & ! horizontal indices
           kd, kr, kdr   , & ! more indices
           n_nd          , & ! number of layers in donor
           n_nr              ! number of layers in receiver

      real (kind=dbl_kind), dimension (nblyr*nlyrn) :: &
           trdr          , & ! trcrn over combined grid
           trr               ! trcrn on receiving grid

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do j = 1, ny_block
      do i = 1, nx_block

         if (nr0 == 0) then ! restart_n  to nblyr
            n_nd = nlyrn
            n_nr = nblyr
         else               ! nblyr to restart_n
            n_nd = nblyr
            n_nr = nlyrn
         endif

         if (n_nd /= n_nr) then ! remap
            do kd = 1, n_nd
               do kr = 1, n_nr
                  kdr = kr + (kd - 1) * n_nr
                  trdr(kdr) = read_S(i,j,it+kd-1)
               enddo
            enddo
            do kr = 1, n_nr
               trr(kr) = c0
               do kd = 1, n_nd
                  kdr = kd + (kr - 1) * n_nd
                  trr(kr) = trr(kr) + trdr(kdr)
               enddo
               dummy_S(i,j,it+kr-1) = trr(kr)/n_nd
            enddo
         else     ! fill trtmp with original trcrn values
            do kr = 1, n_nr
               dummy_S(i,j,it+kr-1) = read_S(i,j,it+kr-1)
            enddo
         endif

      enddo
      enddo

      end subroutine interp_restart

!=======================================================================

      end module ice_zsalinity

!=======================================================================
