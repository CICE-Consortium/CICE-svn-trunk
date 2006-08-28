!=======================================================================
!BOP
!
! !MODULE: ice_dyn_evp - elastic-viscous-plastic sea ice dynamics model
!
! !DESCRIPTION:
!
! Elastic-viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. {\em J. Phys. Oceanogr.}, {\bf 27}, 1849--1867.
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. {\em Journal of Computational Physics}, {\bf 170},
! 18--38.
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere---Incorporation of Metric Terms. {\em Monthly Weather Review},
! {\bf 130}, 1848--1865.
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.
!
!
! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$LastChangedBy$
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
! 2004: Block structure added by William Lipscomb
! 2005: Removed boundary calls for stress arrays (WHL)
! 
! !INTERFACE:
!
      module ice_dyn_evp
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_domain_size
      use ice_constants
!
!EOP
!
      implicit none
      save

      ! namelist parameters

      integer (kind=int_kind) ::
     &   kdyn         ! type of dynamics ( 1 = evp )
     &,  ndte         ! number of subcycles:  ndte=dt/dte

      logical (kind=log_kind) ::
     &   evp_damping  ! if true, use evp damping procedure

      ! other EVP parameters

      character (len=char_len) :: 
     &   yield_curve  ! 'ellipse' or 'teardrop'
                      ! Note: teardrop still needs debugging
                                                                      ! 
      real (kind=dbl_kind), parameter ::
     &   dragw = 0.00536_dbl_kind * rhow
                      ! drag coefficient for water on ice *rhow (kg/m^3)
     &,  eyc = 0.36_dbl_kind
                      ! coefficient for calculating the parameter E
     &,  cosw = c1    ! cos(ocean turning angle)  ! turning angle = 0
     &,  sinw = c0    ! sin(ocean turning angle)  ! turning angle = 0
     &,  a_min = p001 ! minimum ice area
     &,  m_min = p01  ! minimum ice mass

      real (kind=dbl_kind) ::
     &   ecci         ! 1/e^2
     &,  dtei         ! 1/dte, where dte is subcycling timestep (1/s)
     &,  dte2T        ! dte/2T
     &,  dteT         ! dte/T
     &,  denom1       ! constants for stress equation
     &,  denom2       !
     &,  rcon         ! for damping criterion (kg/s)

      real (kind=dbl_kind), allocatable :: 
     &   fcor_blk(:,:,:)   ! Coriolis parameter (1/s)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: evp - elastic-viscous-plastic dynamics driver
!
! !INTERFACE:
!
      subroutine evp (dt)
!
! !DESCRIPTION:
!
! Elastic-viscous-plastic dynamics driver
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_boundary
      use ice_blocks
      use ice_domain
      use ice_state
      use ice_flux
      use ice_grid
      use ice_timers
      use ice_mechred, only: ice_strength
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::
     &   dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: 
     &   ksub            ! subcycle step
     &,  iblk            ! block index

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) ::
     &   tmass        ! total mass of ice and snow
     &,  waterx       ! for ocean stress calculation, x (m/s)
     &,  watery       ! for ocean stress calculation, y (m/s)
     &,  forcex       ! work array: combined atm stress and ocn tilt, x
     &,  forcey       ! work array: combined atm stress and ocn tilt, y
     &,  aiu          ! ice fraction on u-grid
     &,  umass        ! total mass of ice and snow (u grid)
     &,  umassdtei    ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8)::
     &   str          ! stress combinations for momentum equation

      integer (kind=int_kind), 
     &   dimension (nx_block,ny_block,max_blocks) ::
     &   icetmask   ! ice extent mask (T-cell)

      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      str(:,:,:) = c0

       ! This call is needed only if dt changes during runtime.
!      call set_evp_parameters (dt)

      !-----------------------------------------------------------------
      ! boundary updates
      !-----------------------------------------------------------------

      call ice_timer_start(timer_bound)
      call update_ghost_cells(aice,              bndy_info, 
     &                        field_loc_center,  field_type_scalar)
      call update_ghost_cells(vice,              bndy_info, 
     &                        field_loc_center,  field_type_scalar)
      call update_ghost_cells(vsno,              bndy_info, 
     &                        field_loc_center,  field_type_scalar)
      call ice_timer_stop(timer_bound)

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! ice strength
      !-----------------------------------------------------------------

         call ice_strength (nx_block, ny_block,
     &                      nghost,
     &                      aice    (:,:,  iblk),
     &                      vice    (:,:,  iblk),
     &                      aice0   (:,:,  iblk),
     &                      aicen   (:,:,:,iblk),    
     &                      vicen   (:,:,:,iblk),
     &                      strength(:,:,  iblk) )

      !-----------------------------------------------------------------
      ! preparation for dynamics
      !-----------------------------------------------------------------

         call evp_prep1 (nx_block,           ny_block,
     &                   nghost,
     &                   aice    (:,:,iblk), vice    (:,:,iblk),
     &                   vsno    (:,:,iblk), tmask   (:,:,iblk),
     &                   strairxT(:,:,iblk), strairyT(:,:,iblk),
     &                   strairx (:,:,iblk), strairy (:,:,iblk), 
     &                   tmass   (:,:,iblk), icetmask(:,:,iblk))

      enddo                     ! iblk

      !-----------------------------------------------------------------
      ! boundary updates
      !-----------------------------------------------------------------
      call ice_timer_start(timer_bound)

      call update_ghost_cells(strength,         bndy_info,
     &                        field_loc_center, field_type_scalar)

      call update_ghost_cells(icetmask,         bndy_info,
     &                        field_loc_center, field_type_scalar)

      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! convert fields from T to U grid
      !-----------------------------------------------------------------

      call t2ugrid_vector(strairx)
      call t2ugrid_vector(strairy)
      call to_ugrid(tmass,umass)
      call to_ugrid(aice, aiu)

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! more preparation for dynamics
      !-----------------------------------------------------------------

         call evp_prep2 (nx_block,             ny_block,
     &                   nghost,
     &                   aice      (:,:,iblk), aiu       (:,:,iblk),
     &                   tmass     (:,:,iblk), umass     (:,:,iblk),
     &                   umassdtei (:,:,iblk), fcor_blk  (:,:,iblk),
     &                   tmask     (:,:,iblk), umask     (:,:,iblk),      
     &                   uocn      (:,:,iblk), vocn      (:,:,iblk),       
     &                   strairx   (:,:,iblk), strairy   (:,:,iblk),
     &                   ss_tltx   (:,:,iblk), ss_tlty   (:,:,iblk),    
     &                   icetmask  (:,:,iblk), iceumask  (:,:,iblk),
     &                   fm        (:,:,iblk),
     &                   strtltx   (:,:,iblk), strtlty   (:,:,iblk),
     &                   waterx    (:,:,iblk), watery    (:,:,iblk),
     &                   forcex    (:,:,iblk), forcey    (:,:,iblk),
     &                   stressp_1 (:,:,iblk), stressp_2 (:,:,iblk),  
     &                   stressp_3 (:,:,iblk), stressp_4 (:,:,iblk),
     &                   stressm_1 (:,:,iblk), stressm_2 (:,:,iblk),  
     &                   stressm_3 (:,:,iblk), stressm_4 (:,:,iblk),
     &                   stress12_1(:,:,iblk), stress12_2(:,:,iblk), 
     &                   stress12_3(:,:,iblk), stress12_4(:,:,iblk),
     &                   uvel      (:,:,iblk), vvel      (:,:,iblk))

      enddo

      do ksub = 1,ndte        ! subcycling

      !-----------------------------------------------------------------
      ! stress tensor equation, total surface stress
      !-----------------------------------------------------------------

         do iblk = 1, nblocks

            if (trim(yield_curve) == 'teardrop') then
               ! need to compute prs_sig for principal stresses?

               call stress_teardrop
     &                     (nx_block,             ny_block,
     &                      nghost,               ksub,
     &                      icetmask  (:,:,iblk),
     &                      uvel      (:,:,iblk), vvel      (:,:,iblk),       
     &                      dxt       (:,:,iblk), dyt       (:,:,iblk),
     &                      dxhy      (:,:,iblk), dyhx      (:,:,iblk),
     &                      cxp       (:,:,iblk), cyp       (:,:,iblk),
     &                      cxm       (:,:,iblk), cym       (:,:,iblk),
     &                      tarear    (:,:,iblk), tinyarea  (:,:,iblk),
     &                      strength  (:,:,iblk),
     &                      stressp_1 (:,:,iblk), stressp_2 (:,:,iblk),  
     &                      stressp_3 (:,:,iblk), stressp_4 (:,:,iblk),
     &                      stressm_1 (:,:,iblk), stressm_2 (:,:,iblk),
     &                      stressm_3 (:,:,iblk), stressm_4 (:,:,iblk),
     &                      stress12_1(:,:,iblk), stress12_2(:,:,iblk),
     &                      stress12_3(:,:,iblk), stress12_4(:,:,iblk),
     &                      shear     (:,:,iblk), divu      (:,:,iblk),
     &                      rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk),
     &                      str       (:,:,:) )

            else                       ! elliptical yield curve

               call stress (nx_block,             ny_block,
     &                      nghost,               ksub,
     &                      icetmask  (:,:,iblk),
     &                      uvel      (:,:,iblk), vvel      (:,:,iblk),       
     &                      dxt       (:,:,iblk), dyt       (:,:,iblk),
     &                      dxhy      (:,:,iblk), dyhx      (:,:,iblk),
     &                      cxp       (:,:,iblk), cyp       (:,:,iblk),
     &                      cxm       (:,:,iblk), cym       (:,:,iblk),
     &                      tarear    (:,:,iblk), tinyarea  (:,:,iblk),
     &                      strength  (:,:,iblk),
     &                      stressp_1 (:,:,iblk), stressp_2 (:,:,iblk),  
     &                      stressp_3 (:,:,iblk), stressp_4 (:,:,iblk),
     &                      stressm_1 (:,:,iblk), stressm_2 (:,:,iblk),
     &                      stressm_3 (:,:,iblk), stressm_4 (:,:,iblk),
     &                      stress12_1(:,:,iblk), stress12_2(:,:,iblk),
     &                      stress12_3(:,:,iblk), stress12_4(:,:,iblk),
     &                      shear     (:,:,iblk), divu      (:,:,iblk),
     &                      Delt      (:,:,iblk), prs_sig   (:,:,iblk),
     &                      rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk),
     &                      str       (:,:,:) )

            endif               ! yield_curve
      !-----------------------------------------------------------------
      ! momentum equation
      !-----------------------------------------------------------------

            call stepu (nx_block,            ny_block,
     &                  nghost,              iceumask(:,:,iblk),
     &                  aiu      (:,:,iblk), str     (:,:,:),
     &                  uocn     (:,:,iblk), vocn    (:,:,iblk),    
     &                  waterx   (:,:,iblk), watery  (:,:,iblk),
     &                  forcex   (:,:,iblk), forcey  (:,:,iblk),
     &                  umassdtei(:,:,iblk), fm      (:,:,iblk),  
     &                  uarear   (:,:,iblk),
     &                  strocnx  (:,:,iblk), strocny (:,:,iblk), 
     &                  strintx  (:,:,iblk), strinty (:,:,iblk),
     &                  uvel     (:,:,iblk), vvel    (:,:,iblk))

         enddo

         call ice_timer_start(timer_bound)
         call update_ghost_cells (uvel,               bndy_info,
     &                            field_loc_NEcorner, field_type_vector)
         call update_ghost_cells (vvel,               bndy_info,
     &                            field_loc_NEcorner, field_type_vector)
         call ice_timer_stop(timer_bound)

      enddo                     ! subcycling

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      do iblk = 1, nblocks

         call evp_finish
     &        (nx_block,           ny_block,
     &         nghost,
     &         uvel    (:,:,iblk), vvel    (:,:,iblk),    
     &         uocn    (:,:,iblk), vocn    (:,:,iblk),       
     &         strocnx (:,:,iblk), strocny (:,:,iblk), 
     &         strocnxT(:,:,iblk), strocnyT(:,:,iblk))

      enddo

      call u2tgrid_vector(strocnxT)    ! shift
      call u2tgrid_vector(strocnyT)

      call ice_timer_stop(timer_dynamics)    ! dynamics

      end subroutine evp

!=======================================================================
!BOP
!
! !IROUTINE: init_evp - initialize parameters needed for evp dynamics
!
! !INTERFACE:
!
      subroutine init_evp (dt)
!
! !DESCRIPTION:
!
! Initialize parameters and variables needed for the evp dynamics
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_boundary
      use ice_blocks
      use ice_domain
      use ice_state
      use ice_flux
      use ice_grid
      use ice_fileunits
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::
     &   dt      ! time step
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, k
     &,  iblk            ! block index

      real (kind=dbl_kind) ::
     &   dte             ! subcycling timestep for EVP dynamics, s
     &,  ecc             ! (ratio of major to minor ellipse axes)^2
     &,  tdamp2          ! 2(wave damping time scale T)

      call set_evp_parameters (dt)

      if (my_task == master_task) then
         write(nu_diag,*) 'dt  = ',dt
         write(nu_diag,*) 'dte = ',dt/real(ndte,kind=dbl_kind)
         write(nu_diag,*) 'tdamp =', eyc*dt
      endif

      allocate(fcor_blk(nx_block,ny_block,max_blocks))

      do iblk = 1, nblocks

         ! velocity
         uvel(:,:,iblk) = c0    ! m/s
         vvel(:,:,iblk) = c0    ! m/s

         ! strain rates
         divu (:,:,iblk) = c0
         shear(:,:,iblk) = c0
         Delt (:,:,iblk) = c0
         rdg_conv (:,:,iblk) = c0
         rdg_shear(:,:,iblk) = c0

         ! Coriolis parameter
!!         fcor_blk(:,:,iblk) = 1.46e-4_dbl_kind ! Hibler 1979, N. Hem; 1/s
         fcor_blk(:,:,iblk) = c2*omega*sin(ULAT(:,:,iblk)) ! 1/s

         ! stress tensor,  kg/s^2
         stressp_1 (:,:,iblk) = c0
         stressp_2 (:,:,iblk) = c0
         stressp_3 (:,:,iblk) = c0
         stressp_4 (:,:,iblk) = c0
         stressm_1 (:,:,iblk) = c0
         stressm_2 (:,:,iblk) = c0
         stressm_3 (:,:,iblk) = c0
         stressm_4 (:,:,iblk) = c0
         stress12_1(:,:,iblk) = c0
         stress12_2(:,:,iblk) = c0
         stress12_3(:,:,iblk) = c0
         stress12_4(:,:,iblk) = c0

      enddo                     ! iblk

      end subroutine init_evp

!=======================================================================
!BOP
!
! !IROUTINE: set_evp_parameters - set parameters for evp dynamics
!
! !INTERFACE:
!
      subroutine set_evp_parameters (dt)
!
! !DESCRIPTION:
!
! Set parameters needed for the evp dynamics.
! Note: This subroutine is currently called only during initialization.
!       If the dynamics time step can vary during runtime, it should
!        be called whenever the time step changes.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) ::
     &   dt      ! time step
!
!EOP
!
      real (kind=dbl_kind) ::
     &   dte             ! subcycling timestep for EVP dynamics, s
     &,  ecc             ! (ratio of major to minor ellipse axes)^2
     &,  tdamp2          ! 2*(wave damping time scale T)

      ! elastic time step
      dte = dt/real(ndte,kind=dbl_kind)        ! s
      dtei = c1/dte              ! 1/s

      ! major/minor axis length ratio, squared
      ecc  = c4
      ecci = p25                  ! 1/ecc

      ! constants for stress equation
      tdamp2 = c2*eyc*dt                    ! s
      dte2T = dte/tdamp2                    ! ellipse (unitless)
!!!      dteT = dte/(eyc*dt)   ! teardrop (unitless)
      dteT = dte/(c2*eyc*dt)   !lipscomb - larger tdamp seems to work better
      denom1 = c1/(c1+dte2T)
      denom2 = c1/(c1+dte2T*ecc)
      rcon = 1230._dbl_kind*eyc*dt*dtei**2  ! kg/s

      end subroutine set_evp_parameters

!=======================================================================
!BOP
!
! !IROUTINE: evp_prep1 - compute quantities needed for stress tensor and mom eqns
!
! !INTERFACE:
!
      subroutine evp_prep1 (nx_block,  ny_block,
     &                      nghost,
     &                      aice,      vice,    
     &                      vsno,      tmask,
     &                      strairxT,  strairyT,
     &                      strairx,   strairy, 
     &                      tmass,     icetmask)

! !DESCRIPTION:
!
! Computes quantities needed in the stress tensor (sigma)
! and momentum (u) equations, but which do not change during
! the thermodynamics/transport time step:
! --wind stress shift to U grid,
! --ice mass and ice extent masks,
! --pressure (strength), and part of the forcing stresses
! initializes ice velocity for new points to ocean sfc current
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   aice        ! concentration of ice
     &,  vice        ! volume per unit area of ice          (m)
     &,  vsno        ! volume per unit area of snow         (m)
     &,  strairxT    ! stress on ice by air, x-direction
     &,  strairyT    ! stress on ice by air, y-direction

      logical (kind=log_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   tmask       ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(out) ::
     &   strairx     ! stress on ice by air, x-direction
     &,  strairy     ! stress on ice by air, y-direction
     &,  tmass       ! total mass of ice and snow

      integer (kind=int_kind), dimension (nx_block,ny_block),
     &   intent(out) ::
     &   icetmask    ! ice extent mask (T-cell)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j
     &,  ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      logical (kind=log_kind), dimension(nx_block,ny_block) ::
     &   tmphm               ! temporary mask

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      icetmask (:,:) = 0

      do j = 1, ny_block
      do i = 1, nx_block

      !-----------------------------------------------------------------
      ! total mass of ice and snow, centered in T-cell
      ! NOTE: vice and vsno must be up to date in all grid cells,
      !       including ghost cells
      !-----------------------------------------------------------------
         if (tmask(i,j)) then
            tmass(i,j) = (rhoi*vice(i,j) + rhos*vsno(i,j)) ! kg/m^2
         else
            tmass(i,j) = c0
         endif

      !-----------------------------------------------------------------
      ! ice extent mask (T-cells)
      !-----------------------------------------------------------------
         tmphm(i,j) = tmask(i,j) .and. (aice (i,j) > a_min)
     &                           .and. (tmass(i,j) > m_min)

      !-----------------------------------------------------------------
      ! prep to convert to U grid
      !-----------------------------------------------------------------
          ! Factor of aice needed for correct treatment of free drift
          strairx(i,j) = strairxT(i,j)*aice(i,j)
          strairy(i,j) = strairyT(i,j)*aice(i,j)

      enddo
      enddo

      !-----------------------------------------------------------------
      ! augmented mask (land + open ocean)
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi

         ! extend ice extent mask (T-cells) to points around pack
         if (tmphm(i-1,j+1) .or. tmphm(i,j+1) .or. tmphm(i+1,j+1) .or.
     &       tmphm(i-1,j)   .or. tmphm(i,j)   .or. tmphm(i+1,j)   .or.
     &       tmphm(i-1,j-1) .or. tmphm(i,j-1) .or. tmphm(i+1,j-1) ) then
            icetmask(i,j) = 1
         endif

         if (.not.tmask(i,j)) icetmask(i,j) = 0

      enddo
      enddo

      end subroutine evp_prep1

!=======================================================================
!BOP
!
! !IROUTINE: evp_prep2 - compute quantities needed for stress tensor and mom eqns
!
! !INTERFACE:
!
      subroutine evp_prep2 (nx_block,   ny_block,
     &                      nghost,
     &                      aice,       aiu,
     &                      tmass,      umass,
     &                      umassdtei,  fcor,
     &                      tmask,      umask,
     &                      uocn,       vocn,       
     &                      strairx,    strairy,
     &                      ss_tltx,    ss_tlty,    
     &                      icetmask,   iceumask,
     &                      fm,
     &                      strtltx,    strtlty,
     &                      waterx,     watery,
     &                      forcex,     forcey,
     &                      stressp_1,  stressp_2,  
     &                      stressp_3,  stressp_4,
     &                      stressm_1,  stressm_2,  
     &                      stressm_3,  stressm_4,
     &                      stress12_1, stress12_2, 
     &                      stress12_3, stress12_4,
     &                      uvel,       vvel)

! !DESCRIPTION:
!
! Computes quantities needed in the stress tensor (sigma)
! and momentum (u) equations, but which do not change during
! the thermodynamics/transport time step:
! --wind stress shift to U grid,
! --ice mass and ice extent masks,
! --pressure (strength), and part of the forcing stresses
! initializes ice velocity for new points to ocean sfc current
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   aice        ! concentration of ice
     &,  tmass       ! total mass of ice and snow

      logical (kind=log_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   tmask       ! land/boundary mask, thickness (T-cell)
     &,  umask       ! land/boundary mask, thickness (U-cell)

      integer (kind=int_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   icetmask    ! ice extent mask (T-cell)

      logical (kind=log_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   iceumask    ! ice extent mask (U-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aiu         ! ice fraction on u-grid
     &,  umass       ! total mass of ice and snow (u grid)
     &,  fcor        ! Coriolis parameter (1/s)
     &,  strairx     ! stress on ice by air, x-direction
     &,  strairy     ! stress on ice by air, y-direction
     &,  uocn        ! ocean current, x-direction (m/s)
     &,  vocn        ! ocean current, y-direction (m/s)
     &,  ss_tltx     ! sea surface slope, x-direction (m/m)
     &,  ss_tlty     ! sea surface slope, y-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(out) ::
     &   umassdtei   ! mass of U-cell/dte (kg/m^2 s)
     &,  waterx      ! for ocean stress calculation, x (m/s)
     &,  watery      ! for ocean stress calculation, y (m/s)
     &,  forcex      ! work array: combined atm stress and ocn tilt, x
     &,  forcey      ! work array: combined atm stress and ocn tilt, y

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   fm          ! Coriolis param. * mass in U-cell (kg/s)
     &,  stressp_1, stressp_2, stressp_3, stressp_4   ! sigma11+sigma22
     &,  stressm_1, stressm_2, stressm_3, stressm_4   ! sigma11-sigma22
     &,  stress12_1,stress12_2,stress12_3,stress12_4  ! sigma12
     &,  uvel        ! x-component of velocity (m/s)
     &,  vvel        ! y-component of velocity (m/s)
     &,  strtltx     ! stress due to sea surface slope, x-direction
     &,  strtlty     ! stress due to sea surface slope, y-direction
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, k, n
     &,  ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      logical (kind=log_kind), dimension(nx_block,ny_block) ::
     &   iceumask_old      ! old-time iceumask

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      waterx   (:,:) = c0
      watery   (:,:) = c0
      forcex   (:,:) = c0
      forcey   (:,:) = c0
      umassdtei(:,:) = c0

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do j = jlo, jhi
      do i = ilo, ihi

         ! convenient variable for evp
         umassdtei(i,j) = umass(i,j)*dtei ! m/dte, kg/m^2 s

         ! ice extent mask (U-cells)
         iceumask_old(i,j) = iceumask(i,j) ! save
         iceumask(i,j) = (umask(i,j)) .and. (aiu  (i,j) > a_min)
     &                                .and. (umass(i,j) > m_min)

         ! pressure and forcing terms

         if (icetmask(i,j) == 1) then
            fm(i,j) = fcor(i,j)*umass(i,j)   ! Coriolis * mass

            if (umask(i,j)) then
           
               ! for ocean stress
               waterx(i,j) = uocn(i,j)*cosw - vocn(i,j)*sinw
               watery(i,j) = vocn(i,j)*cosw + uocn(i,j)*sinw

               ! combine tilt with wind stress
#ifndef coupled
               ! calculate tilt from geostrophic currents if needed
               strtltx(i,j) = -fm(i,j)*vocn(i,j)
               strtlty(i,j) =  fm(i,j)*uocn(i,j)
#else
               strtltx(i,j) = -gravit*umass(i,j)*ss_tltx(i,j)
               strtlty(i,j) = -gravit*umass(i,j)*ss_tlty(i,j)
#endif
               forcex(i,j) = strairx(i,j) + strtltx(i,j)
               forcey(i,j) = strairy(i,j) + strtlty(i,j)
               
            endif               ! umask
         endif                  ! icetmask

         ! initialize velocity for new ice points to ocean sfc current
         if( iceumask(i,j) .and. (.not. iceumask_old(i,j))) then
            uvel(i,j) = uocn(i,j)
            vvel(i,j) = vocn(i,j)
         endif

         ! set velocity to zero for masked-out points
         if (.not.iceumask(i,j)) then           
            uvel(i,j) = c0
            vvel(i,j) = c0
         endif

      enddo
      enddo

      !-----------------------------------------------------------------
      ! Zero out stress terms where icetmask = 0
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         if (icetmask(i,j)==0) then
            stressp_1 (i,j) = c0
            stressp_2 (i,j) = c0
            stressp_3 (i,j) = c0
            stressp_4 (i,j) = c0
            stressm_1 (i,j) = c0
            stressm_2 (i,j) = c0
            stressm_3 (i,j) = c0
            stressm_4 (i,j) = c0
            stress12_1(i,j) = c0
            stress12_2(i,j) = c0
            stress12_3(i,j) = c0
            stress12_4(i,j) = c0
         endif                  ! icetmask
      enddo                     ! i
      enddo                     ! j

      end subroutine evp_prep2

!=======================================================================
!BOP
!
! !IROUTINE: stress - computes strain rates and internal stress components
!
! !INTERFACE:
!
      subroutine stress (nx_block,   ny_block,
     &                   nghost,     ksub,
     &                   icetmask,
     &                   uvel,       vvel,
     &                   dxt,        dyt,
     &                   dxhy,       dyhx,
     &                   cxp,        cyp,
     &                   cxm,        cym,
     &                   tarear,     tinyarea,
     &                   strength,
     &                   stressp_1,  stressp_2,
     &                   stressp_3,  stressp_4,
     &                   stressm_1,  stressm_2,
     &                   stressm_3,  stressm_4,
     &                   stress12_1, stress12_2, 
     &                   stress12_3, stress12_4,
     &                   shear,      divu,
     &                   Delt,       prs_sig,
     &                   rdg_conv,   rdg_shear,
     &                   str )
!
! !DESCRIPTION:
!
! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: 
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells
     &,  ksub                ! subcycling step

      integer (kind=log_kind), dimension (nx_block,ny_block),
     &   intent(in)::
     &   icetmask    ! ice extent mask (T-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   strength     ! ice strength (N/m)
     &,  uvel         ! x-component of velocity (m/s)
     &,  vvel         ! y-component of velocity (m/s)
     &,  dxt          ! width of T-cell through the middle (m)
     &,  dyt          ! height of T-cell through the middle (m)
     &,  dxhy         ! 0.5*(HTE - HTE)
     &,  dyhx         ! 0.5*(HTN - HTN)
     &,  cyp          ! 1.5*HTE - 0.5*HTE
     &,  cxp          ! 1.5*HTN - 0.5*HTN
     &,  cym          ! 0.5*HTE - 1.5*HTE
     &,  cxm          ! 0.5*HTN - 1.5*HTN
     &,  tarear       ! 1/tarea
     &,  tinyarea     ! puny*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   stressp_1, stressp_2, stressp_3, stressp_4   ! sigma11+sigma22
     &,  stressm_1, stressm_2, stressm_3, stressm_4   ! sigma11-sigma22
     &,  stress12_1,stress12_2,stress12_3,stress12_4  ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   prs_sig      ! replacement pressure, for stress calc
     &,  shear        ! strain rate II component (1/s)
     &,  divu         ! strain rate I component, velocity divergence (1/s)
     &,  Delt         ! function of strain rates (1/s)
     &,  rdg_conv     ! convergence term for ridging (1/s)
     &,  rdg_shear    ! shear term for ridging (1/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8),
     &   intent(inout) ::
     &   str          ! stress combinations
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j
     &,  icellt       ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &   indxti       ! compressed index in i-direction
     &,  indxtj       ! compressed index in j-direction

      real (kind=dbl_kind) ::
     &  divune, divunw, divuse, divusw             ! divergence
     &, tensionne, tensionnw, tensionse, tensionsw ! tension
     &, shearne, shearnw, shearse, shearsw         ! shearing
     &, Deltane, Deltanw, Deltase, Deltasw         ! Delt
     &, c0ne, c0nw, c0se, c0sw                     ! useful combinations
     &, c1ne, c1nw, c1se, c1sw
     &, ssigpn, ssigps, ssigpe, ssigpw
     &, ssigmn, ssigms, ssigme, ssigmw
     &, ssig12n, ssig12s, ssig12e, ssig12w
     &, ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122
     &, csigpne, csigpnw, csigpse, csigpsw
     &, csigmne, csigmnw, csigmse, csigmsw
     &, csig12ne, csig12nw, csig12se, csig12sw
     &, str12ew, str12we, str12ns, str12sn
     &, strp_tmp, strm_tmp

      integer (kind=int_kind) ::
     &  ij      ! loop index, combination of i and j loops
     &, ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      rdg_conv (:,:) = c0 
      rdg_shear(:,:) = c0 
      divu (:,:) = c0 
      Delt (:,:) = c0 
      shear(:,:) = c0 
      prs_sig(:,:) = c0 

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      !-----------------------------------------------------------------
      ! Identify cells where icetmask = 1
      ! Note: The icellt mask includes north and east ghost cells
      !       where stresses are needed.
      !-----------------------------------------------------------------

      icellt = 0
      do j = jlo, jhi+1
      do i = ilo, ihi+1
         if (icetmask(i,j) == 1) then
            icellt = icellt + 1
            indxti(icellt) = i
            indxtj(icellt) = j
         endif
      enddo
      enddo

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         ! divergence  =  e_11 + e_22
         divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  )
     &             + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
         divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  )
     &             + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
         divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1)
     &             + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
         divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1)
     &             + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )

         ! tension strain rate  =  e_11 - e_22
         tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  )
     &             +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
         tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  )
     &             +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
         tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1)
     &             +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
         tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1)
     &             +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )

         ! shearing strain rate  =  e_12
         shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  )
     &           -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
         shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  )
     &           -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
         shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1)
     &           -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
         shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1)
     &           -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )
         
         ! Delta (in the denominator of zeta, eta)
         Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
         Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
         Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
         Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))

      !-----------------------------------------------------------------
      ! on last subcycle, save quantities for mechanical redistribution
      !-----------------------------------------------------------------
         if (ksub == ndte) then
            divu(i,j) = p25*(divune + divunw + divuse + divusw)
     &                * tarear(i,j)
            Delt(i,j) = p25*(Deltane + Deltanw + Deltase + Deltasw)
     &                * tarear(i,j)

! The following terms are more general and can be defined for both
! ellipse and teardrop yield curves.

            rdg_conv(i,j)  = -min(divu(i,j),c0)
            rdg_shear(i,j) = p5*(Delt(i,j)-abs(divu(i,j))) 

            ! diagnostic only
            ! shear = sqrt(tension**2 + shearing**2)
            shear(i,j) = p25*tarear(i,j)*sqrt(
     &           (tensionne + tensionnw + tensionse + tensionsw)**2
     &          +  (shearne +   shearnw +   shearse +   shearsw)**2)
         endif

      !-----------------------------------------------------------------
      ! replacement pressure/Delta                   ! kg/s
      ! save replacement pressure for principal stress calculation
      !-----------------------------------------------------------------
         if (evp_damping) then
        ! enforce damping criterion
            c0ne = min(strength(i,j)/max(Deltane,c4*tinyarea(i,j)),rcon)
            c0nw = min(strength(i,j)/max(Deltanw,c4*tinyarea(i,j)),rcon)
            c0sw = min(strength(i,j)/max(Deltasw,c4*tinyarea(i,j)),rcon)
            c0se = min(strength(i,j)/max(Deltase,c4*tinyarea(i,j)),rcon)
            prs_sig(i,j) = strength(i,j)*
     &                     Deltane/max(Deltane,c4*tinyarea(i,j)) ! ne
         else
            ! original version
            c0ne = strength(i,j)/max(Deltane,tinyarea(i,j))
            c0nw = strength(i,j)/max(Deltanw,tinyarea(i,j))
            c0sw = strength(i,j)/max(Deltasw,tinyarea(i,j))
            c0se = strength(i,j)/max(Deltase,tinyarea(i,j))
            prs_sig(i,j) = c0ne*Deltane ! northeast
         endif

         c1ne = c0ne*dte2T
         c1nw = c0nw*dte2T
         c1sw = c0sw*dte2T
         c1se = c0se*dte2T

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

         stressp_1(i,j) = (stressp_1(i,j) + c1ne*(divune - Deltane))
     &                    * denom1
         stressp_2(i,j) = (stressp_2(i,j) + c1nw*(divunw - Deltanw))
     &                    * denom1
         stressp_3(i,j) = (stressp_3(i,j) + c1sw*(divusw - Deltasw))
     &                    * denom1
         stressp_4(i,j) = (stressp_4(i,j) + c1se*(divuse - Deltase))
     &                    * denom1

         stressm_1(i,j) = (stressm_1(i,j) + c1ne*tensionne) * denom2
         stressm_2(i,j) = (stressm_2(i,j) + c1nw*tensionnw) * denom2
         stressm_3(i,j) = (stressm_3(i,j) + c1sw*tensionsw) * denom2
         stressm_4(i,j) = (stressm_4(i,j) + c1se*tensionse) * denom2

         stress12_1(i,j) = (stress12_1(i,j) + c1ne*shearne*p5) * denom2
         stress12_2(i,j) = (stress12_2(i,j) + c1nw*shearnw*p5) * denom2
         stress12_3(i,j) = (stress12_3(i,j) + c1sw*shearsw*p5) * denom2
         stress12_4(i,j) = (stress12_4(i,j) + c1se*shearse*p5) * denom2

      ! eliminate underflows
!      stressp_1(i,j) = sign(max(abs(stressp_1(i,j)),puny),stressp_1(i,j))
!      stressp_2(i,j) = sign(max(abs(stressp_2(i,j)),puny),stressp_2(i,j))
!      stressp_3(i,j) = sign(max(abs(stressp_3(i,j)),puny),stressp_3(i,j))
!      stressp_4(i,j) = sign(max(abs(stressp_4(i,j)),puny),stressp_4(i,j))

!      stressm_1(i,j) = sign(max(abs(stressm_1(i,j)),puny),stressm_1(i,j))
!      stressm_2(i,j) = sign(max(abs(stressm_2(i,j)),puny),stressm_2(i,j))
!      stressm_3(i,j) = sign(max(abs(stressm_3(i,j)),puny),stressm_3(i,j))
!      stressm_4(i,j) = sign(max(abs(stressm_4(i,j)),puny),stressm_4(i,j))

!      stress12_1(i,j) = sign(max(abs(stress12_1(i,j)),puny),stress12_1(i,j))
!      stress12_2(i,j) = sign(max(abs(stress12_2(i,j)),puny),stress12_2(i,j))
!      stress12_3(i,j) = sign(max(abs(stress12_3(i,j)),puny),stress12_3(i,j))
!      stress12_4(i,j) = sign(max(abs(stress12_4(i,j)),puny),stress12_4(i,j))

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

         ssigpn  = stressp_1(i,j) + stressp_2(i,j)
         ssigps  = stressp_3(i,j) + stressp_4(i,j)
         ssigpe  = stressp_1(i,j) + stressp_4(i,j)
         ssigpw  = stressp_2(i,j) + stressp_3(i,j)
         ssigp1  =(stressp_1(i,j) + stressp_3(i,j))*p055
         ssigp2  =(stressp_2(i,j) + stressp_4(i,j))*p055

         ssigmn  = stressm_1(i,j) + stressm_2(i,j)
         ssigms  = stressm_3(i,j) + stressm_4(i,j)
         ssigme  = stressm_1(i,j) + stressm_4(i,j)
         ssigmw  = stressm_2(i,j) + stressm_3(i,j)
         ssigm1  =(stressm_1(i,j) + stressm_3(i,j))*p055
         ssigm2  =(stressm_2(i,j) + stressm_4(i,j))*p055

         ssig12n = stress12_1(i,j) + stress12_2(i,j)
         ssig12s = stress12_3(i,j) + stress12_4(i,j)
         ssig12e = stress12_1(i,j) + stress12_4(i,j)
         ssig12w = stress12_2(i,j) + stress12_3(i,j)
         ssig121 =(stress12_1(i,j) + stress12_3(i,j))*p111
         ssig122 =(stress12_2(i,j) + stress12_4(i,j))*p111

         csigpne = p111*stressp_1(i,j) + ssigp2 + p027*stressp_3(i,j)
         csigpnw = p111*stressp_2(i,j) + ssigp1 + p027*stressp_4(i,j)
         csigpsw = p111*stressp_3(i,j) + ssigp2 + p027*stressp_1(i,j)
         csigpse = p111*stressp_4(i,j) + ssigp1 + p027*stressp_2(i,j)
         
         csigmne = p111*stressm_1(i,j) + ssigm2 + p027*stressm_3(i,j)
         csigmnw = p111*stressm_2(i,j) + ssigm1 + p027*stressm_4(i,j)
         csigmsw = p111*stressm_3(i,j) + ssigm2 + p027*stressm_1(i,j)
         csigmse = p111*stressm_4(i,j) + ssigm1 + p027*stressm_2(i,j)
         
         csig12ne = p222*stress12_1(i,j) + ssig122
     &            + p055*stress12_3(i,j)
         csig12nw = p222*stress12_2(i,j) + ssig121
     &            + p055*stress12_4(i,j)
         csig12sw = p222*stress12_3(i,j) + ssig122
     &            + p055*stress12_1(i,j)
         csig12se = p222*stress12_4(i,j) + ssig121
     &            + p055*stress12_2(i,j)

         str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         str(i,j,1) = -strp_tmp - strm_tmp - str12ew
     &        + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         str(i,j,2) = strp_tmp + strm_tmp - str12we
     &        + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str(i,j,3) = -strp_tmp - strm_tmp + str12ew
     &        + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         str(i,j,4) = strp_tmp + strm_tmp + str12we
     &        + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str(i,j,5) = -strp_tmp + strm_tmp - str12ns
     &        - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         str(i,j,6) = strp_tmp - strm_tmp - str12sn
     &        - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str(i,j,7) = -strp_tmp + strm_tmp + str12ns
     &        - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         str(i,j,8) = strp_tmp - strm_tmp + str12sn
     &        - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

      enddo                     ! ij


      end subroutine stress

!=======================================================================
!BOP
!
! !IROUTINE: stress_teardrop - strain rates and internal stress (teardrop)
!
! !INTERFACE:
!
      subroutine stress_teardrop
     &                  (nx_block,   ny_block,
     &                   nghost,     ksub,
     &                   icetmask,
     &                   uvel,       vvel,
     &                   dxt,        dyt,
     &                   dxhy,       dyhx,
     &                   cxp,        cyp,
     &                   cxm,        cym,
     &                   tarear,     tinyarea,
     &                   strength,
     &                   stressp_1,  stressp_2,
     &                   stressp_3,  stressp_4,
     &                   stressm_1,  stressm_2,
     &                   stressm_3,  stressm_4,
     &                   stress12_1, stress12_2, 
     &                   stress12_3, stress12_4,
     &                   shear,      divu,
     &                   rdg_conv,   rdg_shear,
     &                   str )
!
! !DESCRIPTION:
!
! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell, assuming a teardrop 
! yield curve.
! 
! This subroutine is based on subroutine stress_ellipse, which
! assumes an elliptical yield curve.
!
! NOTE: This subroutine still needs to be tested and validated!  
!       (As of May 2005)
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!
! !USES
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: 
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells
     &,  ksub                ! subcycling step

      integer (kind=log_kind), dimension (nx_block,ny_block),
     &   intent(in)::
     &   icetmask    ! ice extent mask (T-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   strength     ! ice strength (N/m)
     &,  uvel         ! x-component of velocity (m/s)
     &,  vvel         ! y-component of velocity (m/s)
     &,  dxt          ! width of T-cell through the middle (m)
     &,  dyt          ! height of T-cell through the middle (m)
     &,  dxhy         ! 0.5*(HTE - HTE)
     &,  dyhx         ! 0.5*(HTN - HTN)
     &,  cyp          ! 1.5*HTE - 0.5*HTE
     &,  cxp          ! 1.5*HTN - 0.5*HTN
     &,  cym          ! 0.5*HTE - 1.5*HTE
     &,  cxm          ! 0.5*HTN - 1.5*HTN
     &,  tarear       ! 1/tarea
     &,  tinyarea     ! puny*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   stressp_1, stressp_2, stressp_3, stressp_4   ! sigma11+sigma22
     &,  stressm_1, stressm_2, stressm_3, stressm_4   ! sigma11-sigma22
     &,  stress12_1,stress12_2,stress12_3,stress12_4  ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(out) ::
     &   shear        ! strain rate II component (1/s)
     &,  divu         ! strain rate I component, velocity divergence (1/s)
     &,  rdg_conv     ! convergence term for ridging (1/s)
     &,  rdg_shear    ! shear term for ridging (1/s)
!     &   prs_sig      ! replacement pressure, for stress calc
                       ! not sure if this is needed

      real (kind=dbl_kind), dimension(nx_block,ny_block,8),
     &   intent(inout) ::
     &   str          ! stress combinations
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j
     &,  icellt       ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &   indxti       ! compressed index in i-direction
     &,  indxtj       ! compressed index in j-direction

      real (kind=dbl_kind) ::
     &  divune, divunw, divuse, divusw             ! divergence
     &, tensionne, tensionnw, tensionse, tensionsw ! tension
     &, shearne, shearnw, shearse, shearsw         ! shearing

      real (kind=dbl_kind) ::
     &  totshearne, totshearnw      ! total shearing = 
     &, totshearse, totshearsw      !    sqrt(tension^2 + 4*shear^2) 
     &, kne, knw, kse, ksw          ! k = e_I/e_II in ZR05 
     &, gammaDD_ne, gammaDD_nw      ! gamma*divu 
     &, gammaDD_se, gammaDD_sw 
     &, gammaDT_ne, gammaDT_nw      ! gamma*tension 
     &, gammaDT_se, gammaDT_sw 
     &, gammaDS_ne, gammaDS_nw      ! gamma*shear 
     &, gammaDS_se, gammaDS_sw 
     &, mindiv                      ! minimum divergence rate 
     &, minshear                    ! minimum shear rate 
     &, str2i                       ! inverse strength 
     &, xne, xnw, xse, xsw          ! temporary variables
     &, yne, ynw, yse, ysw          ! temporary variables
     &, Mne,  Mnw,  Mse,  Msw       ! normalized energy dissipation rate 
     &, M1ne, M1nw, M1se, M1sw      ! compressive component of M
     &, M2ne, M2nw, M2se, M2sw      ! shear component of M
     &, denom, denom1, denom2       ! temporary variables 
     &, xstar        ! value of x for which gamma*div = -4*x - 6*x^2
     &, kmax         ! used to choose correct root for xstar 
 
      real (kind=dbl_kind), parameter :: 
     &  atd = c0         ! tensile stress paramater 'a' in ZR05 
     &, adenom = c1/((c1+atd)*(c1+atd))
     &, m827 = -c8/27._dbl_kind         ! in formula for gamma 
     &, p707 = 0.707106781186547_dbl_kind        ! 1/sqrt(2) 

      integer (kind=int_kind) ::
     &  ij      ! loop index, combination of i and j loops
     &, ilo,ihi,jlo,jhi     ! beginning and end of physical domain

!-----------------------------------------------------------------------
! 
! Note: The teardrop yield curve is given by 
! 
!     sigma_II/P = -(sigma_I/P - a) * [1 + a + (sigma_I/P - a)]^(1/2), 
! 
!     where a is the tensile stress parameter (no tensile stress for a = 0), 
!           sigma_I  = (1/2) * (sigma_11 + sigma_22), 
!           sigma_II = (1/2) * sqrt[sigma_11-sigma_22)^2 + 4*sigma_12^2], 
!           sigma_11, sigma_22, and sigma_12 are the 2D stress tensor elements. 
! 
! Squaring both sides, this can be rewritten as 
! 
!     F = y^2 - (c1+a)*u^2 - u^3 = 0, 
! 
!     where x = sigma_I/P, y = sigma_II/P, and u = x - a 
! 
! For more details, see the following: 
! Zhang, J., and D. A. Rothrock, 2005: Effect of sea ice rheology
! in numerical investigations of climate, J. Geophys. Res., 110, C08014,
! doi:10.1029/2004JC002599. 
!-----------------------------------------------------------------------

      !----------------------------------------------------------------- 
      ! Initialize 
      !----------------------------------------------------------------- 
 
      rdg_conv (:,:) = c0 
      rdg_shear(:,:) = c0 
      divu (:,:) = c0 
      shear(:,:) = c0 

      !-----------------------------------------------------------------
      ! Identify cells where icetmask = 1
      ! Note: This mask includes ghost cells for which stresses are needed.
      !-----------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      icellt = 0
      do j = jlo, jhi+1
      do i = ilo, ihi+1
         if (icetmask(i,j) == 1) then
            icellt = icellt + 1
            indxti(icellt) = i
            indxtj(icellt) = j
         endif
      enddo
      enddo

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         ! divergence  =  e_11 + e_22
         divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  )
     &             + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
         divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  )
     &             + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
         divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1)
     &             + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
         divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1)
     &             + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )

         ! tension strain rate  =  e_11 - e_22
         tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  )
     &             +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
         tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  )
     &             +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
         tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1)
     &             +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
         tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1)
     &             +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )

         ! shearing strain rate  =  e_12
         shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  )
     &           -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
         shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  )
     &           -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
         shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1)
     &           -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
         shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1)
     &           -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )
         
      ! minimum deformation rates 
      !lipscomb - might want to preserve sign of small values?
 
         mindiv = tinyarea(i,j)  ! deformation rates have units of area 
         minshear = p707*mindiv  ! so k =1 if all strain rates are tiny 
 
         if (abs(divune) < mindiv) divune = mindiv 
         if (abs(divunw) < mindiv) divunw = mindiv 
         if (abs(divuse) < mindiv) divuse = mindiv 
         if (abs(divusw) < mindiv) divusw = mindiv 
 
         if (abs(tensionne) < minshear) tensionne = minshear 
         if (abs(tensionnw) < minshear) tensionnw = minshear 
         if (abs(tensionse) < minshear) tensionse = minshear 
         if (abs(tensionsw) < minshear) tensionsw = minshear 
 
         if (abs(shearne) < minshear) shearne = minshear 
         if (abs(shearnw) < minshear) shearnw = minshear 
         if (abs(shearse) < minshear) shearse = minshear 
         if (abs(shearsw) < minshear) shearsw = minshear 

!lipscomb - debug
         if (my_task==mtest .and. i==itest .and. j==jtest
     &                      .and. ksub==ndte) then 
            print*, '' 
            print*, 'Strain rates, ksub =', ksub 
            print*, '' 
            print*, 'divu:' 
            print*, divunw, divune 
            print*, divusw, divuse 
 
            print*, '' 
            print*, 'tension:' 
            print*, tensionnw, tensionne 
            print*, tensionsw, tensionse 
 
            print*, '' 
            print*, 'shear:' 
            print*, shearnw, shearne
            print*, shearsw, shearse
 
         endif 
 
!lipscomb - EVP code has replacement pressure at this point.
!           Not sure how this should be handled.

      ! total shear, e_II = sqrt(tension^2 + shear^2)
         totshearne = sqrt(tensionne**2 + shearne**2)
         totshearnw = sqrt(tensionnw**2 + shearnw**2)
         totshearsw = sqrt(tensionsw**2 + shearsw**2)
         totshearse = sqrt(tensionse**2 + shearse**2)

      ! k = e_I/e_II, where e_I  = e_11 + e_22 = divu 
      !                     e_II = sqrt[(e_11-e_22)^2 + 4*e_12^2]
 
         kne = divune / totshearne 
         knw = divunw / totshearnw
         ksw = divusw / totshearsw 
         kse = divuse / totshearse 
 
         ! Require k <=1.  (Stresses go to zero for k > 1) 
         kne = min(kne, c1) 
         knw = min(knw, c1) 
         ksw = min(ksw, c1) 
         kse = min(kse, c1) 
 
         ! Require k >= -1000.  (To avoid roundoff errors in gamma) 
         kne = max(kne, -c1000) 
         knw = max(knw, -c1000) 
         ksw = max(ksw, -c1000) 
         kse = max(kse, -c1000) 

         ! gamma*divu 
         ! Note: gamma*divu has the same sign as divu: 
         !       positive for divergence (k > 0) and negative 
         !       for convergence (k < 0) 
 
         gammaDD_ne = m827 * kne 
     &              * (c2*kne**3 + (c2*kne**2 - c3*(c1-atd)) * 
     &                            sqrt(kne**2 + c3*(c1+atd)) ) 
 
         gammaDD_nw = m827 * knw 
     &              * (c2*knw**3 + (c2*knw**2 - c3*(c1-atd)) * 
     &                            sqrt(knw**2 + c3*(c1+atd)) ) 
 
         gammaDD_sw = m827 * ksw 
     &              * (c2*ksw**3 + (c2*ksw**2 - c3*(c1-atd)) * 
     &                            sqrt(ksw**2 + c3*(c1+atd)) ) 
 
         gammaDD_se = m827 * kse 
     &              * (c2*kse**3 + (c2*kse**2 - c3*(c1-atd)) * 
     &                            sqrt(kse**2 + c3*(c1+atd)) ) 

      !----------------------------------------------------------------- 
      ! Advance the stresses                            ! kg/s^2 
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast 
      !----------------------------------------------------------------- 

         str2i = c1 / (c2*strength(i,j) + puny) 
 
!lipscomb - Not sure that stressp is treated in the best way.
!
! From the flow rule we have 
! 
! gamma*divu = f(x) where f(x) = -4*(1+a)*x - 6*x^2 
!                     and    x = sigma_1/(2P) - a 
! 
! We want to move x in the direction that will result in 
! f(x) closer to gamma*divu.  This may be either to the right 
! or to the left, depending on x.  The graph of gamma*divu versus x
! has a maximum at x = -1/3.  Between x=-2/3 and x=0, a given value
! of gamma*divu corresponds to two different values of x.
!
! This means that if we take dx/dt proportional to the difference
! between gamma*divu and f(x), we can move x in the wrong direction.
! 
! To avoid this problem, try dx/dt = (1/T)*(xstar-x), 
! where xstar is the value of x corresponding to gamma*divu
! (chosen correctly based on k). 
! Not sure if this equation is scaled correctly.

      ! stressp = sig11 + sig22 

         kmax = sqrt(c3*(c1+atd)/c8) ! value of k with max value of gamma*divu 
         denom = c1 / (c1 + dteT) 
 
         if (kne >= kmax) then 
            xstar = -p333*(c1 - sqrt(c1-c1p5*gammaDD_ne*adenom)) 
         else 
            xstar = -p333*(c1 + sqrt(c1-c1p5*gammaDD_ne*adenom)) 
         endif 
         xne = stressp_1(i,j)*str2i - atd                ! sigma_1 to x
         xne = (xne + dteT*xstar)*denom                  ! advance x
         stressp_1(i,j) = c2*strength(i,j)*(xne + atd)   ! x to sigma_1 
 
!lipscomb - debug
         if (my_task==999 .and. i==itest .and. j==jtest
     &                      .and. ksub==ndte) then 
            print*, '' 
            print*,'i, j', i, j 
            print*, 'kne =', kne 
            print*, '' 
            print*, 'gamma*DD =', gammaDD_ne 
            print*, '-4x - 6x^2 =', -c4*xne - c6*xne**2 
            print*, '' 
            print*, 'xstar =', xstar 
            print*, 'x1_ne =', xne 
         endif 

         if (knw >= kmax) then 
            xstar = -p333*(c1 - sqrt(c1-c1p5*gammaDD_nw*adenom)) 
         else 
            xstar = -p333*(c1 + sqrt(c1-c1p5*gammaDD_nw*adenom)) 
         endif 
         xnw = stressp_2(i,j)*str2i - atd 
         xnw = (xnw + dteT*xstar)*denom 
         stressp_2(i,j) = c2*strength(i,j)*(xnw + atd) 
 
         if (ksw >= kmax) then 
            xstar = -p333*(c1 - sqrt(c1-c1p5*gammaDD_sw*adenom)) 
         else  ! kne >= c0 
            xstar = -p333*(c1 + sqrt(c1-c1p5*gammaDD_sw*adenom)) 
         endif 
         xsw = stressp_3(i,j)*str2i - atd 
         xsw = (xsw + dteT*xstar)*denom 
         stressp_3(i,j) = c2*strength(i,j)*(xsw + atd) 
 
         if (kse >= kmax) then 
            xstar = -p333*(c1 - sqrt(c1-c1p5*gammaDD_se*adenom)) 
         else  ! kne >= c0 
            xstar = -p333*(c1 + sqrt(c1-c1p5*gammaDD_se*adenom)) 
         endif 
         xse = stressp_4(i,j)*str2i - atd 
         xse = (xse + dteT*xstar)*denom 
         stressp_4(i,j) = c2*strength(i,j)*(xse + atd) 

!lipscomb - bug check 
 
         if (xne > atd) then 
            print*, 'my_task, i, j, xne:', my_task, i, j, xne 
!            stop 
         endif 
 
         if (xnw > atd) then 
            print*, 'my_task, i, j, xnw:', my_task, i, j, xnw 
!            stop 
         endif 
 
         if (xse > atd) then 
            print*, 'my_task, i, j, xse:', my_task, i, j, xse 
!            stop 
         endif 
 
         if (xsw > atd) then
            print*, 'my_task, i, j, xsw:', my_task, i, j, xsw
!            stop
         endif 

      ! stressm = sig11 - sig22 = sigma_II*cos(theta) 
 
         denom = c1/(c1+c4*dteT) 
 
         gammaDT_ne = gammaDD_ne * tensionne/divune 
         gammaDT_nw = gammaDD_nw * tensionnw/divunw 
         gammaDT_sw = gammaDD_sw * tensionsw/divusw 
         gammaDT_se = gammaDD_se * tensionse/divuse 
 
         xne = stressm_1(i,j)*str2i 
         xnw = stressm_2(i,j)*str2i 
         xsw = stressm_3(i,j)*str2i 
         xse = stressm_4(i,j)*str2i 
 
         xne = (xne + dteT*gammaDT_ne) * denom 
         xnw = (xnw + dteT*gammaDT_nw) * denom 
         xsw = (xsw + dteT*gammaDT_sw) * denom 
         xse = (xse + dteT*gammaDT_se) * denom 
 
         stressm_1(i,j) = c2*strength(i,j)*xne 
         stressm_2(i,j) = c2*strength(i,j)*xnw 
         stressm_3(i,j) = c2*strength(i,j)*xsw 
         stressm_4(i,j) = c2*strength(i,j)*xse 
 
         if (my_task==999 .and. i==itest .and. j==jtest 
     &                      .and. ksub==ndte) then 
            print*, '' 
            print*, 'gamma*DT/4 =', gammaDT_ne/c4 
            print*, 'x2_ne =', xne 
         endif 

      ! stress12 = sig12 = sigma_II*sin(theta)/2 
 
         denom = c1/(c1+c8*dteT) 
 
         gammaDS_ne = gammaDD_ne * shearne/divune 
         gammaDS_nw = gammaDD_nw * shearnw/divunw 
         gammaDS_sw = gammaDD_sw * shearsw/divusw 
         gammaDS_se = gammaDD_se * shearse/divuse 
 
         xne = stress12_1(i,j)*str2i 
         xnw = stress12_2(i,j)*str2i 
         xsw = stress12_3(i,j)*str2i 
         xse = stress12_4(i,j)*str2i 
 
         xne = (xne + dteT*gammaDS_ne) * denom 
         xnw = (xnw + dteT*gammaDS_nw) * denom 
         xsw = (xsw + dteT*gammaDS_sw) * denom 
         xse = (xse + dteT*gammaDS_se) * denom 
 
         stress12_1(i,j) = c2*strength(i,j)*xne 
         stress12_2(i,j) = c2*strength(i,j)*xnw 
         stress12_3(i,j) = c2*strength(i,j)*xsw 
         stress12_4(i,j) = c2*strength(i,j)*xse 

         if (my_task==999 .and. i==itest .and. j==jtest 
     &                      .and. ksub==ndte) then 
            print*, '' 
            print*, 'gamma*DS/8 =', gammaDS_ne/c8 
            print*, 'x12_ne =', xne 
            print*, '' 
            print*, 'gamma*eII/4 =', gammaDD_ne / (c4*kne) 
            yne = sqrt(stressm_1(i,j)**2 +c4* stress12_1(i,j)**2) 
     &          / (c2*strength(i,j)) 
            print*, 'yne =', yne 
         endif 
 
         if (my_task==mtest .and. i==itest .and. j==jtest 
     &                      .and. ksub==ndte) then 
            print*, '' 
            print*, 'New stresses:' 
            print*, 'stressp:' 
            print*, stressp_2(i,j), stressp_1(i,j) 
            print*, stressp_3(i,j), stressp_4(i,j) 
            print*, '' 
            print*, 'stressm:' 
            print*, stressm_2(i,j), stressm_1(i,j) 
            print*, stressm_3(i,j), stressm_4(i,j) 
            print*, '' 
            print*, 'stress12:' 
            print*, stress12_2(i,j), stress12_1(i,j) 
            print*, stress12_3(i,j), stress12_4(i,j) 
         endif 

!lipscomb - bug check 
         if (stressp_1(i,j) > c2*strength(i,j)*atd) 
     &        print*,'my_task, i, j, stressp_1',
     &                my_task, i, j, stressp_1(i,j)
         if (stressp_2(i,j) > c2*strength(i,j)*atd) 
     &        print*,'my_task, i, j, stressp_2', 
     &                my_task, i, j, stressp_2(i,j)
         if (stressp_3(i,j) > c2*strength(i,j)*atd)
     &        print*,'my_task, i, j, stressp_3',
     &                my_task, i, j, stressp_3(i,j)
         if (stressp_4(i,j) > c2*strength(i,j)*atd)
     &        print*,'my_task, i, j, stressp_4', 
     &                my_task, i, j, stressp_4(i,j) 
 
      !----------------------------------------------------------------- 
      ! energy dissipation terms for mechanical redistribution
      !----------------------------------------------------------------- 
 
         if (ksub==ndte) then 
 
            str2i = c1 / (strength(i,j) + puny) 
 
            ! sigma_I, compressive stress component
            ! Cf ZR05, eqn A8 
            ! Note: x = 0 for k = 1 provided atd = 0

            xne = p111*(-c6*(c1+atd) + c2*kne*kne 
     &                + c2*kne*sqrt(kne*kne+c3*(c1+atd))) 
            xnw = p111*(-c6*(c1+atd) + c2*knw*knw 
     &                + c2*knw*sqrt(knw*knw+c3*(c1+atd))) 
            xse = p111*(-c6*(c1+atd) + c2*kse*kse 
     &                + c2*kse*sqrt(kse*kse+c3*(c1+atd))) 
            xsw = p111*(-c6*(c1+atd) + c2*ksw*ksw 
     &                + c2*ksw*sqrt(ksw*ksw+c3*(c1+atd))) 
 
            ! sigma_II, shear stress component 
            ! Note: y = 0 for k = 1 provided atd = 0
 
            yne = -xne * sqrt(c1+atd+xne) 
            ynw = -xnw * sqrt(c1+atd+xnw) 
            yse = -xse * sqrt(c1+atd+xse) 
            ysw = -xsw * sqrt(c1+atd+xsw) 

!lipscomb - make sure xne <= yne, etc.

         ! normalized energy dissipation, M = sigma_I*e_I + sigma_II*e_II >= 0 
         ! cf. ellipse (Flato & Hibler):  M = (delt - e_I)/2
         ! Note: Rheology is dissipative (M > 0) provided atd = 0 

!            Mne = (xne+atd)*divune + yne*totshearne 
!            Mnw = (xnw+atd)*divunw + ynw*totshearnw 
!            Mse = (xse+atd)*divuse + yse*totshearse 
!            Msw = (xsw+atd)*divusw + ysw*totshearsw 
 
         ! As in FH95, break up M into a compressive term and a shear term.
         ! The shear term is multiplied by Cs in the ridging scheme.

!lipscomb - Two ways to do this:
!           (1) Compute M1 and M2 at each corner and average.
!           (2) Compute x, y, divu, and totshear at each corner and average;
!               combine averages to form M1 and M2.
!           Not sure which is better.
!           Divergence and shear for teardrop are larger than for ellipse;
!            not sure why.

         ! compressive part of M, >= 0 
            M1ne = (xne+atd) * min(divune, c0) 
            M1nw = (xnw+atd) * min(divunw, c0) 
            M1se = (xse+atd) * min(divuse, c0) 
            M1sw = (xsw+atd) * min(divusw, c0) 
 
         ! shear part of M, >= 0
            M2ne = yne*totshearne + (xne+atd)*max(divune,c0)
            M2nw = ynw*totshearnw + (xnw+atd)*max(divunw,c0)
            M2se = yse*totshearse + (xse+atd)*max(divuse,c0)
            M2sw = ysw*totshearsw + (xsw+atd)*max(divusw,c0)
 
            rdg_shear(i,j) = p25*tarear(i,j)*(M1ne + M1nw + M1sw + M1se) 
            rdg_conv(i,j)  = p25*tarear(i,j)*(M2ne + M2nw + M2sw + M2se) 

            if (M1ne < c0) print*, 'my_task, i, j, nweshear:', 
     &                              my_task, i, j, M1ne 
            if (M1nw < c0) print*, 'my_task, i, j, nw shear:', 
     &                              my_task, i, j, M1nw 
            if (M1se < c0) print*, 'my_task, i, j, se shear:', 
     &                              my_task, i, j, M1se 
            if (M1sw < c0) print*, 'my_task, i, j, sw shear:', 
     &                              my_task, i, j, M1sw 

            if (M2ne < c0) print*, 'my_task, i, j, nweshear:', 
     &                              my_task, i, j, M2ne 
            if (M2nw < c0) print*, 'my_task, i, j, nw shear:', 
     &                              my_task, i, j, M2nw 
            if (M2se < c0) print*, 'my_task, i, j, se shear:', 
     &                              my_task, i, j, M2se 
            if (M2sw < c0) print*, 'my_task, i, j, sw shear:', 
     &                              my_task, i, j, M2sw 

         ! shear and divergence for diagnostics 
            divu(i,j)  = p25 * tarear(i,j) 
     &                 * (divune + divunw + divuse + divusw) 
            shear(i,j) = p25 * tarear(i,j) 
     &                 * (totshearne + totshearnw 
     &                  + totshearsw + totshearsw) 
 
            if (my_task==mtest .and. i==itest .and. j==jtest 
     &                      .and. ksub==ndte) then 
               print*, '' 
               print*, 'divu =', divu(i,j)
               print*, 'shear =', shear(i,j)
               print*, 'rdg_conv =', rdg_conv(i,j) 
               print*, 'rdg_shear =', rdg_shear(i,j) 
            endif 

         endif                  ! ksub = ndte

      enddo                     ! ij

      end subroutine stress_teardrop

!=======================================================================
!BOP
!
! !IROUTINE: stepu - integrates mom eqn for u,v
!
! !INTERFACE:
!
      subroutine stepu (nx_block,   ny_block,
     &                  nghost,     iceumask,
     &                  aiu,        str,
     &                  uocn,       vocn,
     &                  waterx,     watery,
     &                  forcex,     forcey,
     &                  umassdtei,  fm,
     &                  uarear,
     &                  strocnx,    strocny,
     &                  strintx,    strinty,
     &                  uvel,       vvel)
!
! !DESCRIPTION:
!
! Calculation of the surface stresses
! Integration of the momentum equation to find velocity (u,v)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_work, only: worka, workb
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      logical (kind=log_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   iceumask    ! ice extent mask (U-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aiu         ! ice fraction on u-grid
     &,  waterx      ! for ocean stress calculation, x (m/s)
     &,  watery      ! for ocean stress calculation, y (m/s)
     &,  forcex      ! work array: combined atm stress and ocn tilt, x
     &,  forcey      ! work array: combined atm stress and ocn tilt, y
     &,  umassdtei   ! mass of U-cell/dte (kg/m^2 s)
     &,  uocn        ! ocean current, x-direction (m/s)
     &,  vocn        ! ocean current, y-direction (m/s)
     &,  fm          ! Coriolis param. * mass in U-cell (kg/s)
     &,  uarear      ! 1/uarea

      real (kind=dbl_kind), dimension(nx_block,ny_block,8),
     &   intent(in) ::
     &   str

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   uvel        ! x-component of velocity (m/s)
     &,  vvel        ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   strocnx     ! ice-ocean stress, x-direction
     &,  strocny     ! ice-ocean stress, y-direction
     &,  strintx     ! divergence of internal ice stress, x (N/m^2)
     &,  strinty     ! divergence of internal ice stress, y (N/m^2)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, ij
     &,  ilo,ihi,jlo,jhi     ! beginning and end of physical domain
     &,  icellu              ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &   indxui      ! compressed index in i-direction
     &,  indxuj      ! compressed index in j-direction

      real (kind=dbl_kind) ::
     &   uold, vold            ! old-time uvel, vvel
     &,  vrel                  ! relative ice-ocean velocity
     &,  cca,ccb,ab2,cc1,cc2   ! intermediate variables
     &,  taux, tauy            ! part of ocean stress term          

      !-----------------------------------------------------------------
      ! Identify cells where iceumask is true
      !-----------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      icellu = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (iceumask(i,j)) then
            icellu = icellu + 1
            indxui(icellu) = i
            indxuj(icellu) = j
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         uold = uvel(i,j)
         vold = vvel(i,j)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiu(i,j)*dragw*sqrt((uocn(i,j) - uold)**2 +
     &                              (vocn(i,j) - vold)**2)  ! m/s
         ! ice/ocean stress
         taux = vrel*waterx(i,j) ! NOTE this is not the entire
         tauy = vrel*watery(i,j) ! ocn stress term

         ! alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
         cca = umassdtei(i,j) + vrel * cosw      ! alpha, kg/m^2 s
         ccb = fm(i,j)        + vrel * sinw      ! beta,  kg/m^2 s
         ab2 = cca**2 + ccb**2

         ! divergence of the internal stress tensor
         strintx(i,j) = uarear(i,j)*
     &       (str(i,j,1) + str(i+1,j,2) + str(i,j+1,3) + str(i+1,j+1,4))
         strinty(i,j) = uarear(i,j)*
     &       (str(i,j,5) + str(i,j+1,6) + str(i+1,j,7) + str(i+1,j+1,8))

         ! finally, the velocity components
         cc1 = strintx(i,j) + forcex(i,j) + taux
     &       + umassdtei(i,j)*uold
         cc2 = strinty(i,j) + forcey(i,j) + tauy
     &       + umassdtei(i,j)*vold

         uvel(i,j) = (cca*cc1 + ccb*cc2) / ab2 ! m/s
         vvel(i,j) = (cca*cc2 - ccb*cc1) / ab2

      !-----------------------------------------------------------------
      ! ocean-ice stress for coupling
      ! scale to full grid cell
      !-----------------------------------------------------------------
         strocnx(i,j) = taux / aiu(i,j)
         strocny(i,j) = tauy / aiu(i,j)

      enddo                     ! ij

      end subroutine stepu

!=======================================================================
!BOP
!
! !IROUTINE: evp_finish - calculates ice-ocean stress
!
! !INTERFACE:
!
      subroutine evp_finish (nx_block, ny_block,
     &                       nghost,
     &                       uvel,     vvel,    
     &                       uocn,     vocn,   
     &                       strocnx,  strocny, 
     &                       strocnxT, strocnyT)
!
! !DESCRIPTION:
!
! Calculation of the ice-ocean stress.
! ...the sign will be reversed later...
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   uvel        ! x-component of velocity (m/s)
     &,  vvel        ! y-component of velocity (m/s)
     &,  uocn        ! ocean current, x-direction (m/s)
     &,  vocn        ! ocean current, y-direction (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(out) ::
     &   strocnx     ! ice-ocean stress, x-direction
     &,  strocny     ! ice-ocean stress, y-direction
     &,  strocnxT    ! ice-ocean stress, x-direction
     &,  strocnyT    ! ice-ocean stress, y-direction
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, ij
     &,  ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      real (kind=dbl_kind) :: vrel

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      strocnx(:,:) = c0
      strocny(:,:) = c0

      ! ocean-ice stress for coupling
      do j = jlo, jhi
      do i = ilo, ihi
         vrel = dragw*sqrt((uocn(i,j) - uvel(i,j))**2 +
     &                     (vocn(i,j) - vvel(i,j))**2)  ! m/s
         strocnx(i,j) = strocnx(i,j)
     &                - vrel*(uvel(i,j)*cosw - vvel(i,j)*sinw)
         strocny(i,j) = strocny(i,j)
     &                - vrel*(vvel(i,j)*cosw + uvel(i,j)*sinw)
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Prepare to convert to T grid
      !-----------------------------------------------------------------
      strocnxT(:,:) = strocnx(:,:)
      strocnyT(:,:) = strocny(:,:)

      end subroutine evp_finish

!=======================================================================
!BOP
!
! !IROUTINE: principal_stress - computes principal stress for yield curve
!
! !INTERFACE:
!
      subroutine principal_stress(nx_block,   ny_block,
     &                            stressp_1,  stressm_1,
     &                            stress12_1, prs_sig,
     &                            sig1,       sig2)
!
! !DESCRIPTION:
!
! Computes principal stresses for comparison with the theoretical
! yield curve; northeast values
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   stressp_1   ! sigma11 + sigma22
     &,  stressm_1   ! sigma11 - sigma22
     &,  stress12_1  ! sigma12
     &,  prs_sig     ! replacement pressure, for stress calc

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out)::
     &   sig1        ! principal stress component
     &,  sig2        ! principal stress component
!
!EOP
!
      integer (kind=int_kind) :: i, j

      do j = 1, ny_block
      do i = 1, nx_block
         if (prs_sig(i,j) > puny) then
            sig1(i,j) = (p5*(stressp_1(i,j)
     &                + sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2)))
     &                / prs_sig(i,j)
            sig2(i,j) = (p5*(stressp_1(i,j)
     &                - sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2)))
     &                / prs_sig(i,j)
         else
            sig1(i,j) = spval
            sig2(i,j) = spval
         endif
      enddo
      enddo

      end subroutine principal_stress

!=======================================================================

      end module ice_dyn_evp

!=======================================================================
