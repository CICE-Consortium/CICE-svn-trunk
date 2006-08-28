!=======================================================================
!BOP
!
! !MODULE: ice_grid - spatial grids, masks and boundary conditions
!
! !DESCRIPTION:
!
! Spatial grids, masks, and boundary conditions
!
! !REVISION HISTORY:
!  SVN:$Id$
!  SVN:$LastChangedBy$
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!          Tony Craig, NCAR
!
! 2004: Block structure added by William Lipscomb
!       init_grid split into two parts as in POP 2.0
!       Boundary update routines replaced by POP versions
!
! !INTERFACE:
!
      module ice_grid
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants
      use ice_blocks
      use ice_domain_size
      use ice_domain
      use ice_fileunits
      use ice_read_write
      use ice_timers
!
!EOP
!
      implicit none
      save

      character (len=char_len) ::
     &   grid_file  !  input file for POP grid info
     &,  kmt_file   !  input file for POP grid info
     &,  grid_type  !  current options are rectangular (default),
                    !  displaced_pole, tripole, panarctic, and column

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks)::
     &   dxt        ! width of T-cell through the middle (m)
     &,  dyt        ! height of T-cell through the middle (m)
     &,  dxu        ! width of U-cell through the middle (m)
     &,  dyu        ! height of U-cell through the middle (m)
     &,  HTE        ! length of eastern edge of T-cell (m)
     &,  HTN        ! length of northern edge of T-cell (m)
     &,  tarea      ! area of T-cell (m^2)
     &,  uarea      ! area of U-cell (m^2)
     &,  tarear     ! 1/tarea
     &,  uarear     ! 1/uarea
     &,  tinyarea   ! puny*tarea
     &,  tarean     ! area of NH T-cells
     &,  tareas     ! area of SH T-cells
     &,  ULON       ! longitude of velocity pts (radians)
     &,  ULAT       ! latitude of velocity pts (radians)
     &,  TLON       ! longitude of temp pts (radians)
     &,  TLAT       ! latitude of temp pts (radians)
     &,  ANGLE      ! for conversions between POP grid and lat/lon
     &,  ANGLET     ! ANGLE converted to T-cells

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks)::
     &   cyp        ! 1.5*HTE - 0.5*HTE
     &,  cxp        ! 1.5*HTN - 0.5*HTN
     &,  cym        ! 0.5*HTE - 1.5*HTE
     &,  cxm        ! 0.5*HTN - 1.5*HTN
     &,  dxhy       ! 0.5*(HTE - HTE)
     &,  dyhx       ! 0.5*(HTN - HTN)

      ! geometric quantities used for remapping transport
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks)::
     &   xav    ! mean T-cell value of x
     &,  yav    ! mean T-cell value of y
     &,  xxav   ! mean T-cell value of xx
     &,  xyav   ! mean T-cell value of xy
     &,  yyav   ! mean T-cell value of yy
     &,  xxxav  ! mean T-cell value of xxx
     &,  xxyav  ! mean T-cell value of xxy
     &,  xyyav  ! mean T-cell value of xyy
     &,  yyyav  ! mean T-cell value of yyy

      real (kind=dbl_kind),
     &   dimension (2,2,nx_block,ny_block,max_blocks) ::
     &   mne    ! matrices used for coordinate transformations in remapping
     &,  mnw    ! ne = northeast corner, nw = northwest, etc.
     &,  mse 
     &,  msw

      ! masks
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks)::
     &   hm         ! land/boundary mask, thickness (T-cell)
     &,  uvm        ! land/boundary mask, velocity (U-cell)

      logical (kind=log_kind), 
     &   dimension (nx_block,ny_block,max_blocks) ::
     &   tmask      ! land/boundary mask, thickness (T-cell)
     &,  umask      ! land/boundary mask, velocity (U-cell)
     &,  lmask_n    ! northern hemisphere mask
     &,  lmask_s    ! southern hemisphere mask

!lipscomb - not sure rndex_global is stil needed
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) ::
     &   rndex_global           ! global index for local subdomain (dbl)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_grid1 - - distribute blocks across processors
!
! !INTERFACE:
!
      subroutine init_grid1 
!
! !DESCRIPTION:
!
! Distribute blocks across processors.  The distribution is optimized
! based on latitude and topography, contained in the ULAT and KMT arrays. 
!
! !REVISION HISTORY:
!
! authors: William Lipscomb and Phil Jones, LANL
!
! !USES:
!
      use ice_broadcast
      use ice_work, only: work_g1, work_g2
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, iblk

      !-----------------------------------------------------------------
      ! Get global ULAT and KMT arrays used for block decomposition.
      !-----------------------------------------------------------------

      allocate(work_g1(nx_global,ny_global))
      allocate(work_g2(nx_global,ny_global))

      if (trim(grid_type) == 'displaced_pole' .or.
     &    trim(grid_type) == 'tripole'      ) then

         call ice_open(nu_grid,grid_file,64) ! ULAT
         call ice_open(nu_kmt, kmt_file, 32) ! KMT

         call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)  ! ULAT
         call ice_read_global(nu_kmt, 1,work_g2,'ida4',.true.)  ! KMT

         if (my_task == master_task) then
            close (nu_grid)
            close (nu_kmt)
         endif

      elseif (trim(grid_type) == 'panarctic') then

         call ice_open(nu_grid,grid_file,64) ! ULAT, KMT

         call ice_read_global(nu_grid,1,work_g2,'ida8',.true.)  ! KMT
         call ice_read_global(nu_grid,2,work_g1,'rda8',.true.)  ! ULAT

         if (my_task == master_task) close (nu_grid)

      elseif (trim(grid_type) == 'column') then

         work_g1(:,:) = 75./rad_to_deg  ! arbitrary polar latitude
         work_g2(:,:) = c1

      else   ! rectangular grid

         work_g1(:,:) = 75./rad_to_deg  ! arbitrary polar latitude
         work_g2(:,:) = c1

      endif

      call broadcast_array(work_g1, master_task)   ! ULAT
      call broadcast_array(work_g2, master_task)   ! KMT

      !-----------------------------------------------------------------
      ! distribute blocks among processors
      !-----------------------------------------------------------------

      call init_domain_distribution(work_g2, work_g1)  ! KMT, ULAT

      deallocate(work_g1)
      deallocate(work_g2)

      end subroutine init_grid1

!=======================================================================
!BOP
!
! !IROUTINE: init_grid2 - horizontal grid initialization
!
! !INTERFACE:
!
      subroutine init_grid2
!
! !DESCRIPTION:
!
! Horizontal grid initialization:
!
!     U{LAT,LONG} = true {latitude,longitude} of U points
!     HT{N,E} = cell widths on {N,E} sides of T cell
!     ANGLE = angle between local x direction and true east
!     hm = land mask (c1 for ocean points, c0 for land points)
!     D{X,Y}{T,U} = {x,y} spacing centered at {T,U} points
!     T-grid and ghost cell values
!     Various grid quantities needed for dynamics and transport
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_gather_scatter
      use ice_boundary
      use ice_work, only: work_g1
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, iblk
     &,  ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) ::
     &   angle_0, angle_w, angle_s, angle_sw

      logical (kind=log_kind), dimension(nx_block,ny_block,max_blocks)::
     &   out_of_range

      type (block) ::
     &   this_block           ! block information for current block
      
      !-----------------------------------------------------------------
      ! lat, lon, cell widths, angle, land mask
      !-----------------------------------------------------------------

      if (trim(grid_type) == 'displaced_pole' .or.
     &    trim(grid_type) == 'tripole'      ) then
         call popgrid           ! read POP grid lengths directly
      elseif (trim(grid_type) == 'panarctic') then
         call panarctic_grid    ! pan-Arctic grid
      elseif (trim(grid_type) == 'column') then
         call columngrid        ! column model grid
      else
         call rectgrid          ! regular rectangular grid
      endif

      call ice_timer_start(timer_bound)
      call update_ghost_cells (HTN,                bndy_info,
     &                         field_loc_Nface,    field_type_scalar)
      call update_ghost_cells (HTE,                bndy_info,
     &                         field_loc_Eface,    field_type_scalar)
      call update_ghost_cells (ULAT,               bndy_info,
     &                         field_loc_NEcorner, field_type_scalar)
      call update_ghost_cells (ULON,               bndy_info,
     &                         field_loc_NEcorner, field_type_scalar)
      call update_ghost_cells (ANGLE,              bndy_info,
     &                         field_loc_NEcorner, field_type_angle)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! T-grid cell and U-grid cell quantities
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         this_block = get_block(blocks(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            dxt(i,j,iblk) = p5*(HTN(i,j,iblk) + HTN(i,j-1,iblk))
            dyt(i,j,iblk) = p5*(HTE(i,j,iblk) + HTE(i-1,j,iblk))
            
            tarea(i,j,iblk) = dxt(i,j,iblk)*dyt(i,j,iblk)

            dxu(i,j,iblk) = p5*(HTN(i,j,iblk) + HTN(i+1,j,iblk))
            dyu(i,j,iblk) = p5*(HTE(i,j,iblk) + HTE(i,j+1,iblk))

         enddo
         enddo

      enddo                     ! iblk

      call update_ghost_cells (dxt,                bndy_info,
     &                         field_loc_center,   field_type_scalar)
      call update_ghost_cells (dyt,                bndy_info,
     &                         field_loc_center,   field_type_scalar)
      call update_ghost_cells (dxu,                bndy_info,
     &                         field_loc_NEcorner, field_type_scalar)
      call update_ghost_cells (dyu,                bndy_info,
     &                         field_loc_NEcorner, field_type_scalar)
      call update_ghost_cells (tarea,              bndy_info,
     &                         field_loc_center,   field_type_scalar)

      do iblk = 1, nblocks
         this_block = get_block(blocks(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            uarea(i,j,iblk) = p25*
     &                       (tarea(i,j,  iblk) + tarea(i+1,j,  iblk)
     &                      + tarea(i,j+1,iblk) + tarea(i+1,j+1,iblk))

            tarear(i,j,iblk) = c1/tarea(i,j,iblk)
            uarear(i,j,iblk) = c1/uarea(i,j,iblk)
            tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

            dxhy(i,j,iblk) = p5*(HTE(i,j,iblk) - HTE(i-1,j,iblk))
            dyhx(i,j,iblk) = p5*(HTN(i,j,iblk) - HTN(i,j-1,iblk))

         enddo                  ! i
         enddo                  ! j 

         do j = jlo, jhi+1
         do i = ilo, ihi+1
            cyp(i,j,iblk) = (c1p5*HTE(i,j,iblk) - p5*HTE(i-1,j,iblk))
            cxp(i,j,iblk) = (c1p5*HTN(i,j,iblk) - p5*HTN(i,j-1,iblk))
            cym(i,j,iblk) = (p5*HTE(i,j,iblk) - c1p5*HTE(i-1,j,iblk))
            cxm(i,j,iblk) = (p5*HTN(i,j,iblk) - c1p5*HTN(i,j-1,iblk))
         enddo
         enddo

      enddo                     ! iblk

      !-----------------------------------------------------------------
      ! Ghost cell updates
      ! On the tripole grid, one must be careful with updates of
      !  quantities that involve a difference of cell lengths.
      ! For example, dyhx and dxhy are cell-centered vector components.
      ! Also note that on the tripole grid, cxp and cxm would swap places,
      !  as would cyp and cym.  These quantities are computed only
      !  in north and east ghost cells (above), not south and west.
      !-----------------------------------------------------------------

      call ice_timer_start(timer_bound)
      call update_ghost_cells (uarea,            bndy_info,
     &                         field_loc_NEcorner, field_type_scalar)

      call update_ghost_cells (uarear,           bndy_info,
     &                         field_loc_NEcorner, field_type_scalar)

      call update_ghost_cells (tarear,           bndy_info,
     &                         field_loc_center, field_type_scalar)

      call update_ghost_cells (tinyarea,         bndy_info,
     &                         field_loc_center, field_type_scalar)

      call update_ghost_cells (dxhy,             bndy_info,
     &                         field_loc_center, field_type_vector)

      call update_ghost_cells (dyhx,             bndy_info,
     &                         field_loc_center, field_type_vector)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! Calculate ANGLET to be compatible with POP ocean model
      ! First, ensure that -pi <= ANGLE <= pi
      !-----------------------------------------------------------------

      out_of_range = .false.
      where (ANGLE < -pi .or. ANGLE > pi) out_of_range = .true.
      if (count(out_of_range) > 0) then
         call abort_ice ('ice: init_grid: ANGLE out of expected range')
      endif


      !-----------------------------------------------------------------
      ! Compute ANGLE on T-grid
      !-----------------------------------------------------------------
      ANGLET = c0

      do iblk = 1, nblocks
         this_block = get_block(blocks(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            angle_0  = ANGLE(i  ,j  ,iblk) !   w----0
            angle_w  = ANGLE(i-1,j  ,iblk) !   |    |
            angle_s  = ANGLE(i,  j-1,iblk) !   |    |
            angle_sw = ANGLE(i-1,j-1,iblk) !   sw---s

            if ( angle_0 < c0 ) then
               if ( abs(angle_w - angle_0) > pi)
     &                  angle_w = angle_w  - pi2
               if ( abs(angle_s - angle_0) > pi)
     &                  angle_s = angle_s  - pi2
               if ( abs(angle_sw - angle_0) > pi)
     &                  angle_sw = angle_sw - pi2
            endif

            ANGLET(i,j,iblk) = angle_0 * p25 + angle_w * p25
     &                       + angle_s * p25 + angle_sw* p25
         enddo
         enddo
      enddo
      
      call ice_timer_start(timer_bound)
      call update_ghost_cells (ANGLET,             bndy_info,
     &                         field_loc_NEcorner, field_type_angle)
      call ice_timer_stop(timer_bound)

      call makemask          ! velocity mask, hemisphere masks
      call Tlatlon           ! get lat, lon on the T grid

      !-----------------------------------------------------------------
      ! Compute global index (used for unpacking messages from coupler)
      !-----------------------------------------------------------------

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         do j=1,ny_global
         do i=1,nx_global
            work_g1(i,j) = real((j-1)*nx_global + i,kind=dbl_kind)
         enddo
         enddo
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call scatter_global(rndex_global, work_g1, 
     &                    master_task,  distrb_info,
     &                    field_loc_center, field_type_scalar)

      deallocate(work_g1)

      end subroutine init_grid2

!=======================================================================
!BOP
!
! !IROUTINE: popgrid - read and set POP displaced pole (or tripole)
!                      grid and land mask
!
! !INTERFACE:
!
      subroutine popgrid
!
! !DESCRIPTION:
!
! POP displaced pole grid and land mask. \\
! Grid record number, field and units are: \\
! (1) ULAT  (radians)    \\
! (2) ULON  (radians)    \\
! (3) HTN   (cm)         \\
! (4) HTE   (cm)         \\
! (5) HUS   (cm)         \\
! (6) HUW   (cm)         \\
! (7) ANGLE (radians)    \\
!
! Land mask record number and field is (1) KMT.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_work, only: work1
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, iblk
     &,  ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      type (block) ::
     &   this_block           ! block information for current block

      call ice_open(nu_grid,grid_file,64)
      call ice_open(nu_kmt,kmt_file,32)

      diag = .true.       ! write diagnostic info

      ! lat, lon, cell dimensions, angles
      call ice_read(nu_grid,1,ULAT, 'rda8',diag)
      call ice_read(nu_grid,2,ULON, 'rda8',diag)
      call ice_read(nu_grid,3,HTN,  'rda8',diag)
      call ice_read(nu_grid,4,HTE,  'rda8',diag)
      call ice_read(nu_grid,7,ANGLE,'rda8',diag)

      ! fix units
      HTN(:,:,:) = HTN(:,:,:) * cm_to_m
      HTE(:,:,:) = HTE(:,:,:) * cm_to_m

      ! topography
      call ice_read(nu_kmt,1,work1,'ida4',diag)

      if (my_task == master_task) then
         close (nu_grid)
         close (nu_kmt)
      endif

      hm(:,:,:) = c0
      do iblk = 1, nblocks
         this_block = get_block(blocks(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            hm(i,j,iblk) = work1(i,j,iblk)
            if (hm(i,j,iblk) >= c1) hm(i,j,iblk) = c1

         ! uncomment to mask out tropics
         ! Do this only if running uncoupled
!!!             if (ULAT(i,j,iblk) > shlat/rad_to_deg .and.
!!!     &           ULAT(i,j,iblk) < nhlat/rad_to_deg) 
!!!     &           hm(i,j,iblk) = c0
         enddo
         enddo
      enddo

      end subroutine popgrid

!=======================================================================
!BOP
!
! !IROUTINE: panarctic_grid - read and set Pan-Arctic grid and land mask
!
! !INTERFACE:
!
      subroutine panarctic_grid
!
! !DESCRIPTION:
!
! Pan-Arctic grid and mask developed by Wieslaw Maslowski
!
! !REVISION HISTORY:
!
! authors: Wieslaw Maslowki, Naval Postgraduate School (based on popgrid)
!          William H. Lipscomb, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_gather_scatter
      use ice_work, only: work1
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      !-----------------------------------------------------------------
      ! 
      ! PIPS rotated spherical grid and land mask
      !      rec no.         field         units
      !      -------         -----         -----
      !   land mask
      !         1             KMT         
      !   grid
      !         2            ULAT         radians
      !         3            ULON         radians
      !         4             HTN           cm
      !         5             HTE           cm
      !         6             HUS           cm
      !         7             HUW           cm
      !         8            ANGLE        radians
      !
      ! NOTE: There is no separate kmt file.  Land mask is part of grid file.
      !-----------------------------------------------------------------

      integer (kind=int_kind) ::
     &   i, j, iblk
     &,  ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      type (block) ::
     &   this_block           ! block information for current block

      call ice_open(nu_grid,grid_file,64)

      diag = .true.       ! write diagnostic info

      if (my_task == master_task)
     &     write (nu_diag,*) '** Reading pan-Arctic grid **'

      ! read topography
      call ice_read(nu_grid,1,work1,'ida8',diag)

      hm(:,:,:) = c0
      do iblk = 1, nblocks
         this_block = get_block(blocks(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            hm(i,j,iblk) = work1(i,j,iblk)
            if (hm(i,j,iblk) >= c1) hm(i,j,iblk) = c1
         enddo
         enddo
      enddo                     ! iblk

      ! read other grid quantities

      call ice_read(nu_grid,2,ULAT, 'rda8',diag)
      call ice_read(nu_grid,3,ULON, 'rda8',diag)
      call ice_read(nu_grid,4,HTN,  'rda8',diag)
      call ice_read(nu_grid,5,HTE,  'rda8',diag)
      call ice_read(nu_grid,8,ANGLE,'rda8',diag)

      ! fix units
      HTN(:,:,:) = HTN(:,:,:)*cm_to_m
      HTE(:,:,:) = HTE(:,:,:)*cm_to_m

      if (my_task == master_task) close (nu_grid)

      end subroutine panarctic_grid

!=======================================================================
!BOP
!
! !IROUTINE: columngrid - column grid and mask
!
! !INTERFACE:
!
      subroutine columngrid
!
! !DESCRIPTION:
!
! Column grid and mask
!
! !REVISION HISTORY:
!
! author: C. M. Bitz UW, (based on rectgrid by Hunke)
!
! modified Nov. 2003 by William H. Lipscomb, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_gather_scatter
      use ice_work, only: work_g1
      use ice_exit      
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, iblk

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block

            HTN  (i,j,iblk) = 1.6e4_dbl_kind
            HTE  (i,j,iblk) = 1.6e4_dbl_kind

            ! used to find hemisphere and init_state, need not be exact
            ULAT (i,j,iblk) = 75.0_dbl_kind/rad_to_deg  
            ULON (i,j,iblk) = 170.0_dbl_kind/rad_to_deg 

            ANGLE(i,j,iblk) = c0               ! "square with the world"

         enddo                  ! i
         enddo                  ! j
      enddo                     ! iblk

      !-----------------------------------------------------------------
      ! Verify that nx_global and ny_global are 1
      !-----------------------------------------------------------------

      if ((nx_global /= 1).or. (ny_global /= 1)) then
         write(nu_diag,*)
     &        'Because you have selected the column model flag'
         write(nu_diag,*) 'Please set nx_global=ny_global=1 in file'
         write(nu_diag,*) 'ice_domain_size.F and recompile'
         call abort_ice ('ice: columngrid: check nx_global, ny_global')
      endif

      !-----------------------------------------------------------------
      ! Construct T-cell land mask
      !-----------------------------------------------------------------

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         do j = 1, ny_global
         do i = 1, nx_global
            work_g1(i,j) = c1
         enddo
         enddo
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call scatter_global(hm, work_g1, master_task, distrb_info,
     &                    field_loc_center, field_type_scalar)

      deallocate(work_g1)

      end subroutine columngrid

!=======================================================================
!BOP
!
! !IROUTINE: rectgrid - regular rectangular grid and mask
!
! !INTERFACE:
!
      subroutine rectgrid
!
! !DESCRIPTION:
!
! Regular rectangular grid and mask
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_gather_scatter
      use ice_work, only: work_g1
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, iblk
     &,  ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) ::
     &   this_block           ! block information for current block

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            HTN  (i,j,iblk) = 1.6e4_dbl_kind  ! constant longitude spacing =
                                              ! POP <2/3> min, m
            HTE  (i,j,iblk) = 1.6e4_dbl_kind  ! constant latitude  spacing =
                                              ! POP <2/3> min, m
            ULAT (i,j,iblk) = c0              ! remember to set Coriolis !
            ULON (i,j,iblk) = c0
            ANGLE(i,j,iblk) = c0              ! "square with the world"

         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Construct T-cell land mask
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
         work_g1(:,:) = c0      ! initialize hm as land

!!!      do j=1,ny_global        ! open
!!!      do i=1,nx_global        ! open
         do j=3,ny_global-2     ! closed: NOTE ny_global > 5
            do i=3,nx_global-2  ! closed: NOTE nx_global > 5
               work_g1(i,j) = c1
            enddo
         enddo
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call scatter_global(hm, work_g1, master_task, distrb_info,
     &                    field_loc_center, field_type_scalar)

      deallocate(work_g1)

      end subroutine rectgrid

!=======================================================================
!BOP
!
! !IROUTINE: makemask - makes logical land masks (T,U) and hemispheric masks
!
! !INTERFACE:
!
      subroutine makemask
!
! !DESCRIPTION:
!
! Sets the boundary values for the T cell land mask (hm) and
! makes the logical land masks for T and U cells (tmask, umask).
! Also creates hemisphere masks (mask-n northern, mask-s southern)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_boundary
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, iblk
     &,  ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) ::
     &   this_block           ! block information for current block

      call ice_timer_start(timer_bound)
      call update_ghost_cells (hm,               bndy_info,
     &                         field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! construct T-cell and U-cell masks
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         this_block = get_block(blocks(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            uvm(i,j,iblk) = min (hm(i,j,  iblk), hm(i+1,j,  iblk),
     &                           hm(i,j+1,iblk), hm(i+1,j+1,iblk))
         enddo
         enddo
      enddo

      call ice_timer_start(timer_bound)
      call update_ghost_cells (uvm,                bndy_info,
     &                         field_loc_NEcorner, field_type_scalar)
      call ice_timer_stop(timer_bound)

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            tmask(i,j,iblk) = .false.
            umask(i,j,iblk) = .false.
            if ( hm(i,j,iblk) > p5) tmask(i,j,iblk) = .true.
            if (uvm(i,j,iblk) > p5) umask(i,j,iblk) = .true.
         enddo
         enddo

      !-----------------------------------------------------------------
      ! create hemisphere masks
      !-----------------------------------------------------------------

         lmask_n(:,:,iblk) = .false.
         lmask_s(:,:,iblk) = .false.

         tarean(:,:,iblk) = c0
         tareas(:,:,iblk) = c0

         do j = 1, ny_block
         do i = 1, nx_block

            if (ULAT(i,j,iblk) >= -puny) lmask_n(i,j,iblk) = .true. ! N. Hem.
            if (ULAT(i,j,iblk) <  -puny) lmask_s(i,j,iblk) = .true. ! S. Hem.

            ! N hemisphere area mask (m^2)
            if (lmask_n(i,j,iblk)) tarean(i,j,iblk) = tarea(i,j,iblk)

            ! S hemisphere area mask (m^2)
            if (lmask_s(i,j,iblk)) tareas(i,j,iblk) = tarea(i,j,iblk)

         enddo
         enddo


      enddo  ! iblk

      end subroutine makemask

!=======================================================================
!BOP
!
! !IROUTINE: Tlatlon - initializes latitude and longitudes on T grid
!
! !INTERFACE:
!
      subroutine Tlatlon
!
! !DESCRIPTION:
!
! Initializes latitude and longitude on T grid
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL; code originally based on POP grid
! generation routine
!
! !USES:
!
      use ice_domain_size
      use ice_boundary
      use ice_gather_scatter
      use ice_global_reductions, only: global_minval, global_maxval
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &     i, j, iblk           ! horizontal indices
     &,    ig, jg               ! global horizontal indices
     &,    im1                  ! ig - 1
     &,    ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) ::
     &     z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da

      type (block) ::
     &     this_block           ! block information for current block

      TLAT(:,:,:) = c0
      TLON(:,:,:) = c0

      do iblk = 1, nblocks
         this_block = get_block(blocks(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            z1 = cos(ULAT(i-1,j-1,iblk))
            x1 = cos(ULON(i-1,j-1,iblk))*z1
            y1 = sin(ULON(i-1,j-1,iblk))*z1
            z1 = sin(ULAT(i-1,j-1,iblk))

            z2 = cos(ULAT(i,j-1,iblk))
            x2 = cos(ULON(i,j-1,iblk))*z2
            y2 = sin(ULON(i,j-1,iblk))*z2
            z2 = sin(ULAT(i,j-1,iblk))

            z3 = cos(ULAT(i-1,j,iblk))
            x3 = cos(ULON(i-1,j,iblk))*z3
            y3 = sin(ULON(i-1,j,iblk))*z3
            z3 = sin(ULAT(i-1,j,iblk))

            z4 = cos(ULAT(i,j,iblk))
            x4 = cos(ULON(i,j,iblk))*z4
            y4 = sin(ULON(i,j,iblk))*z4
            z4 = sin(ULAT(i,j,iblk))

            tx = (x1+x2+x3+x4)/c4
            ty = (y1+y2+y3+y4)/c4
            tz = (z1+z2+z3+z4)/c4
            da = sqrt(tx**2+ty**2+tz**2)

            tz = tz/da

            ! TLON in radians East
            TLON(i,j,iblk) = c0
            if (tx /= c0 .or. ty /= c0)
     &           TLON(i,j,iblk) = atan2(ty,tx)

            ! TLAT in radians North
            TLAT(i,j,iblk) = asin(tz)
            
         enddo                  ! i
         enddo                  ! j         
      enddo                     ! iblk

      call ice_timer_start(timer_bound)
      call update_ghost_cells (TLON,             bndy_info,
     &                         field_loc_center, field_type_scalar)
      call update_ghost_cells (TLAT,             bndy_info,
     &                         field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      x1 = global_minval(TLON, distrb_info, field_loc_center, tmask)
      x2 = global_maxval(TLON, distrb_info, field_loc_center, tmask)
      x3 = global_minval(TLAT, distrb_info, field_loc_center, tmask)
      x4 = global_maxval(TLAT, distrb_info, field_loc_center, tmask)

      y1 = global_minval(ULON, distrb_info, field_loc_NEcorner, umask)
      y2 = global_maxval(ULON, distrb_info, field_loc_NEcorner, umask)
      y3 = global_minval(ULAT, distrb_info, field_loc_NEcorner, umask)
      y4 = global_maxval(ULAT, distrb_info, field_loc_NEcorner, umask)

      if (my_task==master_task) then
         write(nu_diag,*) ' '
         write(nu_diag,*) 'min/max ULON:', y1*rad_to_deg, y2*rad_to_deg
         write(nu_diag,*) 'min/max TLON:', x1*rad_to_deg, x2*rad_to_deg
         write(nu_diag,*) 'min/max ULAT:', y3*rad_to_deg, y4*rad_to_deg
         write(nu_diag,*) 'min/max TLAT:', x3*rad_to_deg, x4*rad_to_deg
      endif                     ! my_task

      end subroutine Tlatlon

!=======================================================================
!BOP
!
! !IROUTINE: t2ugrid_vector - transfer vector from T-cells to U-cells
!
! !INTERFACE:
!
      subroutine t2ugrid_vector (work)
!
! !DESCRIPTION:
!
! Transfer vector component from T-cell centers to U-cell centers.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_boundary
      use ice_work, only: work1
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks),
     &     intent(inout) :: 
     &     work
!
!EOP
!
      work1(:,:,:) = work(:,:,:)

      call ice_timer_start(timer_bound)
      call update_ghost_cells(work1,            bndy_info,
     &                        field_loc_center, field_type_vector)
      call ice_timer_stop(timer_bound)

      call to_ugrid(work1,work)

      end subroutine t2ugrid_vector

!=======================================================================
!BOP
!
! !IROUTINE: to_ugrid - shift from T-cell to U-cell midpoints
!
! !INTERFACE:
!
      subroutine to_ugrid(work1,work2)
!
! !DESCRIPTION:
!
! Shifts quantities from the T-cell midpoint (work1) to the U-cell
! midpoint (work2)
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
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
     &   work1(nx_block,ny_block,max_blocks)

      real (kind=dbl_kind), intent(out) ::
     &   work2(nx_block,ny_block,max_blocks)

      type (block) ::
     &   this_block           ! block information for current block
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, iblk
     &,  ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      work2(:,:,:) = c0

      do iblk = 1, nblocks
         this_block = get_block(blocks(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            work2(i,j,iblk) = p25 *
     &                        (work1(i,  j,  iblk)*tarea(i,  j,  iblk)
     &                       + work1(i+1,j,  iblk)*tarea(i+1,j,  iblk)
     &                       + work1(i,  j+1,iblk)*tarea(i,  j+1,iblk)
     &                       + work1(i+1,j+1,iblk)*tarea(i+1,j+1,iblk))
     &                       / uarea(i,  j,  iblk)
         enddo
         enddo
      enddo

      end subroutine to_ugrid

!=======================================================================
!BOP
!
! !IROUTINE: u2tgrid_vector - transfer vector from U-cells to T-cells
!
! !INTERFACE:
!
      subroutine u2tgrid_vector (work)
!
! !DESCRIPTION:
!
! Transfer from U-cell centers to T-cell centers. Writes work into
! another array that has ghost cells
! NOTE: Input array is dimensioned only over physical cells.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_boundary
      use ice_work, only: work1
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks),
     &   intent(inout) ::
     &   work
!
!EOP
!
      work1(:,:,:) = work(:,:,:)

      call ice_timer_start(timer_bound)
      call update_ghost_cells(work1,            bndy_info,
     &                        field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      call to_tgrid(work1,work)

      end subroutine u2tgrid_vector

!=======================================================================
!BOP
!
! !IROUTINE: to_tgrid - shifts array from U-cell to T-cell midpoints
!
! !INTERFACE:
!
      subroutine to_tgrid(work1, work2)
!
! !DESCRIPTION:
!
! Shifts quantities from the U-cell midpoint (work1) to the T-cell
! midpoint (work2)
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind) :: work1(nx_block,ny_block,max_blocks)
     &,                       work2(nx_block,ny_block,max_blocks)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, iblk
     &,  ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) ::
     &   this_block           ! block information for current block
      
      do iblk = 1, nblocks
         this_block = get_block(blocks(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            work2(i,j,iblk) = p25 *
     &                       (work1(i,  j  ,iblk) * uarea(i,  j,  iblk)
     &                      + work1(i-1,j  ,iblk) * uarea(i-1,j,  iblk)
     &                      + work1(i,  j-1,iblk) * uarea(i,  j-1,iblk)
     &                      + work1(i-1,j-1,iblk) * uarea(i-1,j-1,iblk))
     &                      / tarea(i,  j,  iblk)
         enddo
         enddo
      enddo

      end subroutine to_tgrid

!=======================================================================

      end module ice_grid

!=======================================================================
