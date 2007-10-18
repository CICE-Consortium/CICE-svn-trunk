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
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!          Tony Craig, NCAR
!
! 2004: Block structure added by William Lipscomb
!       init_grid split into two parts as in POP 2.0
!       Boundary update routines replaced by POP versions
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2007: Neumann boundary condition option added by E. Hunke
!       Option to read from netcdf files (A. Keen, Met Office)
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
      use ice_exit
!
!EOP
!
      implicit none
      save

      character (len=char_len_long) :: &
         grid_format  , & ! file format ('bin'=binary or 'nc'=netcdf)
         grid_file    , & !  input file for POP grid info
         kmt_file     , & !  input file for POP grid info
         grid_type        !  current options are rectangular (default),
                          !  displaced_pole, tripole, panarctic, latlon

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks):: &
         dxt    , & ! width of T-cell through the middle (m)
         dyt    , & ! height of T-cell through the middle (m)
         dxu    , & ! width of U-cell through the middle (m)
         dyu    , & ! height of U-cell through the middle (m)
         HTE    , & ! length of eastern edge of T-cell (m)
         HTN    , & ! length of northern edge of T-cell (m)
         tarea  , & ! area of T-cell (m^2)
         uarea  , & ! area of U-cell (m^2)
         tarear , & ! 1/tarea
         uarear , & ! 1/uarea
         tinyarea,& ! puny*tarea
         tarean , & ! area of NH T-cells
         tareas , & ! area of SH T-cells
         ULON   , & ! longitude of velocity pts (radians)
         ULAT   , & ! latitude of velocity pts (radians)
         TLON   , & ! longitude of temp pts (radians)
         TLAT   , & ! latitude of temp pts (radians)
         ANGLE  , & ! for conversions between POP grid and lat/lon
         ANGLET     ! ANGLE converted to T-cells

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks):: &
         cyp    , & ! 1.5*HTE - 0.5*HTE
         cxp    , & ! 1.5*HTN - 0.5*HTN
         cym    , & ! 0.5*HTE - 1.5*HTE
         cxm    , & ! 0.5*HTN - 1.5*HTN
         dxhy   , & ! 0.5*(HTE - HTE)
         dyhx       ! 0.5*(HTN - HTN)

      ! geometric quantities used for remapping transport
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks):: &
         xav  , & ! mean T-cell value of x
         yav  , & ! mean T-cell value of y
         xxav , & ! mean T-cell value of xx
         xyav , & ! mean T-cell value of xy
         yyav , & ! mean T-cell value of yy
         xxxav, & ! mean T-cell value of xxx
         xxyav, & ! mean T-cell value of xxy
         xyyav, & ! mean T-cell value of xyy
         yyyav    ! mean T-cell value of yyy

      real (kind=dbl_kind), &
         dimension (2,2,nx_block,ny_block,max_blocks) :: &
         mne, & ! matrices used for coordinate transformations in remapping
         mnw, & ! ne = northeast corner, nw = northwest, etc.
         mse, & 
         msw

      ! masks
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks):: &
         hm     , & ! land/boundary mask, thickness (T-cell)
         uvm        ! land/boundary mask, velocity (U-cell)

      logical (kind=log_kind), &
         dimension (nx_block,ny_block,max_blocks) :: &
         tmask  , & ! land/boundary mask, thickness (T-cell)
         umask  , & ! land/boundary mask, velocity (U-cell)
         lmask_n, & ! northern hemisphere mask
         lmask_s    ! southern hemisphere mask

      ! grid dimensions for rectangular grid
      real (kind=dbl_kind), parameter ::  &
         dxrect = 1.6e4_dbl_kind   ,&! uniform HTN (m)
         dyrect = 1.6e4_dbl_kind     ! uniform HTE (m)

!lipscomb - not sure rndex_global is stil needed
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         rndex_global       ! global index for local subdomain (dbl)


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
      integer (kind=int_kind) :: &
         i, j, iblk, &
         fid_grid, &     ! file id for netCDF grid file
         fid_kmt         ! file id for netCDF kmt file


      character (char_len) :: &
         fieldname       ! field name in netCDF file

      !-----------------------------------------------------------------
      ! Get global ULAT and KMT arrays used for block decomposition.
      !-----------------------------------------------------------------

      allocate(work_g1(nx_global,ny_global))
      allocate(work_g2(nx_global,ny_global))

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole'      ) then

         if (trim(grid_format) == 'nc') then

            call ice_open_nc(grid_file,fid_grid)
            call ice_open_nc(kmt_file,fid_kmt)

            fieldname='ulat'
            call ice_read_global_nc(fid_grid,1,fieldname,work_g1,.true.)
            fieldname='kmt'
            call ice_read_global_nc(fid_kmt,1,fieldname,work_g2,.true.)

            if (my_task == master_task) then
               call ice_close_nc(fid_grid)
               call ice_close_nc(fid_kmt)
            endif

         else

            call ice_open(nu_grid,grid_file,64) ! ULAT
            call ice_open(nu_kmt, kmt_file, 32) ! KMT

            call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)  ! ULAT
            call ice_read_global(nu_kmt, 1,work_g2,'ida4',.true.)  ! KMT

            if (my_task == master_task) then
               close (nu_grid)
               close (nu_kmt)
            endif

         endif

      elseif (trim(grid_type) == 'panarctic') then

         call ice_open(nu_grid,grid_file,64) ! ULAT, KMT

         call ice_read_global(nu_grid,1,work_g2,'ida8',.true.)  ! KMT
         call ice_read_global(nu_grid,2,work_g1,'rda8',.true.)  ! ULAT

         if (my_task == master_task) close (nu_grid)

      elseif (trim(grid_type) == 'latlon') then

         work_g1(:,:) = 75._dbl_kind/rad_to_deg  ! arbitrary polar latitude
         work_g2(:,:) = c1

      else   ! rectangular grid

         work_g1(:,:) = 75._dbl_kind/rad_to_deg  ! arbitrary polar latitude
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
      use ice_blocks, only: get_block, block
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
         angle_0, angle_w, angle_s, angle_sw

      logical (kind=log_kind), dimension(nx_block,ny_block,max_blocks):: &
         out_of_range

      character (char_len) :: &
         bc                   ! boundary condition type (Dirichlet, Neumann)

      type (block) :: &
         this_block           ! block information for current block
      
      !-----------------------------------------------------------------
      ! lat, lon, cell widths, angle, land mask
      !-----------------------------------------------------------------

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole'      ) then
         if (trim(grid_format) == 'nc') then
            call popgrid_nc     ! read POP grid lengths from nc file
         else
            call popgrid        ! read POP grid lengths directly
         endif 
      elseif (trim(grid_type) == 'panarctic') then
         call panarctic_grid    ! pan-Arctic grid
#ifdef SEQ_MCT
      elseif (trim(grid_type) == 'latlon') then
         call latlongrid        ! lat lon grid for sequential CCSM (CAM mode)
         return
#endif
      else
         call rectgrid          ! regular rectangular grid
      endif

      !-----------------------------------------------------------------
      ! T-grid cell and U-grid cell quantities
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
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

            uarea(i,j,iblk) = dxu(i,j,iblk)*dyu(i,j,iblk)

            if (tarea(i,j,iblk) > c0) then
               tarear(i,j,iblk) = c1/tarea(i,j,iblk)
            else
               tarear(i,j,iblk) = c0 ! possible on boundaries
            endif
            if (uarea(i,j,iblk) > c0) then
               uarear(i,j,iblk) = c1/uarea(i,j,iblk)
            else
               uarear(i,j,iblk) = c0 ! possible on boundaries
            endif
            tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

            dxhy(i,j,iblk) = p5*(HTE(i,j,iblk) - HTE(i-1,j,iblk))
            dyhx(i,j,iblk) = p5*(HTN(i,j,iblk) - HTN(i,j-1,iblk))
         enddo
         enddo

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
      bc = 'Neumann'
      call update_ghost_cells (dxt,                bndy_info, &
                               field_loc_center,   field_type_scalar, bc)
      call update_ghost_cells (dyt,                bndy_info, &
                               field_loc_center,   field_type_scalar, bc)
      call update_ghost_cells (tarea,              bndy_info, &
                               field_loc_center,   field_type_scalar, bc)
      call update_ghost_cells (dxu,                bndy_info, &
                               field_loc_NEcorner, field_type_scalar, bc)
      call update_ghost_cells (dyu,                bndy_info, &
                               field_loc_NEcorner, field_type_scalar, bc)
      call update_ghost_cells (uarea,              bndy_info, &
                               field_loc_NEcorner, field_type_scalar, bc)
      call update_ghost_cells (uarear,             bndy_info, &
                               field_loc_NEcorner, field_type_scalar, bc)
      call update_ghost_cells (tarear,             bndy_info, &
                               field_loc_center,   field_type_scalar, bc)
      call update_ghost_cells (tinyarea,           bndy_info, &
                               field_loc_center,   field_type_scalar, bc)
      call update_ghost_cells (dxhy,               bndy_info, &
                               field_loc_center,   field_type_vector, bc)
      call update_ghost_cells (dyhx,               bndy_info, &
                               field_loc_center,   field_type_vector, bc)
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
         this_block = get_block(blocks_ice(iblk),iblk)         
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
               if ( abs(angle_w - angle_0) > pi) &
                        angle_w = angle_w  - pi2
               if ( abs(angle_s - angle_0) > pi) &
                        angle_s = angle_s  - pi2
               if ( abs(angle_sw - angle_0) > pi) &
                        angle_sw = angle_sw - pi2
            endif

            ANGLET(i,j,iblk) = angle_0 * p25 + angle_w * p25 &
                             + angle_s * p25 + angle_sw* p25
         enddo
         enddo
      enddo
      
      call ice_timer_start(timer_bound)
      bc = 'Neumann'
      call update_ghost_cells (ANGLET,             bndy_info, &
                               field_loc_center, field_type_angle, bc)
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

      call scatter_global(rndex_global, work_g1,  &
                          master_task,  distrb_info, &
                          field_loc_center, field_type_scalar)

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
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      type (block) :: &
         this_block           ! block information for current block

      call ice_open(nu_grid,grid_file,64)
      call ice_open(nu_kmt,kmt_file,32)

      diag = .true.       ! write diagnostic info

      ! lat, lon, cell dimensions, angles
      call ice_read(nu_grid,1,ULAT, 'rda8',diag, &
                    field_loc_NEcorner, field_type_scalar)
      call ice_read(nu_grid,2,ULON, 'rda8',diag, &
                    field_loc_NEcorner, field_type_scalar)
      call ice_read(nu_grid,3,HTN,  'rda8',diag, &
                    field_loc_Nface, field_type_scalar)
      call ice_read(nu_grid,4,HTE,  'rda8',diag, &
                    field_loc_Eface, field_type_scalar)
      call ice_read(nu_grid,7,ANGLE,'rda8',diag, &
                    field_loc_NEcorner, field_type_scalar)

      ! fix units
      HTN(:,:,:) = HTN(:,:,:) * cm_to_m
      HTE(:,:,:) = HTE(:,:,:) * cm_to_m

      ! topography
      call ice_read(nu_kmt,1,work1,'ida4',diag, &
                    field_loc_center, field_type_scalar)

      if (my_task == master_task) then
         close (nu_grid)
         close (nu_kmt)
      endif

      hm(:,:,:) = c0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
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
!!!             if (ULAT(i,j,iblk) > shlat/rad_to_deg .and. &
!!!                 ULAT(i,j,iblk) < nhlat/rad_to_deg) &
!!!                 hm(i,j,iblk) = c0
         enddo
         enddo
      enddo

      end subroutine popgrid

!=======================================================================
!BOP
!
! !IROUTINE: popgrid_nc - read and set POP tripole
!                      grid and land mask from netCDF file 
!
! !INTERFACE:
!
      subroutine popgrid_nc

#ifdef ncdf
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
! Revised for netcdf input: Ann Keen, Met Office, May 2007
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
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi, &     ! beginning and end of physical domain
	 fid_grid, &		! file id for netCDF grid file
	 fid_kmt		! file id for netCDF kmt file

      logical (kind=log_kind) :: diag

      character (char_len) :: &
         fieldname		! field name in netCDF file

      type (block) :: &
         this_block           ! block information for current block

      call ice_open_nc(grid_file,fid_grid)
      call ice_open_nc(kmt_file,fid_kmt)

      diag = .true.       ! write diagnostic info

      ! lat, lon, cell dimensions, angles
      fieldname='ulat'
      call ice_read_nc(fid_grid,1,fieldname,ULAT,diag, &
                       field_loc_NEcorner, field_type_scalar)
      fieldname='ulon'
      call ice_read_nc(fid_grid,2,fieldname,ULON,diag, &
                       field_loc_NEcorner, field_type_scalar)
      fieldname='htn'
      call ice_read_nc(fid_grid,3,fieldname,HTN,diag, &
                       field_loc_Nface, field_type_scalar)
      fieldname='hte'
      call ice_read_nc(fid_grid,4,fieldname,HTE,diag, &
                       field_loc_Eface, field_type_scalar)
      fieldname='angle'
      call ice_read_nc(fid_grid,7,fieldname,ANGLE,diag, &
                       field_loc_NEcorner, field_type_scalar)

      ! fix units
      HTN(:,:,:) = HTN(:,:,:) * cm_to_m
      HTE(:,:,:) = HTE(:,:,:) * cm_to_m

      ! fix ANGLE: roundoff error due to single precision
      where (ANGLE >  pi) ANGLE =  pi
      where (ANGLE < -pi) ANGLE = -pi

      ! topography
      fieldname='kmt'
      call ice_read_nc(fid_kmt,1,fieldname,work1,diag, &
                       field_loc_center, field_type_scalar)

      if (my_task == master_task) then
         call ice_close_nc(fid_grid)
         call ice_close_nc(fid_kmt)
      endif

      hm(:,:,:) = c0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
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
!!!             if (ULAT(i,j,iblk) > shlat/rad_to_deg .and. &
!!!                 ULAT(i,j,iblk) < nhlat/rad_to_deg) &
!!!                 hm(i,j,iblk) = c0
         enddo
         enddo
      enddo

#endif
      end subroutine popgrid_nc

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

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      type (block) :: &
         this_block           ! block information for current block

      call ice_open(nu_grid,grid_file,64)

      diag = .true.       ! write diagnostic info

      if (my_task == master_task) &
           write (nu_diag,*) '** Reading pan-Arctic grid **'

      ! read topography
      call ice_read(nu_grid,1,work1,'ida8',diag, &
                    field_loc_center, field_type_scalar)

      hm(:,:,:) = c0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
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

      call ice_read(nu_grid,2,ULAT, 'rda8',diag, &
                    field_loc_NEcorner, field_type_scalar)
      call ice_read(nu_grid,3,ULON, 'rda8',diag, &
                    field_loc_NEcorner, field_type_scalar)
      call ice_read(nu_grid,4,HTN,  'rda8',diag, &
                    field_loc_Nface, field_type_scalar)
      call ice_read(nu_grid,5,HTE,  'rda8',diag, &
                    field_loc_Eface, field_type_scalar)
      call ice_read(nu_grid,8,ANGLE,'rda8',diag, &
                    field_loc_NEcorner, field_type_scalar)

      ! fix units
      HTN(:,:,:) = HTN(:,:,:)*cm_to_m
      HTE(:,:,:) = HTE(:,:,:)*cm_to_m

      if (my_task == master_task) close (nu_grid)

      end subroutine panarctic_grid
#ifdef SEQ_MCT
!=======================================================================
!BOP
!
! !IROUTINE: latlongrid- lat and lon grid for coupling to standalone CAM
!
! !INTERFACE:
!
      subroutine latlongrid
!
! !DESCRIPTION:
!
! Read in kmt file that matches CAM lat-lon grid and has single column 
! functionality
!
! !REVISION HISTORY:
!
! author: Mariana Vertenstein
! 2007: Elizabeth Hunke upgraded to netcdf90 and cice ncdf calls
!
! !USES:
!
!     use ice_boundary
      use ice_domain_size
      use ice_gather_scatter
      use ice_scam, only : scmlat, scmlon, single_column
      use netcdf
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk    
      
      integer (kind=int_kind) :: &
         ni, nj, ncid, dimid, varid, ier

      character (len=char_len) :: &
         subname='latlongrid' ! subroutine name

      type (block) :: &
         this_block           ! block information for current block
      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
           closelat, &        ! Single-column latitude value
           closelon, &        ! Single-column longitude value
           closelatidx, &     ! Single-column latitude index to retrieve
           closelonidx        ! Single-column longitude index to retrieve

      integer (kind=int_kind) :: &
           start(2), &        ! Start index to read in
           count(2)           ! Number of points to read in

      real (kind=dbl_kind), allocatable :: &
           lats(:),lons(:),pos_lons(:), glob_grid(:,:)  ! temporaries 

      real (kind=dbl_kind) :: &
         pos_scmlon           ! temporary

      !-----------------------------------------------------------------
      ! - kmt file is actually clm fractional land file
      ! - Determine consistency of dimensions
      ! - Read in lon/lat centers in degrees from kmt file
      ! - Read in ocean from "kmt" file (1 for ocean, 0 for land)
      !-----------------------------------------------------------------

      ! Determine dimension of domain file and check for consistency

         call ice_open_nc(kmt_file, ncid))

         call check_ret(nf90_inq_dimid (ncid, 'ni', dimid), subname)
         call check_ret(nf90_inq_dimension(ncid, dimid, len=ni), subname)
         call check_ret(nf90_inq_dimid (ncid, 'nj', dimid), subname)
         call check_ret(nf90_inq_dimension(ncid, dimid, len=nj), subname)

         if (single_column) then
            if ((nx_global /= 1).or. (ny_global /= 1)) then
               write(nu_diag,*) 'Because you have selected the column model flag'
               write(nu_diag,*) 'Please set nx_global=ny_global=1 in file'
               write(nu_diag,*) 'ice_domain_size.F and recompile'
               call abort_ice ('latlongrid: check nx_global, ny_global')
            endif
         else
            if (nx_global /= ni .and. ny_global /= nj) then
               call abort_ice ('latlongrid: ni,ny not equal to nx_global,ny_global')
            end if
         end if
         
      ! Determine start/count to read in for either single column or global lat-lon grid

      if (single_column) then
         allocate(lats(nj))
         allocate(lons(ni))
         allocate(pos_lons(ni))
         allocate(glob_grid(ni,nj))

         call ice_read_global_nc(ncid, 1, 'xc', glob_grid, dbug)
         do i = 1,ni
            lons(i) = glob_grid(i,1)
         end do
         call ice_read_global_nc(ncid, 1, 'yc', glob_grid, dbug)
         do j = 1,nj
            lats(j) = glob_grid(1,j) 
         end do
         
         ! convert lons array and scmlon to 0,360 and find index of value closest to 0
         ! and obtain single-column longitude/latitude indices to retrieve
         
         pos_lons(:)= mod(lons(:) + 360._r8,360._r8)
         pos_scmlon = mod(scmlon  + 360._r8,360._r8)
         start(1) = (MINLOC(abs(pos_lons-pos_scmlon),dim=1))
         start(2) = (MINLOC(abs(lats    -scmlat    ),dim=1))

         deallocate(lats)
         deallocate(lons)
         deallocate(pos_lons)
         deallocate(glob_grid)
      else
          start(1)=1
          start(2)=1
      endif
      count(1)=nx_global
      count(2)=ny_global
      
      ! Read in domain file for either single column or for global lat-lon grid

         call ice_read_nc(ncid, 1, 'xc', TLON, dbug)
         do j = 1, ny_global
         do i = 1, nx_global
            ! Convert from degrees to radians
            TLON(i,j) = pi*TLON(i,j)/180._dbl_kind 
         end do
         end do 
         call ice_read_nc(ncid, 1, 'yc', TLAT, dbug)
         do j = 1, ny_global
         do i = 1, nx_global
            ! Convert from degrees to radians
            TLAT(i,j) = pi*TLAT(i,j)/180._dbl_kind 
         end do
         end do 
         ! Note that area read in from domain file is in km^2 - must first convert to m^2 and 
         ! then convert to radians^2
         call ice_read_nc(ncid, 1, 'area', tarea, dbug)
         do j = 1, ny_global
         do i = 1, nx_global
            ! Convert from km^2 to m^2
            tarea(i,j) = tarea(i,j) * 1.e6_dbl_kind
         end do
         end do 
         call ice_read_nc(ncid, 1, 'mask', hm, dbug)
         call ice_close(ncid)

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      ! The U grid (velocity) is not used when run with sequential CAM
      ! because we only use thermodynamic sea ice.  However, ULAT is used
      ! in the default initialization of CICE so we calculate it here as 
      ! a "dummy" so that CICE will initialize with ice.  If a no ice
      ! initialization is OK (or desired) this can be commented out and
      ! ULAT will remain 0 as specified above.  ULAT is located at the
      ! NE corner of the grid cell, TLAT at the center, so here ULAT is
      ! hacked by adding half the latitudinal spacing (in radians) to
      ! TLAT.
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            uarea(i,j,iblk)     = p25*  &
                                 (tarea(i,j,  iblk) + tarea(i+1,j,  iblk) &
                                + tarea(i,j+1,iblk) + tarea(i+1,j+1,iblk))
            tarear(i,j,iblk)   = c1/tarea(i,j,iblk)
            uarear(i,j,iblk)   = c1/uarea(i,j,iblk)
            tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

            ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/ny_global)  
            ULON  (i,j,iblk) = c0
            ANGLE (i,j,iblk) = c0                             

            ANGLET(i,j,iblk) = c0                             
            HTN   (i,j,iblk) = 1.e36_dbl_kind
            HTE   (i,j,iblk) = 1.e36_dbl_kind
            dxt   (i,j,iblk) = 1.e36_dbl_kind
            dyt   (i,j,iblk) = 1.e36_dbl_kind
            dxu   (i,j,iblk) = 1.e36_dbl_kind
            dyu   (i,j,iblk) = 1.e36_dbl_kind
            dxhy  (i,j,iblk) = 1.e36_dbl_kind
            dyhx  (i,j,iblk) = 1.e36_dbl_kind
            cyp   (i,j,iblk) = 1.e36_dbl_kind
            cxp   (i,j,iblk) = 1.e36_dbl_kind
            cym   (i,j,iblk) = 1.e36_dbl_kind
            cxm   (i,j,iblk) = 1.e36_dbl_kind
         enddo
         enddo
      enddo

      call makemask

      end subroutine latlongrid

!=======================================================================
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
      subroutine check_ret(ret, calling)
!
! !DESCRIPTION:
!     Check return status from netcdf call
!
! !USES:
        use netcdf
! !ARGUMENTS:
        implicit none
        integer, intent(in) :: ret
        character(len=*) :: calling
!
! !REVISION HISTORY:
! author: Mariana Vertenstein
! 2007: Elizabeth Hunke upgraded to netcdf90
!
!EOP
!
      if (ret /= NF90_NOERR) then
         write(nu_diag,*)'netcdf error from ',trim(calling)
         write(nu_diag,*)'netcdf strerror = ',trim(NF90_STRERROR(ret))
         call abort_ice('ice ice_grid: netcdf check_ret error')
      end if
        
      end subroutine check_ret

#endif      
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
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            HTN  (i,j,iblk) = dxrect
            HTE  (i,j,iblk) = dyrect
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

      call scatter_global(hm, work_g1, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

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
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      character (char_len) :: &
         bc                   ! boundary condition type (Dirichlet, Neumann)

      type (block) :: &
         this_block           ! block information for current block

      call ice_timer_start(timer_bound)
      bc = 'Dirichlet'
      call update_ghost_cells (hm,               bndy_info, &
                               field_loc_center, field_type_scalar, bc)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! construct T-cell and U-cell masks
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            uvm(i,j,iblk) = min (hm(i,j,  iblk), hm(i+1,j,  iblk), &
                                 hm(i,j+1,iblk), hm(i+1,j+1,iblk))
         enddo
         enddo
      enddo

      call ice_timer_start(timer_bound)
      call update_ghost_cells (uvm,                bndy_info, &
                               field_loc_NEcorner, field_type_scalar, bc)
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
            if (lmask_n(i,j,iblk)) tarean(i,j,iblk) = tarea(i,j,iblk) &
                                                    * hm(i,j,iblk)

            ! S hemisphere area mask (m^2)
            if (lmask_s(i,j,iblk)) tareas(i,j,iblk) = tarea(i,j,iblk) &
                                                    * hm(i,j,iblk)

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
      integer (kind=int_kind) :: &
           i, j, iblk       , & ! horizontal indices
           ig, jg           , & ! global horizontal indices
           im1              , & ! ig - 1
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
           z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da

      character (char_len) :: &
         bc                   ! boundary condition type (Dirichlet, Neumann)

      type (block) :: &
           this_block           ! block information for current block

      TLAT(:,:,:) = c0
      TLON(:,:,:) = c0

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
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
            if (tx /= c0 .or. ty /= c0) TLON(i,j,iblk) = atan2(ty,tx)

            ! TLAT in radians North
            TLAT(i,j,iblk) = asin(tz)
            
         enddo                  ! i
         enddo                  ! j         
      enddo                     ! iblk

      call ice_timer_start(timer_bound)
      bc = 'Neumann'
      call update_ghost_cells (TLON,             bndy_info, &
                               field_loc_center, field_type_scalar, bc)
      call update_ghost_cells (TLAT,             bndy_info, &
                               field_loc_center, field_type_scalar, bc)
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
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
           intent(inout) :: & 
           work
!
!EOP
!
      work1(:,:,:) = work(:,:,:)

      call ice_timer_start(timer_bound)
      call update_ghost_cells(work1,            bndy_info, &
                              field_loc_center, field_type_vector)
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
      real (kind=dbl_kind), intent(in) :: &
         work1(nx_block,ny_block,max_blocks)

      real (kind=dbl_kind), intent(out) :: &
         work2(nx_block,ny_block,max_blocks)

      type (block) :: &
         this_block           ! block information for current block
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      work2(:,:,:) = c0

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            work2(i,j,iblk) = p25 * &
                              (work1(i,  j,  iblk)*tarea(i,  j,  iblk)  &
                             + work1(i+1,j,  iblk)*tarea(i+1,j,  iblk)  &
                             + work1(i,  j+1,iblk)*tarea(i,  j+1,iblk)  &
                             + work1(i+1,j+1,iblk)*tarea(i+1,j+1,iblk)) &
                             / uarea(i,  j,  iblk)
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
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), &
         intent(inout) :: &
         work
!
!EOP
!
      work1(:,:,:) = work(:,:,:)

      call ice_timer_start(timer_bound)
      call update_ghost_cells(work1,            bndy_info, &
                              field_loc_center, field_type_scalar)
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
      real (kind=dbl_kind) :: work1(nx_block,ny_block,max_blocks), &
                              work2(nx_block,ny_block,max_blocks)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block
      
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            work2(i,j,iblk) = p25 *  &
                             (work1(i,  j  ,iblk) * uarea(i,  j,  iblk)  &
                            + work1(i-1,j  ,iblk) * uarea(i-1,j,  iblk)  &
                            + work1(i,  j-1,iblk) * uarea(i,  j-1,iblk)  & 
                            + work1(i-1,j-1,iblk) * uarea(i-1,j-1,iblk)) &
                            / tarea(i,  j,  iblk)
         enddo
         enddo
      enddo

      end subroutine to_tgrid

!=======================================================================
!BOP
!
! !IROUTINE: sinlat - calculates sin of latitudes
!
! !INTERFACE:
!
      subroutine sinlat(a, k)
!
! !DESCRIPTION:
!
! Calculates the sin of latitudes on the grid based on ny_global and
! using code from CAM (/control/gauaw_mod.F90)  In CAM, gauaw_mod.F90 is
! only used to calculate the sin of latitudes (and thence latitude) if
! the dynamical core is set to eul.  If using one of the other dynamical
! cores and coupling to stand alone CAM the latitudes calculated in this
! way may not match the grid from CAM.
!
! !REVISION HISTORY:
!
! author: Jacob Sewall
!
! !USES:
      use ice_exit

!
! !INPUT/OUTPUT PARAMETERS:
!

      real (kind=dbl_kind), dimension(ny_global), intent(out) :: & 
           a            ! sin of latitudes

      integer (kind=int_kind), intent(in) :: &
           k            ! number of latitudes (ny_global)
!
!EOP
!
      real (kind=dbl_kind), dimension(k) :: &
           sinlats      ! sine of latitudes

      real (kind=dbl_kind) :: &
           eps      , & ! convergence criterion
           c        , & ! constant combination
           fk       , & ! real k
           xz       , & ! abscissa estimate
           pkm1     , & ! |
           pkm2     , & ! |-polynomials
           pkmrk    , & ! |
           pk       , & ! |
           sp       , & ! current iteration latitude increment
           avsp     , & ! |sp|
           fn       , & ! real n
           avsp_prev, &
           testdiff

      real (kind=dbl_kind), parameter :: &
           eps27 = 1.e-27_dbl_kind

      integer (kind=int_kind) :: &
           n, l     , &
           iter     , & ! iteration counter
           is       , & ! latitude index
           kk           ! k/2 (number of latitudes in a hemisphere)

!
!----------------------------------------------------------------------
!

      c = (c1-(c2/pi)**2)*p25
      fk = k
      kk = k/2
      call bsslzr(sinlats,kk)
      do is=1,kk
        xz = cos(sinlats(is)/sqrt((fk+p5)**2+c))
!
! This is the first approximation to xz
!
        iter = 0
        avsp = c100     !initialize avsp to a very large number
10      continue
        avsp_prev = avsp
        pkm2 = c1
        pkm1 = xz
        iter = iter + 1
        if (iter > 100) then  ! Error exit
           call abort_ice('SINLAT: no convergence in 100 iterations')
        end if
!
! Computation of the legendre polynomial
!
        do n=2,k
          fn = n
          pk = ((c2*fn-1._r8)*xz*pkm1-(fn-c1)*pkm2)/fn
          pkm2 = pkm1
          pkm1 = pk
        enddo
        pkm1 = pkm2
        pkmrk = (fk*(pkm1-xz*pk))/(c1-xz**2)
        sp = pk/pkmrk
        xz = xz - sp
        avsp = abs(sp)
        testdiff = avsp_prev - avsp
        if (testdiff > eps27) go to 10
        sinlats(is) = xz
      end do
!
! Complete the sets of abscissas and weights, using the symmetry.
! Also note truncation from real(r8) to real*8
!
      do n=1,kk
        l = k + 1 - n
        a(n) = sinlats(n)
        a(l) = -sinlats(n)
      end do
 
      end subroutine sinlat

!=======================================================================
!BOP
!
! !IROUTINE: bsslzr - Return n zeros (if n<50) of the Bessel function
!
! !INTERFACE:
!
      subroutine bsslzr(bes, n)
!
! !DESCRIPTION:
!
!
! Return n zeros (or if n>50, approximate zeros), of the Bessel function
! j0,in the array bes. The first 50 zeros will be given exactly, and the
! remaining zeros are computed by extrapolation,and therefore not exact.
!
! Modified 1/23/97 by Jim Rosinski to use real*16 arithmetic
! placed in ice_grid.F 4/26/2006 by Jacob Sewall

! !REVISION HISTORY:
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, D. Williamson, August 1992
! Reviewed:          J. Hack, D. Williamson, April 1996
!
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!

      integer (kind=int_kind), intent(in) :: &
           n          ! number of latitudes in hemisphere (ny_global/2)

      real (kind=dbl_kind), dimension(n), intent(inout) :: & 
           bes        ! sin of latitudes
!
!EOP
!

!----------------------------------------
! Local Variables
!----------------------------------------

      integer (kind=int_kind) :: &
         nn, j
     
      real (kind=dbl_kind), dimension(50) :: bz
      
      save bz
!
!-----------------------------------------
! Local Workspace
!-----------------------------------------

      data bz/ 2.4048255577_r8, 5.5200781103_r8, 8.6537279129_r8 ,&
         11.7915344391_r8, 14.9309177086_r8, 18.0710639679_r8    ,&
         21.2116366299_r8, 24.3524715308_r8, 27.4934791320_r8    ,&
         30.6346064684_r8, 33.7758202136_r8, 36.9170983537_r8    ,&
         40.0584257646_r8, 43.1997917132_r8, 46.3411883717_r8    ,&
         49.4826098974_r8, 52.6240518411_r8, 55.7655107550_r8    ,&
         58.9069839261_r8, 62.0484691902_r8, 65.1899648002_r8    ,&
         68.3314693299_r8, 71.4729816036_r8, 74.6145006437_r8    ,&
         77.7560256304_r8, 80.8975558711_r8, 84.0390907769_r8    ,&
         87.1806298436_r8, 90.3221726372_r8, 93.4637187819_r8    ,&
         96.6052679510_r8, 99.7468198587_r8, 102.8883742542_r8   ,&
         106.0299309165_r8, 109.1714896498_r8, 112.3130502805_r8 ,&
         115.4546126537_r8, 118.5961766309_r8, 121.7377420880_r8 ,&
         124.8793089132_r8, 128.0208770059_r8, 131.1624462752_r8 ,&
         134.3040166383_r8, 137.4455880203_r8, 140.5871603528_r8 ,&
         143.7287335737_r8, 146.8703076258_r8, 150.0118824570_r8 ,&
         153.1534580192_r8, 156.2950342685_r8/  
!
      nn = n 
      if (n > 50) then 
         bes(50) = bz(50) 
         do j = 51, n 
            bes(j) = bes(j-1) + pi 
         end do 
         nn = 49 
      endif 
      bes(:nn) = bz(:nn) 

      end subroutine bsslzr 

!=======================================================================

      end module ice_grid

!=======================================================================
