c $Id: $
c=======================================================================
!---! spatial grids, masks, and boundary conditions
!---!
!---! authors Elizabeth C. Hunke, LANL
!---!         Tony Craig, NCAR
c=======================================================================

      module ice_grid

      use ice_kinds_mod
      use ice_constants
      use ice_domain
      use ice_fileunits

      implicit none

      character (len=char_len) :: 
     &   grid_file  !  input file for POP grid info
     &,  kmt_file   !  input file for POP grid info
     &,  grid_type  !  rectangular (default) or displaced_pole 

      real (kind=dbl_kind), dimension (imt_global,jmt_global) ::
     &   TLAT_G     ! latitude  of cell center
     &,  TLON_G     ! longitude of cell center

      real (kind=dbl_kind), dimension (imt_local,jmt_local) ::
     &   dxt        ! width of T-cell through the middle (m)
     &,  dyt        ! height of T-cell through the middle (m)
     &,  dxu        ! width of U-cell through the middle (m)
     &,  dyu        ! height of U-cell through the middle (m)
     &,  HTE        ! length of eastern edge of T-cell (m)
     &,  HTN        ! length of northern edge of T-cell (m)
     &,  HTS        ! length of southern edge of T-cell
     &,  HTW        ! length of western edge of T-cell
     &,  tarea      ! area of T-cell (m^2)
     &,  uarea      ! area of U-cell (m^2)
     &,  ULON       ! longitude of velocity pts (radians)
     &,  ULAT       ! latitude of velocity pts (radians)
     &,  dxt2r      ! 1. / (dxt(i+1,j) + dxt(i,j))
     &,  dyt2r      ! 1. / (dyt(i,j+1) + dyt(i,j))
     &,  cx, cy, cxx, cxy, cyy  ! mean T-cell values of x, y, xx, xy, yy  
     &,  cxxx, cxxy, cxyy, cyyy ! mean T-cell values of xxx, xxy, xyy, yyy 

      real (kind=dbl_kind), dimension (2,2,imt_local,jmt_local) ::
     &   mne, mnw   ! matrices used for coordinate transformations
     &,  msw, mse   ! ne = northeast corner, nw = northwest corner, etc.

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   cyp        ! 1.5*HTE - 0.5*HTE
     &,  cxp        ! 1.5*HTN - 0.5*HTN
     &,  cym        ! 0.5*HTE - 1.5*HTE
     &,  cxm        ! 0.5*HTN - 1.5*HTN
     &,  dxt2       ! 0.5*dxt
     &,  dyt2       ! 0.5*dyt
     &,  dxt4       ! 0.25*dxt
     &,  dyt4       ! 0.25*dyt
     &,  dxhy       ! 0.5*(HTE - HTE)
     &,  dyhx       ! 0.5*(HTN - HTN)
     &,  tarear     ! 1/tarea
     &,  uarear     ! 1/uarea
     &,  tinyarea   ! puny*tarea
     &,  ANGLE      ! for conversions between POP grid and lat/lon
     &,  ANGLET     ! ANGLE converted to T-cells
     &,  tarean     ! area of NH T-cells
     &,  tareas     ! area of SH T-cells

      ! Masks
      real (kind=dbl_kind), dimension (imt_local,jmt_local) ::
     &   hm         ! land/boundary mask, thickness (T-cell)
     &,  uvm        ! land/boundary mask, velocity (U-cell)
     &,  mask_n     ! northern hemisphere
     &,  mask_s     ! southern hemisphere

      logical (kind=log_kind) ::
     &   tmask(imt_local,jmt_local) ! land/boundary mask, thickness (T-cell)
     &,  umask(imt_local,jmt_local) ! land/boundary mask, velocity (U-cell)
     &,  icetmask(ilo:ihi,jlo:jhi)  ! ice extent mask (T-cell)
     &,  iceumask(ilo:ihi,jlo:jhi)  ! ice extent mask (U-cell)

c=======================================================================

      contains

c=======================================================================

      subroutine init_grid

!---!-------------------------------------------------------------------
!---! Horizontal grid initialization
!---!     HT{N,E} = cell widths on {N,E} sides of T cell
!---!     U{LAT,LONG} = true {latitude,longitude} of U points
!---!     D{X,Y}{T,U} = {x,y} spacing centered at {T,U} points
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) ::
     &     tmpg(imt_global,jmt_global)    ! temporary array
     &,    tmpg4(4,imt_global,jmt_global) ! temporary array
     &,    work1(imt_local,jmt_local)     ! temporary array
     &,    rangle

      if (grid_type .eq. 'displaced_pole') then
        call popgrid              ! read POP grid lengths directly
      elseif (grid_type .eq. 'rectangular') then
        call rectgrid             ! regular rectangular grid
      endif

      call bound(HTN)
      call bound(HTE)
      call bound(ULAT)
      call bound(ULON)

      !-----------------------------------------------------------------
      ! construct T-grid cell and U-grid cell widths
      !-----------------------------------------------------------------

      do j=jlo,jhi
       do i=ilo,ihi
        dxt(i,j) = p5*(HTN(i,j) + HTN(i,j-1))
        dyt(i,j) = p5*(HTE(i,j) + HTE(i-1,j))

        tarea(i,j) = dxt(i,j)*dyt(i,j)

        dxu(i,j) = p5*(HTN(i,j) + HTN(i+1,j))
        dyu(i,j) = p5*(HTE(i,j) + HTE(i,j+1))
       enddo
      enddo

      call bound(dxt)
      call bound(dyt)
      call bound(dxu)
      call bound(dyu)
      call bound(tarea)

      do j=jlo,jhi
       do i=ilo,ihi
        uarea(i,j) = p25*(tarea(i,j) + tarea(i+1,j)
     &         + tarea(i,j+1) + tarea(i+1,j+1))        ! m^2
       enddo
      enddo
      call bound(uarea)

      ! grid length combinations
      do j=jlo,jhi
       do i=ilo,ihi
        dxt2(i,j) = 0.5*dxt(i,j)
        dyt2(i,j) = 0.5*dyt(i,j)
        dxt4(i,j) = 0.25*dxt(i,j)
        dyt4(i,j) = 0.25*dyt(i,j)
        tarear(i,j) = 1./tarea(i,j)
        uarear(i,j) = 1./uarea(i,j)
        tinyarea(i,j) = puny*tarea(i,j)

        cyp(i,j) = (1.5*HTE(i,j) - 0.5*HTE(i-1,j))
        cxp(i,j) = (1.5*HTN(i,j) - 0.5*HTN(i,j-1))
        cym(i,j) = (0.5*HTE(i,j) - 1.5*HTE(i-1,j))
        cxm(i,j) = (0.5*HTN(i,j) - 1.5*HTN(i,j-1))

        dxhy(i,j) = 0.5*(HTE(i,j) - HTE(i-1,j))
        dyhx(i,j) = 0.5*(HTN(i,j) - HTN(i,j-1))
       enddo
      enddo

      call u2tgrid(ANGLET)   ! ANGLE on the T grid
      ! correct anglet calculation for branch point problem
      work1(ilo:ihi,jlo:jhi) = ANGLE(ilo:ihi,jlo:jhi)
      call bound(work1)
      do j=jlo,jhi
      do i=ilo,ihi
        if (abs(work1(i-1,j)-work1(i,j)) .gt. pih)
     &     work1(i-1,j)=work1(i-1,j)
     &     - 2*pi*abs(work1(i-1,j)-work1(i,j))/(work1(i-1,j)-work1(i,j))
        if (abs(work1(i,j-1)-work1(i,j)) .gt. pih)
     &     work1(i,j-1)=work1(i,j-1)
     &     - 2*pi*abs(work1(i,j-1)-work1(i,j))/(work1(i,j-1)-work1(i,j))
        if (abs(work1(i-1,j-1)-work1(i,j)) .gt. pih)
     &     work1(i-1,j-1)=work1(i-1,j-1)
     &     - 2*pi*abs(work1(i-1,j-1)-work1(i,j))
     &              /(work1(i-1,j-1)-work1(i,j))
        rangle = p25*(work1(i,j)  *uarea(i,j)
     &              + work1(i-1,j)*uarea(i-1,j)
     &              + work1(i,j-1)*uarea(i,j-1)
     &              + work1(i-1,j-1)*uarea(i-1,j-1))/tarea(i,j)
        if (rangle.ne.ANGLET(i,j)) then
c          write(nu_diag,*) ' changing ANGLET ',
c     &                     my_task,i,j,ANGLET(i,j),rangle
          ANGLET(i,j)=rangle
        endif
      enddo
      enddo

      call Tlatlon           ! lat, lon on the T grid
      call makemask          ! velocity mask, hemisphere masks

      end subroutine init_grid

c=======================================================================

      subroutine popgrid

!---!-------------------------------------------------------------------
!---! POP displaced pole grid and land mask
!---!   grid
!---!      rec no.         field         units
!---!      -------         -----         -----
!---!         1            ULAT         radians
!---!         2            ULON         radians
!---!         3             HTN           cm
!---!         4             HTE           cm
!---!         5             HUS           cm
!---!         6             HUW           cm
!---!         7            ANGLE        radians
!---!   land mask
!---!      rec no.         field         units
!---!      -------         -----         -----
!---!         1             KMT
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_read_write

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: work(ilo:ihi,jlo:jhi)
      logical (kind=log_kind) :: scatter, diag

      call ice_open(nu_grid,grid_file,64)
      call ice_open(nu_kmt,kmt_file,32)

      scatter = .true.    ! scatter data to all processors
      diag = .true.       ! write diagnostic info

      call ice_read(nu_grid,1,work,'rda8',scatter,diag)
      ULAT(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(nu_grid,2,work,'rda8',scatter,diag)
      ULON(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(nu_grid,3,work,'rda8',scatter,diag)
      HTN(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)*cm_to_m
      call ice_read(nu_grid,4,work,'rda8',scatter,diag)
      HTE(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)*cm_to_m
      call ice_read(nu_grid,7,work,'rda8',scatter,diag)
      ANGLE(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(nu_kmt,1,work,'ida4',scatter,diag)

      if (my_task.eq.master_task) then
         close (nu_grid) 
         close (nu_kmt) 
      endif

      do j=jlo,jhi
       do i=ilo,ihi
         hm(i,j) = work(i,j)
         if (hm(i,j).ge.c1) hm(i,j) = c1 
         ANGLET(i,j) = ANGLE(i,j)   
       enddo
      enddo

      end subroutine popgrid

c=======================================================================

      subroutine rectgrid

!---!-------------------------------------------------------------------
!---! Regular rectangular grid and mask
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_mpi_internal

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: 
     &   hmg(imt_global,jmt_global)
     &,  work(ilo:ihi,jlo:jhi)

c.. calculate various geometric 2d arrays
      do j=jlo,jhi
       do i=ilo,ihi
c         HTN  (i,j) = 3.1e4_dbl_kind  ! constant longitude spacing =  
                                       ! POP <4/3> min, m
c         HTE  (i,j) = 3.1e4_dbl_kind  ! constant latitude  spacing = 
                                       ! POP <4/3> min, m
         HTN  (i,j) = 1.6e4_dbl_kind   ! constant longitude spacing = 
                                       ! POP <2/3> min, m
         HTE  (i,j) = 1.6e4_dbl_kind   ! constant latitude  spacing = 
                                       ! POP <2/3> min, m
         ULAT (i,j) = c0               ! remember to set Coriolis !
         ULON (i,j) = c0
         ANGLE(i,j) = c0               ! "square with the world"
         ANGLET(i,j) = c0 
       enddo
      enddo

      !-----------------------------------------------------------------
      ! construct T-cell land mask
      !-----------------------------------------------------------------
      do j=1,jmt_global         ! initialize hm
       do i=1,imt_global        ! open
        hmg(i,j) = c0
       enddo
      enddo

c      do j=1,jmt_global        ! open
c       do i=1,imt_global       ! open
      do j=3,jmt_global-2       ! closed: NOTE jmt_global > 5
       do i=3,imt_global-2      ! closed: NOTE imt_global > 5
        hmg(i,j) = c1
       enddo
      enddo

      call global_scatter(hmg,work)

      do j=jlo,jhi
       do i=ilo,ihi
        hm(i,j) = work(i,j)
       enddo
      enddo

      end subroutine rectgrid

c=======================================================================

      subroutine makemask

!---!-------------------------------------------------------------------
!---! Sets the boundary values for the T cell land mask (hm) and
!---! makes the logical land masks for T and U cells (tmask, umask).
!---! Also creates hemisphere masks (mask_n northern, mask_s southern)
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j

      call bound(hm)  !!! use real arrays to get boundary conditions

      !-----------------------------------------------------------------
      ! construct T-cell and U-cell masks
      !-----------------------------------------------------------------

      do j=jlo,jhi
       do i=ilo,ihi
        uvm(i,j) = min(hm(i,j),hm(i+1,j),hm(i,j+1),hm(i+1,j+1))
       enddo
      enddo
      call bound(uvm)  !!! use real arrays to get boundary conditions

      do j=1,jmt_local
       do i=1,imt_local
        tmask(i,j) = .false.
        umask(i,j) = .false.
        if ( hm(i,j).gt.p5) tmask(i,j) = .true. 
        if (uvm(i,j).gt.p5) umask(i,j) = .true. 
       enddo
      enddo

      !-----------------------------------------------------------------
      ! create hemisphere masks
      !-----------------------------------------------------------------

      do j=1,jmt_local
       do i=1,imt_local
        mask_n(i,j) = c0
        mask_s(i,j) = c0
       enddo
      enddo
      do j=jlo,jhi
       do i=ilo,ihi
        if (ULAT(i,j).ge.-puny) mask_n(i,j) = c1  ! northern hemisphere
        if (ULAT(i,j).lt.-puny) mask_s(i,j) = c1  ! southern hemisphere

        tarean(i,j) = tarea(i,j)*mask_n(i,j)  ! N hemisphere area mask (m^2)
        tareas(i,j) = tarea(i,j)*mask_s(i,j)  ! S hemisphere area mask (m^2)
       enddo
      enddo

      end subroutine makemask

c=======================================================================

      subroutine Tlatlon

!---!-------------------------------------------------------------------
!---! initializes latitude and longitude on T grid
!---!
!---! author Elizabeth C. Hunke, LANL
!---! code originally based on POP grid generation routine
!---!-------------------------------------------------------------------
      use ice_model_size
      use ice_mpi_internal

      integer (kind=int_kind) ::
     &   i, j                           ! generic indices

      integer (kind=int_kind) :: im1
      real (kind=dbl_kind) ::
     &   z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da
     &,  ULATG(imt_global,jmt_global)   ! latitude  of NE cell corner
     &,  ULONG(imt_global,jmt_global)   ! longitude of NE cell corner

      call global_gather(ULONG,ULON(ilo:ihi,jlo:jhi))
      call global_gather(ULATG,ULAT(ilo:ihi,jlo:jhi))

      if (my_task == master_task) then

      do j=2,jmt_global
        do i=1,imt_global

            if (i.eq.1) then
               im1=imt_global
            else
               im1=i-1
            endif

            z1 = cos(ULATG(im1,j-1))
            x1 = cos(ULONG(im1,j-1))*z1
            y1 = sin(ULONG(im1,j-1))*z1
            z1 = sin(ULATG(im1,j-1))

            z2 = cos(ULATG(i,j-1))
            x2 = cos(ULONG(i,j-1))*z2
            y2 = sin(ULONG(i,j-1))*z2
            z2 = sin(ULATG(i,j-1))

            z3 = cos(ULATG(im1,j))
            x3 = cos(ULONG(im1,j))*z3
            y3 = sin(ULONG(im1,j))*z3
            z3 = sin(ULATG(im1,j))

            z4 = cos(ULATG(i,j))
            x4 = cos(ULONG(i,j))*z4
            y4 = sin(ULONG(i,j))*z4
            z4 = sin(ULATG(i,j))

            tx = (x1+x2+x3+x4)/c4
            ty = (y1+y2+y3+y4)/c4
            tz = (z1+z2+z3+z4)/c4
            da = sqrt(tx**2+ty**2+tz**2)

            tz = tz/da
            TLON_G(i,j) = c0
            if (tx.ne.c0.or.ty.ne.c0) TLON_G(i,j) = atan2(ty,tx)
            TLAT_G(i,j) = asin(tz)

        end do
      end do

      ! j=1: linear approximation
      do i=1,imt_global
         TLON_G(i,1) = TLON_G(i,2)
         TLAT_G(i,1) = c2*TLAT_G(i,2) - TLAT_G(i,3)
      end do

      write (nu_diag,*) 'min/max ULONG',minval(ULONG),maxval(ULONG)
      write (nu_diag,*) 'min/max TLON_G',minval(TLON_G),maxval(TLON_G)
      write (nu_diag,*) 'min/max ULATG',minval(ULATG),maxval(ULATG)
      write (nu_diag,*) 'min/max TLAT_G',minval(TLAT_G),maxval(TLAT_G)
      endif ! master_task

      end subroutine Tlatlon

c=======================================================================

      subroutine init_remap

!---!-------------------------------------------------------------------
!---! Grid quantities used by the remapping transport scheme
!---!     
!---! author William H. Lipscomb, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) ::
     &     i, j, m                         ! counting indices

      real (kind=dbl_kind), dimension (4) ::  
     &     ar,                             ! triangle area
     &     xv1, xv2, xv3, yv1, yv2, yv3    ! vertices

      real (kind=dbl_kind) ::
     &     x0, x1, x2, x3, y0, y1, y2, y3  ! interior points
     &,    theta1, theta2, theta3, theta4  ! angles

      !-----------------------------------------------------------------
      ! Compute lengths of south and west cell faces.
      ! Also compute cell widths used in limited_gradient subroutine.
      !-----------------------------------------------------------------

      do j=jlo,jhi
         do i=ilo,ihi
            HTS(i,j) = HTN(i,j-1)
            HTW(i,j) = HTE(i-1,j)
            dxt2r(i,j) = c1 / (dxt(i+1,j) + dxt(i,j))
            dyt2r(i,j) = c1 / (dyt(i,j+1) + dyt(i,j))
         enddo
      enddo
      call bound(HTS)
      call bound(HTW)  
      call bound(dxt2r)
      call bound(dyt2r)

      !-----------------------------------------------------------------
      ! Construct mean values of geometic quantities over the T-cells,
      ! relative to the origin (0,0) at the geometric center of the cell.
      ! (The geometric center is located at the intersection of the
      ! line joining the midpoints of the north and south faces with
      ! the line joining the midpoints of the east and west faces.
      ! The intersection is assumed to form a right angle.)
      ! These mean values are used to compute the cell centroid, center of 
      ! ice area, and centers of ice and snow volume.
      ! The calculation is done by summing contributions from each of
      ! four triangles, labeled N, E, S, and W.  These triangles are
      ! formed by connecting the geometric center to the four cell 
      ! corners.  Integrals are computed using the method described 
      ! in module ice_transport.F, subroutine sum_fluxes.
      !
      ! cx  = mean value of x
      ! cy  = mean value of y
      ! cxx = mean value of x^2
      ! cxy = mean value of xy
      ! cyy = mean value of y^2
      ! cxxx = mean value of x^3
      ! cxxy = mean value of x^2 * y
      ! cxyy = mean value of x * y^2
      ! cyyy = mean value of y^3
      !------------------------------------------------------------------


      do j = jlo,jhi
         do i = ilo,ihi
            ! Initialize
            cx(i,j) = c0
            cy(i,j) = c0
            cxx(i,j) = c0
            cxy(i,j) = c0
            cyy(i,j) = c0
            cxxx(i,j) = c0
            cxxy(i,j) = c0
            cxyy(i,j) = c0
            cyyy(i,j) = c0

            ! Compute vertices and area of each triangle.
            xv1(1) = c0         ! East triangle
            yv1(1) = c0
            xv2(1) = p5*HTN(i,j)
            yv2(1) = p5*HTE(i,j)
            xv3(1) = p5*HTS(i,j)
            yv3(1) = -p5*HTE(i,j)
            ar(1)  = p25*HTE(i,j)*dxt(i,j)

            xv1(2) = c0         ! West triangle
            yv1(2) = c0
            xv2(2) = -p5*HTN(i,j)
            yv2(2) = p5*HTW(i,j)
            xv3(2) = -p5*HTS(i,j)
            yv3(2) = -p5*HTW(i,j)
            ar(2)  = p25*HTW(i,j)*dxt(i,j)

            xv1(3) = c0         ! North triangle
            yv1(3) = c0
            xv2(3) = -p5*HTN(i,j)
            yv2(3) = p5*HTW(i,j)
            xv3(3) = p5*HTN(i,j)
            yv3(3) = p5*HTE(i,j)
            ar(3)  = p25*HTN(i,j)*dyt(i,j)

            xv1(4) = c0         ! South triangle
            yv1(4) = c0
            xv2(4) = -p5*HTS(i,j)
            yv2(4) = -p5*HTW(i,j)
            xv3(4) = p5*HTS(i,j)
            yv3(4) = -p5*HTE(i,j)
            ar(4)  = p25*HTS(i,j)*dyt(i,j)

            ! Compute contribution to means from each of 
            ! 4 triangles (E, W, N, S)
            do m = 1, 4

               x0 = (xv1(m) + xv2(m) + xv3(m)) / c3 ! midpoint
               y0 = (yv1(m) + yv2(m) + yv3(m)) / c3 
               x1 = p6*x0 + p4*xv1(m) ! other 3 points needed for integral 
               x2 = p6*x0 + p4*xv2(m)
               x3 = p6*x0 + p4*xv3(m)
               y1 = p6*y0 + p4*yv1(m)
               y2 = p6*y0 + p4*yv2(m)
               y3 = p6*y0 + p4*yv3(m)

               cx(i,j) = cx(i,j) + ar(m)*x0
               cy(i,j) = cy(i,j) + ar(m)*y0
               cxx(i,j) = cxx(i,j) + ar(m) * 
     &              ( p5625m * x0*x0 
     &              + p52083 * (x1*x1 + x2*x2 + x3*x3) )
               cxy(i,j) = cxy(i,j) + ar(m) * 
     &              ( p5625m * x0*y0 
     &              + p52083 * (x1*y1 + x2*y2 + x3*y3) )
               cyy(i,j) = cyy(i,j) + ar(m) * 
     &              ( p5625m * y0*y0 
     &              + p52083 * (y1*y1 + y2*y2 + y3*y3) )
               cxxx(i,j) = cxxx(i,j) + ar(m) * 
     &              ( p5625m * x0*x0*x0
     &              + p52083 * (x1*x1*x1 + x2*x2*x2 + x3*x3*x3) )
               cxxy(i,j) = cxxy(i,j) + ar(m) * 
     &              ( p5625m * x0*x0*y0
     &              + p52083 * (x1*x1*y1 + x2*x2*y2 + x3*x3*y3) )
               cxyy(i,j) = cxyy(i,j) + ar(m) * 
     &              ( p5625m * x0*y0*y0
     &              + p52083 * (x1*y1*y1 + x2*y2*y2 + x3*y3*y3) )
               cyyy(i,j) = cyyy(i,j) + ar(m) * 
     &              ( p5625m * y0*y0*y0
     &              + p52083 * (y1*y1*y1 + y2*y2*y2 + y3*y3*y3) )

            enddo               ! m (loop over 4 triangles)
            ! Divide by grid cell area
            cx(i,j) = cx(i,j) * tarear(i,j) 
            cy(i,j) = cy(i,j) * tarear(i,j)
            cxx(i,j) = cxx(i,j) * tarear(i,j)
            cxy(i,j) = cxy(i,j) * tarear(i,j)
            cyy(i,j) = cyy(i,j) * tarear(i,j)
            cxxx(i,j) = cxxx(i,j) * tarear(i,j)
            cxxy(i,j) = cxxy(i,j) * tarear(i,j)
            cxyy(i,j) = cxyy(i,j) * tarear(i,j)
            cyyy(i,j) = cyyy(i,j) * tarear(i,j)
            
         enddo                  ! i
      enddo                     ! j

      ! Compute ghost cell values
      call bound(cx)
      call bound(cy)
      call bound(cxx)
      call bound(cxy)
      call bound(cyy)
      call bound(cxxx)
      call bound(cxxy)
      call bound(cxyy)
      call bound(cyyy)

      !-----------------------------------------------------------------
      ! Compute matrix elements needed to transform from a reference frame
      ! whose origin is at a cell corner (a U-cell frame) to a reference 
      ! frame whose origin is at a cell center (a T-cell frame).
      ! Transformations are needed because the axes of these two reference
      ! frames are not parallel on a curved grid.
      !----------------------------------------------------------------- 

      do j = 1, jmt_local
         do i = 1, imt_local

      !----------------------------------------------------------------- 
      ! Compute angles between the U-cell coordinate axes and T-cell
      ! coordinate axes.  The U-cell axes lie along the cell faces,
      ! which do not meet at right angles in the T-cell reference frame.
      ! The angles are defined as positive if the U-cell axes are obtained 
      ! by a counterclockwise rotation from the T-cell axes.  
      ! Four independent angles are needed.
      !-----------------------------------------------------------------

            theta1 = atan(p5*(HTE(i,j)-HTW(i,j))/HTN(i,j))
            theta2 = atan(p5*(HTN(i,j)-HTS(i,j))/HTW(i,j))
            theta3 = atan(p5*(HTW(i,j)-HTE(i,j))/HTS(i,j))
            theta4 = atan(p5*(HTS(i,j)-HTN(i,j))/HTE(i,j))

      !-----------------------------------------------------------------
      ! Compute matrix elements.
      !-----------------------------------------------------------------

            mne(1,1,i,j) =  cos(theta1)
            mne(2,1,i,j) =  sin(theta1)
            mne(1,2,i,j) = -sin(theta4)
            mne(2,2,i,j) =  cos(theta4)

            mnw(1,1,i,j) =  cos(theta1)
            mnw(2,1,i,j) =  sin(theta1)
            mnw(1,2,i,j) = -sin(theta2)
            mnw(2,2,i,j) =  cos(theta2)

            msw(1,1,i,j) =  cos(theta3)
            msw(2,1,i,j) =  sin(theta3)
            msw(1,2,i,j) = -sin(theta2)
            msw(2,2,i,j) =  cos(theta2)

            mse(1,1,i,j) =  cos(theta3)
            mse(2,1,i,j) =  sin(theta3)
            mse(1,2,i,j) = -sin(theta4)
            mse(2,2,i,j) =  cos(theta4)

         enddo
      enddo

      end subroutine init_remap

c=======================================================================

      subroutine t2ugrid(work)

!---!-------------------------------------------------------------------
!---! transfer from T-cell centers to U-cell centers
!---! writes work into another array that has ghost cells
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: work (ilo:ihi,jlo:jhi)
     &,   work1(imt_local,jmt_local)

      do j=jlo,jhi
       do i=ilo,ihi
        work1(i,j) = work(i,j)
       enddo
      enddo
      call bound(work1)
      call to_ugrid(work1,work)

      end subroutine t2ugrid

c=======================================================================

      subroutine to_ugrid(work1,work2)

!---!-------------------------------------------------------------------
!---! shifts quantities from the T-cell midpoint (work1) to the U-cell 
!---! midpoint (work2)
!---! NOTE results may be wrong near branch cuts
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: work1(imt_local,jmt_local)
     &,   work2(ilo:ihi,jlo:jhi)

      do j=jlo,jhi
       do i=ilo,ihi
       work2(i,j) = p25*(work1(i,j)*tarea(i,j) 
     &                   + work1(i+1,j)*tarea(i+1,j)
     &                   + work1(i,j+1)*tarea(i,j+1) 
     &                   + work1(i+1,j+1)*tarea(i+1,j+1))/uarea(i,j)
       enddo
      enddo

      end subroutine to_ugrid

c=======================================================================

      subroutine u2tgrid(work)

!---!-------------------------------------------------------------------
!---! transfer from U-cell centers to T-cell centers
!---! writes work into another array that has ghost cells
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: work (ilo:ihi,jlo:jhi)
     &,   work1(imt_local,jmt_local)

      do j=jlo,jhi
       do i=ilo,ihi
        work1(i,j) = work(i,j)
       enddo
      enddo
      call bound(work1)
      call to_tgrid(work1,work)

      end subroutine u2tgrid

c=======================================================================

      subroutine to_tgrid(work1,work2)

!---!-------------------------------------------------------------------
!---! shifts quantities from the U-cell midpoint (work1) to the T-cell 
!---! midpoint (work2)
!---! NOTE results may be wrong near branch cuts
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: work1(imt_local,jmt_local)
     &,   work2(ilo:ihi,jlo:jhi)

      do j=jlo,jhi
       do i=ilo,ihi
       work2(i,j) = p25*(work1(i,j)*uarea(i,j) 
     &                   + work1(i-1,j)*uarea(i-1,j)
     &                   + work1(i,j-1)*uarea(i,j-1) 
     &                   + work1(i-1,j-1)*uarea(i-1,j-1))/tarea(i,j)
       enddo
      enddo

      end subroutine to_tgrid

c=======================================================================

      subroutine bound(work1)

!---!-------------------------------------------------------------------
!---! fills ghost cells with boundary information
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      real (kind=dbl_kind) :: work1(1)
      call bound_ijn(1,work1,.true.,.true.,.true.,.true.)

      end subroutine bound

c=======================================================================

      subroutine bound_sw(work1)

!---!-------------------------------------------------------------------
!---! fills south and west ghost cells with boundary information
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      real (kind=dbl_kind) :: work1(1)
      call bound_ijn(1,work1,.false.,.true.,.false.,.true.)

      end subroutine bound_sw

c=======================================================================

      subroutine bound_narr_ne(narrays,work1)

!---!-------------------------------------------------------------------
!---! fills north and east ghost cells with boundary information
!---!
!---! NOTE: work1 array has form (number_arrays,i_index,j_index)
!---!       for evp dynamics performance
!---!
!---! authors Tony Craig, NCAR
!---!         Elizabeth Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: narrays
      real (kind=dbl_kind) :: work1(1)
      call bound_nij(narrays,work1,.true.,.false.,.true.,.false.)

      end subroutine bound_narr_ne

c=======================================================================

      subroutine bound_narr(narrays,work1)

!---!-------------------------------------------------------------------
!---! fills north and east ghost cells with boundary information
!---! narr arrays at once (for performance)
!---!
!---! authors Tony Craig, NCAR
!---!         Elizabeth Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: narrays
      real (kind=dbl_kind) :: work1(1)
      call bound_ijn(narrays,work1,.true.,.true.,.true.,.true.)

      end subroutine bound_narr

c=======================================================================

      subroutine bound_ijn(nd,work1,north,south,east,west)

!---!-------------------------------------------------------------------
!---! Periodic/Neumann conditions for global domain boundaries
!---! Assumptions:  a *single* row of ghost cells (num_ghost_cells=1)
!---!
!---! work1 array has form (i_index,j_index,number_arrays)
!---!
!---! authors Tony Craig, NCAR
!---!         Elizabeth Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_timers
      use ice_mpi_internal

      integer (kind=int_kind) :: i, j, nd, n
      real (kind=dbl_kind) :: work1(imt_local,jmt_local,nd)
      logical north,south,east,west
#ifdef _MPI
      integer (kind=int_kind) :: icnt,jcnt
     &,   status(MPI_STATUS_SIZE),request(4)
      real (kind=dbl_kind) :: workw(jlo:jhi,nd),worke(jlo:jhi,nd)
     &,  workn(ilo-1:ihi+1,nd),works(ilo-1:ihi+1,nd)
#endif

      call ice_timer_start(10) ! bound

#ifdef _MPI
      jcnt=(jhi-jlo+1)*nd
      icnt=(ihi-ilo+1+2*num_ghost_cells)*nd

      !-----------------------------------------------------------------
      ! west data to east data, west shift
      !-----------------------------------------------------------------
      if (east) then

      do n=1,nd
      do j=jlo,jhi
      workw(j,n)=work1(ilo,j,n)
      enddo
      enddo

      call MPI_SENDRECV(workw,jcnt,MPI_REAL8,nbr_west,my_task,
     &                  worke,jcnt,MPI_REAL8,nbr_east,nbr_east,
     &                  MPI_COMM_ICE,status,ierr)

      do n=1,nd
      do j=jlo,jhi
      work1(ihi+1,j,n)=worke(j,n)
      enddo
      enddo

      endif

      !-----------------------------------------------------------------
      ! east data to west data, east shift
      !-----------------------------------------------------------------
      if (west) then

      do n=1,nd
      do j=jlo,jhi
      worke(j,n)=work1(ihi,j,n)
      enddo
      enddo

      call MPI_SENDRECV(worke,jcnt,MPI_REAL8,nbr_east,my_task,
     &                  workw,jcnt,MPI_REAL8,nbr_west,nbr_west,
     &                  MPI_COMM_ICE,status,ierr)

      do n=1,nd
      do j=jlo,jhi
      work1(ilo-1,j,n)=workw(j,n)
      enddo
      enddo

      endif

      !-----------------------------------------------------------------
      ! north data to south data, north shift
      !-----------------------------------------------------------------
      if (south) then
      if (nbr_south.ne.-1) then
        call MPI_IRECV(works,
     &                 icnt,MPI_REAL8,nbr_south,nbr_south,
     &                 MPI_COMM_ICE,request(1),ierr)
      else

        do n=1,nd
        do i=ilo-1,ihi+1
        work1(i,jlo-1,n)=work1(i,jlo,n)
        enddo
        enddo

      endif
      if (nbr_north.ne.-1) then

        do n=1,nd
        do i=ilo-1,ihi+1
        workn(i,n)=work1(i,jhi,n)
        enddo
        enddo

        call MPI_ISEND (workn,
     &                 icnt,MPI_REAL8,nbr_north,my_task,
     &                 MPI_COMM_ICE,request(2),ierr)
      endif
      if (nbr_north.ne.-1) then
        call MPI_WAIT(request(2), status, ierr)
      endif
      if (nbr_south.ne.-1) then
        call MPI_WAIT(request(1), status, ierr)

        do n=1,nd
        do i=ilo-1,ihi+1
        work1(i,jlo-1,n)=works(i,n)
        enddo
        enddo

      endif
      endif

      !-----------------------------------------------------------------
      ! south data to north data, south shift
      !-----------------------------------------------------------------
      if (north) then
      if (nbr_north.ne.-1) then
        call MPI_IRECV(workn,
     &                 icnt,MPI_REAL8,nbr_north,nbr_north,
     &                 MPI_COMM_ICE,request(3),ierr)
      else

        do n=1,nd
        do i=ilo-1,ihi+1
        work1(i,jhi+1,n)=work1(i,jhi,n)
        enddo
        enddo

      endif
      if (nbr_south.ne.-1) then

        do n=1,nd
        do i=ilo-1,ihi+1
        works(i,n)=work1(i,jlo,n)
        enddo
        enddo

        call MPI_ISEND (works,
     &                 icnt,MPI_REAL8,nbr_south,my_task,
     &                 MPI_COMM_ICE,request(4),ierr)
      endif
      if (nbr_south.ne.-1) then
        call MPI_WAIT(request(4), status, ierr)
      endif
      if (nbr_north.ne.-1) then
        call MPI_WAIT(request(3), status, ierr)

        do n=1,nd
        do i=ilo-1,ihi+1
        work1(i,jhi+1,n)=workn(i,n)
        enddo
        enddo

      endif
      endif

#else
      !-----------------------------------------------------------------
      ! single domain
      !-----------------------------------------------------------------
 
      do n=1,nd

      ! Periodic conditions
      do j=jlo,jhi
       work1(ilo-1,j,n) = work1(ihi,j,n)
       work1(ihi+1,j,n) = work1(ilo,j,n)
      enddo

      ! Neumann conditions (POP grid land points)
      do i=ilo-1,ihi+1
        work1(i,jlo-1,n) = work1(i,jlo,n)
        work1(i,jhi+1,n) = work1(i,jhi,n)
      enddo

      enddo  ! n
#endif
      call ice_timer_stop(10)  ! bound

      end subroutine bound_ijn

c=======================================================================

      subroutine bound_nij(nd,work1,north,south,east,west)

!---!-------------------------------------------------------------------
!---! Periodic/Neumann conditions for global domain boundaries
!---! Assumptions:  a *single* row of ghost cells (num_ghost_cells=1)
!---!
!---! work1 array has form (number_arrays,i_index,j_index)
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      use ice_timers
      use ice_mpi_internal

      integer (kind=int_kind) :: i, j, nd, n
      real (kind=dbl_kind) :: work1(nd,imt_local,jmt_local)
      logical north,south,east,west
#ifdef _MPI
      integer (kind=int_kind) :: icnt,jcnt
     &,   status(MPI_STATUS_SIZE),request(4)
      real (kind=dbl_kind) :: workw(nd,jlo:jhi),worke(nd,jlo:jhi)
     &,  workn(nd,ilo-1:ihi+1),works(nd,ilo-1:ihi+1)
#endif

      call ice_timer_start(10) ! bound

#ifdef _MPI
      jcnt=(jhi-jlo+1)*nd
      icnt=(ihi-ilo+1+2*num_ghost_cells)*nd

      !-----------------------------------------------------------------
      ! west data to east data, west shift
      !-----------------------------------------------------------------
      if (east) then

      do j=jlo,jhi
      do n=1,nd
      workw(n,j)=work1(n,ilo,j)
      enddo
      enddo

      call MPI_SENDRECV(workw,jcnt,MPI_REAL8,nbr_west,my_task,
     &                  worke,jcnt,MPI_REAL8,nbr_east,nbr_east,
     &                  MPI_COMM_ICE,status,ierr)

      do j=jlo,jhi
      do n=1,nd
      work1(n,ihi+1,j)=worke(n,j)
      enddo
      enddo

      endif

      !-----------------------------------------------------------------
      ! east data to west data, east shift
      !-----------------------------------------------------------------
      if (west) then

      do j=jlo,jhi
      do n=1,nd
      worke(n,j)=work1(n,ihi,j)
      enddo
      enddo

      call MPI_SENDRECV(worke,jcnt,MPI_REAL8,nbr_east,my_task,
     &                  workw,jcnt,MPI_REAL8,nbr_west,nbr_west,
     &                  MPI_COMM_ICE,status,ierr)

      do j=jlo,jhi
      do n=1,nd
      work1(n,ilo-1,j)=workw(n,j)
      enddo
      enddo

      endif

      !-----------------------------------------------------------------
      ! north data to south data, north shift
      !-----------------------------------------------------------------
      if (south) then
      if (nbr_south.ne.-1) then
        call MPI_IRECV(works,
     &                 icnt,MPI_REAL8,nbr_south,nbr_south,
     &                 MPI_COMM_ICE,request(1),ierr)
      else

        do i=ilo-1,ihi+1
        do n=1,nd
        work1(n,i,jlo-1)=work1(n,i,jlo)
        enddo
        enddo

      endif
      if (nbr_north.ne.-1) then

        do i=ilo-1,ihi+1
        do n=1,nd
        workn(n,i)=work1(n,i,jhi)
        enddo
        enddo

        call MPI_ISEND (workn,
     &                 icnt,MPI_REAL8,nbr_north,my_task,
     &                 MPI_COMM_ICE,request(2),ierr)
      endif
      if (nbr_north.ne.-1) then
        call MPI_WAIT(request(2), status, ierr)
      endif
      if (nbr_south.ne.-1) then
        call MPI_WAIT(request(1), status, ierr)

        do i=ilo-1,ihi+1
        do n=1,nd
        work1(n,i,jlo-1)=works(n,i)
        enddo
        enddo

      endif
      endif

      !-----------------------------------------------------------------
      ! south data to north data, south shift
      !-----------------------------------------------------------------
      if (north) then
      if (nbr_north.ne.-1) then
        call MPI_IRECV(workn,
     &                 icnt,MPI_REAL8,nbr_north,nbr_north,
     &                 MPI_COMM_ICE,request(3),ierr)
      else

        do i=ilo-1,ihi+1
        do n=1,nd
        work1(n,i,jhi+1)=work1(n,i,jhi)
        enddo
        enddo

      endif
      if (nbr_south.ne.-1) then

        do i=ilo-1,ihi+1
        do n=1,nd
        works(n,i)=work1(n,i,jlo)
        enddo
        enddo

        call MPI_ISEND (works,
     &                 icnt,MPI_REAL8,nbr_south,my_task,
     &                 MPI_COMM_ICE,request(4),ierr)
      endif
      if (nbr_south.ne.-1) then
        call MPI_WAIT(request(4), status, ierr)
      endif
      if (nbr_north.ne.-1) then
        call MPI_WAIT(request(3), status, ierr)

        do i=ilo-1,ihi+1
        do n=1,nd
        work1(n,i,jhi+1)=workn(n,i)
        enddo
        enddo

      endif
      endif

#else
      !-----------------------------------------------------------------
      ! single domain
      !-----------------------------------------------------------------

      ! Periodic conditions
      do j=jlo,jhi
      do n=1,nd
       work1(n,ilo-1,j) = work1(n,ihi,j)
       work1(n,ihi+1,j) = work1(n,ilo,j)
      enddo
      enddo

      ! Neumann conditions (POP grid land points)
      do i=ilo-1,ihi+1
      do n=1,nd
        work1(n,i,jlo-1) = work1(n,i,jlo)
        work1(n,i,jhi+1) = work1(n,i,jhi)
      enddo
      enddo

#endif
      call ice_timer_stop(10)  ! bound

      end subroutine bound_nij

c=======================================================================

      end module ice_grid

c=======================================================================
