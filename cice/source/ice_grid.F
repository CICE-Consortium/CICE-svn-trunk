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

      implicit none

      character (len=char_len) :: 
     &   grid_file         !  input file for POP grid info
     &,  kmt_file          !  input file for POP grid info
     &,  grid_type         !  rectangular (default) or displaced_pole 

      real (kind=dbl_kind), dimension (imt_global,jmt_global) ::
     &   TLAT_G   ! latitude  of cell center
     &,  TLON_G   ! longitude of cell center

      real (kind=dbl_kind), dimension (imt_local,jmt_local) ::
     &   dxt  ! width of T-cell through the middle (m)
     &,  dyt  ! height of T-cell through the middle (m)
     &,  dxu  ! width of U-cell through the middle (m)
     &,  dyu  ! height of U-cell through the middle (m)
     &,  HTE  ! length of eastern edge of T-cell (m)
     &,  HTN  ! length of northern edge of T-cell (m)
     &,  tarea! area of T-cell (m^2)
     &,  uarea! area of U-cell (m^2)
     &,  ULON ! longitude of velocity pts (radians)
     &,  ULAT ! latitude of velocity pts (radians)

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::
     &   cyp   ! 1.5*HTE - 0.5*HTE
     &,  cxp   ! 1.5*HTN - 0.5*HTN
     &,  cym   ! 0.5*HTE - 1.5*HTE
     &,  cxm   ! 0.5*HTN - 1.5*HTN
     &,  dxt2   ! 0.5*dxt
     &,  dyt2   ! 0.5*dyt
     &,  dxt4   ! 0.25*dxt
     &,  dyt4   ! 0.25*dyt
     &,  dxhy   ! 0.5*(HTE - HTE)
     &,  dyhx   ! 0.5*(HTN - HTN)
     &,  tarear   ! 1/tarea
     &,  uarear   ! 1/uarea
     &,  tinyarea ! puny*tarea
     &,  ANGLE    ! for conversions between POP grid and lat/lon
     &,  ANGLET   ! ANGLE converted to T-cells
     &,  tarean      ! area of NH cells
     &,  tareas      ! area of SH cells

      ! Masks
      real (kind=dbl_kind), dimension (imt_local,jmt_local) ::
     &   hm   ! land/boundary mask, thickness (T-cell)
     &,  uvm  ! land/boundary mask, velocity (U-cell)
     &,  mask_n   ! northern hemisphere
     &,  mask_s   ! southern hemisphere

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
     &     tmpg(imt_global,jmt_global)   ! temporary array
     &,    tmpg4(4,imt_global,jmt_global)   ! temporary array

      if (grid_type .eq. 'displaced_pole') then
        call popgrid              ! read POP grid lengths directly
      else
        call rectgrid             ! regular rectangular grid
      endif

      call bound(HTN)
      call bound(HTE)
      call bound(ULAT)
      call bound(ULON)

      !-----------------------------------------------------------------
      ! construct T-grid cell and U-grid cell widths
      !-----------------------------------------------------------------

c$OMP PARALLEL DO PRIVATE(i,j)
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

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        uarea(i,j) = p25*(tarea(i,j) + tarea(i+1,j)
     &         + tarea(i,j+1) + tarea(i+1,j+1))        ! m^2
       enddo
      enddo
      call bound(uarea)

      ! grid length combinations
c$OMP PARALLEL DO PRIVATE(i,j)
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

      call ice_open(11,grid_file,64)
      call ice_open(12,kmt_file,32)

      scatter = .true.
      diag = .true.

      call ice_read(11,1,work,'rda8',scatter,diag)
      ULAT(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(11,2,work,'rda8',scatter,diag)
      ULON(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(11,3,work,'rda8',scatter,diag)
      HTN(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)*cm_to_m
      call ice_read(11,4,work,'rda8',scatter,diag)
      HTE(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)*cm_to_m
      call ice_read(11,7,work,'rda8',scatter,diag)
      ANGLE(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(12,1,work,'ida4',scatter,diag)

      if (my_task.eq.master_task) then
         close (11) 
         close (12) 
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
c$OMP PARALLEL DO PRIVATE(i,j)
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
         ULAT (i,j) = c0      ! remember to set Coriolis !
         ULON (i,j) = c0
         ANGLE(i,j) = c0      ! "square with the world"
         ANGLET(i,j) = c0 
       enddo
      enddo

      !-----------------------------------------------------------------
      ! construct T-cell land mask
      !-----------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(i,j)
        do j=1,jmt_global        ! initialize hm
         do i=1,imt_global       ! open
          hmg(i,j) = c0
         enddo
        enddo

c        do j=1,jmt_global        ! open
c         do i=1,imt_global       ! open
c$OMP PARALLEL DO PRIVATE(i,j)
        do j=3,jmt_global-2       ! closed: NOTE jmt_global > 5
         do i=3,imt_global-2      ! closed: NOTE imt_global > 5
          hmg(i,j) = c1
         enddo
        enddo

      call global_scatter(hmg,work)

c$OMP PARALLEL DO PRIVATE(i,j)
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

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        uvm(i,j) = min(hm(i,j),hm(i+1,j),hm(i,j+1),hm(i+1,j+1))
       enddo
      enddo
      call bound(uvm)  !!! use real arrays to get boundary conditions

c$OMP PARALLEL DO PRIVATE(i,j)
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

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        mask_n(i,j) = c0
        mask_s(i,j) = c0
       enddo
      enddo
c$OMP PARALLEL DO PRIVATE(i,j)
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

      use ice_read_write ! if reading ULAT, ULON directly from file

      integer (kind=int_kind) ::
     &   i, j                          ! generic indices

      integer (kind=int_kind) :: im1
      real (kind=dbl_kind) ::
     &     z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da
     &,    ULATG(imt_global,jmt_global)   ! latitude  of NE cell corner
     &,    ULONG(imt_global,jmt_global)   ! longitude of NE cell corner

c      call global_gather(ULONG,ULON(ilo:ihi,jlo:jhi))
c      call global_gather(ULATG,ULAT(ilo:ihi,jlo:jhi))

      if (my_task == master_task) then
      call ice_open(11,grid_file,64)
      read(11,rec=1) ULATG
      read(11,rec=2) ULONG
      close (11) 

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

      write (6,*) 'min/max ULONG',minval(ULONG),maxval(ULONG)
      write (6,*) 'min/max TLON_G',minval(TLON_G),maxval(TLON_G)
      write (6,*) 'min/max ULATG',minval(ULATG),maxval(ULATG)
      write (6,*) 'min/max TLAT_G',minval(TLAT_G),maxval(TLAT_G)
      endif ! master_task

      end subroutine Tlatlon

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

c$OMP PARALLEL DO PRIVATE(i,j)
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
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: work1(imt_local,jmt_local)
     &,   work2(ilo:ihi,jlo:jhi)

c$OMP PARALLEL DO PRIVATE(i,j)
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

c$OMP PARALLEL DO PRIVATE(i,j)
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
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) :: work1(imt_local,jmt_local)
     &,   work2(ilo:ihi,jlo:jhi)

c$OMP PARALLEL DO PRIVATE(i,j)
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
c$OMP PARALLEL DO PRIVATE(j)
      do j=jlo,jhi
       work1(ilo-1,j,n) = work1(ihi,j,n)
       work1(ihi+1,j,n) = work1(ilo,j,n)
      enddo

      ! Neumann conditions (POP grid land points)
c$OMP PARALLEL DO PRIVATE(i)
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
c$OMP PARALLEL DO PRIVATE(j)
      do j=jlo,jhi
      do n=1,nd
       work1(n,ilo-1,j) = work1(n,ihi,j)
       work1(n,ihi+1,j) = work1(n,ilo,j)
      enddo
      enddo

      ! Neumann conditions (POP grid land points)
c$OMP PARALLEL DO PRIVATE(i)
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
