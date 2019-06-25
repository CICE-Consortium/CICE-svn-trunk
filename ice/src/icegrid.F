c $Id: icegrid.F,v 2.2 2000/02/03 21:34:50 eclare Exp $

c-----------------------------------------------------------------------
c   spatial grids, masks, boundary conditions and calendar
c-----------------------------------------------------------------------
      subroutine grid
c.. Horizontal grid initialization
c..     HT{N,E} = cell widths on {N,E} sides of T cell
c..     U{LAT,LONG} = true {latitude,longitude} of U points
c..     D{X,Y}{T,U} = {x,y} spacing centered at {T,U} points
      implicit none
      include "ice.h"            ! ice code common blocks
      integer i, j

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        hm(i,j) = 0.             ! initialize T-cell mask
       enddo
      enddo

      call popgrid               ! read POP grid lengths directly
c      call rectgrid             ! regular rectangular grid

      call bound(HTN)
      call bound(HTE)
      call bound(ULAT)
      call bound(ULONG)

      call makemask              ! velocity mask, hemisphere masks

c.. construct T-grid cell and U-grid cell widths

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        dxt(i,j) = 0.5*(HTN(i,j) + HTN(i,j-1))
        dyt(i,j) = 0.5*(HTE(i,j) + HTE(i-1,j))

        dxtr(i,j) = 1./dxt(i,j)
        dytr(i,j) = 1./dyt(i,j)
        dxtr4(i,j) = 0.25/dxt(i,j)
        dytr4(i,j) = 0.25/dyt(i,j)

        tarea(i,j) = dxt(i,j)*dyt(i,j)

        dxu(i,j) = 0.5*(HTN(i,j) + HTN(i+1,j))
        dyu(i,j) = 0.5*(HTE(i,j) + HTE(i,j+1))
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
        uarea(i,j) = 0.25*(tarea(i,j) + tarea(i+1,j)
     &         + tarea(i,j+1) + tarea(i+1,j+1))        ! m^2
       enddo
      enddo
      call bound(uarea)

      call grid_netcdf   ! write latitude, longitude, land masks

      return
      end
c-----------------------------------------------------------------------
      subroutine popgrid
c.. POP displaced pole grid and land mask
c..   grid
c..      rec no.         field         units
c..      -------         -----         -----
c..         1            ULAT         radians
c..         2            ULON         radians
c..         3             HTN           cm
c..         4             HTE           cm
c..         5             HUS           cm
c..         6             HUW           cm
c..         7            ANGLE        radians
c..   land mask
c..      rec no.         field         units
c..      -------         -----         -----
c..         1             KMT
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real work(ilo:ihi,jlo:jhi)

      call ice_open(11,grid_file,64)
      call ice_open(12,kmt_file,32)

      call ice_read(11,1,work,'rda8')
      ULAT(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(11,2,work,'rda8')
      ULONG(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(11,3,work,'rda8')
      HTN(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)/100.
      call ice_read(11,4,work,'rda8')
      HTE(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)/100.
      call ice_read(11,7,work,'rda8')
      ANGLE(ilo:ihi,jlo:jhi)=work(ilo:ihi,jlo:jhi)
      call ice_read(12,1,work,'ida4')
      do j=jlo,jhi
      do i=ilo,ihi
        hm(i,j)  = work(i,j)
        if (hm(i,j).ge.1.) hm(i,j) = 1.
      enddo
      enddo

      if (my_task.eq.master_task) then
         close (11) 
         close (12) 
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rectgrid
c.. Regular rectangular grid and mask
      implicit none
      include "ice.h"            ! ice code common blocks
      integer i, j
      real hmg(imt_global,jmt_global)

c.. calculate various geometric 2d arrays
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
         HTN  (i,j) = 1.8e4      ! constant longitude spacing, m
         HTE  (i,j) = 1.8e4      ! constant latitude  spacing, m
         ULAT (i,j) = 0.         ! remember to set Coriolis !
         ULONG(i,j) = 0.
         ANGLE(i,j) = 0.         ! "square with the world"

c         HTN  (i,j) = 3.1e4     ! POP <4/3> min, m
c         HTE  (i,j) = 3.1e4     ! POP <4/3> min, m
c         HTN  (i,j) = 1.6e4     ! POP <2/3> min, m
c         HTE  (i,j) = 1.6e4     ! POP <2/3> min, m
       enddo
      enddo

c.. construct T-cell land mask
c$OMP PARALLEL DO PRIVATE(i,j)
        do j=1,jmt_global        ! initialize hm
         do i=1,imt_global       ! open
          hmg(i,j) = 0.
         enddo
        enddo

c        do j=1,jmt_global        ! open
c         do i=1,imt_global       ! open
c$OMP PARALLEL DO PRIVATE(i,j)
        do j=3,jmt_global-2       ! closed: NOTE jmt_global > 5
         do i=3,imt_global-2      ! closed: NOTE imt_global > 5
          hmg(i,j) = 1.
         enddo
        enddo

      call global_scatter(hmg,hm)

      return
      end
c-----------------------------------------------------------------------
      subroutine makemask
c.. Sets the boundary values for the T cell land mask (hm) and
c.. makes the land mask for U cells (uvm) from hm
c.. Also creates hemisphere masks (mask_n northern, mask_s southern)
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j

      call bound(hm)

c.. construct velocity mask uvm

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        uvm(i,j) = 0.
       enddo
      enddo

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        uvm(i,j) = min(hm(i,j),hm(i+1,j),hm(i,j+1),hm(i+1,j+1))
       enddo
      enddo
      call bound(uvm)

c.. create hemisphere masks

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
       do i=1,imt_local
        mask_n(i,j) = 0.
        mask_s(i,j) = 0.
       enddo
      enddo
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        if (ULAT(i,j).gt.0.) mask_n(i,j) = 1.  ! northern hemisphere
        if (ULAT(i,j).lt.0.) mask_s(i,j) = 1.  ! southern hemisphere
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine t2ugrid(work)
c.. transfer from T-cell centers to U-cell centers
c.. writes work into another array that has ghost cells
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real work (ilo:ihi,jlo:jhi)
      real work1(imt_local,jmt_local)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        work1(i,j) = work(i,j)
       enddo
      enddo
      call bound(work1)
      call to_ugrid(work1,work)

      return
      end
c-----------------------------------------------------------------------
      subroutine to_ugrid(work1,work2)
c.. shifts quantities from the T-cell midpoint (work1) to the U-cell 
c.. midpoint (work2)
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real work1(imt_local,jmt_local)
      real work2(ilo:ihi,jlo:jhi)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
       work2(i,j) = 0.25*(work1(i,j)*tarea(i,j) 
     &                   + work1(i+1,j)*tarea(i+1,j)
     &                   + work1(i,j+1)*tarea(i,j+1) 
     &                   + work1(i+1,j+1)*tarea(i+1,j+1))/uarea(i,j)
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine u2tgrid(work)
c.. transfer from U-cell centers to T-cell centers
c.. writes work into another array that has ghost cells
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real work (ilo:ihi,jlo:jhi)
      real work1(imt_local,jmt_local)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        work1(i,j) = work(i,j)
       enddo
      enddo
      call bound(work1)
      call to_tgrid(work1,work)

      return
      end
c-----------------------------------------------------------------------
      subroutine to_tgrid(work1,work2)
c.. shifts quantities from the U-cell midpoint (work1) to the T-cell 
c.. midpoint (work2)
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j
      real work1(imt_local,jmt_local)
      real work2(ilo:ihi,jlo:jhi)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
       work2(i,j) = 0.25*(work1(i,j)*uarea(i,j) 
     &                   + work1(i-1,j)*uarea(i-1,j)
     &                   + work1(i,j-1)*uarea(i,j-1) 
     &                   + work1(i-1,j-1)*uarea(i-1,j-1))/tarea(i,j)
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine bound(work1)
      implicit none
      real work1(1)
      call bound_dir(1,work1,1,1,1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine bound_sw(work1)
      implicit none
      real work1(1)
      call bound_dir(1,work1,0,1,0,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine bound_ne(work1)
      implicit none
      real work1(1)
      call bound_dir(1,work1,1,0,1,0)
      return
      end
c-----------------------------------------------------------------------
      subroutine bound8_ne(work1)
      implicit none
      real work1(1)
      call bound_dir(8,work1,1,0,1,0)
      return
      end
c-----------------------------------------------------------------------
      subroutine bound_dir(nd,work1,north,south,east,west)
c.. Periodic/Neumann conditions for global domain boundaries
c.. Assumptions:  a *single* local domain
c..               a *single* row of ghost cells (nghost=1)
      implicit none
      include "ice.h"          ! ice code common blocks
      include "ice_mpi.h"      ! common blocks for mpi
      integer i, j, nd
      real work1(nd,imt_local,jmt_local)
      logical north,south,east,west
#ifdef _MPI
      include "mpif.h"         ! MPI library definitions
      integer icnt,jcnt
      real workw(nd,jlo:jhi),worke(nd,jlo:jhi)
      real workn(nd,ilo-1:ihi+1),works(nd,ilo-1:ihi+1)
      integer status(MPI_STATUS_SIZE),request(4)
#endif

      call timer_start(8)

#ifdef _MPI
      jcnt=(jhi-jlo+1)*nd
      icnt=(ihi-ilo+1+2*nghost)*nd

c west data to east data, west shift
      if (east) then
      workw=work1(:,ilo,jlo:jhi)
      call MPI_SENDRECV(workw,jcnt,MPI_REAL8,nbr_west,my_task,
     &                  worke,jcnt,MPI_REAL8,nbr_east,nbr_east,
     &                  MPI_COMM_ICE,status,ierr)
      work1(:,ihi+1,jlo:jhi)=worke
      endif

c east data to west data, east shift
      if (west) then
      worke=work1(:,ihi,jlo:jhi)
      call MPI_SENDRECV(worke,jcnt,MPI_REAL8,nbr_east,my_task,
     &                  workw,jcnt,MPI_REAL8,nbr_west,nbr_west,
     &                  MPI_COMM_ICE,status,ierr)
      work1(:,ilo-1,jlo:jhi)=workw
      endif

c north data to south data, north shift
      if (south) then
      if (nbr_south.ne.-1) then
        call MPI_IRECV(works,
     &                 icnt,MPI_REAL8,nbr_south,nbr_south,
     &                 MPI_COMM_ICE,request(1),ierr)
      else
        work1(:,:,jlo-1)=work1(:,:,jlo)
      endif
      if (nbr_north.ne.-1) then
        workn=work1(:,ilo-1:ihi+1,jhi)
        call MPI_ISEND (workn,
     &                 icnt,MPI_REAL8,nbr_north,my_task,
     &                 MPI_COMM_ICE,request(2),ierr)
      endif
      endif

c south data to north data, south shift
      if (north) then
      if (nbr_north.ne.-1) then
        call MPI_IRECV(workn,
     &                 icnt,MPI_REAL8,nbr_north,nbr_north,
     &                 MPI_COMM_ICE,request(3),ierr)
      else
        work1(:,:,jhi+1)=work1(:,:,jhi)
      endif
      if (nbr_south.ne.-1) then
        works=work1(:,ilo-1:ihi+1,jlo)
        call MPI_ISEND (works,
     &                 icnt,MPI_REAL8,nbr_south,my_task,
     &                 MPI_COMM_ICE,request(4),ierr)
      endif
      endif

      if (south.and.nbr_north.ne.-1) then
        call MPI_WAIT(request(2), status, ierr)
      endif
      if (north.and.nbr_south.ne.-1) then
        call MPI_WAIT(request(4), status, ierr)
      endif
      if (south.and.nbr_south.ne.-1) then
        call MPI_WAIT(request(1), status, ierr)
        work1(:,ilo-1:ihi+1,jlo-1)=works
      endif
      if (north.and.nbr_north.ne.-1) then
        call MPI_WAIT(request(3), status, ierr)
        work1(:,ilo-1:ihi+1,jhi+1)=workn
      endif

c      call MPI_barrier(MPI_COMM_ICE,ierr)

#else
c.. Periodic conditions
c$OMP PARALLEL DO PRIVATE(j)
      do j=jlo,jhi
       work1(:,ilo-1,j) = work1(:,ihi,j)
       work1(:,ihi+1,j) = work1(:,ilo,j)
      enddo

c.. Neumann conditions (POP grid land points)
c$OMP PARALLEL DO PRIVATE(i)
      do i=ilo-1,ihi+1
        work1(:,i,jlo-1) = work1(:,i,jlo)
        work1(:,i,jhi+1) = work1(:,i,jhi)
      enddo

#endif
      call timer_stop(8)

      return
      end
c--------------------------------------------------------------
      subroutine global_scatter(workg,work)
      implicit none
      include "ice.h"          ! ice code common blocks
      include "ice_mpi.h"      ! common blocks for mpi
      integer i,j
      real workg(imt_global,jmt_global)
      real work(ilo:ihi,jlo:jhi)

#ifdef _MPI
      include "mpif.h"         ! MPI library definitions
      integer request(2),n
      integer status(MPI_STATUS_SIZE)
      work=0.
      call MPI_IRECV(work(ilo,jlo), 1, mpi_interior_real(my_task+1), 
     &         master_task, my_task, MPI_COMM_ICE, request(2), ierr)
      if (my_task .eq. master_task) then
        do n=1,nproc_s
          call MPI_ISEND(workg(local_start(1,n),local_start(2,n)), 
     &                  1, mpi_interior_real_global(n), n-1,
     &                  n-1, MPI_COMM_ICE, request(1), ierr)
          call MPI_WAIT(request(1), status, ierr)
        end do
      endif
      call MPI_WAIT(request(2), status, ierr)
c      call MPI_barrier(MPI_COMM_ICE,ierr)
#else
      ! shift indices from global domain to local domain
c      work(ilo:ihi,jlo:jhi)=workg(1:imt_global,1:jmt_global)
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
        work(i,j)=workg(i-nghost,j-nghost)
      enddo
      enddo
#endif

      return
      end
c--------------------------------------------------------------
      subroutine global_gather(workg,work)
      implicit none
      include "ice.h"          ! ice code common blocks
      include "ice_mpi.h"      ! common blocks for mpi
      integer i,j
      real workg(imt_global,jmt_global)
      real work(ilo:ihi,jlo:jhi)

#ifdef _MPI
      include "mpif.h"         ! MPI library definitions
      integer request(2),recv_req(nproc_s),n
      integer status(MPI_STATUS_SIZE)
      integer recv_status(MPI_STATUS_SIZE,nproc_s)

      if (my_task .eq. master_task) then
        workg=0.
        do n=1,nproc_s
          call MPI_IRECV(workg(local_start(1,n),local_start(2,n)), 
     &                  1, mpi_interior_real_global(n), n-1,
     &                  n-1, MPI_COMM_ICE, recv_req(n), ierr)
        end do
      endif
      call MPI_ISEND(work(ilo,jlo), 1, mpi_interior_real(my_task+1), 
     &         master_task, my_task, MPI_COMM_ICE, request(2), ierr)
      call MPI_WAIT(request(2), status, ierr)

      if (my_task .eq. master_task) then
          call MPI_WAITALL(nproc_s, recv_req, recv_status, ierr)
      endif
c      call MPI_barrier(MPI_COMM_ICE,ierr)
#else
      ! shift indices from local domain to global domain
c      workg(1:imt_global,1:jmt_global)=work(ilo:ihi,jlo:jhi)
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
        workg(i-nghost,j-nghost)=work(i,j)
      enddo
      enddo
#endif

      return
      end
c-----------------------------------------------------------------------
