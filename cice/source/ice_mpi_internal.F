c $Id: ice_mpi.H,v 1.3 2000/01/23 02:51:47 tcraig Exp $
c=======================================================================
!---! parameters and commons blocks for MPI parallelization 
!---! internal to ice model
!---!
!---! authors Tony Craig, NCAR
!---!         Elizabeth C. Hunke, LANL
c=======================================================================

      module ice_mpi_internal

      use ice_kinds_mod
      use ice_domain

      implicit none

      integer (kind=int_kind), parameter ::
     &   nproc_x = NPROC_X   ! number of processors assigned in x direction
     &,  nproc_y = NPROC_Y   ! number of processors assigned in y direction
     &,  nproc_s = nproc_x*nproc_y  ! total number of processors assigned, 
                                    ! both directions

      integer (kind=int_kind) ::
     &   nb_tasks               ! total number of tasks
     &,  ierr                   ! general-use error flag
     &,  MPI_COMM_ICE           ! communicator for internal ice comms

#ifdef _MPI

      integer (kind=int_kind) ::
     &   nbr_west, nbr_east,        ! processor numbers for the 
     &   nbr_north, nbr_south       ! n,s,e,w neighbor processors

      integer (kind=int_kind), dimension(2,nproc_s) ::
     &   local_start                ! global i,j start for each processor

      integer (kind=int_kind), dimension(nproc_s) ::
     &   mpi_interior_int,          ! MPI data type definitions for
     &   mpi_interior_real,         ! physical domain blocks
     &   mpi_interior_int_global,   ! in case of padding
     &   mpi_interior_real_global

#endif

#ifdef coupled  ! this option is specifically for routine setup_mpi
      include "mpif.h"              ! MPI library definitions
#else
#ifdef _MPI
      include "mpif.h"              ! MPI library definitions
#endif
#endif

c=======================================================================

      contains

c=======================================================================

      subroutine ice_global_real_minmax(nc,work,string)

!---!-------------------------------------------------------------------
!---! determines and writes both minimum and maximum over global grid 
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: nc
      real (kind=dbl_kind), intent(in) :: work(nc)
      character (len=8), intent(in) :: string
      real (kind=dbl_kind) :: amin,amax

      amin = ice_global_real_minval(nc,work)
      amax = ice_global_real_maxval(nc,work)
      if (my_task.eq.master_task) 
     &    write (6,*) ' min/max ',string,amin,amax

      end subroutine ice_global_real_minmax

c=======================================================================

      real function ice_global_real_minval(nc,work)

!---!-------------------------------------------------------------------
!---! computes minimum over the global grid
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: nc
      integer (kind=int_kind) :: n
      real (kind=dbl_kind), intent(in) :: work(nc)
      real (kind=dbl_kind) :: local_val

      local_val = minval(work)
#ifdef _MPI
      call MPI_ALLREDUCE(local_val, ice_global_real_minval, 1,
     &                   MPI_REAL8, MPI_MIN, MPI_COMM_ICE, ierr)
#else
      ice_global_real_minval=local_val
#endif

      end function ice_global_real_minval

c=======================================================================

      real function ice_global_real_maxval(nc,work)

!---!-------------------------------------------------------------------
!---! computes maximum over the global grid
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: nc
      integer (kind=int_kind) :: n
      real (kind=dbl_kind), intent(in) :: work(nc)
      real (kind=dbl_kind) :: local_val

      local_val = maxval(work)
#ifdef _MPI
      call MPI_ALLREDUCE(local_val, ice_global_real_maxval, 1,
     &                   MPI_REAL8, MPI_MAX, MPI_COMM_ICE, ierr)
#else
      ice_global_real_maxval=local_val
#endif

      end function ice_global_real_maxval

c=======================================================================

      real function ice_global_real_sum(nc,work)

!---!-------------------------------------------------------------------
!---! sums given array over the global grid
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: nc
      real (kind=dbl_kind), intent(in) :: work(nc)
      real (kind=dbl_kind) :: local_val

      local_val = sum(work)
#ifdef _MPI
      call MPI_ALLREDUCE(local_val, ice_global_real_sum, 1,
     &                   MPI_REAL8, MPI_SUM, MPI_COMM_ICE, ierr)
#else
      ice_global_real_sum=local_val
#endif

      end function ice_global_real_sum

c=======================================================================

      subroutine ice_bcast_rscalar(val)

!---!-------------------------------------------------------------------
!---! broadcasts real scalar value to all processors
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) :: val
#ifdef _MPI
      call MPI_BCAST(val,1,MPI_REAL8,master_task,MPI_COMM_ICE,ierr)
#endif

      end subroutine ice_bcast_rscalar

c=======================================================================

      subroutine ice_bcast_iscalar(ival)

!---!-------------------------------------------------------------------
!---! broadcasts integer scalar value to all processors
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: ival

#ifdef _MPI
      call MPI_BCAST(ival,1,MPI_INTEGER,master_task,MPI_COMM_ICE,ierr)
#endif

      end subroutine ice_bcast_iscalar

c=======================================================================

      subroutine global_scatter(workg,work)

!---!-------------------------------------------------------------------
!---! scatters a global array and adds ghost cells to create a 
!---! distributed array
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_constants

      integer (kind=int_kind) :: i,j
      real (kind=dbl_kind) :: workg(imt_global,jmt_global)
     &,   work(ilo:ihi,jlo:jhi)

#ifdef _MPI

      integer (kind=int_kind) :: request(2),n
     &,   status(MPI_STATUS_SIZE)
      work=c0
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
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
        work(i,j)=workg(i-num_ghost_cells,j-num_ghost_cells)
      enddo
      enddo
#endif

      end subroutine global_scatter

c=======================================================================

      subroutine global_gather(workg,work)

!---!-------------------------------------------------------------------
!---! gathers a distributed array and strips off ghost
!---! cells to create a local array with global dimensions
!---!
!---! author Tony Craig, NCAR
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_constants

      integer (kind=int_kind) :: i,j
      real (kind=dbl_kind) :: workg(imt_global,jmt_global)
     &,   work(ilo:ihi,jlo:jhi)

#ifdef _MPI

      integer (kind=int_kind) :: request(2),recv_req(nproc_s),n
     &, status(MPI_STATUS_SIZE),recv_status(MPI_STATUS_SIZE,nproc_s)

      if (my_task .eq. master_task) then
        workg=c0
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
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
        workg(i-num_ghost_cells,j-num_ghost_cells)=work(i,j)
      enddo
      enddo
#endif

      end subroutine global_gather

c=======================================================================

      subroutine get_sum(flag,work1,work2,work3,gsum)

!---!-------------------------------------------------------------------
!---! computes a (weighted) sum over the global grid
!---! if flag = 1 then work1 is weighted by work2 before being
!---! added to work3
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: flag

      real (kind=dbl_kind), intent(in) ::
     &   work1(ilo:ihi,jlo:jhi)
     &,  work2(ilo:ihi,jlo:jhi)

      real (kind=dbl_kind), intent(in), dimension(:,:) ::
     &  work3

      real (kind=dbl_kind), intent(out) ::
     &   gsum

      real (kind=dbl_kind) ::
     &   work (ilo:ihi,jlo:jhi)

      integer (kind=int_kind) :: i,j
 
      do j=jlo,jhi
       do i=ilo,ihi
        if (flag.eq.1) then
          work(i,j) = work1(i,j) * work2(i,j)
        else
          work(i,j) = work1(i,j)
        endif
        work(i,j) = work(i,j) * work3(i,j)
       enddo
      enddo
      gsum = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),work)

      end subroutine get_sum

c=======================================================================

      subroutine end_run

#ifdef _MPI
      call MPI_FINALIZE(ierr)
#endif

      end subroutine end_run

c=======================================================================

      end module ice_mpi_internal

c=======================================================================
