! CVS: $Id: ice_mpi_internal.F,v 1.7 2004/02/25 17:35:14 eclare Exp $
! CVS: $Source: /home/climate/CVS-COSIM/cice/source/ice_mpi_internal.F,v $
! CVS: $Name:  $
!=======================================================================
!BOP
!
! !MODULE: ice_mpi_internal - parameters and common blocks for MPI internal to ice model
!
! !DESCRIPTION:
!
! Parameters and commons blocks for MPI parallelization 
! internal to ice model
!
! !REVISION HISTORY:
!
! authors: Tony Craig, NCAR
!          Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      module ice_mpi_internal
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain
!
!EOP
!
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

      integer (kind=int_kind), dimension(2,nproc_s) ::
     &   local_start                ! global i,j start for each processor

#ifdef _MPI

      integer (kind=int_kind) ::
     &   nbr_west, nbr_east,        ! processor numbers for the 
     &   nbr_north, nbr_south       ! n,s,e,w neighbor processors

      integer (kind=int_kind), dimension(nproc_s) ::
     &   mpi_interior_int,          ! MPI data type definitions for
     &   mpi_interior_real,         ! physical domain blocks
     &   mpi_interior_int_global,   ! in case of padding
     &   mpi_interior_real_global

#endif

#ifdef coupled  
      ! this option is specifically for routine setup_mpi
      include "mpif.h"              ! MPI library definitions
#else
#ifdef _MPI
      include "mpif.h"              ! MPI library definitions
#endif
#endif

!=======================================================================
 
      contains

!=======================================================================
!BOP
!
! !IROUTINE: ice_global_real_minmax - computes global min/max and prints 
!
! !INTERFACE:
!
      subroutine ice_global_real_minmax(nc,work,string)
!
! !DESCRIPTION:
!
! Determines and writes both minimum and maximum over global grid 
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
      use ice_fileunits
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: nc
      real (kind=dbl_kind), intent(in) :: work(nc)
      character (len=8), intent(in) :: string
!
!EOP
!
      real (kind=dbl_kind) :: amin,amax

      amin = ice_global_real_minval(nc,work)
      amax = ice_global_real_maxval(nc,work)
      if (my_task.eq.master_task) 
     &    write (nu_diag,*) ' min/max ',string,amin,amax

      end subroutine ice_global_real_minmax

!=======================================================================
!BOP
!
! !IROUTINE: real function ice_global_real_minval - computes global min
!
! !INTERFACE:
!
      real (kind=dbl_kind) function ice_global_real_minval(nc,work)
!
! !DESCRIPTION:
!
! Computes minimum over the global grid
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: nc
      real (kind=dbl_kind), intent(in) :: work(nc)
!
!EOP
!
      integer (kind=int_kind) :: n

      real (kind=dbl_kind) :: local_val

      local_val = minval(work)
#ifdef _MPI
      call MPI_ALLREDUCE(local_val, ice_global_real_minval, 1,
     &                   MPI_REAL8, MPI_MIN, MPI_COMM_ICE, ierr)
#else
      ice_global_real_minval=local_val
#endif

      end function ice_global_real_minval
 
!=======================================================================
!BOP
!
! !IROUTINE: real function ice_global_real_maxval - computes global max 
!
! !INTERFACE:
!
      real (kind=dbl_kind) function ice_global_real_maxval(nc,work)
!
! !DESCRIPTION:
!
! Computes maximum over the global grid
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: nc
      real (kind=dbl_kind), intent(in) :: work(nc)
!
!EOP
!
      integer (kind=int_kind) :: n
      real (kind=dbl_kind) :: local_val

      local_val = maxval(work)
#ifdef _MPI
      call MPI_ALLREDUCE(local_val, ice_global_real_maxval, 1,
     &                   MPI_REAL8, MPI_MAX, MPI_COMM_ICE, ierr)
#else
      ice_global_real_maxval=local_val
#endif

      end function ice_global_real_maxval

!=======================================================================
!BOP
!
! !IROUTINE: real function ice_global_real_sum - sums input array over global grid
!
! !INTERFACE:
!
      real (kind=dbl_kind) function ice_global_real_sum(nc,work)
!
! !DESCRIPTION:
!
! Sums given array over the global grid
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: nc
      real (kind=dbl_kind), intent(in) :: work(nc)
!
!EOP
!
      real (kind=dbl_kind) :: local_val

      local_val = sum(work)
#ifdef _MPI
      call MPI_ALLREDUCE(local_val, ice_global_real_sum, 1,
     &                   MPI_REAL8, MPI_SUM, MPI_COMM_ICE, ierr)
#else
      ice_global_real_sum=local_val
#endif

      end function ice_global_real_sum

!=======================================================================
!BOP
!
! !IROUTINE: ice_bcast_logical - broadcasts scalar logical to all processors
!
! !INTERFACE:
!
      subroutine ice_bcast_logical(logval)
!
! !DESCRIPTION:
!
! Broadcasts a scalar logical value to all processors
!
! !REVISION HISTORY:
!
! author: Julie Schramm, NCAR
!         Added 22 Oct 2001 
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      implicit none

      logical (kind=log_kind), intent(inout) :: logval
!
!EOP
!
      integer (kind=int_kind) :: reduce_int   ! Local temporary

      if (logval) then
        reduce_int = 1
      else
        reduce_int = 0
      endif

#ifdef _MPI
      call MPI_BCAST(reduce_int,1,MPI_INTEGER,master_task,
     &               MPI_COMM_ICE,ierr)
#endif
      if (reduce_int .eq. 1) then
        logval = .true.
      else
        logval = .false.
      endif

      end subroutine ice_bcast_logical

!=======================================================================
!BOP
!
! !IROUTINE: ice_bcast_char - broadcasts scalar character to all processors
!
! !INTERFACE:
!
      subroutine ice_bcast_char(charval)
!
! !DESCRIPTION:
!
! Broadcasts a scalar character value to all processors
!
! !REVISION HISTORY:
!
! author: Julie Schramm, NCAR
!         Added 22 Oct 2001 
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      implicit none
      character (*), intent(inout) :: charval
!
!EOP
!
      integer (kind=int_kind) :: clength    ! length of character string

      clength = len(charval)

#ifdef _MPI
      call MPI_BCAST(charval,clength,MPI_CHARACTER,master_task,
     &                       MPI_COMM_ICE,ierr)
#endif

      end subroutine ice_bcast_char

!=======================================================================
!BOP
!
! !IROUTINE: ice_bcast_rscalar - broadcasts real scalar all processors
!
! !INTERFACE:
!
      subroutine ice_bcast_rscalar(val)
!
! !DESCRIPTION:
!
! Broadcasts a real scalar character value to all processors
!
! !REVISION HISTORY:
!
! author: Julie Schramm, NCAR
!         Added 22 Oct 2001 
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(inout) :: val
!
!EOP
!
#ifdef _MPI
      call MPI_BCAST(val,1,MPI_REAL8,master_task,MPI_COMM_ICE,ierr)
#endif

      end subroutine ice_bcast_rscalar

!=======================================================================
!BOP
!
! !IROUTINE: ice_bcast_iscalar - broadcasts integer scalar all processors
!
! !INTERFACE:
!
      subroutine ice_bcast_iscalar(ival)
!
! !DESCRIPTION:
!
! Broadcasts an integer scalar character value to all processors
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(inout) :: ival
!
!EOP
!
#ifdef _MPI
      call MPI_BCAST(ival,1,MPI_INTEGER,master_task,MPI_COMM_ICE,ierr)
#endif

      end subroutine ice_bcast_iscalar

!=======================================================================
!BOP
!
! !IROUTINE: global_scatter - scatters global to distributed array, adding ghost cells
!
! !INTERFACE: 
!
      subroutine global_scatter(workg,work)
!
! !DESCRIPTION:
!
! Scatters a global array and adds ghost cells to create a 
! distributed array
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
      use ice_model_size
      use ice_constants
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind) :: workg(imt_global,jmt_global)
     &,                       work(ilo:ihi,jlo:jhi)
!
!EOP
!
      integer (kind=int_kind) :: i,j

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
      do j=jlo,jhi
      do i=ilo,ihi
        work(i,j)=workg(i-num_ghost_cells,j-num_ghost_cells)
      enddo
      enddo
#endif

      end subroutine global_scatter

!=======================================================================
!BOP
!
! !IROUTINE: global_gather - gathers distributed array into global array
!
! !INTERFACE: 
!
      subroutine global_gather(workg,work)
!
! !DESCRIPTION:
!
! Gathers a distributed array and strips off ghost
! cells to create a local array with global dimensions
!
! !REVISION HISTORY:
!
! author: Tony Craig, NCAR
!
! !USES:
!
      use ice_model_size
      use ice_constants
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind) :: workg(imt_global,jmt_global)
     &,                       work(ilo:ihi,jlo:jhi)
!
!EOP
!
      integer (kind=int_kind) :: i,j
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

! Barrier added by JLS 3 May 2001 to avoid exhausting memory when
! attempting to run on the gx1v2 grid on 16 CPUs
      call MPI_barrier(MPI_COMM_ICE,ierr)

      call MPI_ISEND(work(ilo,jlo), 1, mpi_interior_real(my_task+1), 
     &         master_task, my_task, MPI_COMM_ICE, request(2), ierr)
      call MPI_WAIT(request(2), status, ierr)

      if (my_task .eq. master_task) then
          call MPI_WAITALL(nproc_s, recv_req, recv_status, ierr)
      endif
#else
      ! shift indices from local domain to global domain
      do j=jlo,jhi
      do i=ilo,ihi
        workg(i-num_ghost_cells,j-num_ghost_cells)=work(i,j)
      enddo
      enddo
#endif

      end subroutine global_gather

!=======================================================================
!BOP
!
! !IROUTINE: get_sum - computes weighted sum over global grid
!
! !INTERFACE: 
!
      subroutine get_sum(flag,work1,work2,work3,gsum)
!
! !DESCRIPTION:
!
! Computes a (weighted) sum over the global grid;
! if flag = 1 then work1 is weighted by work2 before being
! added to work3
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: flag

      real (kind=dbl_kind), intent(in) ::
     &   work1(ilo:ihi,jlo:jhi)
     &,  work2(ilo:ihi,jlo:jhi)

      real (kind=dbl_kind), intent(in) ::
     &   work3(imt_local,jmt_local)

      real (kind=dbl_kind), intent(out) ::
     &   gsum
!
!EOP
!
      real (kind=dbl_kind) ::
     &   worka (ilo:ihi,jlo:jhi)

      integer (kind=int_kind) :: i,j

      do j=jlo,jhi
       do i=ilo,ihi
        if (flag.eq.1) then
          worka(i,j) = work1(i,j) * work2(i,j)
        else
          worka(i,j) = work1(i,j)
        endif
        worka(i,j) = worka(i,j) * work3(i,j)
       enddo
      enddo
      gsum = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),worka)

      end subroutine get_sum

!=======================================================================
!BOP
!
! !IROUTINE: end_run - ends run
!
! !INTERFACE:
!
      subroutine end_run
!
! !DESCRIPTION:
!
! Ends run by calling MPI_FINALIZE.
!
! !REVISION HISTORY:
! 
! author: ?
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
#ifdef _MPI
      call MPI_FINALIZE(ierr)
#endif
!
!EOP
!
      end subroutine end_run
 
!=======================================================================

      end module ice_mpi_internal

!=======================================================================
