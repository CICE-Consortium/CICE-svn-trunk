!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!BOP
! !MODULE: ice_global_reductions

 module ice_global_reductions

! !DESCRIPTION:
!  This module contains all the routines for performing global
!  reductions like global sums, minvals, maxvals, etc.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
!
! !USES:

   use ice_kinds_mod
   use ice_communicate
   use ice_constants
   use ice_blocks
   use ice_distribution
   use ice_domain_size
   !use ice_domain   ! commented out because it gives circular dependence
   !use timers

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: global_sum,      &
             global_sum_prod, &
             global_maxval,   &
             global_minval,   &
             init_global_reductions

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  generic interfaces for module procedures
!
!-----------------------------------------------------------------------

   interface global_sum
     module procedure global_sum_dbl,              &
                      global_sum_real,             &
                      global_sum_int,              &
                      global_sum_scalar_dbl,       &
                      global_sum_scalar_real,      &
                      global_sum_scalar_int
   end interface

   interface global_sum_prod
     module procedure global_sum_prod_dbl,         &
                      global_sum_prod_real,        &
                      global_sum_prod_int
   end interface

   interface global_maxval
     module procedure global_maxval_dbl,           &
                      global_maxval_real,          &
                      global_maxval_int,           &
                      global_maxval_scalar_dbl,    &
                      global_maxval_scalar_real,   &
                      global_maxval_scalar_int
   end interface

   interface global_minval
     module procedure global_minval_dbl,           &
                      global_minval_real,          &
                      global_minval_int,           &
                      global_minval_scalar_dbl,    &
                      global_minval_scalar_real,   &
                      global_minval_scalar_int
   end interface

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   !integer (int_kind) :: timer_local, timer_mpi
   logical(log_kind) :: ltripole_grid  ! in lieu of use domain

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_global_reductions
! !INTERFACE:

 subroutine init_global_reductions(tripole_flag)

! !DESCRIPTION:
!  Initializes necessary buffers for global reductions.
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:
!
   logical(log_kind), intent(in) :: tripole_flag
!
!EOP
!BOC

   ltripole_grid = tripole_flag

   !call get_timer(timer_local, 'SUM_LOCAL')
   !call get_timer(timer_mpi  , 'SUM_MPI')

!EOC

 end subroutine init_global_reductions

!***********************************************************************
!BOP
! !IROUTINE: global_sum
! !INTERFACE:

 function global_sum_dbl(X, dist, field_loc, MASK)

! !DESCRIPTION:
!  computes the global sum of the _physical domain_ of a 2-d
!  array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_sum
!  function corresponding to double precision arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

! !USES:

   include 'mpif.h'  ! MPI Fortran include file

! !INPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      X                    ! array to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), dimension(size(X,dim=1), &
                        size(X,dim=2), &
                        size(X,dim=3)), intent(in), optional :: &
      MASK                 ! real multiplicative mask

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
      global_sum_dbl       ! resulting global sum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!   real (dbl_kind), dimension(:), allocatable :: &
!      local_block_sum,    &! sum of local blocks
!      global_block_sum     ! sum of all blocks

   real (dbl_kind) ::          &
      local_sum           ! sum of all local blocks

   integer (int_kind) :: &
      i,j,n,             &! local counters
      ib,ie,jb,je,       &! beg,end of physical domain
      bid,               &! block location
      ierr                ! MPI error flag

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------
!
!  use this code for sums that are not reproducible for best performance
!
!-----------------------------------------------------------------------

   local_sum = c0

   !call timer_start(timer_local)
   if (ltripole_grid .and. (field_loc == field_loc_Nface .or. &
                            field_loc == field_loc_NEcorner)) then
      !*** must remove redundant points from sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            this_block = get_block(n,bid)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            if (this_block%jblock == nblocks_y) then
               !*** for topmost row in tripole only sum
               !*** 1st half of domain - others are redundant
               if (present(MASK)) then
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + X(i,je,bid)*MASK(i,je,bid)
                  end do
               else ! no mask
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + X(i,je,bid)
                  end do
               endif
               je = je - 1
            endif
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   else ! regular global sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   endif
   !call timer_stop(timer_local)

   !call timer_start(timer_mpi)
   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_sum, global_sum_dbl, 1, &
                            mpiR8, MPI_SUM, dist%communicator, ierr)
      else
         global_sum_dbl = c0
      endif
   else
      global_sum_dbl = local_sum
   endif
   !call timer_stop(timer_mpi)

!-----------------------------------------------------------------------
!
!  use this code for sums that are more reproducible - performance
!  of this code does not scale well for large numbers of blocks
!  NOTE: THIS CODE NOT MODIFIED FOR TRIPOLE GRIDS YET
!
!-----------------------------------------------------------------------
!
!   allocate (local_block_sum(nblocks_tot), &
!             global_block_sum(nblocks_tot))
!
!   local_block_sum = c0
!
!   !call timer_start(timer_local)
!   if (present(MASK)) then
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)*MASK(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   else
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   endif
!   !call timer_stop(timer_local)
!
!   !call timer_start(timer_mpi)
!   call MPI_ALLREDUCE(local_block_sum, global_block_sum, nblocks_tot, &
!                      mpiR8, MPI_SUM, dist%communicator, ierr)
!   !call timer_stop(timer_mpi)
!
!   global_sum_dbl = c0
!   do n=1,nblocks_tot
!     global_sum_dbl = global_sum_dbl + global_block_sum(n)
!   end do
!
!   deallocate ( local_block_sum, global_block_sum)
!
!-----------------------------------------------------------------------

 end function global_sum_dbl

!***********************************************************************

 function global_sum_real(X, dist, field_loc, MASK)

!-----------------------------------------------------------------------
!
!  computes the global sum of the _physical domain_ of a 2-d
!  array.
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

!-----------------------------------------------------------------------
!
!  input vars
!
!-----------------------------------------------------------------------

   real (real_kind), dimension(:,:,:), intent(in) :: &
      X                    ! array to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), dimension(size(X,dim=1), &
                              size(X,dim=2), &
                              size(X,dim=3)), intent(in), optional :: &
      MASK                 ! real multiplicative mask

!-----------------------------------------------------------------------
!
!  output result
!
!-----------------------------------------------------------------------

   real (real_kind) :: &
      global_sum_real       ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!   real (real_kind), dimension(:), allocatable :: &
!      local_block_sum,    &! sum of local blocks
!      global_block_sum     ! sum of all blocks

   real (real_kind) ::          &
      local_sum           ! sum of local blocks

   integer (int_kind) :: &
      i,j,n,             &! local counters
      ib,ie,jb,je,       &! beg,end of physical domain
      bid,               &! block location
      ierr                ! MPI error flag

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------
!
!  use this code for sums that perform better but are not reproducible
!
!-----------------------------------------------------------------------

   local_sum = c0

   !call timer_start(timer_local)
   if (ltripole_grid .and. (field_loc == field_loc_Nface .or. &
                            field_loc == field_loc_NEcorner)) then
      !*** must remove redundant points from sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            this_block = get_block(n,bid)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            if (this_block%jblock == nblocks_y) then
               !*** for topmost row in tripole only sum
               !*** 1st half of domain - others are redundant
               if (present(MASK)) then
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + X(i,je,bid)*MASK(i,je,bid)
                  end do
               else ! no mask
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + X(i,je,bid)
                  end do
               endif
               je = je - 1
            endif
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   else ! regular global sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   endif
   !call timer_stop(timer_local)

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_sum, global_sum_real, 1, &
                            mpiR4, MPI_SUM, dist%communicator, ierr)
      else
         global_sum_real = c0
      endif
   else
      global_sum_real = local_sum
   endif

!-----------------------------------------------------------------------
!
!  use this code for sums that are more reproducible - performance
!  of this code does not scale well for large numbers of blocks
!  NOTE: THIS CODE NOT MODIFIED FOR TRIPOLE GRIDS YET
!
!-----------------------------------------------------------------------
!
!   allocate (local_block_sum(nblocks_tot), &
!             global_block_sum(nblocks_tot))
!
!   local_block_sum = c0
!
!   if (present(MASK)) then
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)*MASK(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   else
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   endif
!
!   call MPI_ALLREDUCE(local_block_sum, global_block_sum, nblocks_tot, &
!                      mpiR4, MPI_SUM, dist%communicator, ierr)
!
!   global_sum_real = c0
!   do n=1,nblocks_tot
!     global_sum_real = global_sum_real + global_block_sum(n)
!   end do
!
!   deallocate ( local_block_sum, global_block_sum)
!
!-----------------------------------------------------------------------

 end function global_sum_real

!***********************************************************************

 function global_sum_int(X, dist, field_loc, MASK)

!-----------------------------------------------------------------------
!
!  computes the global sum of the _physical domain_ of a 2-d
!  array.
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

!-----------------------------------------------------------------------
!
!  input vars
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      X                    ! array to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), dimension(size(X,dim=1), &
                              size(X,dim=2), &
                              size(X,dim=3)), intent(in), optional :: &
      MASK                 ! real multiplicative mask

!-----------------------------------------------------------------------
!
!  output result
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      global_sum_int       ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!   integer (int_kind), dimension(:), allocatable :: &
!      local_block_sum,    &! sum of local blocks
!      global_block_sum     ! sum of all blocks

   integer (int_kind) :: &
      local_sum           ! sum of local blocks

   integer (int_kind) :: &
      i,j,n,             &! local counters
      ib,ie,jb,je,       &! beg,end of physical domain
      bid,               &! block location
      ierr                ! MPI error flag

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------
!
!  use this code for sums that are not reproducible but perform better
!
!-----------------------------------------------------------------------

   local_sum = c0

   !call timer_start(timer_local)
   if (ltripole_grid .and. (field_loc == field_loc_Nface .or. &
                            field_loc == field_loc_NEcorner)) then
      !*** must remove redundant points from sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            this_block = get_block(n,bid)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            if (this_block%jblock == nblocks_y) then
               !*** for topmost row in tripole only sum
               !*** 1st half of domain - others are redundant
               if (present(MASK)) then
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + X(i,je,bid)*MASK(i,je,bid)
                  end do
               else ! no mask
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + X(i,je,bid)
                  end do
               endif
               je = je - 1
            endif
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   else ! regular global sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   endif
   !call timer_stop(timer_local)

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_sum, global_sum_int, 1, &
                       mpi_integer, MPI_SUM, dist%communicator, ierr)
      else
         global_sum_int = 0
      endif
   else
      global_sum_int = local_sum
   endif

!-----------------------------------------------------------------------
!
!  use this code for sums that are more reproducible - performance
!  of this code does not scale well for large numbers of blocks
!  NOTE: THIS CODE NOT YET CHANGED FOR TRIPOLE GRID
!
!-----------------------------------------------------------------------
!
!   allocate (local_block_sum(nblocks_tot), &
!             global_block_sum(nblocks_tot))
!
!   local_block_sum = c0
!
!   if (present(MASK)) then
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)*MASK(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   else
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   endif
!
!   call MPI_ALLREDUCE(local_block_sum, global_block_sum, nblocks_tot, &
!                      mpi_integer, MPI_SUM, dist%communicator, ierr)
!
!   global_sum_int = 0
!   do n=1,nblocks_tot
!     global_sum_int = global_sum_int + global_block_sum(n)
!   end do
!
!   deallocate ( local_block_sum, global_block_sum)
!
!-----------------------------------------------------------------------

 end function global_sum_int

!***********************************************************************

 function global_sum_scalar_dbl(local_scalar, dist)

!-----------------------------------------------------------------------
!
!  this function returns the sum of scalar value across processors
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

   type (distrb), intent(in) :: &
      dist                 ! distribution from which this is called

   real (dbl_kind), intent(inout) :: &
      local_scalar                ! local scalar to be compared

   real (dbl_kind) :: &
      global_sum_scalar_dbl   ! resulting global sum

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_scalar, global_sum_scalar_dbl, 1, &
                            mpiR8, MPI_SUM, dist%communicator, ierr)
      else
         global_sum_scalar_dbl = c0
      endif
   else
      global_sum_scalar_dbl = local_scalar
   endif

!-----------------------------------------------------------------------

 end function global_sum_scalar_dbl

!***********************************************************************

 function global_sum_scalar_real(local_scalar, dist)

!-----------------------------------------------------------------------
!
!  this function returns the sum of scalar value across processors
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

   real (real_kind), intent(inout) :: &
      local_scalar                ! local scalar to be compared

   type (distrb), intent(in) :: &
      dist                 ! distribution from which this is called

   real (real_kind) :: &
      global_sum_scalar_real   ! resulting global sum

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_scalar, global_sum_scalar_real, 1, &
                            mpiR4, MPI_SUM, dist%communicator, ierr)
      else
         global_sum_scalar_real = c0
      endif
   else
      global_sum_scalar_real = local_scalar
   endif

!-----------------------------------------------------------------------

 end function global_sum_scalar_real

!***********************************************************************

 function global_sum_scalar_int(local_scalar, dist)

!-----------------------------------------------------------------------
!
!  this function returns the sum of scalar value across processors
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(inout) :: &
      local_scalar                ! local scalar to be compared

   type (distrb), intent(in) :: &
      dist                 ! distribution from which this is called

   integer (int_kind) :: &
      global_sum_scalar_int   ! resulting global sum

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_scalar, global_sum_scalar_int, 1, &
                         MPI_INTEGER, MPI_SUM, dist%communicator, ierr)
      else
         global_sum_scalar_int = 0
      endif
   else
      global_sum_scalar_int = local_scalar
   endif

!-----------------------------------------------------------------------

 end function global_sum_scalar_int

!EOC
!***********************************************************************
!BOP
! !IROUTINE: global_sum_prod
! !INTERFACE:

 function global_sum_prod_dbl (X,Y,dist,field_loc, MASK)

! !DESCRIPTION:
!  this routine performs a global sum over the physical domain
!  of a product of two 2-d arrays.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic 
!  global_sum_prod function corresponding to double precision arrays.
!  The generic interface is identical but will handle real and integer 
!  2-d slabs.

! !USES:

   include 'mpif.h'  ! MPI Fortran include file

! !INPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
     X,                &! first array in product to be summed
     Y                  ! second array in product to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X,Y

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), &
     dimension(size(X,dim=1),size(X,dim=2),size(X,dim=3)), &
     intent(in), optional :: &
     MASK               ! real multiplicative mask

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
     global_sum_prod_dbl ! resulting global sum of X*Y

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!   real (dbl_kind), dimension(:), allocatable :: &
!     local_block_sum,  &! sum of each block
!     global_block_sum   ! global sum each block

   real (dbl_kind) ::         & 
     local_sum           ! sum of each block

   integer (int_kind) :: &
      i,j,n,            &! dummy counters
      bid,              &! local block location
      ib,ie,jb,je,      &! beg,end indices of physical domain
      ierr               ! MPI error flag

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------
!
!  use this code for sums that are not reproducible but perform better
!
!-----------------------------------------------------------------------

   local_sum = c0

   !call timer_start(timer_local)
   if (ltripole_grid .and. (field_loc == field_loc_Nface .or. &
                            field_loc == field_loc_NEcorner)) then
      !*** must remove redundant points from sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            this_block = get_block(n,bid)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            if (this_block%jblock == nblocks_y) then
               !*** for topmost row in tripole only sum
               !*** 1st half of domain - others are redundant
               if (present(MASK)) then
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + &
                                 X(i,je,bid)*Y(i,je,bid)*MASK(i,je,bid)
                  end do
               else ! no mask
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + X(i,je,bid)*Y(i,je,bid)
                  end do
               endif
               je = je - 1
            endif
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + &
                              X(i,j,bid)*Y(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*Y(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   else ! regular global sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + &
                              X(i,j,bid)*Y(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*Y(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   endif
   !call timer_stop(timer_local)

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_sum, global_sum_prod_dbl, 1, &
                      mpiR8, MPI_SUM, dist%communicator, ierr)
      else
         global_sum_prod_dbl = 0
      endif
   else
      global_sum_prod_dbl = local_sum
   endif

!-----------------------------------------------------------------------
!
!  use this code for sums that are more reproducible - performance
!  of this code does not scale well for large numbers of blocks
!  NOTE: THIS CODE NOT YET CHANGED FOR TRIPOLE GRIDS
!
!-----------------------------------------------------------------------
!
!   allocate (local_block_sum(nblocks_tot), &
!             global_block_sum(nblocks_tot))
!
!   local_block_sum = c0
!
!   if (present(MASK)) then
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)*Y(i,j,bid)*MASK(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   else
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)*Y(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   endif
!
!   call MPI_ALLREDUCE(local_block_sum, global_block_sum, nblocks_tot, &
!                      mpiR8, MPI_SUM, dist%communicator, ierr)
!
!   global_sum_prod_dbl = c0
!   do n=1,nblocks_tot
!     global_sum_prod_dbl = global_sum_prod_dbl + global_block_sum(n)
!   end do
!
!   deallocate ( local_block_sum, global_block_sum)
!
!-----------------------------------------------------------------------

 end function global_sum_prod_dbl

!***********************************************************************

 function global_sum_prod_real (X, Y, dist, field_loc, MASK)

!-----------------------------------------------------------------------
!
!  this routine performs a global sum over the physical domain
!  of a product of two 2-d arrays.
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   real (real_kind), dimension(:,:,:), intent(in) :: &
     X,                &! first array in product to be summed
     Y                  ! second array in product to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X,Y

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), &
     dimension(size(X,dim=1),size(X,dim=2),size(X,dim=3)), &
     intent(in), optional :: &
     MASK               ! real multiplicative mask

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (real_kind) :: &
     global_sum_prod_real ! resulting global sum of X*Y

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!   real (dbl_kind), dimension(:), allocatable :: &
!     local_block_sum,  &! sum of each block
!     global_block_sum   ! global sum each block

   real (dbl_kind) ::         &
     local_sum,         &! sum of local blocks
     global_sum_prod_tmp ! sum of global blocks

   integer (int_kind) :: &
      i,j,n,            &! dummy counters
      bid,              &! local block location
      ib,ie,jb,je,      &! beg,end indices of physical domain
      ierr               ! MPI error flag

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------
!
!  use this code for sums that are not reproducible but perform better
!
!-----------------------------------------------------------------------

   local_sum = c0

   !call timer_start(timer_local)
   if (ltripole_grid .and. (field_loc == field_loc_Nface .or. &
                            field_loc == field_loc_NEcorner)) then
      !*** must remove redundant points from sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            this_block = get_block(n,bid)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            if (this_block%jblock == nblocks_y) then
               !*** for topmost row in tripole only sum
               !*** 1st half of domain - others are redundant
               if (present(MASK)) then
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + &
                                 X(i,je,bid)*Y(i,je,bid)*MASK(i,je,bid)
                  end do
               else ! no mask
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + X(i,je,bid)*Y(i,je,bid)
                  end do
               endif
               je = je - 1
            endif
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + &
                              X(i,j,bid)*Y(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*Y(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   else ! regular global sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + &
                              X(i,j,bid)*Y(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*Y(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   endif
   !call timer_stop(timer_local)

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_sum, global_sum_prod_tmp, 1, &
                      mpiR8, MPI_SUM, dist%communicator, ierr)
         global_sum_prod_real = global_sum_prod_tmp
      else
         global_sum_prod_real = 0
      endif
   else
      global_sum_prod_real = local_sum
   endif

!-----------------------------------------------------------------------
!
!  use this code for sums that are more reproducible - performance
!  of this code does not scale well for large numbers of blocks
!  NOTE: THIS CODE NOT YET CHANGED FOR TRIPOLE GRIDS
!
!-----------------------------------------------------------------------
!
!
!   allocate (local_block_sum(nblocks_tot), &
!             global_block_sum(nblocks_tot))
!
!   local_block_sum = c0
!
!   if (present(MASK)) then
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)*Y(i,j,bid)*MASK(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   else
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)*Y(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   endif
!
!   call MPI_ALLREDUCE(local_block_sum, global_block_sum, nblocks_tot, &
!                      mpiR8, MPI_SUM, dist%communicator, ierr)
!
!   global_sum_prod_real = c0
!   do n=1,nblocks_tot
!     global_sum_prod_real = global_sum_prod_real + global_block_sum(n)
!   end do
!
!   deallocate ( local_block_sum, global_block_sum)
!
!-----------------------------------------------------------------------

 end function global_sum_prod_real

!***********************************************************************

 function global_sum_prod_int (X, Y, dist, field_loc, MASK)

!-----------------------------------------------------------------------
!
!  this routine performs a global sum over the physical domain
!  of a product of two 2-d arrays.
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), intent(in) :: &
     X,                &! first array in product to be summed
     Y                  ! second array in product to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X,Y

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), &
     dimension(size(X,dim=1),size(X,dim=2),size(X,dim=3)), &
     intent(in), optional :: &
     MASK               ! real multiplicative mask

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     global_sum_prod_int ! resulting global sum of X*Y

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!   integer (int_kind), dimension(:), allocatable :: &
!     local_block_sum,  &! sum of each block
!     global_block_sum   ! global sum each block

   integer (int_kind) :: &
     local_sum           ! sum of local blocks

   integer (int_kind) :: &
      i,j,n,            &! dummy counters
      bid,              &! local block location
      ib,ie,jb,je,      &! beg,end indices of physical domain
      ierr               ! MPI error flag

   type (block) :: &
      this_block          ! holds local block information

!-----------------------------------------------------------------------
!
!  use this code for sums that are not reproducible but perform better
!
!-----------------------------------------------------------------------

   local_sum = c0

   !call timer_start(timer_local)
   if (ltripole_grid .and. (field_loc == field_loc_Nface .or. &
                            field_loc == field_loc_NEcorner)) then
      !*** must remove redundant points from sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            this_block = get_block(n,bid)
            ib = this_block%ilo
            ie = this_block%ihi
            jb = this_block%jlo
            je = this_block%jhi
            if (this_block%jblock == nblocks_y) then
               !*** for topmost row in tripole only sum
               !*** 1st half of domain - others are redundant
               if (present(MASK)) then
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + &
                                 X(i,je,bid)*Y(i,je,bid)*MASK(i,je,bid)
                  end do
               else ! no mask
                  do i=ib,ie
                     if (this_block%i_glob(i) <= nx_global/2) &
                     local_sum = local_sum + X(i,je,bid)*Y(i,je,bid)
                  end do
               endif
               je = je - 1
            endif
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + &
                              X(i,j,bid)*Y(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*Y(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   else ! regular global sum
      do n=1,nblocks_tot
         if (dist%blockLocation(n) == my_task+1) then
            bid = dist%blockLocalID(n)
            call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
            if (present(MASK)) then
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + &
                              X(i,j,bid)*Y(i,j,bid)*MASK(i,j,bid)
               end do
               end do
            else ! no mask
               do j=jb,je
               do i=ib,ie
                  local_sum = local_sum + X(i,j,bid)*Y(i,j,bid)
               end do
               end do
            endif
         endif
      end do !block loop
   endif
   !call timer_stop(timer_local)

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_sum, global_sum_prod_int, 1, &
                      mpi_integer, MPI_SUM, dist%communicator, ierr)
      else
         global_sum_prod_int = 0
      endif
   else
      global_sum_prod_int = local_sum
   endif

!-----------------------------------------------------------------------
!
!  use this code for sums that are more reproducible - performance
!  of this code does not scale well for large numbers of blocks
!  NOTE: THIS CODE NOT YET CHANGED FOR TRIPOLE GRIDS
!
!-----------------------------------------------------------------------
!
!   allocate (local_block_sum(nblocks_tot), &
!             global_block_sum(nblocks_tot))
!
!   local_block_sum = c0
!
!   if (present(MASK)) then
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)*Y(i,j,bid)*MASK(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   else
!     do n=1,nblocks_tot
!       if (dist%blockLocation(n) == my_task+1) then
!         bid = dist%blockLocalID(n)
!         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
!         do j=jb,je
!         do i=ib,ie
!           local_block_sum(n) = &
!           local_block_sum(n) + X(i,j,bid)*Y(i,j,bid)
!         end do
!         end do
!       endif
!     end do
!   endif
!
!   call MPI_ALLREDUCE(local_block_sum, global_block_sum, nblocks_tot, &
!                      mpi_integer, MPI_SUM, dist%communicator, ierr)
!
!   global_sum_prod_int = 0
!   do n=1,nblocks_tot
!     global_sum_prod_int = global_sum_prod_int + global_block_sum(n)
!   end do
!
!   deallocate ( local_block_sum, global_block_sum)
!
!-----------------------------------------------------------------------

 end function global_sum_prod_int

!EOC
!***********************************************************************
!BOP
! !IROUTINE: global_maxval
! !INTERFACE:

 function global_maxval_dbl (X, dist, field_loc, LMASK)

! !DESCRIPTION:
!  This function computes the global maxval of the physical domain
!  of a 2-d field
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to double precision arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs.

! !USES:

   include 'mpif.h'  ! MPI Fortran include file

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  input vars
!
!-----------------------------------------------------------------------

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      X            ! array containing field for which max required

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   logical (log_kind), &
     dimension(size(X,dim=1),size(X,dim=2),size(X,dim=3)), &
     intent(in), optional :: &
      LMASK                ! mask for excluding parts of domain

!-----------------------------------------------------------------------
!
!  output result
!
!-----------------------------------------------------------------------

   real (dbl_kind) :: &
      global_maxval_dbl   ! resulting max val of global domain

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind) :: &
      local_maxval         ! max value of local subdomain

   integer (int_kind) :: &
      i,j,n,bid, &
      ib, ie, jb, je,     &! start,end of physical domain
      ierr                 ! MPI error flag

!-----------------------------------------------------------------------

   local_maxval = -bignum
   do n=1,nblocks_tot
      if (dist%blockLocation(n) == my_task+1) then
         bid = dist%blockLocalID(n)
         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
         if (present(LMASK)) then
            do j=jb,je
            do i=ib,ie
               if (LMASK(i,j,bid)) &
                  local_maxval = max(local_maxval,X(i,j,bid))
            end do
            end do
         else ! no mask
            do j=jb,je
            do i=ib,ie
               local_maxval = max(local_maxval,X(i,j,bid))
            end do
            end do
         endif
      endif
   end do ! block loop

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_maxval, global_maxval_dbl, 1, &
                            mpiR8, MPI_MAX, dist%communicator, ierr)
      else
         global_maxval_dbl = c0
      endif
   else
      global_maxval_dbl = local_maxval
   endif

!-----------------------------------------------------------------------

 end function global_maxval_dbl

!***********************************************************************

 function global_maxval_real (X, dist, field_loc, LMASK)

!-----------------------------------------------------------------------
!
!  this function computes the global maxval of the physical domain
!  of a 2-d field
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

!-----------------------------------------------------------------------
!
!  input vars
!
!-----------------------------------------------------------------------

   real (real_kind), dimension(:,:,:), intent(in) :: &
      X            ! array containing field for which max required

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   logical (log_kind), &
     dimension(size(X,dim=1),size(X,dim=2),size(X,dim=3)), &
     intent(in), optional :: &
      LMASK                ! mask for excluding parts of domain

!-----------------------------------------------------------------------
!
!  output result
!
!-----------------------------------------------------------------------

   real (real_kind) :: &
      global_maxval_real   ! resulting max val of global domain

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (real_kind) :: &
      local_maxval         ! max value of local subdomain

   integer (int_kind) :: &
      i,j,n,bid, &
      ib, ie, jb, je,     &! start,end of physical domain
      ierr                 ! MPI error flag

!-----------------------------------------------------------------------

   local_maxval = -bignum
   do n=1,nblocks_tot
      if (dist%blockLocation(n) == my_task+1) then
         bid = dist%blockLocalID(n)
         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
         if (present(LMASK)) then
            do j=jb,je
            do i=ib,ie
               if (LMASK(i,j,bid)) &
                  local_maxval = max(local_maxval,X(i,j,bid))
            end do
            end do
         else ! no mask
            do j=jb,je
            do i=ib,ie
               local_maxval = max(local_maxval,X(i,j,bid))
            end do
            end do
         endif
      endif
   end do ! block loop

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_maxval, global_maxval_real, 1, &
                            mpiR4, MPI_MAX, dist%communicator, ierr)
      else
         global_maxval_real = c0
      endif
   else
      global_maxval_real = local_maxval
   endif

!-----------------------------------------------------------------------

 end function global_maxval_real

!***********************************************************************

 function global_maxval_int (X, dist, field_loc, LMASK)

!-----------------------------------------------------------------------
!
!  this function computes the global maxval of the physical domain
!  of a 2-d field
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

!-----------------------------------------------------------------------
!
!  input vars
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      X            ! array containing field for which max required

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   logical (log_kind), &
     dimension(size(X,dim=1),size(X,dim=2),size(X,dim=3)), &
     intent(in), optional :: &
      LMASK                ! mask for excluding parts of domain

!-----------------------------------------------------------------------
!
!  output result
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      global_maxval_int    ! resulting max val of global domain

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      local_maxval         ! max value of local subdomain

   integer (int_kind) :: &
      i,j,n,bid, &
      ib, ie, jb, je,     &! start,end of physical domain
      ierr                 ! MPI error flag

!-----------------------------------------------------------------------

   local_maxval = -1000000
   do n=1,nblocks_tot
      if (dist%blockLocation(n) == my_task+1) then
         bid = dist%blockLocalID(n)
         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
         if (present(LMASK)) then
            do j=jb,je
            do i=ib,ie
               if (LMASK(i,j,bid)) &
                  local_maxval = max(local_maxval,X(i,j,bid))
            end do
            end do
         else ! no mask
            do j=jb,je
            do i=ib,ie
               local_maxval = max(local_maxval,X(i,j,bid))
            end do
            end do
         endif
      endif
   end do ! block loop

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_maxval, global_maxval_int, 1, &
                         mpi_integer, MPI_MAX, dist%communicator, ierr)
      else
         global_maxval_int = 0
      endif
   else
      global_maxval_int = local_maxval
   endif

!-----------------------------------------------------------------------
!EOC

 end function global_maxval_int

!***********************************************************************
!BOP
! !IROUTINE: global_minval
! !INTERFACE:

 function global_minval_dbl (X, dist, field_loc, LMASK)

! !DESCRIPTION:
!  This function computes the global minval of the physical domain
!  of a 2-d field
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is actually the specific interface for the generic global_minval
!  function corresponding to double precision arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs.

! !USES:

   include 'mpif.h'  ! MPI Fortran include file

! !INPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      X            ! array containing field for which min required

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   logical (log_kind), &
     dimension(size(X,dim=1),size(X,dim=2),size(X,dim=3)), &
     intent(in), optional :: &
      LMASK                ! mask for excluding parts of domain

! !OUTPUT PARAMETERS:

   real (dbl_kind) :: &
      global_minval_dbl   ! resulting min val of global domain

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind) :: &
      local_minval         ! min value of local subdomain

   integer (int_kind) :: &
      i,j,n,bid, &
      ib, ie, jb, je,     &! start,end of physical domain
      ierr                 ! MPI error flag

!-----------------------------------------------------------------------

   local_minval = bignum
   do n=1,nblocks_tot
      if (dist%blockLocation(n) == my_task+1) then
         bid = dist%blockLocalID(n)
         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
         if (present(LMASK)) then
            do j=jb,je
            do i=ib,ie
               if (LMASK(i,j,bid)) &
                  local_minval = min(local_minval,X(i,j,bid))
            end do
            end do
         else ! no mask
            do j=jb,je
            do i=ib,ie
               local_minval = min(local_minval,X(i,j,bid))
            end do
            end do
         endif
      endif
   end do

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_minval, global_minval_dbl, 1, &
                            mpiR8, MPI_MIN, dist%communicator, ierr)
      else
         global_minval_dbl = c0
      endif
   else
      global_minval_dbl = local_minval
   endif

!-----------------------------------------------------------------------

 end function global_minval_dbl

!***********************************************************************

 function global_minval_real (X, dist, field_loc, LMASK)

!-----------------------------------------------------------------------
!
!  this function computes the global minval of the physical domain
!  of a 2-d field
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

!-----------------------------------------------------------------------
!
!  input vars
!
!-----------------------------------------------------------------------

   real (real_kind), dimension(:,:,:), intent(in) :: &
      X            ! array containing field for which min required

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   logical (log_kind), &
     dimension(size(X,dim=1),size(X,dim=2),size(X,dim=3)), &
     intent(in), optional :: &
      LMASK                ! mask for excluding parts of domain

!-----------------------------------------------------------------------
!
!  output result
!
!-----------------------------------------------------------------------

   real (real_kind) :: &
      global_minval_real   ! resulting min val of global domain

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (real_kind) :: &
      local_minval         ! min value of local subdomain

   integer (int_kind) :: &
      i,j,n,bid, &
      ib, ie, jb, je,     &! start,end of physical domain
      ierr                 ! MPI error flag

!-----------------------------------------------------------------------

   local_minval = bignum
   do n=1,nblocks_tot
      if (dist%blockLocation(n) == my_task+1) then
         bid = dist%blockLocalID(n)
         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
         if (present(LMASK)) then
            do j=jb,je
            do i=ib,ie
               if (LMASK(i,j,bid)) &
                  local_minval = min(local_minval,X(i,j,bid))
            end do
            end do
         else ! no mask
            do j=jb,je
            do i=ib,ie
               local_minval = min(local_minval,X(i,j,bid))
            end do
            end do
         endif
      endif
   end do ! block loop

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_minval, global_minval_real, 1, &
                            mpiR4, MPI_MIN, dist%communicator, ierr)
      else
         global_minval_real = c0
      endif
   else
      global_minval_real = local_minval
   endif

!-----------------------------------------------------------------------

 end function global_minval_real

!***********************************************************************

 function global_minval_int (X, dist, field_loc, LMASK)

!-----------------------------------------------------------------------
!
!  this function computes the global minval of the physical domain
!  of a 2-d field
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

!-----------------------------------------------------------------------
!
!  input vars
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      X            ! array containing field for which min required

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   logical (log_kind), &
     dimension(size(X,dim=1),size(X,dim=2),size(X,dim=3)), &
     intent(in), optional :: &
      LMASK                ! mask for excluding parts of domain

!-----------------------------------------------------------------------
!
!  output result
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      global_minval_int    ! resulting min val of global domain

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      local_minval         ! min value of local subdomain

   integer (int_kind) :: &
      i,j,n,bid, &
      ib, ie, jb, je,     &! start,end of physical domain
      ierr                 ! MPI error flag

!-----------------------------------------------------------------------

   local_minval = 1000000
   do n=1,nblocks_tot
      if (dist%blockLocation(n) == my_task+1) then
         bid = dist%blockLocalID(n)
         call get_block_parameter(n,ilo=ib,ihi=ie,jlo=jb,jhi=je)
         if (present(LMASK)) then
            do j=jb,je
            do i=ib,ie
               if (LMASK(i,j,bid)) &
                  local_minval = min(local_minval,X(i,j,bid))
            end do
            end do
         else ! no mask
            do j=jb,je
            do i=ib,ie
               local_minval = min(local_minval,X(i,j,bid))
            end do
            end do
         endif
      endif
   end do

   if (dist%nprocs > 1) then
      if (my_task < dist%nprocs) then
         call MPI_ALLREDUCE(local_minval, global_minval_int, 1, &
                         mpi_integer, MPI_MIN, dist%communicator, ierr)
      else
         global_minval_int = 0
      endif
   else
      global_minval_int = local_minval
   endif

!-----------------------------------------------------------------------

 end function global_minval_int

!***********************************************************************

 function global_maxval_scalar_dbl (local_scalar)

!-----------------------------------------------------------------------
!
!  this function returns the maximum scalar value across processors
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

   real (dbl_kind), intent(inout) :: &
      local_scalar                ! local scalar to be compared

   real (dbl_kind) :: &
      global_maxval_scalar_dbl   ! resulting global max

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------

   call MPI_ALLREDUCE(local_scalar, global_maxval_scalar_dbl, 1, &
                      mpiR8, MPI_MAX, MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end function global_maxval_scalar_dbl

!***********************************************************************

 function global_maxval_scalar_real (local_scalar)

!-----------------------------------------------------------------------
!
!  this function returns the maximum scalar value across processors
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

   real (real_kind), intent(inout) :: &
      local_scalar                ! local scalar to be compared

   real (real_kind) :: &
      global_maxval_scalar_real   ! resulting global max

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------

   call MPI_ALLREDUCE(local_scalar, global_maxval_scalar_real, 1, &
                      mpiR4, MPI_MAX, MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end function global_maxval_scalar_real

!***********************************************************************

 function global_maxval_scalar_int (local_scalar)

!-----------------------------------------------------------------------
!
!  this function returns the maximum scalar value across processors
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(inout) :: &
      local_scalar                ! local scalar to be compared

   integer (int_kind) :: &
      global_maxval_scalar_int   ! resulting global max

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------

   call MPI_ALLREDUCE(local_scalar, global_maxval_scalar_int, 1, &
                      MPI_INTEGER, MPI_MAX, MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end function global_maxval_scalar_int

!***********************************************************************

 function global_minval_scalar_dbl (local_scalar)

!-----------------------------------------------------------------------
!
!  this function returns the minimum scalar value across processors
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

   real (dbl_kind), intent(inout) :: &
      local_scalar                ! local scalar to be compared

   real (dbl_kind) :: &
      global_minval_scalar_dbl   ! resulting global min

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------

   call MPI_ALLREDUCE(local_scalar, global_minval_scalar_dbl, 1, &
                      mpiR8, MPI_MIN, MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end function global_minval_scalar_dbl

!***********************************************************************

 function global_minval_scalar_real (local_scalar)

!-----------------------------------------------------------------------
!
!  this function returns the minimum scalar value across processors
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

   real (real_kind), intent(inout) :: &
      local_scalar                ! local scalar to be compared

   real (real_kind) :: &
      global_minval_scalar_real   ! resulting global min

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------

   call MPI_ALLREDUCE(local_scalar, global_minval_scalar_real, 1, &
                      mpiR4, MPI_MIN, MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end function global_minval_scalar_real

!***********************************************************************

 function global_minval_scalar_int (local_scalar)

!-----------------------------------------------------------------------
!
!  this function returns the minimum scalar value across processors
!
!-----------------------------------------------------------------------

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(inout) :: &
      local_scalar                ! local scalar to be compared

   integer (int_kind) :: &
      global_minval_scalar_int   ! resulting global min

   integer (int_kind) :: ierr ! MPI error flag

!-----------------------------------------------------------------------

   call MPI_ALLREDUCE(local_scalar, global_minval_scalar_int, 1, &
                      MPI_INTEGER, MPI_MIN, MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------
!EOC

 end function global_minval_scalar_int

!***********************************************************************

 end module ice_global_reductions

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
