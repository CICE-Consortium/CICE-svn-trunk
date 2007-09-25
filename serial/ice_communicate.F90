!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP

 module ice_communicate

! !MODULE: ice_communicate
! !DESCRIPTION:
!  This module contains the necessary routines and variables for
!  communicating between processors.  This instance of the module
!  is for serial execution so not much is done.
!
! !REVISION HISTORY:
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
!
! !USES:

   use ice_kinds_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public  :: init_communicate,          &
              get_num_procs,             &
              create_communicator

! !PUBLIC DATA MEMBERS:

   integer (int_kind), public :: &
      MPI_COMM_ICE,             &! MPI communicator for ice comms
      mpi_dbl,                  &! MPI type for dbl_kind
      my_task,                  &! MPI task number for this task
      master_task                ! task number of master task

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_communicate
! !INTERFACE:

 subroutine init_communicate

! !DESCRIPTION:
!  This routine sets up MPI environment and defines ice communicator.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

#ifdef coupled
   include 'mpif.h'   ! MPI Fortran include file

   integer (int_kind) :: ierr  ! MPI error flag
#endif

!-----------------------------------------------------------------------
!
!  initiate mpi environment and create communicator for internal
!  ice communications
!
!-----------------------------------------------------------------------

#ifdef coupled
   call MPI_INIT(ierr)
   call create_ice_communicator
   call MPI_COMM_RANK  (MPI_COMM_ICE, my_task, ierr)
#else
   my_task = 0
#endif

   master_task = 0

#ifdef coupled
!-----------------------------------------------------------------------
!
!  On some 64-bit machines where real_kind and dbl_kind are
!  identical, the MPI implementation uses MPI_REAL for both.
!  In these cases, set MPI_DBL to MPI_REAL.
!
!-----------------------------------------------------------------------

   MPI_DBL = MPI_DOUBLE_PRECISION

#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine init_communicate

!***********************************************************************
!BOP
! !IROUTINE: get_num_procs
! !INTERFACE:

 function get_num_procs()

! !DESCRIPTION:
!  This function returns the number of processors assigned to
!  the ice model.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (int_kind) :: get_num_procs

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  serial execution, must be only 1
!
!-----------------------------------------------------------------------

   get_num_procs = 1

!-----------------------------------------------------------------------
!EOC

 end function get_num_procs

!***********************************************************************
!BOP
! !IROUTINE: create_ice_communicator
! !INTERFACE:

 subroutine create_ice_communicator

! !DESCRIPTION:
!  This routine queries all the tasks in MPI_COMM_WORLD to see
!  which belong to the ice.  In standalone mode, this should
!  be all tasks, but in coupled mode CICE needs to determine
!  which tasks are assigned to the ice component.
!
!  this routine should be called after mpi_init, but before
!  setting up any internal mpi setups (since these will require
!  the internal communicators returned by this routine)
!
! !REVISION HISTORY:
!  same as module

#ifdef coupled
! !INCLUDES:

   include 'mpif.h'

#endif
!EOP
!BOC
#ifdef coupled
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (3) :: cmodel   ! model name temporary

   integer (int_kind) :: &
     MPI_GROUP_WORLD,       &! group id for MPI_COMM_WORLD
     MPI_GROUP_ATM,         &! group of processors assigned to atm
     MPI_GROUP_OCN,         &! group of processors assigned to ocn
     MPI_GROUP_ICE,         &! group of processors assigned to ice
     MPI_GROUP_LND,         &! group of processors assigned to lnd
     MPI_GROUP_CPL,         &! group of processors assigned to cpl
     MPI_COMM_ATM,          &! group of processors assigned to atm
     MPI_COMM_ICE,          &! group of processors assigned to ice
     MPI_COMM_LND,          &! group of processors assigned to lnd
     MPI_COMM_CPL            ! group of processors assigned to cpl

   integer (int_kind) :: &
     n,                     &! dummy loop counter
     ierr,                  &! error flag for MPI comms
     nprocs_all,            &! total processor count
     my_task_all,           &! rank of process in coupled domain
     ntasks_atm,            &! num tasks assigned to atm
     ntasks_ocn,            &! num tasks assigned to ocn
     ntasks_ice,            &! num tasks assigned to ice
     ntasks_lnd,            &! num tasks assigned to lnd
     ntasks_cpl              ! num tasks assigned to cpl

   integer (int_kind), dimension(3) :: &
     range_ocn,             &! range of tasks assigned to ocean
     range_atm,             &! range of tasks assigned to atmos
     range_ice,             &! range of tasks assigned to ice
     range_lnd,             &! range of tasks assigned to land
     range_cpl               ! range of tasks assigned to coupler

!-----------------------------------------------------------------------
!
!  determine processor rank in full (coupled) domain
!
!-----------------------------------------------------------------------

   call MPI_COMM_RANK (MPI_COMM_WORLD, my_task_all, ierr)

!-----------------------------------------------------------------------
!
!  determine which group of processes assigned to each model
!  assume the first processor assigned to a model is the task that
!  will communicate coupled model messages
!
!-----------------------------------------------------------------------

   call MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs_all, ierr)

   ntasks_atm = 0
   ntasks_ocn = 0
   ntasks_ice = 0
   ntasks_lnd = 0
   ntasks_cpl = 0
   range_ocn(1) = nprocs_all
   range_atm(1) = nprocs_all
   range_ice(1) = nprocs_all
   range_lnd(1) = nprocs_all
   range_cpl(1) = nprocs_all
   range_ocn(2) = 0
   range_atm(2) = 0
   range_ice(2) = 0
   range_lnd(2) = 0
   range_cpl(2) = 0
   range_ocn(3) = 1
   range_atm(3) = 1
   range_ice(3) = 1
   range_lnd(3) = 1
   range_cpl(3) = 1

   !***
   !*** each processor broadcasts its model to all the processors
   !*** in the coupled domain
   !***

   do n=0,nprocs_all-1
     if (n == my_task_all) then
       cmodel = 'ice'
     else
       cmodel = 'unk'
     endif

     call MPI_BCAST(cmodel, 3, MPI_CHARACTER, n, MPI_COMM_WORLD, ierr)

     select case(cmodel)
     case ('ocn')
       ntasks_ocn = ntasks_ocn + 1
       range_ocn(1) = min(n,range_ocn(1))
       range_ocn(2) = max(n,range_ocn(2))
     case ('atm')
       ntasks_atm = ntasks_atm + 1
       range_atm(1) = min(n,range_atm(1))
       range_atm(2) = max(n,range_atm(2))
     case ('ice')
       ntasks_ice = ntasks_ice + 1
       range_ice(1) = min(n,range_ice(1))
       range_ice(2) = max(n,range_ice(2))
     case ('lnd')
       ntasks_lnd = ntasks_lnd + 1
       range_lnd(1) = min(n,range_lnd(1))
       range_lnd(2) = max(n,range_lnd(2))
     case ('cpl')
       ntasks_cpl = ntasks_cpl + 1
       range_cpl(1) = min(n,range_cpl(1))
       range_cpl(2) = max(n,range_cpl(2))
     case default
       stop 'Unknown model name in comm setup'
     end select

   end do

!-----------------------------------------------------------------------
!
!  create subroup and communicator for each models internal
!  communciations, note that MPI_COMM_CREATE must be called by
!  all processes in MPI_COMM_WORLD so this must be done by all
!  models consistently and in the same order.
!
!-----------------------------------------------------------------------

   call MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr)

   if (ntasks_atm > 0) &
     call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range_atm, &
                               MPI_GROUP_ATM, ierr)

   if (ntasks_ocn > 0) &
     call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range_ocn,  &
                               MPI_GROUP_OCN, ierr)

   if (ntasks_ice > 0) &
     call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range_ice,  &
                               MPI_GROUP_ICE, ierr)

   if (ntasks_lnd > 0) &
     call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range_lnd,  &
                               MPI_GROUP_LND, ierr)

   if (ntasks_cpl > 0) &
     call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range_cpl,  &
                               MPI_GROUP_CPL, ierr)

   if (ntasks_atm > 0) &
     call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_ATM,  &
                           MPI_COMM_ATM, ierr)

   if (ntasks_ocn > 0) &
     call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_OCN,  &
                           MPI_COMM_OCN, ierr)

   if (ntasks_ice > 0) &
     call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_ICE,  &
                           MPI_COMM_ICE, ierr)

   if (ntasks_lnd > 0) &
     call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_LND,  &
                           MPI_COMM_LND, ierr)

   if (ntasks_cpl > 0) &
     call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_CPL,  &
                           MPI_COMM_CPL, ierr)

#else
   MPI_COMM_ICE = 0
#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine create_ice_communicator

!***********************************************************************
!BOP
! !IROUTINE: create_communicator
! !INTERFACE:

 subroutine create_communicator(new_comm, num_procs)

! !DESCRIPTION:
!  This routine creates a separate communicator for a subset of
!  processors under default ice communicator.
!
!  this routine should be called from init_domain1 when the
!  domain configuration (e.g. nprocs_btrop) has been determined
!
! !REVISION HISTORY:
!  same as module

#ifdef coupled
! !INCLUDES:

   include 'mpif.h'

#endif
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      num_procs         ! num of procs in new distribution

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      new_comm          ! new communicator for this distribution

!EOP
!BOC
#ifdef coupled
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     MPI_GROUP_ICE,         &! group of processors assigned to ice
     MPI_GROUP_NEW           ! group of processors assigned to new dist

   integer (int_kind) :: &
     ierr                    ! error flag for MPI comms

   integer (int_kind), dimension(3) :: &
     range                   ! range of tasks assigned to new dist
                             !  (assumed 0,num_procs-1)

!-----------------------------------------------------------------------
!
!  determine group of processes assigned to distribution
!
!-----------------------------------------------------------------------

   call MPI_COMM_GROUP (MPI_COMM_ICE, MPI_GROUP_ICE, ierr)

   range(1) = 0
   range(2) = num_procs-1
   range(3) = 1

!-----------------------------------------------------------------------
!
!  create subroup and communicator for new distribution
!  note: MPI_COMM_CREATE must be called by all procs in MPI_COMM_ICE
!
!-----------------------------------------------------------------------

   call MPI_GROUP_RANGE_INCL(MPI_GROUP_ICE, 1, range, &
                             MPI_GROUP_NEW, ierr)

   call MPI_COMM_CREATE (MPI_COMM_ICE, MPI_GROUP_NEW,  &
                         new_comm, ierr)

#else
   new_comm = MPI_COMM_ICE
#endif
!-----------------------------------------------------------------------
!EOC

 end subroutine create_communicator

!***********************************************************************

 end module ice_communicate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
