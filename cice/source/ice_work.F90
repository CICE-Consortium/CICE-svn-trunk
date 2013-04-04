!=======================================================================
!
!BOP
!
! !MODULE: ice_work - globally accessible, temporary work arrays
!
! !DESCRIPTION:
!
! Declare globally accessible, temporary work arrays to conserve memory.
! These arrays should be used only within a single subroutine!
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_work
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks
!
!EOP
!
      implicit none

      private
      public :: init_work

      ! global

      real (kind=dbl_kind), dimension(:,:), allocatable, public :: &
         work_g1, &
         work_g2, &
         work_g3

      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
         work_g4

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
         work_g5

      real (kind=real_kind), dimension(:,:), allocatable, public :: &
         work_gr

      real (kind=real_kind), dimension(:,:,:), allocatable, public :: &
         work_gr3
 
      real (kind=real_kind), dimension(:,:,:,:), allocatable, public :: &
         work_gr4

      integer(kind=int_kind), dimension(:,:), allocatable, public :: &
         work_gi4

      integer(selected_int_kind(13)), dimension(:,:), allocatable, public :: &
         work_gi8
 
      integer(kind=int_kind), dimension(:,:,:), allocatable, public :: &
         work_gi5

      integer(selected_int_kind(13)), dimension(:,:,:), allocatable, public :: &
         work_gi9

      ! all blocks
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         work1, &
         work2

      ! local (single block)
      real (kind=dbl_kind), dimension (nx_block,ny_block), public :: &
         worka, &
         workb, &
         workc, &
         workd

      ! vertical
      real (kind=dbl_kind), dimension(:), allocatable, public :: &
         work_z

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_work - initialize work arrays
!
! !DESCRIPTION:
!
! Initialize work arrays
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine init_work
!
! !USES:
!
      use ice_constants, only: c0
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      work1(:,:,:) = c0
      work2(:,:,:) = c0

      worka(:,:) = c0
      workb(:,:) = c0
      workc(:,:) = c0
      workd(:,:) = c0

      end subroutine init_work

!=======================================================================

      end module ice_work

!=======================================================================
