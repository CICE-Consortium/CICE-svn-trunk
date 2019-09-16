!  SVN:$Id: $
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_extrapolate

!  This module contains routines for updating halo regions
!  (ghost cells) via local extrapolation, for regional grids

   use ice_kinds_mod
   use ice_blocks, only: nx_block, ny_block, nghost, &
           nblocks_tot, ice_blocksNorth, &
           ice_blocksSouth, ice_blocksEast, ice_blocksWest, &
           ice_blocksEast2, ice_blocksWest2, &
           ice_blocksNorthEast, ice_blocksNorthWest, &
           ice_blocksEastNorthEast, ice_blocksWestNorthWest, &
           ice_blocksSouthEast, ice_blocksSouthWest, &
           ice_blocksGetNbrID, get_block_parameter
   use ice_distribution, only: distrb, ice_distributionGet

   implicit none
   private
   save

   public :: ice_HaloExtrapolate, ice_HaloNeumann

   interface ice_HaloExtrapolate  ! generic interface
      module procedure ice_HaloExtrapolate2DR8, &
                       ice_HaloExtrapolate2DR4, &
                       ice_HaloExtrapolate2DI4, &
                       ice_HaloExtrapolate3DR8, &
                       ice_HaloExtrapolate4DR8
   end interface

   interface ice_HaloNeumann  ! generic interface
      module procedure ice_HaloNeumann2DR8, &
                       ice_HaloNeumann2DR4, &
                       ice_HaloNeumann2DI4, &
                       ice_HaloNeumann3DR8, &
                       ice_HaloNeumann4DR8
   end interface

!***********************************************************************

contains

!***********************************************************************

 subroutine ice_HaloExtrapolate2DR8(ARRAY,dist,ew_bndy_type,ns_bndy_type)

!  This subroutine extrapolates ARRAY values into the first row or column 
!  of ghost cells, and is intended for grid variables whose ghost cells 
!  would otherwise be set using the default boundary conditions (Dirichlet 
!  or Neumann).
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface ice_HaloExtrapolate

   use ice_blocks, only: block, nblocks_x, nblocks_y, get_block
   use ice_constants, only: c2
   use ice_distribution, only: ice_distributionGetBlockID

   character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc                  ! ghost cell column or row

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,iblk) = c2*ARRAY(2,j,iblk) - ARRAY(3,j,iblk)
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, nghost + 1, -1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,iblk) = c2*ARRAY(ibc-1,j,iblk) - ARRAY(ibc-2,j,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,iblk) = c2*ARRAY(i,2,iblk) - ARRAY(i,3,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, nghost + 1, -1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,iblk) = c2*ARRAY(i,ibc-1,iblk) - ARRAY(i,ibc-2,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloExtrapolate2DR8

!***********************************************************************

 subroutine ice_HaloExtrapolate2DR4(ARRAY,dist,ew_bndy_type,ns_bndy_type)

!  This subroutine extrapolates ARRAY values into the first row or column 
!  of ghost cells, and is intended for grid variables whose ghost cells 
!  would otherwise be set using the default boundary conditions (Dirichlet 
!  or Neumann).
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
!  This is the specific interface for single precision arrays 
!  corresponding to the generic interface ice_HaloExtrapolate

   use ice_blocks, only: block, nblocks_x, nblocks_y, get_block
   use ice_constants, only: c2
   use ice_distribution, only: ice_distributionGetBlockID

   character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   real (real_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc                  ! ghost cell column or row

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,iblk) = c2*ARRAY(2,j,iblk) - ARRAY(3,j,iblk)
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, nghost + 1, -1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,iblk) = c2*ARRAY(ibc-1,j,iblk) - ARRAY(ibc-2,j,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,iblk) = c2*ARRAY(i,2,iblk) - ARRAY(i,3,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, nghost + 1, -1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,iblk) = c2*ARRAY(i,ibc-1,iblk) - ARRAY(i,ibc-2,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloExtrapolate2DR4

!***********************************************************************

 subroutine ice_HaloExtrapolate2DI4(ARRAY,dist,ew_bndy_type,ns_bndy_type)

!  This subroutine extrapolates ARRAY values into the first row or column 
!  of ghost cells, and is intended for grid variables whose ghost cells 
!  would otherwise be set using the default boundary conditions (Dirichlet 
!  or Neumann).
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
!  This is the specific interface for integer arrays 
!  corresponding to the generic interface ice_HaloExtrapolate

   use ice_blocks, only: block, nblocks_x, nblocks_y, get_block
   use ice_distribution, only: ice_distributionGetBlockID

   character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc                  ! ghost cell column or row

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,iblk) = 2*ARRAY(2,j,iblk) - ARRAY(3,j,iblk)
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, nghost + 1, -1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,iblk) = 2*ARRAY(ibc-1,j,iblk) - ARRAY(ibc-2,j,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,iblk) = 2*ARRAY(i,2,iblk) - ARRAY(i,3,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, nghost + 1, -1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,iblk) = 2*ARRAY(i,ibc-1,iblk) - ARRAY(i,ibc-2,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloExtrapolate2DI4

!***********************************************************************

 subroutine ice_HaloExtrapolate3DR8(ARRAY,dist,ew_bndy_type,ns_bndy_type)

!  This subroutine extrapolates ARRAY values into the first row or column 
!  of ghost cells, and is intended for grid variables whose ghost cells 
!  would otherwise be set using the default boundary conditions (Dirichlet 
!  or Neumann).
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface ice_HaloExtrapolate

   use ice_blocks, only: block, nblocks_x, nblocks_y, get_block
   use ice_constants, only: c2
   use ice_distribution, only: ice_distributionGetBlockID

   character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   real (dbl_kind), dimension(:,:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc                  ! ghost cell column or row

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,:,iblk) = c2*ARRAY(2,j,:,iblk) - ARRAY(3,j,:,iblk)
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, nghost + 1, -1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,:,iblk) = c2*ARRAY(ibc-1,j,:,iblk) - ARRAY(ibc-2,j,:,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,:,iblk) = c2*ARRAY(i,2,:,iblk) - ARRAY(i,3,:,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, nghost + 1, -1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,:,iblk) = c2*ARRAY(i,ibc-1,:,iblk) - ARRAY(i,ibc-2,:,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloExtrapolate3DR8

!***********************************************************************

 subroutine ice_HaloExtrapolate4DR8(ARRAY,dist,ew_bndy_type,ns_bndy_type)

!  This subroutine extrapolates ARRAY values into the first row or column 
!  of ghost cells, and is intended for grid variables whose ghost cells 
!  would otherwise be set using the default boundary conditions (Dirichlet 
!  or Neumann).
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface ice_HaloExtrapolate

   use ice_blocks, only: block, nblocks_x, nblocks_y, get_block
   use ice_constants, only: c2
   use ice_distribution, only: ice_distributionGetBlockID

   character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   real (dbl_kind), dimension(:,:,:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc                  ! ghost cell column or row

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,:,:,iblk) = c2*ARRAY(2,j,:,:,iblk) - ARRAY(3,j,:,:,iblk)
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, nghost + 1, -1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,:,:,iblk) = c2*ARRAY(ibc-1,j,:,:,iblk) - ARRAY(ibc-2,j,:,:,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,:,:,iblk) = c2*ARRAY(i,2,:,:,iblk) - ARRAY(i,3,:,:,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, nghost + 1, -1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,:,:,iblk) = c2*ARRAY(i,ibc-1,:,:,iblk) - ARRAY(i,ibc-2,:,:,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloExtrapolate4DR8

!***********************************************************************

 subroutine ice_HaloNeumann2DR8(ARRAY,dist,ew_bndy_type,ns_bndy_type)

!  This subroutine sets ARRAY values in the first row or column 
!  of ghost cells equal to the values just inside the physical domain.
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface ice_HaloNeumann

   use ice_blocks, only: block, nblocks_x, nblocks_y, get_block
   use ice_distribution, only: ice_distributionGetBlockID

   character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc                  ! ghost cell column or row

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,iblk) = ARRAY(2,j,iblk)
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, nghost + 1, -1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,iblk) = ARRAY(ibc-1,j,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,iblk) = ARRAY(i,2,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, nghost + 1, -1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,iblk) = ARRAY(i,ibc-1,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloNeumann2DR8

!***********************************************************************

 subroutine ice_HaloNeumann2DR4(ARRAY,dist,ew_bndy_type,ns_bndy_type)

!  This subroutine sets ARRAY values in the first row or column 
!  of ghost cells equal to the values just inside the physical domain.
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
!  This is the specific interface for single precision arrays 
!  corresponding to the generic interface ice_HaloNeumann

   use ice_blocks, only: block, nblocks_x, nblocks_y, get_block
   use ice_distribution, only: ice_distributionGetBlockID

   character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   real (real_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc                  ! ghost cell column or row

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,iblk) = ARRAY(2,j,iblk)
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, nghost + 1, -1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,iblk) = ARRAY(ibc-1,j,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,iblk) = ARRAY(i,2,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, nghost + 1, -1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,iblk) = ARRAY(i,ibc-1,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloNeumann2DR4

!***********************************************************************

 subroutine ice_HaloNeumann2DI4(ARRAY,dist,ew_bndy_type,ns_bndy_type)

!  This subroutine sets ARRAY values in the first row or column 
!  of ghost cells equal to the values just inside the physical domain.
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
!  This is the specific interface for integer arrays 
!  corresponding to the generic interface ice_HaloNeumann

   use ice_blocks, only: block, nblocks_x, nblocks_y, get_block
   use ice_distribution, only: ice_distributionGetBlockID

   character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc                  ! ghost cell column or row

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,iblk) = 2*ARRAY(2,j,iblk) - ARRAY(3,j,iblk)
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, nghost + 1, -1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,iblk) = 2*ARRAY(ibc-1,j,iblk) - ARRAY(ibc-2,j,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,iblk) = 2*ARRAY(i,2,iblk) - ARRAY(i,3,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, nghost + 1, -1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,iblk) = 2*ARRAY(i,ibc-1,iblk) - ARRAY(i,ibc-2,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloNeumann2DI4

!***********************************************************************

 subroutine ice_HaloNeumann3DR8(ARRAY,dist,ew_bndy_type,ns_bndy_type)

!  This subroutine sets ARRAY values in the first row or column 
!  of ghost cells equal to the values just inside the physical domain.
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface ice_HaloNeumann

   use ice_blocks, only: block, nblocks_x, nblocks_y, get_block
   use ice_distribution, only: ice_distributionGetBlockID

   character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   real (dbl_kind), dimension(:,:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc                  ! ghost cell column or row

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,:,iblk) = ARRAY(2,j,:,iblk)
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, nghost + 1, -1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,:,iblk) = ARRAY(ibc-1,j,:,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,:,iblk) = ARRAY(i,2,:,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, nghost + 1, -1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,:,iblk) = ARRAY(i,ibc-1,:,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloNeumann3DR8

!***********************************************************************

 subroutine ice_HaloNeumann4DR8(ARRAY,dist,ew_bndy_type,ns_bndy_type)

!  This subroutine sets ARRAY values in the first row or column 
!  of ghost cells equal to the values just inside the physical domain.
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.
!
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface ice_HaloNeumann

   use ice_blocks, only: block, nblocks_x, nblocks_y, get_block
   use ice_distribution, only: ice_distributionGetBlockID

   character (char_len) :: &
       ew_bndy_type,    &! type of domain bndy in each logical
       ns_bndy_type      !    direction (ew is i, ns is j)

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   real (dbl_kind), dimension(:,:,:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,           &! dummy loop indices
     numBlocks,       &! number of local blocks
     blockID,            &! block location
     ibc                  ! ghost cell column or row

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  Linear extrapolation
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblk = 1, numBlocks
      call ice_distributionGetBlockID(dist, iblk, blockID)
      this_block = get_block(blockID, blockID)

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            do j = 1, ny_block
               ARRAY(1,j,:,:,iblk) = ARRAY(2,j,:,:,iblk)
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_bndy_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, nghost + 1, -1
               if (this_block%i_glob(i) == 0) ibc = ibc - 1
            enddo
            do j = 1, ny_block
               ARRAY(ibc,j,:,:,iblk) = ARRAY(ibc-1,j,:,:,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_bndy_type) /= 'cyclic') then
            do i = 1, nx_block
               ARRAY(i,1,:,:,iblk) = ARRAY(i,2,:,:,iblk)
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_bndy_type) /= 'cyclic' .and. &
             trim(ns_bndy_type) /= 'tripole' .and. &
             trim(ns_bndy_type) /= 'tripoleT' ) then
            ! locate ghost cell column (avoid padding)
            ibc = ny_block
            do j = ny_block, nghost + 1, -1
               if (this_block%j_glob(j) == 0) ibc = ibc - 1
            enddo
            do i = 1, nx_block
               ARRAY(i,ibc,:,:,iblk) = ARRAY(i,ibc-1,:,:,iblk)
            enddo
         endif
      endif

   enddo ! iblk

!-----------------------------------------------------------------------

 end subroutine ice_HaloNeumann4DR8

!***********************************************************************

 end module ice_extrapolate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
