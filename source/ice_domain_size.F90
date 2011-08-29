!=======================================================================
!BOP
!
! !MODULE: ice_domain_size
!
! !DESCRIPTION:
!
! Defines the global domain size and number of categories and layers.
! Code originally based on domain_size.F in POP
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! author Elizabeth C. Hunke, LANL
! 2004: Block structure and snow parameters added by William Lipscomb
!       Renamed (used to be ice_model_size)
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!       Removed hardwired sizes (NX...can now be set in compile scripts)
!
! !INTERFACE:
!
      module ice_domain_size
!
! !USES:
!
      use ice_kinds_mod
!
!EOP
!=======================================================================

      implicit none
      save

      integer (kind=int_kind), parameter :: &
        nx_global = NXGLOB    , & ! i-axis size
        ny_global = NYGLOB    , & ! j-axis size
        ncat      =   5       , & ! number of categories
        nilyr     =   4       , & ! number of ice layers per category
        ntilyr    = ncat*nilyr, & ! number of ice layers in all categories
        nslyr     =   1       , & ! number of snow layers per category
        ntslyr    = ncat*nslyr, & ! number of snow layers in all categories
        max_aero  =   6       , & ! maximum number of aerosols 
        n_aero    =   3       , & ! number of aerosols in use
        max_ntrcr =             & ! number of tracers (defined in ice_init)
                      1         & ! 1 = surface temperature
                  +   2         & ! number of ice, snow tracers excl. aero
                  + n_aero*4  , & ! number of aerosols * 4 aero layers
        max_nstrm =   5           ! max number of history output streams

      integer (kind=int_kind), parameter :: &
        block_size_x = BLCKX  , & ! size of block in first horiz dimension
        block_size_y = BLCKY      ! size of block in second horiz dimension

   !*** The model will inform the user of the correct
   !*** values for the parameter below.  A value higher than
   !*** necessary will not cause the code to fail, but will
   !*** allocate more memory than is necessary.  A value that
   !*** is too low will cause the code to exit.  
   !*** A good initial guess is found using
   !*** max_blocks = (nx_global/block_size_x)*(ny_global/block_size_y)/
   !***               num_procs
 
      integer (kind=int_kind), parameter :: &
        max_blocks = MXBLCKS      ! max number of blocks per processor

!=======================================================================

      end module ice_domain_size

!=======================================================================
