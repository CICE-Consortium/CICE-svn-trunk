!  SVN:$Id$
!=======================================================================

! Defines the global domain size and number of categories and layers.
! Code originally based on domain_size.F in POP
!
! author Elizabeth C. Hunke, LANL
! 2004: Block structure and snow parameters added by William Lipscomb
!       Renamed (used to be ice_model_size)
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!       Removed hardwired sizes (NX...can now be set in compile scripts)

      module ice_domain_size

      use ice_kinds_mod

!=======================================================================

      implicit none
      private
      save

      integer (kind=int_kind), parameter, public :: &
        nx_global = NXGLOB    , & ! i-axis size
        ny_global = NYGLOB    , & ! j-axis size
        ncat      = NICECAT   , & ! number of categories
        nilyr     = NICELYR   , & ! number of ice layers per category
        nslyr     = NSNWLYR   , & ! number of snow layers per category
        max_algae =   3       , & ! maximum number of algal types 
        max_dic   =   1       , & ! maximum number of dissolved inorganic carbon types 
        max_doc   =   3       , & ! maximum number of dissolved organic carbon types
        max_don   =   1       , & ! maximum number of dissolved organic nitrogen types
        max_fe    =   2       , & ! maximum number of iron types
        max_aero  =   6       , & ! maximum number of aerosols 
        n_aero    = NTRAERO   , & ! number of aerosols in use 
        n_zaero   = TRZAERO   , & ! number of z aerosols in use 
        n_algae   = TRALG     , & ! number of algae in use 
        n_doc     = TRDOC     , & ! number of DOC pools in use
        n_dic     = TRDIC     , & ! number of DIC pools in use
        n_don     = TRDON     , & ! number of DON pools in use
        n_fed     = TRFED     , & ! number of Fe  pools in use dissolved Fe
        n_fep     = TRFEP     , & ! number of Fe  pools in use particulate Fe
        nblyr     = NBGCLYR   , & ! number of bio/brine layers per category 
                                  ! maximum number of biology tracers + aerosols
                                  ! *** add to kscavz in ice_zbgc_shared.F90 
        max_nbtrcr= max_algae*2 & ! algal nitrogen and chlorophyll
                   + max_dic    & ! dissolved inorganic carbon
                   + max_doc    & ! dissolved organic carbon
                   + max_don    & ! dissolved organic nitrogen
                   + 4          & ! nitrate, ammonium, silicate, and PON
                   + 3          & ! DMSPp, DMSPd, DMS
                   + max_fe*2   & ! dissolved Fe and  particulate Fe
                   + max_aero,  & ! aerosols
        nltrcr    = (TRBGCZ+TRZS)*TRBRI, & ! number of zbgc (includes zaero) + 
                                           ! zsalinity tracers 
        max_nsw   = (nilyr+nslyr+2) & ! total chlorophyll plus aerosols
                  * (1+TRZAERO),& ! number of tracers active in shortwave calculation
        max_ntrcr =   1         & ! 1 = surface temperature              
                  + nilyr       & ! ice salinity
                  + nilyr       & ! ice enthalpy
                  + nslyr       & ! snow enthalpy
                              !!!!! optional tracers:
                  + TRAGE       & ! age
                  + TRFY        & ! first-year area
                  + TRLVL*2     & ! level/deformed ice
                  + TRPND*3     & ! ponds
                  + n_aero*4    & ! number of aerosols * 4 aero layers
                  + TRBRI       & ! brine height 
                  + TRBGCS      & ! skeletal layer BGC 
                  + TRZS  *TRBRI* nblyr    & ! zsalinity (off if TRBRI=0)
                  + TRBGCZ*TRBRI*(nblyr+3) & ! zbgc      (off if TRBRI=0) 
                  + TRBGCZ    , & ! mobile/stationary phase tracer 
        max_nstrm =   5           ! max number of history output streams

      integer (kind=int_kind), parameter, public :: &
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
 
      integer (kind=int_kind), parameter, public :: &
        max_blocks = MXBLCKS      ! max number of blocks per processor

!=======================================================================

      end module ice_domain_size

!=======================================================================
