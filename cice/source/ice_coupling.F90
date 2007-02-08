!=======================================================================
!
!BOP
!
! !MODULE: ice_coupling - message passing to and from the coupler
!
! !DESCRIPTION:
! (1) Message passing to and from the coupler
! (2) ESMF import/export state routines
! These are mutually exclusive.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! author: Elizabeth C. Hunke, LANL
!         Tony Craig, NCAR, Dec-30-2002, modified for cpl6
!         Philip W Jones, LANL, modified for ESMF import/export
!
! 2004: Block structure added by William Lipscomb
! 2005: Added ESMF import/export state routines
! 2006: Updated the code for newest CCSM cpl6 - D. Bailey
!
! !INTERFACE:
!
      module ice_coupling
!
! !USES:
!
#ifdef USE_ESMF
      use esmf_mod
#endif
      use ice_kinds_mod
      use ice_blocks
      use ice_boundary
      use ice_domain_size
      use ice_domain
      use ice_constants
      use ice_calendar
      use ice_grid
      use ice_state
      use ice_flux
      use ice_timers
      use ice_fileunits
      use ice_communicate, only: my_task, master_task
#ifdef CCSM
      use shr_sys_mod, only : shr_sys_flush
      use ice_restart, only : runtype
      use cpl_contract_mod
      use cpl_interface_mod
      use cpl_fields_mod
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

#ifdef USE_ESMF
   public :: CICE_CoupledInit,          &
             CICE_CoupledExtractImport, &
             CICE_CoupledFillExport,    &
             CICE_CoupledAccumulateExport
#elif CCSM
   public :: ice_coupling_setup, &
             init_cpl,           &
             to_coupler,         &
             from_coupler,       &
             exit_coupler
#endif

! !PUBLIC DATA MEMBERS:

   logical(kind=log_kind), public :: &
    &  l_CCSMcoupled     ! true if coupled to CCSM

!
!EOP
!BOC
!

#ifdef USE_ESMF
   !--------------------------------------------------------------------
   !  import state fields
   !--------------------------------------------------------------------

   integer (kind=int_kind), parameter :: &
      numImportFields = 21   ! number of fields in import state

   !*** pointers for physical domain only

   real (kind=dbl_kind), dimension(:,:), pointer :: &
      importPtrZlvl   ,   &! altitude of atm fields
      importPtrUatm   ,   &! E-W wind speed
      importPtrVatm   ,   &! N-S wind speed
      importPtrPotT   ,   &! atm potential temperature
      importPtrTair   ,   &! atm air temperature
      importPtrQa     ,   &! humidity
      importPtrRhoa   ,   &! air density
      importPtrSWvdr  ,   &! downward shortwave heat flux v direct
      importPtrSWvdf  ,   &! downward shortwave heat flux v diffuse
      importPtrSWidr  ,   &! downward shortwave heat flux i direct
      importPtrSWidf  ,   &! downward shortwave heat flux i diffuse
      importPtrFlw    ,   &! downward longwave  heat flux
      importPtrFrain  ,   &! downward water flux due to rain
      importPtrFsnow  ,   &! downward water flux due to snow
      importPtrSST    ,   &! sea surface temperature
      importPtrSSS    ,   &! sea surface salinity
      importPtrUocn   ,   &! E-W ocean current speed
      importPtrVocn   ,   &! N-S ocean current speed
      importPtrSSHtltx,   &! E-W sea surface slope
      importPtrSSHtlty,   &! N-S sea surface slope
      importPtrFrzmlt      ! Freeze flux or melt potential

   type (ESMF_Array) ::   &! ESMF array for:
      importArrayZlvl   , &! altitude of atm fields
      importArrayUatm   , &! E-W wind speed
      importArrayVatm   , &! N-S wind speed
      importArrayPotT   , &! atm potential temperature
      importArrayTair   , &! atm air temperature
      importArrayQa     , &! humidity
      importArrayRhoa   , &! air density
      importArraySWvdr  , &! downward shortwave heat flux v direct
      importArraySWvdf  , &! downward shortwave heat flux v diffuse
      importArraySWidr  , &! downward shortwave heat flux i direct
      importArraySWidf  , &! downward shortwave heat flux i diffuse
      importArrayFlw    , &! downward longwave  heat flux
      importArrayFrain  , &! downward water flux due to rain
      importArrayFsnow  , &! downward water flux due to snow
      importArraySST    , &! sea surface temperature
      importArraySSS    , &! sea surface salinity
      importArrayUocn   , &! E-W ocean current speed
      importArrayVocn   , &! N-S ocean current speed
      importArraySSHtltx, &! E-W sea surface slope
      importArraySSHtlty, &! N-S sea surface slope
      importArrayFrzmlt    ! Freeze flux or melt potential

   type (ESMF_Field), dimension(numImportFields) ::   &
      importFieldList       ! List of ESMF Fields for import state

   integer (kind=int_kind) :: &! location in field list for
      importIndxZlvl   ,   &! altitude of atm fields
      importIndxUatm   ,   &! E-W wind speed
      importIndxVatm   ,   &! N-S wind speed
      importIndxPotT   ,   &! atm potential temperature
      importIndxTair   ,   &! atm air temperature
      importIndxQa     ,   &! humidity
      importIndxRhoa   ,   &! air density
      importIndxSWvdr  ,   &! downward shortwave heat flux v direct
      importIndxSWvdf  ,   &! downward shortwave heat flux v diffuse
      importIndxSWidr  ,   &! downward shortwave heat flux i direct
      importIndxSWidf  ,   &! downward shortwave heat flux i diffuse
      importIndxFlw    ,   &! downward longwave  heat flux
      importIndxFrain  ,   &! downward water flux due to rain
      importIndxFsnow  ,   &! downward water flux due to snow
      importIndxSST    ,   &! sea surface temperature
      importIndxSSS    ,   &! sea surface salinity
      importIndxUocn   ,   &! E-W ocean current speed
      importIndxVocn   ,   &! N-S ocean current speed
      importIndxSSHtltx,   &! E-W sea surface slope
      importIndxSSHtlty,   &! N-S sea surface slope
      importIndxFrzmlt      ! Freeze flux or melt potential

   !--------------------------------------------------------------------
   !  export state fields
   !--------------------------------------------------------------------

   integer (kind=int_kind), parameter :: &
      numExportFields = 21     ! number of fields in export state

   real (kind=dbl_kind), dimension(:,:), pointer :: &
      exportPtrIfrac,         &! ice fraction
      exportPtrTsrf,          &! ice surface temperature
      exportPtrAlbvdr,        &! albedo, visible, direct
      exportPtrAlbidr,        &! albedo, nearIR , direct
      exportPtrAlbvdf,        &! albedo, visible, diffuse
      exportPtrAlbidf,        &! albedo, nearIR , diffuse
      exportPtrTauxa,         &! atm-ice wind stress e-w
      exportPtrTauya,         &! atm-ice wind stress n-s
      exportPtrFlat,          &! latent heat flx
      exportPtrFsens,         &! sensible heat flx
      exportPtrFLWout,        &! upward longwave heat flx
      exportPtrEvap,          &! evaporation water flx
      exportPtrTref,          &! 2m reference temperature
      exportPtrQref,          &! 2m reference spec humidity
      exportPtrFSWabs,        &! net SW absorb heat flx
      exportPtrFSWthru,       &! solar heat flux thru ice to ocean
      exportPtrFmelt,         &! melt heat flx
      exportPtrFresh,         &! fresh water flx from melt
      exportPtrFsalt,         &! salt flx from melt
      exportPtrTauxo,         &! ice-ocean stress e-w
      exportPtrTauyo           ! ice-ocean stress n-s

   type (ESMF_Array) ::        &! ESMF arrays for
      exportArrayIfrac,         &! ice fraction
      exportArrayTsrf,          &! ice surface temperature
      exportArrayAlbvdr,        &! albedo, visible, direct
      exportArrayAlbidr,        &! albedo, nearIR , direct
      exportArrayAlbvdf,        &! albedo, visible, diffuse
      exportArrayAlbidf,        &! albedo, nearIR , diffuse
      exportArrayTauxa,         &! atm-ice wind stress e-w
      exportArrayTauya,         &! atm-ice wind stress n-s
      exportArrayFlat,          &! latent heat flx
      exportArrayFsens,         &! sensible heat flx
      exportArrayFLWout,        &! upward longwave heat flx
      exportArrayEvap,          &! evaporation water flx
      exportArrayTref,          &! 2m reference temperature
      exportArrayQref,          &! 2m reference spec humidity
      exportArrayFSWabs,        &! net SW absorb heat flx
      exportArrayFSWthru,       &! solar heat flux thru ice to ocean
      exportArrayFmelt,         &! melt heat flx
      exportArrayFresh,         &! fresh water flx from melt
      exportArrayFsalt,         &! salt flx from melt
      exportArrayTauxo,         &! ice-ocean stress e-w
      exportArrayTauyo           ! ice-ocean stress n-s

   type (ESMF_Field), dimension(numImportFields) ::   &
      exportFieldList           ! List of ESMF Fields for export state
   integer (kind=int_kind) ::  &! location in field list, avg buffer for
      exportIndxIfrac,         &! ice fraction
      exportIndxTsrf,          &! ice surface temperature
      exportIndxAlbvdr,        &! albedo, visible, direct
      exportIndxAlbidr,        &! albedo, nearIR , direct
      exportIndxAlbvdf,        &! albedo, visible, diffuse
      exportIndxAlbidf,        &! albedo, nearIR , diffuse
      exportIndxTauxa,         &! atm-ice wind stress e-w
      exportIndxTauya,         &! atm-ice wind stress n-s
      exportIndxFlat,          &! latent heat flx
      exportIndxFsens,         &! sensible heat flx
      exportIndxFLWout,        &! upward longwave heat flx
      exportIndxEvap,          &! evaporation water flx
      exportIndxTref,          &! 2m reference temperature
      exportIndxQref,          &! 2m reference spec humidity
      exportIndxFSWabs,        &! net SW absorb heat flx
      exportIndxFSWthru,       &! solar heat flux thru ice to ocean
      exportIndxFmelt,         &! melt heat flx
      exportIndxFresh,         &! fresh water flx from melt
      exportIndxFsalt,         &! salt flx from melt
      exportIndxTauxo,         &! ice-ocean stress e-w
      exportIndxTauyo           ! ice-ocean stress n-s
   real (kind=dbl_kind) :: &
      timeLastCoupled        ! accumulated time since last coupling

   real (kind=dbl_kind), &
      dimension(nx_block, ny_block, max_blocks, numExportFields) :: &
      exportBufferSum     ! buffer for time avg of export fields
#endif
! ifdef esmf

#ifdef CCSM  

      integer (kind=int_kind), dimension (cpl_fields_ibuf_total) :: &
         isbuf &
      ,  irbuf

      real (kind=dbl_kind), allocatable :: &
         sbuf(:,:)

      integer (kind=int_kind), save :: &
         nrecv, nsend

      real (kind=dbl_kind), allocatable :: &
         buffs(:,:) &
      ,  buffr(:,:)

      type(cpl_contract) :: &
         contractS &
      ,  contractR

      integer(kind=int_kind), save :: &
         nadv_i &
      ,  info_dbug

      integer(kind=int_kind), save :: index_i2c_Si_t        ! temperature
      integer(kind=int_kind), save :: index_i2c_Si_tref     ! 2m reference temperature
      integer(kind=int_kind), save :: index_i2c_Si_qref     ! 2m reference specific humidity
      integer(kind=int_kind), save :: index_i2c_Si_ifrac    ! fractional ice coverage
      integer(kind=int_kind), save :: index_i2c_Si_avsdr    ! albedo: visible, direct
      integer(kind=int_kind), save :: index_i2c_Si_anidr    ! albedo: near ir, direct
      integer(kind=int_kind), save :: index_i2c_Si_avsdf    ! albedo: visible, diffuse
      integer(kind=int_kind), save :: index_i2c_Si_anidf    ! albedo: near ir, diffuse
      integer(kind=int_kind), save :: index_i2c_index       ! global data compr index
      integer(kind=int_kind), save :: index_i2c_Faii_taux   ! wind stress, zonal
      integer(kind=int_kind), save :: index_i2c_Faii_tauy   ! wind stress, meridional
      integer(kind=int_kind), save :: index_i2c_Faii_lat    ! latent          heat flux
      integer(kind=int_kind), save :: index_i2c_Faii_sen    ! sensible        heat flux
      integer(kind=int_kind), save :: index_i2c_Faii_lwup   ! upward longwave heat flux
      integer(kind=int_kind), save :: index_i2c_Faii_evap   ! evaporation    water flux
      integer(kind=int_kind), save :: index_i2c_Faii_swnet  ! shortwave: net absorbed
      integer(kind=int_kind), save :: index_i2c_Fioi_swpen  ! net SW penetrating ice
      integer(kind=int_kind), save :: index_i2c_Fioi_melth  ! heat  flux from melting ice
      integer(kind=int_kind), save :: index_i2c_Fioi_meltw  ! water flux from melting ice
      integer(kind=int_kind), save :: index_i2c_Fioi_salt   ! salt  flux from melting ice
      integer(kind=int_kind), save :: index_i2c_Fioi_taux   ! ice/ocn stress, zonal
      integer(kind=int_kind), save :: index_i2c_Fioi_tauy   ! ice/ocn stress, meridional
      
      integer(kind=int_kind), save :: index_c2i_So_t        ! ocn temp
      integer(kind=int_kind), save :: index_c2i_So_s        ! ocn salinity
      integer(kind=int_kind), save :: index_c2i_So_u        ! ocn u velocity
      integer(kind=int_kind), save :: index_c2i_So_v        ! ocn v velocity
      integer(kind=int_kind), save :: index_c2i_So_dhdx     ! ocn surface slope, zonal
      integer(kind=int_kind), save :: index_c2i_So_dhdy     ! ocn surface slope, merid
      integer(kind=int_kind), save :: index_c2i_Sa_z        ! atm bottom layer height
      integer(kind=int_kind), save :: index_c2i_Sa_u        ! atm u velocity
      integer(kind=int_kind), save :: index_c2i_Sa_v        ! atm v velocity
      integer(kind=int_kind), save :: index_c2i_Sa_ptem     ! atm potential temp
      integer(kind=int_kind), save :: index_c2i_Sa_tbot     ! atm bottom temp
      integer(kind=int_kind), save :: index_c2i_Sa_shum     ! atm specfic humidity
      integer(kind=int_kind), save :: index_c2i_Sa_dens     ! atm air density
      integer(kind=int_kind), save :: index_c2i_Fioo_q      ! ocn freeze or melt heat
      integer(kind=int_kind), save :: index_c2i_Faxa_swndr  ! atm sw near-ir, direct
      integer(kind=int_kind), save :: index_c2i_Faxa_swvdr  ! atm sw visable, direct
      integer(kind=int_kind), save :: index_c2i_Faxa_swndf  ! atm sw near-ir, diffuse
      integer(kind=int_kind), save :: index_c2i_Faxa_swvdf  ! atm sw visable, diffuse
      integer(kind=int_kind), save :: index_c2i_Faxa_lwdn   ! long-wave down
      integer(kind=int_kind), save :: index_c2i_Faxc_rain   ! rain
      integer(kind=int_kind), save :: index_c2i_Faxc_snow   ! snow

      type (block) :: this_block
#endif
! endif CCSM
!
! EOC
!
#ifdef USE_ESMF   
! next four subroutines
!=======================================================================

 contains

!=======================================================================
!
!BOP
!
! !IROUTINE: CICE_CoupledInit - initializes states and ESMF coupling
!
! !INTERFACE:
!
 subroutine CICE_CoupledInit(importState, exportState, errorCode)
!
! !DESCRIPTION:
!
!  This routine sets up everything necessary for coupling with
!  other ESMF components, including setting up the import and
!  export states for the run method and initializing the data
!  in the export state for use in other models.
!
! !REVISION HISTORY:
!
!  author: Philip W Jones, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:

   type (ESMF_State), intent(inout) :: &
      importState,        &! import state for CICE run method
      exportState          ! export state for CICE run method

! !OUTPUT PARAMETERS:

   integer (kind=int_kind), intent(out) :: &
      errorCode        ! error code: ESMF_Success or ESMF_Failure

!
!EOP
!BOC
!
   !--------------------------------------------------------------------
   !  local variables
   !--------------------------------------------------------------------

   integer (kind=int_kind) :: &
      i,j,k,             &! dummy loop index
      iblk,              &! dummy block index
      istat               ! error flag for allocates

   !integer (kind=int_kind), dimension(nibuff) :: &
   !   ibuff               ! integer control buffer

   !--------------------------------------------------------------------
   !  initialize some scalars and perform some checks
   !--------------------------------------------------------------------

   errorCode         = ESMF_SUCCESS
   timeLastCoupled   = c0

   if (max_blocks > 1) then
      write(nu_diag,*) &
         '(ice) CICE_CoupledInit: ESMF only supports 1 block/proc'
      return
   endif

   !if (.not. l_CCSMcoupled) return

   !--------------------------------------------------------------------
   !  create CICE import state with ESMF fields and arrays
   !--------------------------------------------------------------------

   !*** Currently assume 1 block/node

   importIndxZlvl   =  1
   importIndxUatm   =  2
   importIndxVatm   =  3
   importIndxPotT   =  4
   importIndxTair   =  5
   importIndxQa     =  6
   importIndxRhoa   =  7
   importIndxSWvdr  =  8
   importIndxSWvdf  =  9
   importIndxSWidr  = 10
   importIndxSWidf  = 11
   importIndxFlw    = 12
   importIndxFrain  = 13
   importIndxFsnow  = 14
   importIndxSST    = 15
   importIndxSSS    = 16
   importIndxUocn   = 17
   importIndxVocn   = 18
   importIndxSSHtltx= 19
   importIndxSSHtlty= 20
   importIndxFrzmlt = 21

   allocate( importPtrZlvl   (nx_block, ny_block), &
             importPtrUatm   (nx_block, ny_block), &
             importPtrVatm   (nx_block, ny_block), &
             importPtrPotT   (nx_block, ny_block), &
             importPtrTair   (nx_block, ny_block), &
             importPtrQa     (nx_block, ny_block), &
             importPtrRhoa   (nx_block, ny_block), &
             importPtrSWvdr  (nx_block, ny_block), &
             importPtrSWvdf  (nx_block, ny_block), &
             importPtrSWidr  (nx_block, ny_block), &
             importPtrSWidf  (nx_block, ny_block), &
             importPtrFlw    (nx_block, ny_block), &
             importPtrFrain  (nx_block, ny_block), &
             importPtrFsnow  (nx_block, ny_block), &
             importPtrSST    (nx_block, ny_block), &
             importPtrSSS    (nx_block, ny_block), &
             importPtrUocn   (nx_block, ny_block), &
             importPtrVocn   (nx_block, ny_block), &
             importPtrSSHtltx(nx_block, ny_block), &
             importPtrSSHtlty(nx_block, ny_block), &
             importPtrFrzmlt (nx_block, ny_block), &
             stat = istat)

   if (istat > 0) then
      write(nu_diag,*) &
        '(ice) CICE_CoupledInit: error allocating import state pointers'
      return
   endif

   importArrayZlvl    = ESMF_ArrayCreate(importPtrZlvl,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayUatm    = ESMF_ArrayCreate(importPtrUatm,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayVatm    = ESMF_ArrayCreate(importPtrVatm,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayPotT    = ESMF_ArrayCreate(importPtrPotT,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayTair    = ESMF_ArrayCreate(importPtrTair,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayQa      = ESMF_ArrayCreate(importPtrQa,            &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayRhoa    = ESMF_ArrayCreate(importPtrRhoa,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArraySWvdr   = ESMF_ArrayCreate(importPtrSWvdr,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArraySWvdf   = ESMF_ArrayCreate(importPtrSWvdf,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArraySWidr   = ESMF_ArrayCreate(importPtrSWidr,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArraySWidf   = ESMF_ArrayCreate(importPtrSWidf,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayFlw     = ESMF_ArrayCreate(importPtrFlw,           &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayFrain   = ESMF_ArrayCreate(importPtrFrain,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayFsnow   = ESMF_ArrayCreate(importPtrFsnow,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArraySST     = ESMF_ArrayCreate(importPtrSST,           &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArraySSS     = ESMF_ArrayCreate(importPtrSSS,           &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayUocn    = ESMF_ArrayCreate(importPtrUocn,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayVocn    = ESMF_ArrayCreate(importPtrVocn,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArraySSHtltx = ESMF_ArrayCreate(importPtrSSHtltx,       &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArraySSHtlty = ESMF_ArrayCreate(importPtrSSHtlty,       &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   importArrayFrzmlt  = ESMF_ArrayCreate(importPtrFrzmlt,        &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)

   if (errorCode /= ESMF_SUCCESS) then
      write(nu_diag,*) &
         '(ice) CICE_CoupledInit: error creating import state arrays'
      return
   endif

!----------------------------------------------------------------------
! Due to problems with grids, we add arrays to the import state
! rather than Fields, so give the arrays names and then add to states
!----------------------------------------------------------------------

   call ESMF_ArraySet(importArrayZlvl, &
                      name = 'Height', rc = errorCode)
   call ESMF_ArraySet(importArrayUatm, &
                      name = 'E-W wind speed', rc = errorCode)
   call ESMF_ArraySet(importArrayVatm, &
                      name = 'N-S wind speed', rc = errorCode)
   call ESMF_ArraySet(importArrayPotT, &
                      name = 'Potential temperature', rc = errorCode)
   call ESMF_ArraySet(importArrayTair, &
                      name = 'Air temperature', rc = errorCode)
   call ESMF_ArraySet(importArrayQa, &
                      name = 'Humidity', rc = errorCode)
   call ESMF_ArraySet(importArrayRhoa, &
                      name = 'Air density', rc = errorCode)
   call ESMF_ArraySet(importArraySWvdr, &
                      name = 'SW heat visible direct', rc = errorCode)
   call ESMF_ArraySet(importArraySWvdf, &
                      name = 'SW heat visible diffuse', rc = errorCode)
   call ESMF_ArraySet(importArraySWidr, &
                      name = 'SW heat nearIR direct', rc = errorCode)
   call ESMF_ArraySet(importArraySWidf, &
                      name = 'SW heat nearIR diffuse', rc = errorCode)
   call ESMF_ArraySet(importArrayFlw, &
                      name = 'Longwave heat flux', rc = errorCode)
   call ESMF_ArraySet(importArrayFrain, &
                      name = 'Rain water flux', rc = errorCode)
   call ESMF_ArraySet(importArrayFsnow, &
                      name = 'Snow water flux', rc = errorCode)
   call ESMF_ArraySet(importArraySST, &
                      name = 'Sea surface temperature', rc = errorCode)
   call ESMF_ArraySet(importArraySSS, &
                      name = 'Sea surface salinity', rc = errorCode)
   call ESMF_ArraySet(importArrayUocn, &
                      name = 'Ocean current e-w', rc = errorCode)
   call ESMF_ArraySet(importArrayVocn, &
                      name = 'Ocean current n-s', rc = errorCode)
   call ESMF_ArraySet(importArraySSHtltx, &
                      name = 'Sea surface slope e-w', rc = errorCode)
   call ESMF_ArraySet(importArraySSHtlty, &
                      name = 'Sea surface slope n-s',rc = errorCode)
   call ESMF_ArraySet(importArrayFrzmlt, &
                      name = 'Freeze heat flux or melt potential',&
                      rc = errorCode)

   if (errorCode /= ESMF_SUCCESS) then
      write(nu_diag,*) &
         '(ice) CICE_CoupledInit: error naming import state arrays'
      return
   endif

   !
   ! create empty import state
   !
   importState = ESMF_StateCreate(stateName='CICE Import State', &
                                  statetype = ESMF_STATE_IMPORT, &
                                  rc = errorCode)

   if (errorCode /= ESMF_SUCCESS) then
      write(nu_diag,*) &
         '(ice) CICE_CoupledInit: error creating import state'
      return
   endif

   call ESMF_StateAddArray(importState, importArrayZlvl, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayUatm, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayVatm, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayPotT, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayTair, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayQa, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayRhoa, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArraySWvdr, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArraySWvdf, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArraySWidr, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArraySWidf, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayFlw, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayFrain, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayFsnow, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArraySST, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArraySSS, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayUocn, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayVocn, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArraySSHtltx, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArraySSHtlty, &
                           rc = errorCode)
   call ESMF_StateAddArray(importState, importArrayFrzmlt, &
                           rc = errorCode)

   if (errorCode /= ESMF_SUCCESS) then
      write(nu_diag,*) &
         '(ice) CICE_CoupledInit: error adding import state arrays'
      return
   endif

!----------------------------------------------------------------------
!
!  Here is the proper field code once the grid pieces are working
!
!----------------------------------------------------------------------
!   importFieldList(importIndxZlvl   ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayZlvl   ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Height',              &
!                                          rc = errorCode)
!
!   importFieldList(importIndxUatm   ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayUatm   ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'E-W wind speed',&
!                                          rc = errorCode)
!   importFieldList(importIndxVatm   ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayVatm   ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'N-S wind speed',      &
!                                          rc = errorCode)
!   importFieldList(importIndxPotT   ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayPotT   ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Potential temperature',&
!                                          rc = errorCode)
!   importFieldList(importIndxTair   ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayTair   ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Air temperature',     &
!                                          rc = errorCode)
!   importFieldList(importIndxQa     ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayQa     ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Humidity',            &
!                                          rc = errorCode)
!   importFieldList(importIndxRhoa   ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayRhoa   ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Air density',         &
!                                          rc = errorCode)
!   importFieldList(importIndxSWvdr  ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArraySWvdr  ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'SW heat visible direct',&
!                                          rc = errorCode)
!   importFieldList(importIndxSWvdf  ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArraySWvdf  ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'SW heat visible diffuse',&
!                                          rc = errorCode)
!   importFieldList(importIndxSWidr  ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArraySWidr  ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'SW heat nearIR direct',&
!                                          rc = errorCode)
!   importFieldList(importIndxSWidf  ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArraySWidf  ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'SW heat nearIR diffuse',&
!                                          rc = errorCode)
!   importFieldList(importIndxFlw    ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayFlw    ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Longwave heat flux',  &
!                                          rc = errorCode)
!   importFieldList(importIndxFrain  ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayFrain  ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Rain water flux',     &
!                                          rc = errorCode)
!   importFieldList(importIndxFsnow  ) = ESMF_FieldCreate(               &
!!                                          iceGrid, importArrayFsnow  ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Snow water flux',     &
!                                          rc = errorCode)
!   importFieldList(importIndxSST    ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArraySST    ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Sea surface temperature',&
!                                          rc = errorCode)
!   importFieldList(importIndxSSS    ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArraySSS    ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Sea surface salinity',&
!                                          rc = errorCode)
!   importFieldList(importIndxUocn   ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayUocn   ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Ocean current e-w',   &
!                                          rc = errorCode)
!   importFieldList(importIndxVocn   ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayVocn   ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Ocean current n-s',   &
!                                          rc = errorCode)
!   importFieldList(importIndxSSHtltx) = ESMF_FieldCreate(               &
!                                          iceGrid, importArraySSHtltx,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Sea surface slope e-w',&
!                                          rc = errorCode)
!   importFieldList(importIndxSSHtlty) = ESMF_FieldCreate(               &
!                                          iceGrid, importArraySSHtlty,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Sea surface slope n-s',&
!                                          rc = errorCode)
!   importFieldList(importIndxFrzmlt ) = ESMF_FieldCreate(               &
!                                          iceGrid, importArrayFrzmlt ,  &
!                                          copyflag = ESMF_DATA_REF,     &
!                                          horzRelloc = ESMF_CELL_CENTER,&
!                                          haloWidth = nghost,           &
!                                          name = 'Freeze heat flux or melt potential',&
!                                          rc = errorCode)
!
!   if (errorCode /= ESMF_SUCCESS) then
!      write(nu_diag,*) &
!         '(ice) CICE_CoupledInit: error creating import state fields'
!      return
!   endif
!
!   !--------------------------------------------------------------------
!   ! add attributes to import fields
!   !--------------------------------------------------------------------
!
!   call ESMF_FieldSetAttribute(importFieldList(importIndxZlvl   ),  &
!                                               'units', 'm',        &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxUatm   ),  &
!                                               'units', 'm/s',      &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxVatm   ),  &
!                                               'units', 'm/s',      &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxPotT   ),  &
!                                               'units', 'deg C',    &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxTair   ),  &
!                                               'units', 'deg C',    &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxQa     ),  &
!                                               'units', ' ',        &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxRhoa   ),  &
!                                               'units', 'kg/m3',    &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxSWvdr  ),  &
!                                               'units', 'W/m2',     &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxSWvdf  ),  &
!                                               'units', 'W/m2',     &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxSWidr  ),  &
!                                               'units', 'W/m2',     &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxSWidf  ),  &
!                                               'units', 'W/m2',     &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxFlw    ),  &
!                                               'units', 'W/m2',     &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxFrain  ),  &
!                                               'units', 'kg/m2/s',  &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxFsnow  ),  &
!                                               'units', 'kg/m2/s',  &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxSST    ),  &
!                                               'units', 'deg C',    &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxSSS    ),  &
!                                               'units', 'psu',      &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxUocn   ),  &
!                                               'units', 'm/s',      &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxVocn   ),  &
!                                               'units', 'm/s',      &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxSSHtltx),  &
!                                               'units', 'm/m',      &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxSSHtlty),  &
!                                               'units', 'm/m',      &
!                                               rc = errorCode)
!   call ESMF_FieldSetAttribute(importFieldList(importIndxFrzmlt ),  &
!                                               'units', 'W/m2',     &
!                                               rc = errorCode)
!
!   if (errorCode /= ESMF_SUCCESS) then
!      write(nu_diag,*) &
!         '(ice) CICE_CoupledInit: error adding import attributes'
!      return
!   endif
!
   !--------------------------------------------------------------------
   !  create export state with ESMF fields and arrays
   !--------------------------------------------------------------------

   exportIndxIfrac   =  1
   exportIndxTsrf    =  2
   exportIndxAlbvdr  =  3
   exportIndxAlbidr  =  4
   exportIndxAlbvdf  =  5
   exportIndxAlbidf  =  6
   exportIndxTauxa   =  7
   exportIndxTauya   =  8
   exportIndxFlat    =  9
   exportIndxFsens   = 10
   exportIndxFLWout  = 11
   exportIndxEvap    = 12
   exportIndxTref    = 13
   exportIndxQref    = 14
   exportIndxFSWabs  = 15
   exportIndxFSWthru = 16
   exportIndxFmelt   = 17
   exportIndxFresh   = 18
   exportIndxFsalt   = 19
   exportIndxTauxo   = 20
   exportIndxTauyo   = 21

   allocate ( exportPtrIfrac  (nx_block, ny_block), &
              exportPtrTsrf   (nx_block, ny_block), &
              exportPtrAlbvdr (nx_block, ny_block), &
              exportPtrAlbidr (nx_block, ny_block), &
              exportPtrAlbvdf (nx_block, ny_block), &
              exportPtrAlbidf (nx_block, ny_block), &
              exportPtrTauxa  (nx_block, ny_block), &
              exportPtrTauya  (nx_block, ny_block), &
              exportPtrFlat   (nx_block, ny_block), &
              exportPtrFsens  (nx_block, ny_block), &
              exportPtrFLWout (nx_block, ny_block), &
              exportPtrEvap   (nx_block, ny_block), &
              exportPtrTref   (nx_block, ny_block), &
              exportPtrQref   (nx_block, ny_block), &
              exportPtrFSWabs (nx_block, ny_block), &
              exportPtrFSWthru(nx_block, ny_block), &
              exportPtrFmelt  (nx_block, ny_block), &
              exportPtrFresh  (nx_block, ny_block), &
              exportPtrFsalt  (nx_block, ny_block), &
              exportPtrTauxo  (nx_block, ny_block), &
              exportPtrTauyo  (nx_block, ny_block), &
              stat = istat)

   if (istat > 0) then
      write(nu_diag, *) &
         '(ice) CICE_CoupledInit: error allocating export state pointers'
      errorCode = ESMF_Failure
      return
   endif

   exportArrayIfrac   = ESMF_ArrayCreate(exportPtrIfrac,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayTsrf    = ESMF_ArrayCreate(exportPtrTsrf,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayAlbvdr  = ESMF_ArrayCreate(exportPtrAlbvdr,        &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayAlbidr  = ESMF_ArrayCreate(exportPtrAlbidr,        &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayAlbvdf  = ESMF_ArrayCreate(exportPtrAlbvdf,        &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayAlbidf  = ESMF_ArrayCreate(exportPtrAlbidf,        &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayTauxa   = ESMF_ArrayCreate(exportPtrTauxa,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayTauya   = ESMF_ArrayCreate(exportPtrTauya,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayFlat    = ESMF_ArrayCreate(exportPtrFlat,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayFsens   = ESMF_ArrayCreate(exportPtrFsens,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayFLWout  = ESMF_ArrayCreate(exportPtrFLWout,        &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayEvap    = ESMF_ArrayCreate(exportPtrEvap,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayTref    = ESMF_ArrayCreate(exportPtrTref,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayQref    = ESMF_ArrayCreate(exportPtrQref,          &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayFSWabs  = ESMF_ArrayCreate(exportPtrFSWabs,        &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayFSWthru = ESMF_ArrayCreate(exportPtrFSWthru,       &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayFmelt   = ESMF_ArrayCreate(exportPtrFmelt,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayFresh   = ESMF_ArrayCreate(exportPtrFresh,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayFsalt   = ESMF_ArrayCreate(exportPtrFsalt,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayTauxo   = ESMF_ArrayCreate(exportPtrTauxo,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)
   exportArrayTauyo   = ESMF_ArrayCreate(exportPtrTauyo,         &
                                         docopy = ESMF_DATA_REF, &
                                         haloWidth = nghost, rc=errorCode)

   if (errorCode /= ESMF_SUCCESS) then
      write(nu_diag,*) &
         '(ice) CICE_CoupledInit: error creating export state arrays'
      return
   endif

   !-------------------------------------------------------------------
   ! due to problems with grids, add arrays to export state instead
   ! of fields.  add names and then add arrays individually
   !-------------------------------------------------------------------

   call ESMF_ArraySet(exportArrayIfrac, &
                      name = 'Ice fraction', rc = errorCode)
   call ESMF_ArraySet(exportArrayTsrf, &
                      name = 'Surface temperature', rc = errorCode)
   call ESMF_ArraySet(exportArrayAlbvdr, &
                      name = 'Albedo: visible direct', rc = errorCode)
   call ESMF_ArraySet(exportArrayAlbidr, &
                      name = 'Albedo: near-ir direct', rc = errorCode)
   call ESMF_ArraySet(exportArrayAlbvdf, &
                      name = 'Albedo: visible diffuse', rc = errorCode)
   call ESMF_ArraySet(exportArrayAlbidf, &
                      name = 'Albedo: near-ir diffuse', rc = errorCode)
   call ESMF_ArraySet(exportArrayTauxa, &
                      name = 'Wind stress E-W',      rc = errorCode)
   call ESMF_ArraySet(exportArrayTauya, &
                      name = 'Wind stress N-S',      rc = errorCode)
   call ESMF_ArraySet(exportArrayFlat, &
                      name = 'Latent heat flux',     rc = errorCode)
   call ESMF_ArraySet(exportArrayFsens, &
                      name = 'Sensible heat flux',   rc = errorCode)
   call ESMF_ArraySet(exportArrayFLWout, &
                      name = 'Longwave radiative heat flux', rc = errorCode)
   call ESMF_ArraySet(exportArrayEvap, &
                      name = 'Evaporation water flux', rc = errorCode)
   call ESMF_ArraySet(exportArrayTref, &
                      name = '2m reference temperature', rc = errorCode)
   call ESMF_ArraySet(exportArrayQref, &
                      name = '2m reference humidity', rc = errorCode)
   call ESMF_ArraySet(exportArrayFSWabs, &
                      name = 'Short wave heat flux absorbed by ice',&
                      rc = errorCode)
   call ESMF_ArraySet(exportArrayFSWthru, &
                      name = 'Short wave heat flux into ocean', &
                      rc = errorCode)
   call ESMF_ArraySet(exportArrayFmelt, &
                      name = 'Melt heat flux',       rc = errorCode)
   call ESMF_ArraySet(exportArrayFresh, &
                      name = 'Fresh water flux',     rc = errorCode)
   call ESMF_ArraySet(exportArrayFsalt, &
                      name = 'Salt flux',            rc = errorCode)
   call ESMF_ArraySet(exportArrayTauxo, &
                      name = 'Ocean-ice stress E-W', rc = errorCode)
   call ESMF_ArraySet(exportArrayTauyo, &
                      name = 'Ocean-ice stress N-S', rc = errorCode)

   if (errorCode /= ESMF_SUCCESS) then
      write(nu_diag,*) &
         '(ice) CICE_CoupledInit: error naming export state fields'
      return
   endif

   exportState = ESMF_StateCreate(stateName='CICE Export State', &
                                  statetype = ESMF_STATE_EXPORT, &
                                  rc = errorCode)

   if (errorCode /= ESMF_SUCCESS) then
      write(nu_diag,*) &
         '(ice) CICE_CoupledInit: error creating export state'
      return
   endif

   call ESMF_StateAddArray(exportState, exportArrayIfrac, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayTsrf, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayAlbvdr, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayAlbidr, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayAlbvdf, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayAlbidf, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayTauxa, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayTauya, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayFlat, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayFsens, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayFLWout, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayEvap, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayTref, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayQref, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayFSWabs, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayFSWthru, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayFmelt, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayFresh, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayFsalt, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayTauxo, &
                           rc = errorCode)
   call ESMF_StateAddArray(exportState, exportArrayTauyo, &
                           rc = errorCode)

   if (errorCode /= ESMF_SUCCESS) then
      write(nu_diag,*) &
         '(ice) CICE_CoupledInit: error adding export state arrays'
      return
   endif

!----------------------------------------------------------------------
!  here is the proper field code once grids are working
!----------------------------------------------------------------------
!   exportFieldList(exportIndxIfrac)   = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayIfrac,     &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Ice fraction',         &
!                                        rc = errorCode)
!   exportFieldList(exportIndxTsrf)    = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayTsrf,      &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Surface temperature',  &
!                                        rc = errorCode)
!   exportFieldList(exportIndxAlbvdr)  = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayAlbvdr ,   &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Albedo: visible direct',&
!                                        rc = errorCode)
!   exportFieldList(exportIndxAlbidr)  = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayAlbidr,    &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Albedo: near-ir direct',&
!                                        rc = errorCode)
!   exportFieldList(exportIndxAlbvdf)  = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayAlbvdf,    &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Albedo: visible diffuse',&
!                                        rc = errorCode)
!   exportFieldList(exportIndxAlbidf)  = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayAlbidf,    &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Albedo: near-ir diffuse',&
!                                        rc = errorCode)
!   exportFieldList(exportIndxTauxa)   = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayTauxa,     &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Wind stress E-W',      &
!                                        rc = errorCode)
!   exportFieldList(exportIndxTauya)   = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayTauya,     &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Wind stress N-S',      &
!                                        rc = errorCode)
!   exportFieldList(exportIndxFlat)    = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayFlat,      &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Latent heat flux',     &
!                                        rc = errorCode)
!   exportFieldList(exportIndxFsens)   = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayFsens,     &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Sensible heat flux',   &
!                                        rc = errorCode)
!   exportFieldList(exportIndxFLWout)  = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayFLWout,    &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Longwave radiative heat flux',&
!                                        rc = errorCode)
!   exportFieldList(exportIndxEvap)    = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayEvap,      &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Evaporation water flux',&
!                                        rc = errorCode)
!   exportFieldList(exportIndxTref)    = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayTref ,     &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = '2m reference temperature',&
!                                        rc = errorCode)
!   exportFieldList(exportIndxQref)    = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayQref,      &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = '2m reference humidity',&
!                                        rc = errorCode)
!   exportFieldList(exportIndxFSWabs)  = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayFSWabs,    &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Short wave heat flux absorbed by ice',&
!                                        rc = errorCode)
!   exportFieldList(exportIndxFSWthru) = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayFSWthru,   &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Short wave heat flux into ocean',&
!                                        rc = errorCode)
!   exportFieldList(exportIndxFmelt)   = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayFmelt,     &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Melt heat flux',       &
!                                        rc = errorCode)
!   exportFieldList(exportIndxFresh)   = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayFresh,     &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Fresh water flux',     &
!                                        rc = errorCode)
!   exportFieldList(exportIndxFsalt)   = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayFsalt,     &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Salt flux',            &
!                                        rc = errorCode)
!   exportFieldList(exportIndxTauxo)   = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayTauxo,     &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Ocean-ice stress E-W', &
!                                        rc = errorCode)
!   exportFieldList(exportIndxTauyo)   = ESMF_FieldCreate(              &
!                                        iceGrid, exportArrayTauyo,     &
!                                        copyflag = ESMF_DATA_REF,      &
!                                        horzRelloc = ESMF_CELL_CENTER, &
!                                        haloWidth = nghost,            &
!                                        name = 'Ocean-ice stress N-S', &
!                                        rc = errorCode)
!
!   if (errorCode /= ESMF_SUCCESS) then
!      write(nu_diag,*) &
!         '(ice) CICE_CoupledInit: error creating export state fields'
!      return
!   endif
!
!   !--------------------------------------------------------------------
!   ! add units attribute to fields
!   !--------------------------------------------------------------------
!
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxIfrac  ), &
!                               'units', 'unitless', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxTsrf   ), &
!                               'units', 'deg C', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxAlbvdr ), &
!                               'units', 'unitless', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxAlbidr ), &
!                               'units', 'unitless', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxAlbvdf ), &
!                               'units', 'unitless', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxAlbidf ), &
!                               'units', 'unitless', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxTauxa  ), &
!                               'units', 'N/m2', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxTauya  ), &
!                               'units', 'N/m2', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxFlat   ), &
!                               'units', 'W/m2', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxFsens  ), &
!                               'units', 'W/m2', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxFLWout ), &
!                               'units', 'W/m2', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxEvap   ), &
!                               'units', 'kg/m2/s', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxTref   ), &
!                               'units', 'deg C', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxQref   ), &
!                               'units', ' ', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxFSWabs ), &
!                               'units', 'W/m2', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxFSWthru), &
!                               'units', 'W/m2', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxFmelt  ), &
!                               'units', 'W/m2', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxFresh  ), &
!                               'units', 'kg/m2/s', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxFsalt  ), &
!                               'units', 'kg/m2/s', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxTauxo  ), &
!                               'units', 'N/m2', rc = errorCode)
!   call ESMF_FieldSetAttribute(exportFieldList(exportIndxTauyo  ), &
!                               'units', 'N/m2', rc = errorCode)
!
!   if (errorCode /= ESMF_SUCCESS) then
!      write(nu_diag,*) &
!         '(ice) CICE_CoupledInit: error adding units to export fields'
!      return
!   endif
!
   !--------------------------------------------------------------------
   !  now that all fields are defined, set up actual import and export
   !  states
   !--------------------------------------------------------------------
!
!   importState = ESMF_StateCreate(stateName='CICE Import State', &
!                                  statetype = ESMF_STATE_IMPORT, &
!                                  fieldList = importFieldList,   &
!                                  rc = errorCode)
!
!   if (errorCode /= ESMF_SUCCESS) then
!      write(nu_diag,*) &
!         '(ice) CICE_CoupledInit: error creating import state'
!      return
!   endif
!
!   exportState = ESMF_StateCreate(stateName='CICE Export State', &
!                                  statetype = ESMF_STATE_EXPORT, &
!                                  fieldList = exportFieldList,   &
!                                  rc = errorCode)
!
!   if (errorCode /= ESMF_SUCCESS) then
!      write(nu_diag,*) &
!         '(ice) CICE_CoupledInit: error creating export state'
!      return
!   endif
!
   !--------------------------------------------------------------------
   ! fill initial export state
   !--------------------------------------------------------------------

   call CICE_CoupledAccumulateExport(errorCode)

   call CICE_CoupledFillExport(exportState, errorCode)

!-----------------------------------------------------------------------
!EOC
!

 end subroutine CICE_CoupledInit

!=======================================================================
!BOP
!
! !IROUTINE: CICE_CoupledExtractImport - extract fields from import state
!
! !INTERFACE:
!
 subroutine CICE_CoupledExtractImport(importState, errorCode)
!
! !DESCRIPTION:
!
!  Extracts fields from ESMF import state and stores them in
!  appropriate CICE arrays
!
! !REVISION HISTORY:
!
!  author: Philip W Jones, LANL
!
! !USES:
!
! !INPUT PARAMETERS:

   type (ESMF_State), intent(in) :: &
      importState        ! import state from which to extract fields

! !OUTPUT PARAMETERS:

   integer (kind=int_kind), intent(out) :: &
      errorCode        ! error code: ESMF_Success or ESMF_Failure

!
!EOP
!BOC
!
   !--------------------------------------------------------------------
   ! local variables
   !--------------------------------------------------------------------

   integer (kind=int_kind) :: &
      i,j,iblk                 ! local loop indices

   real (kind=dbl_kind) ::  &
      workx, worky,         &! temps for vector rotations
      cosT, sinT             ! rotation cosines and sines

   !--------------------------------------------------------------------
   !  import state is defined in init routine so pointers already
   !  reference correct data
   !--------------------------------------------------------------------

   errorCode = ESMF_SUCCESS

   !if (.not. l_CCSMcoupled) return

   call ice_timer_start(timer_couple)  ! time spent coupling

   !--------------------------------------------------------------------
   ! Zero stuff while waiting, only filling in active cells.
   !--------------------------------------------------------------------

   zlvl   (:,:,:) = c0
   uatm   (:,:,:) = c0
   vatm   (:,:,:) = c0
   potT   (:,:,:) = c0
   Tair   (:,:,:) = c0
   Qa     (:,:,:) = c0
   rhoa   (:,:,:) = c0
   swvdr  (:,:,:) = c0
   swvdf  (:,:,:) = c0
   swidr  (:,:,:) = c0
   swidf  (:,:,:) = c0
   flw    (:,:,:) = c0
   frain  (:,:,:) = c0
   fsnow  (:,:,:) = c0

   sst    (:,:,:) = c0
   sss    (:,:,:) = c0
   uocn   (:,:,:) = c0
   vocn   (:,:,:) = c0
   ss_tltx(:,:,:) = c0
   ss_tlty(:,:,:) = c0
   frzmlt (:,:,:) = c0

   !--------------------------------------------------------------------
   ! unpack import state
   !--------------------------------------------------------------------

   !--- unpack message

   do iblk = 1, nblocks ! nblocks must be one currently
   do j = 1, ny_block
   do i = 1, nx_block

      !--- ocn states--
      sst  (i,j,iblk) = importPtrSST (i,j)
      sss  (i,j,iblk) = importPtrSSS (i,j)
      uocn (i,j,iblk) = importPtrUocn(i,j)
      vocn (i,j,iblk) = importPtrVocn(i,j)

      !--- atm states-
      zlvl (i,j,iblk) = importPtrZlvl(i,j)
      uatm (i,j,iblk) = importPtrUatm(i,j)
      vatm (i,j,iblk) = importPtrVatm(i,j)
      potT (i,j,iblk) = importPtrPotT(i,j)
      Tair (i,j,iblk) = importPtrTair(i,j)
      Qa   (i,j,iblk) = importPtrQa  (i,j)
      rhoa (i,j,iblk) = importPtrRhoa(i,j)

      !--- ocn states--
      ss_tltx(i,j,iblk) = importPtrSSHtltx(i,j)
      ss_tlty(i,j,iblk) = importPtrSSHtlty(i,j)
      frzmlt (i,j,iblk) = importPtrFrzmlt (i,j)

      !--- atm fluxes--
      swvdr(i,j,iblk) = importPtrSWvdr(i,j)
      swidr(i,j,iblk) = importPtrSWidr(i,j)
      swvdf(i,j,iblk) = importPtrSWvdf(i,j)
      swidf(i,j,iblk) = importPtrSWidf(i,j)
      flw  (i,j,iblk) = importPtrFlw  (i,j)
      frain(i,j,iblk) = importPtrFrain(i,j)
      fsnow(i,j,iblk) = importPtrFsnow(i,j)

   enddo                    ! i
   enddo                    ! j
   enddo                       ! iblk

   !--------------------------------------------------------------------
   ! update ghost cells of import states
   !--------------------------------------------------------------------
   !
   !call update_ghost_cells(sst    , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(sss    , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(uocn   , bndy_info, field_loc_center, &
   !                                            field_type_vector)
   !call update_ghost_cells(vocn   , bndy_info, field_loc_center, &
   !                                            field_type_vector)
   !call update_ghost_cells(zlvl   , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(uatm   , bndy_info, field_loc_center, &
   !                                            field_type_vector)
   !call update_ghost_cells(vatm   , bndy_info, field_loc_center, &
   !                                            field_type_vector)
   !call update_ghost_cells(potT   , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(Tair   , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(Qa     , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(rhoa   , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(ss_tltx, bndy_info, field_loc_center, &
   !                                            field_type_vector)
   !call update_ghost_cells(ss_tlty, bndy_info, field_loc_center, &
   !                                            field_type_vector)
   !call update_ghost_cells(frzmlt , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(swvdr  , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(swidr  , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(swvdf  , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(swidf  , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(flw    , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(frain  , bndy_info, field_loc_center, &
   !                                            field_type_scalar)
   !call update_ghost_cells(fsnow  , bndy_info, field_loc_center, &
   !                                            field_type_scalar)

   !--------------------------------------------------------------------
   ! rotate zonal/meridional vectors to local coordinates
   ! Vector fields come in on T grid, but are oriented geographically
   ! need to rotate to pop-grid FIRST using ANGLET
   ! then interpolate to the U-cell centers  (otherwise we
   ! interpolate across the pole)
   ! use ANGLET which is on the T grid
   ! compute additional data derived quantities
   !--------------------------------------------------------------------


   do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block

         ! rotation quantities
         cosT = cos(ANGLET(i,j,iblk))
         sinT = sin(ANGLET(i,j,iblk))

         ! ocean
         workx      = uocn  (i,j,iblk) ! currents, m/s
         worky      = vocn  (i,j,iblk)
         uocn(i,j,iblk) = workx*cosT + worky*sinT
         vocn(i,j,iblk) = worky*cosT - workx*sinT

         workx      = ss_tltx  (i,j,iblk) ! sea sfc tilt, m/m
         worky      = ss_tlty  (i,j,iblk)
         ss_tltx(i,j,iblk) = workx*cosT + worky*sinT
         ss_tlty(i,j,iblk) = worky*cosT - workx*sinT

         ! atmosphere
         workx      = uatm(i,j,iblk) ! wind velocity, m/s
         worky      = vatm(i,j,iblk)
         uatm (i,j,iblk) = workx*cosT + worky*sinT
         vatm (i,j,iblk) = worky*cosT - workx*sinT

         ! compute some derived quantities
         sst(i,j,iblk) = sst(i,j,iblk) - Tffresh    ! sea sfc temp (C)
         Tf (i,j,iblk) = -1.8_dbl_kind              ! hardwired for NCOM
!        Tf (i,j,iblk) = -depressT*sss(i,j,iblk)   ! freezing temp (C)
!        Tf (i,j,iblk) = -depressT*max(sss(i,j,iblk),ice_ref_salinity)

         wind (i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
         fsw  (i,j,iblk) = swvdr(i,j,iblk) + swvdf(i,j,iblk) &
                         + swidr(i,j,iblk) + swidf(i,j,iblk)

      enddo                    ! i
      enddo                    ! j
   enddo                       ! iblk

   !--------------------------------------------------------------------
   ! Interpolate ocean dynamics variables from T-cell centers to
   ! U-cell centers.
   !--------------------------------------------------------------------

   call t2ugrid_vector(uocn)
   call t2ugrid_vector(vocn)
   call t2ugrid_vector(ss_tltx)
   call t2ugrid_vector(ss_tlty)

   time_forc=time !???

   call ice_timer_stop(timer_couple)   ! time spent coupling

!-----------------------------------------------------------------------
!EOC

 end subroutine CICE_CoupledExtractImport

!=======================================================================
!BOP
!
! !IROUTINE: CICE_CoupledFillExport - fills export state with fields
!
! !INTERFACE:
!
 subroutine CICE_CoupledFillExport(exportState, errorCode)
!
! !DESCRIPTION:
!
! Fills ESMF export state with CICE data
!
! !REVISION HISTORY:
!
! author: Philip W Jones, LANL
!
! !USES
!

! !OUTPUT PARAMETERS:

   integer (kind=int_kind), intent(out) :: &
      errorCode          ! output error code (success or fail)

   type (ESMF_state), intent(out) :: &
      exportState        ! export state to be filled

!
!EOP
!BOC
!
   !--------------------------------------------------------------------
   ! local variables
   !--------------------------------------------------------------------

   integer(kind=int_kind) :: i,j,iblk     ! local loop indices

   real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
      Tsrf,      &! surface temperature
      tauxa,     &! atmo/ice stress
      tauya,     &
      tauxo,     &! ice/ocean stress
      tauyo,     &
      ailohi      ! fractional ice area

   real (kind=dbl_kind) :: &
      workx, worky,        &! tmps for converting grid
      cosT, sinT

   !--------------------------------------------------------------------
   ! start timer and set flags
   !--------------------------------------------------------------------

   errorCode = ESMF_Success
   call ice_timer_start(timer_couple)  ! time spent coupling

   !--------------------------------------------------------------------
   ! manipulate some variables into proper form and rotate vectors
   ! from model grid to geographic orientation (E-W, N-S)
   !--------------------------------------------------------------------

   do iblk = 1, nblocks
   do j = 1, ny_block
   do i = 1, nx_block

      ! ice fraction
      ailohi(i,j,iblk) = exportBufferSum(i,j,iblk,exportIndxIfrac)/ &
                         timeLastCoupled

      ! surface temperature convert to Kelvin
      Tsrf(i,j,iblk)  = Tffresh + &
              exportBufferSum(i,j,iblk,exportIndxTsrf)/timeLastCoupled

      ! rotation factors
      cosT = cos(ANGLET(i,j,iblk))
      sinT = sin(ANGLET(i,j,iblk))

      ! rotate wind stress to geographic directions (EW,NS)
      workx = exportBufferSum(i,j,iblk,exportIndxTauxa)/timeLastCoupled
      worky = exportBufferSum(i,j,iblk,exportIndxTauya)/timeLastCoupled
      tauxa(i,j,iblk) = workx*cosT - worky*sinT
      tauya(i,j,iblk) = worky*cosT + workx*sinT

      ! rotate ice/ocean stress to geographic directions (EW,NS)
      workx = -exportBufferSum(i,j,iblk,exportIndxTauxo)/timeLastCoupled
      worky = -exportBufferSum(i,j,iblk,exportIndxTauyo)/timeLastCoupled
      tauxo(i,j,iblk) = workx*cosT - worky*sinT
      tauyo(i,j,iblk) = worky*cosT + workx*sinT

   enddo                    ! i
   enddo                    ! j
   enddo                    ! iblk

   !--------------------------------------------------------------------
   ! check for bad ice fractions
   !--------------------------------------------------------------------

   if (count(tmask .and. ailohi < c0) > 1) then
      write(nu_diag,*) &
         '(ice) CICE_CoupledFillExport: ERROR ailohi < 0.'
      errorCode = ESMF_Failure
      return
   endif

   !--------------------------------------------------------------------
   ! pack data into ESMF export state
   !--------------------------------------------------------------------

   !
   ! initialize export state with special values
   !

   exportPtrIfrac   = spval_dbl
   exportPtrTsrf    = spval_dbl
   exportPtrAlbvdr  = spval_dbl
   exportPtrAlbidr  = spval_dbl
   exportPtrAlbvdf  = spval_dbl
   exportPtrAlbidf  = spval_dbl
   exportPtrTauxa   = spval_dbl
   exportPtrTauya   = spval_dbl
   exportPtrFlat    = spval_dbl
   exportPtrFsens   = spval_dbl
   exportPtrFLWout  = spval_dbl
   exportPtrEvap    = spval_dbl
   exportPtrTref    = spval_dbl
   exportPtrQref    = spval_dbl
   exportPtrFSWabs  = spval_dbl
   exportPtrFSWthru = spval_dbl
   exportPtrFmelt   = spval_dbl
   exportPtrFresh   = spval_dbl
   exportPtrFsalt   = spval_dbl
   exportPtrTauxo   = spval_dbl
   exportPtrTauyo   = spval_dbl

   do iblk = 1, nblocks
   do j = 1, ny_block
   do i = 1, nx_block

      !--- ice states
      exportPtrIfrac(i,j) = ailohi(i,j,iblk)

      !--- only fill the rest where there is ice

      if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then
         exportPtrTsrf(i,j) = Tsrf(i,j,iblk) ! temperature

         exportPtrAlbvdr(i,j) = &
         exportBufferSum(i,j,iblk,exportIndxAlbvdr)/timeLastCoupled
         exportPtrAlbidr(i,j) = &
         exportBufferSum(i,j,iblk,exportIndxAlbidr)/timeLastCoupled
         exportPtrAlbvdf(i,j) = &
         exportBufferSum(i,j,iblk,exportIndxAlbvdf)/timeLastCoupled
         exportPtrAlbidf(i,j) = &
         exportBufferSum(i,j,iblk,exportIndxAlbidf)/timeLastCoupled

         !--- a/i fluxes computed by ice
         exportPtrTauxa(i,j) = tauxa(i,j,iblk)
         exportPtrTauya(i,j) = tauya(i,j,iblk)
         exportPtrFlat  (i,j) = &
         exportBufferSum(i,j,iblk,exportIndxFlat  )/timeLastCoupled
         exportPtrFsens (i,j) = &
         exportBufferSum(i,j,iblk,exportIndxFsens )/timeLastCoupled
         exportPtrFLWout(i,j) = &
         exportBufferSum(i,j,iblk,exportIndxFLWout)/timeLastCoupled
         exportPtrEvap  (i,j) = &
         exportBufferSum(i,j,iblk,exportIndxEvap  )/timeLastCoupled
         exportPtrTref  (i,j) = &
         exportBufferSum(i,j,iblk,exportIndxTref  )/timeLastCoupled
         exportPtrQref  (i,j) = &
         exportBufferSum(i,j,iblk,exportIndxQref  )/timeLastCoupled
         exportPtrFSWabs(i,j) = &
         exportBufferSum(i,j,iblk,exportIndxFSWabs)/timeLastCoupled

         !--- i/o fluxes computed by ice
         exportPtrFSWthru(i,j) = &
         exportBufferSum (i,j,iblk,exportIndxFSWthru)/timeLastCoupled
         exportPtrFmelt  (i,j) = &
         exportBufferSum (i,j,iblk,exportIndxFmelt  )/timeLastCoupled
         exportPtrFresh  (i,j) = &
         exportBufferSum (i,j,iblk,exportIndxFresh  )/timeLastCoupled
         exportPtrFsalt  (i,j) = &
         exportBufferSum (i,j,iblk,exportIndxFsalt  )/timeLastCoupled
         exportPtrTauxo  (i,j) = tauxo(i,j,iblk)
         exportPtrTauyo  (i,j) = tauyo(i,j,iblk)
      endif  ! tmask and ailohi > c0
   enddo                    ! i
   enddo                    ! j
   enddo                     ! iblk

   call ice_timer_stop(timer_couple)    ! time spent coupling

!-----------------------------------------------------------------------
!EOC

 end subroutine CICE_CoupledFillExport

!=======================================================================
!BOP
!
! !IROUTINE: CICE_CoupledAccumulateExport - accumulates time avg fields
!
! !INTERFACE:
!
 subroutine CICE_CoupledAccumulateExport(errorCode)
!
! !DESCRIPTION:
!
! Accumulates time averages of export state fields for the case when
! the coupling frequency is longer than the internal CICE time step.
!
! !REVISION HISTORY:
!
! author: Philip W Jones, LANL
!
! !USES
!
!
! !OUTPUT PARAMETERS:
!
   integer (kind=int_kind), intent(out) :: &
      errorCode              ! output error code

!
!EOP
!BOC
!

   !--------------------------------------------------------------------
   ! local variables
   !--------------------------------------------------------------------

   integer(kind=int_kind) :: i,j,iblk     ! local loop indices

   !--------------------------------------------------------------------
   ! zero buffer if this is the first time after a coupling interval
   !--------------------------------------------------------------------

   if (timeLastCoupled == c0) exportBufferSum = c0

   !--------------------------------------------------------------------
   ! accumulate sums for time averages
   !--------------------------------------------------------------------

   do iblk = 1, nblocks
   do j = 1, ny_block
   do i = 1, nx_block

      exportBufferSum(i,j,iblk,exportIndxIfrac  ) = &
      exportBufferSum(i,j,iblk,exportIndxIfrac  ) + dt*aice(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxTsrf   ) = &
      exportBufferSum(i,j,iblk,exportIndxTsrf   ) + dt*trcr(i,j,1,iblk)
      exportBufferSum(i,j,iblk,exportIndxAlbvdr ) = &
      exportBufferSum(i,j,iblk,exportIndxAlbvdr ) + dt*alvdr(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxAlbidr ) = &
      exportBufferSum(i,j,iblk,exportIndxAlbidr ) + dt*alidr(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxAlbvdf ) = &
      exportBufferSum(i,j,iblk,exportIndxAlbvdf ) + dt*alvdf(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxAlbidf ) = &
      exportBufferSum(i,j,iblk,exportIndxAlbidf ) + dt*alidf(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxTauxa  ) = &
      exportBufferSum(i,j,iblk,exportIndxTauxa  ) + dt*strairxT(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxTauya  ) = &
      exportBufferSum(i,j,iblk,exportIndxTauya  ) + dt*strairyT(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxFlat   ) = &
      exportBufferSum(i,j,iblk,exportIndxFlat   ) + dt*flat(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxFsens  ) = &
      exportBufferSum(i,j,iblk,exportIndxFsens  ) + dt*fsens(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxFLWout ) = &
      exportBufferSum(i,j,iblk,exportIndxFLWout ) + dt*flwout(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxEvap   ) = &
      exportBufferSum(i,j,iblk,exportIndxEvap   ) + dt*evap(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxTref   ) = &
      exportBufferSum(i,j,iblk,exportIndxTref   ) + dt*Tref(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxQref   ) = &
      exportBufferSum(i,j,iblk,exportIndxQref   ) + dt*Qref(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxFSWabs ) = &
      exportBufferSum(i,j,iblk,exportIndxFSWabs ) + dt*fswabs(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxFSWthru) = &
      exportBufferSum(i,j,iblk,exportIndxFSWthru) + dt*fswthru(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxFmelt  ) = &
      exportBufferSum(i,j,iblk,exportIndxFmelt  ) + dt*fhocn(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxFresh  ) = &
      exportBufferSum(i,j,iblk,exportIndxFresh  ) + dt*fresh(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxFsalt  ) = &
      exportBufferSum(i,j,iblk,exportIndxFsalt  ) + dt*fsalt(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxTauxo  ) = &
      exportBufferSum(i,j,iblk,exportIndxTauxo  ) + dt*strocnxT(i,j,iblk)
      exportBufferSum(i,j,iblk,exportIndxTauyo  ) = &
      exportBufferSum(i,j,iblk,exportIndxTauyo  ) + dt*strocnyT(i,j,iblk)

   enddo                    ! i
   enddo                    ! j
   enddo                     ! iblk

   !--------------------------------------------------------------------
   !  update time since last coupling
   !--------------------------------------------------------------------

   timeLastCoupled = timeLastCoupled + dt

!-----------------------------------------------------------------------
!EOC

 end subroutine CICE_CoupledAccumulateExport
!=======================================================================

#endif   
! ifdef esmf (previous four subroutines)

#ifdef CCSM
! next five subroutines
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: ice_coupling_setup - sets mpi communicators and task ids
!
! !INTERFACE:
!
      subroutine ice_coupling_setup(in_model_name,model_comm)
!
! !DESCRIPTION:
!
! This routine uses get the model communicator from ccsm share code
!
! !REVISION HISTORY:
!
! author: T. Craig, NCAR, Dec 30, 2002: for cpl6
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (3), intent(in) :: in_model_name   

      integer, intent(out) ::  model_comm     ! communicator for model
!
!EOP
!
      write(nu_diag,*) 'calling cpl_interface_init for model: ', &
           in_model_name,' ', trim(cpl_fields_icename)

      call cpl_interface_init(cpl_fields_icename,model_comm)

      call shr_sys_flush(nu_diag)

      end subroutine ice_coupling_setup

!=======================================================================
!BOP
!
! !IROUTINE: init_cpl - initializes message passing between ice and coupler
!
! !INTERFACE:
!
      subroutine init_cpl
!
! !DESCRIPTION:
!
! Initializes message passing between ice and coupler
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_broadcast

! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer(kind=int_kind) :: i,j,n,iblk      ! local loop indices
      integer(kind=int_kind) :: ilo,ihi,jlo,jhi ! local loop indices
      integer(kind=int_kind) :: iyear     ! yyyy

      write(nu_diag,*) &
           '(ice_coupling,init_cpl) send initial msg. set contract'
      call shr_sys_flush(nu_diag)

      nadv_i = nint(secday/dt)

      isbuf                          = 0         ! default info-buffer value
      isbuf(cpl_fields_ibuf_cdate  ) = idate     ! initial date (coded: yyyymmdd)
      isbuf(cpl_fields_ibuf_sec    ) = sec       ! elapsed seconds into date
      isbuf(cpl_fields_ibuf_stopnow) = stop_now  ! stop now flag
      isbuf(cpl_fields_ibuf_userest) = 0         ! use model restart data initally
      isbuf(cpl_fields_ibuf_ncpl   ) = nadv_i    ! number of comms per day
      isbuf(cpl_fields_ibuf_lsize  ) &
         = (nx_block-2*nghost)*(ny_block-2*nghost)*max_blocks
      isbuf(cpl_fields_ibuf_lisize ) = nx_block-2*nghost ! local size wrt i-index
      isbuf(cpl_fields_ibuf_ljsize ) = ny_block-2*nghost ! local size wrt j-index
      isbuf(cpl_fields_ibuf_gsize  ) = nx_global*ny_global ! size of global grid
      isbuf(cpl_fields_ibuf_gisize ) = nx_global  ! global size wrt i-index
      isbuf(cpl_fields_ibuf_gjsize ) = ny_global  ! global size wrt j-index
      isbuf(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
      isbuf(cpl_fields_ibuf_dead   ) = 0           ! not a dead model

      allocate(sbuf((nx_block-2*nghost)*(ny_block-2*nghost), &
         cpl_fields_grid_total))
      sbuf = -888.0_dbl_kind
      n=0
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo,jhi
         do i = ilo,ihi
            n=n+1
            sbuf(n,cpl_fields_grid_lon  ) = TLON(i,j,iblk)*rad_to_deg
            sbuf(n,cpl_fields_grid_lat  ) = TLAT(i,j,iblk)*rad_to_deg
            sbuf(n,cpl_fields_grid_area ) = tarea(i,j,iblk)/(radius*radius)
            sbuf(n,cpl_fields_grid_mask ) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
            sbuf(n,cpl_fields_grid_index) = rndex_global(i,j,iblk)
         enddo
         enddo
      enddo

      call cpl_interface_contractInit &
           (contractS, cpl_fields_icename, cpl_fields_cplname, &
            cpl_fields_i2c_fields, isbuf, sbuf)

      call cpl_interface_contractInit &
           (contractR, cpl_fields_icename, cpl_fields_cplname, &
            cpl_fields_c2i_fields, isbuf, sbuf)

      write(nu_diag,*) '(init_cpl) Initialized contracts with coupler'
      call shr_sys_flush(nu_diag)

      !-----------------------------------------------------------------
      ! Allocate memory for send and receive buffers
      !-----------------------------------------------------------------

      nsend = cpl_interface_contractNumatt(contractS)
      allocate(buffs((nx_block-2*nghost)*(ny_block-2*nghost)*max_blocks, &
         nsend))

      nrecv = cpl_interface_contractNumatt(contractR)
      allocate(buffr((nx_block-2*nghost)*(ny_block-2*nghost)*max_blocks, &
         nrecv))

      !-----------------------------------------------------------------
      ! Determine send indices
      !-----------------------------------------------------------------

      index_i2c_index      = cpl_interface_contractIndex(contractS,'index')
      index_i2c_Si_ifrac   = cpl_interface_contractIndex(contractS,'Si_ifrac')
      index_i2c_Si_t       = cpl_interface_contractIndex(contractS,'Si_t')
      index_i2c_Si_avsdr   = cpl_interface_contractIndex(contractS,'Si_avsdr')
      index_i2c_Si_anidr   = cpl_interface_contractIndex(contractS,'Si_anidr')
      index_i2c_Si_avsdf   = cpl_interface_contractIndex(contractS,'Si_avsdf')
      index_i2c_Si_anidf   = cpl_interface_contractIndex(contractS,'Si_anidf')
      index_i2c_Si_tref    = cpl_interface_contractIndex(contractS,'Si_tref')
      index_i2c_Si_qref    = cpl_interface_contractIndex(contractS,'Si_qref')
      index_i2c_Faii_taux  = cpl_interface_contractIndex(contractS,'Faii_taux')
      index_i2c_Faii_tauy  = cpl_interface_contractIndex(contractS,'Faii_tauy')
      index_i2c_Faii_lat   = cpl_interface_contractIndex(contractS,'Faii_lat')
      index_i2c_Faii_sen   = cpl_interface_contractIndex(contractS,'Faii_sen')
      index_i2c_Faii_lwup  = cpl_interface_contractIndex(contractS,'Faii_lwup')
      index_i2c_Faii_evap  = cpl_interface_contractIndex(contractS,'Faii_evap')
      index_i2c_Faii_swnet = cpl_interface_contractIndex(contractS,'Faii_swnet')
      index_i2c_Fioi_swpen = cpl_interface_contractIndex(contractS,'Fioi_swpen')
      index_i2c_Fioi_melth = cpl_interface_contractIndex(contractS,'Fioi_melth')
      index_i2c_Fioi_meltw = cpl_interface_contractIndex(contractS,'Fioi_meltw')
      index_i2c_Fioi_salt  = cpl_interface_contractIndex(contractS,'Fioi_salt')
      index_i2c_Fioi_taux  = cpl_interface_contractIndex(contractS,'Fioi_taux')
      index_i2c_Fioi_tauy  = cpl_interface_contractIndex(contractS,'Fioi_tauy')

      !-----------------------------------------------------------------
      ! Determine receive indices
      !-----------------------------------------------------------------

      index_c2i_So_t       = cpl_interface_contractIndex(contractR,'So_t')
      index_c2i_So_s       = cpl_interface_contractIndex(contractR,'So_s')
      index_c2i_So_u       = cpl_interface_contractIndex(contractR,'So_u')
      index_c2i_So_v       = cpl_interface_contractIndex(contractR,'So_v')
      index_c2i_Sa_z       = cpl_interface_contractIndex(contractR,'Sa_z')
      index_c2i_Sa_u       = cpl_interface_contractIndex(contractR,'Sa_u')
      index_c2i_Sa_v       = cpl_interface_contractIndex(contractR,'Sa_v')
      index_c2i_Sa_ptem    = cpl_interface_contractIndex(contractR,'Sa_ptem')
      index_c2i_Sa_tbot    = cpl_interface_contractIndex(contractR,'Sa_tbot')
      index_c2i_Sa_shum    = cpl_interface_contractIndex(contractR,'Sa_shum')
      index_c2i_Sa_dens    = cpl_interface_contractIndex(contractR,'Sa_dens')
      index_c2i_So_dhdx    = cpl_interface_contractIndex(contractR,'So_dhdx')
      index_c2i_So_dhdy    = cpl_interface_contractIndex(contractR,'So_dhdy')
      index_c2i_Fioo_q     = cpl_interface_contractIndex(contractR,'Fioo_q')
      index_c2i_Faxa_swvdr = cpl_interface_contractIndex(contractR,'Faxa_swvdr')
      index_c2i_Faxa_swndr = cpl_interface_contractIndex(contractR,'Faxa_swndr')
      index_c2i_Faxa_swvdf = cpl_interface_contractIndex(contractR,'Faxa_swvdf')
      index_c2i_Faxa_swndf = cpl_interface_contractIndex(contractR,'Faxa_swndf')
      index_c2i_Faxa_lwdn  = cpl_interface_contractIndex(contractR,'Faxa_lwdn')
      index_c2i_Faxc_rain  = cpl_interface_contractIndex(contractR,'Faxc_rain')
      index_c2i_Faxc_snow  = cpl_interface_contractIndex(contractR,'Faxc_snow')

      !-----------------------------------------------------------------
      ! Receive initial message from coupler.
      !-----------------------------------------------------------------

      call cpl_interface_ibufRecv(cpl_fields_cplname,irbuf)

      if (my_task==master_task) then
         write(nu_diag,*) &
              '(init_cpl) Received control buffer from coupler'
         call shr_sys_flush(nu_diag)

         if (trim(runtype)=='startup' .or. &
             trim(runtype)== 'hybrid') then
            idate = irbuf(cpl_fields_ibuf_cdate)
            write(nu_diag,*) '(init_cpl) idate from coupler = ',idate
            iyear   = (idate/10000)               ! integer year of basedate
            month = (idate-iyear*10000)/100       ! integer month of basedate
            mday  = idate-iyear*10000-month*100-1 ! day of year of basedate
            time  = ((iyear-1)*daycal(13)+daycal(month)+mday)*secday
            call calendar(time)                 ! recompute calendar info
            time_forc = time
            call shr_sys_flush(nu_diag)
         endif

      endif                     ! my_task==master_task

      call broadcast_scalar(idate,    master_task)
      call broadcast_scalar(time,     master_task)
      call broadcast_scalar(time_forc,master_task)

      deallocate(sbuf)

      write(nu_diag,*) '(ice_coupling,init_cpl) done setting contract'

      !-----------------------------------------------------------------
      ! Send initial state info to coupler.
      !-----------------------------------------------------------------

      call to_coupler

      end subroutine init_cpl

!=======================================================================
!BOP
!
! !IROUTINE: from_coupler - input from coupler to sea ice model
!
! !INTERFACE:
!
      subroutine from_coupler
!
! !DESCRIPTION:
!
! Reads input data from coupler to sea ice model
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_global_reductions
      use ice_work, only: work1

! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i,j,n,n2,iblk   ! local loop indices
      integer(kind=int_kind)  :: ilo,ihi,jlo,jhi ! local loop indices
      integer (kind=int_kind) :: index           ! attribute vector property

      real (kind=dbl_kind) :: &
         gsum, workx, worky

      call ice_timer_start(timer_couple)  ! time spent coupling

      !-----------------------------------------------------------------
      ! Zero stuff while waiting, only filling in active cells.
      !-----------------------------------------------------------------

      zlvl   (:,:,:) = c0
      uatm   (:,:,:) = c0
      vatm   (:,:,:) = c0
      potT   (:,:,:) = c0
      Tair   (:,:,:) = c0
      Qa     (:,:,:) = c0
      rhoa   (:,:,:) = c0
      swvdr  (:,:,:) = c0
      swvdf  (:,:,:) = c0
      swidr  (:,:,:) = c0
      swidf  (:,:,:) = c0
      flw    (:,:,:) = c0
      frain  (:,:,:) = c0
      fsnow  (:,:,:) = c0
      sst    (:,:,:) = c0
      sss    (:,:,:) = c0
      uocn   (:,:,:) = c0
      vocn   (:,:,:) = c0
      ss_tltx(:,:,:) = c0
      ss_tlty(:,:,:) = c0
      frzmlt (:,:,:) = c0

      !-----------------------------------------------------------------
      ! recv input field msg
      !-----------------------------------------------------------------
     
      call ice_timer_start(timer_cplrecv)  ! time spent receiving

      call cpl_interface_contractRecv &
           (cpl_fields_cplname, contractR, irbuf, buffr)

      call ice_timer_stop(timer_cplrecv)
!     call ice_timer_start(17)  ! time spent cr-unpacking

      !--- unpack message
      n=0
      do iblk = 1, nblocks

       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
 
       do j = jlo,jhi
       do i = ilo,ihi

          n=n+1

          !--- ocn states--

          sst(i,j,iblk)  = buffr(n,index_c2i_So_t)
          sss(i,j,iblk)  = buffr(n,index_c2i_So_s)
          uocn(i,j,iblk) = buffr(n,index_c2i_So_u)
          vocn(i,j,iblk) = buffr(n,index_c2i_So_v)

          !--- atm states-

          zlvl (i,j,iblk) = buffr(n,index_c2i_Sa_z)
          uatm (i,j,iblk) = buffr(n,index_c2i_Sa_u)
          vatm (i,j,iblk) = buffr(n,index_c2i_Sa_v)
          potT (i,j,iblk) = buffr(n,index_c2i_Sa_ptem)
          Tair (i,j,iblk) = buffr(n,index_c2i_Sa_tbot)
          Qa   (i,j,iblk) = buffr(n,index_c2i_Sa_shum)
          rhoa (i,j,iblk) = buffr(n,index_c2i_Sa_dens)

          !--- ocn states--

          ss_tltx(i,j,iblk) = buffr(n,index_c2i_So_dhdx)
          ss_tlty(i,j,iblk) = buffr(n,index_c2i_So_dhdy)
          frzmlt (i,j,iblk) = buffr(n,index_c2i_Fioo_q)

          !--- atm fluxes--

          swvdr(i,j,iblk) = buffr(n,index_c2i_Faxa_swvdr)
          swidr(i,j,iblk) = buffr(n,index_c2i_Faxa_swndr)
          swvdf(i,j,iblk) = buffr(n,index_c2i_Faxa_swvdf)
          swidf(i,j,iblk) = buffr(n,index_c2i_Faxa_swndf)
          flw  (i,j,iblk) = buffr(n,index_c2i_Faxa_lwdn)
          frain(i,j,iblk) = buffr(n,index_c2i_Faxc_rain)
          fsnow(i,j,iblk) = buffr(n,index_c2i_Faxc_snow)

       end do
       end do
      end do

      call update_ghost_cells(sst    , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(sss    , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(uocn   , bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(vocn   , bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(zlvl   , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(uatm   , bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(vatm   , bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(potT   , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(Tair   , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(Qa     , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(rhoa   , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(ss_tltx, bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(ss_tlty, bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(frzmlt , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(swvdr  , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(swidr  , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(swvdf  , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(swidf  , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(flw    , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(frain  , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(fsnow  , bndy_info, field_loc_center, &
                                                  field_type_scalar)

!     call ice_timer_stop(17)  ! time spent cr-unpacking

      !-----------------------------------------------------------------
      ! broadcast dbug diagnostic level
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_infobug) >= 2 ) then
         if (my_task == master_task) write (nu_diag,*) &
              '(from_coupler) dbug level >= 2'
         info_dbug = 1
      endif

      !-----------------------------------------------------------------
      ! broadcast write_restart flag
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_resteod) == 1 .AND. new_day) then
         if (my_task == master_task) write (nu_diag,*) &
              '(from_coupler) received write restart signal'
         write_restart = 1
      endif

      !-----------------------------------------------------------------
      ! broadcast cpl_write_history flag
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_histeod) == 1 .AND. new_day) then
         if (my_task == master_task) write (nu_diag,*) &
              '(from_coupler) received write history signal'
         cpl_write_history = 1
      endif

      !-----------------------------------------------------------------
      ! broadcast stop_now flag
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_stopnow) == 1) then
         if (my_task==master_task) write (nu_diag,*) &
              '(from_coupler) received terminate signal'
         stop_now = 1
      endif

      if (info_dbug == 1 .AND. stop_now /= 1) then

        do n=1,nrecv
           work1 = c0
           n2   = 0
           do iblk = 1, nblocks

              this_block = get_block(blocks_ice(iblk),iblk)         
              ilo = this_block%ilo
              ihi = this_block%ihi
              jlo = this_block%jlo
              jhi = this_block%jhi

              do j = jlo,jhi
              do i = ilo,ihi
                 n2 = n2 + 1 
                 if (hm(i,j,iblk) > p5) work1(i,j,iblk) = buffr(n2,n)
              enddo
              enddo
           enddo

           call update_ghost_cells (work1,            bndy_info, &
                                    field_loc_center, field_type_scalar)

           gsum = global_sum(work1, distrb_info, field_loc_center, &
                             tarea)

           if (my_task == master_task) then
              write (nu_diag,100) 'ice', 'recv', n, gsum
           endif
        enddo                   ! nrecv
      endif
 100  format ('comm_diag',1x,a3,1x,a4,1x,i3,es26.19)

      !-----------------------------------------------------------------
      ! rotate zonal/meridional vectors to local coordinates
      ! compute data derived quantities
      !-----------------------------------------------------------------

      ! Vector fields come in on T grid, but are oriented geographically
      ! need to rotate to pop-grid FIRST using ANGLET
      ! then interpolate to the U-cell centers  (otherwise we
      ! interpolate across the pole)
      ! use ANGLET which is on the T grid !

      do iblk = 1, nblocks

       do j = 1,ny_block
       do i = 1,nx_block

          ! ocean
          workx      = uocn  (i,j,iblk) ! currents, m/s 
          worky      = vocn  (i,j,iblk)
          uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                         + worky*sin(ANGLET(i,j,iblk))
          vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                         - workx*sin(ANGLET(i,j,iblk))

          workx      = ss_tltx  (i,j,iblk)           ! sea sfc tilt, m/m
          worky      = ss_tlty  (i,j,iblk)
          ss_tltx(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                            + worky*sin(ANGLET(i,j,iblk))
          ss_tlty(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                            - workx*sin(ANGLET(i,j,iblk))

          sst(i,j,iblk) = sst(i,j,iblk) - Tffresh       ! sea sfc temp (C)
          Tf (i,j,iblk) = -1.8_dbl_kind                 ! hardwired for NCOM
!         Tf (i,j,iblk) = -depressT*sss(i,j,iblk)       ! freezing temp (C)
!         Tf (i,j,iblk) = -depressT*max(sss(i,j,iblk),ice_ref_salinity)

       enddo
       enddo
      enddo

      ! Interpolate ocean dynamics variables from T-cell centers to 
      ! U-cell centers.

      call t2ugrid_vector(uocn)
      call t2ugrid_vector(vocn)
      call t2ugrid_vector(ss_tltx)
      call t2ugrid_vector(ss_tlty)

      ! Atmosphere variables are needed in T cell centers in
      ! subroutine stability and are interpolated to the U grid
      ! later as necessary.

      do iblk = 1, nblocks
       do j = 1, ny_block
       do i = 1, nx_block
         ! atmosphere
         workx      = uatm(i,j,iblk) ! wind velocity, m/s
         worky      = vatm(i,j,iblk) 
         uatm (i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                         + worky*sin(ANGLET(i,j,iblk)) ! note uatm, vatm, wind
         vatm (i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & !  are on the T-grid here
                         - workx*sin(ANGLET(i,j,iblk))

         wind (i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
         fsw  (i,j,iblk) = swvdr(i,j,iblk) + swvdf(i,j,iblk) &
                         + swidr(i,j,iblk) + swidf(i,j,iblk)
       enddo
       enddo
      enddo

      time_forc=time

      call ice_timer_stop(timer_couple)   ! time spent coupling

      end subroutine from_coupler

!=======================================================================
!BOP
!
! !IROUTINE: to_coupler - send data from sea ice model to coupler
!
! !INTERFACE:
!
      subroutine to_coupler
!
! !DESCRIPTION:
!
! Sea ice model to coupler
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
      use ice_global_reductions
      use ice_work, only: work1

! !INPUT/OUTPUT PARAMETERS:
!
!EOP

      integer(kind=int_kind) :: i,j,n,n2,iblk    ! local loop indices
      integer(kind=int_kind) :: ilo,ihi,jlo,jhi ! local loop indices
      integer(kind=int_kind) :: index        ! attribute vector property

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         Tsrf  &      ! surface temperature
      ,  tauxa &     ! atmo/ice stress
      ,  tauya &
      ,  tauxo &       ! ice/ocean stress
      ,  tauyo &
      ,  ailohi        ! fractional ice area

      real (kind=dbl_kind) :: &
         gsum, workx, worky           ! tmps for converting grid

      logical :: flag
      flag=.false.

      call ice_timer_start(timer_couple)  ! time spent coupling

      do iblk = 1, nblocks
       do j = 1, ny_block
       do i = 1, nx_block

        ! ice fraction
        ailohi(i,j,iblk) = aice(i,j,iblk)

        ! surface temperature
        Tsrf(i,j,iblk)  = Tffresh + trcr(i,j,1,iblk)             !Kelvin

        ! wind stress  (on POP T-grid:  convert to lat-lon)
        workx = strairxT(i,j,iblk)                             ! N/m^2
        worky = strairyT(i,j,iblk)                             ! N/m^2
        tauxa(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                        - worky*sin(ANGLET(i,j,iblk))
        tauya(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                        + workx*sin(ANGLET(i,j,iblk))

        ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
        workx = -strocnxT(i,j,iblk)                            ! N/m^2
        worky = -strocnyT(i,j,iblk)                            ! N/m^2
        tauxo(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                        - worky*sin(ANGLET(i,j,iblk))
        tauyo(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                        + workx*sin(ANGLET(i,j,iblk))

       enddo
       enddo
      enddo

      !--- set info buffer flags ---
      isbuf                          = 0       ! unused
      isbuf(cpl_fields_ibuf_stopnow) = 0       ! stop flag: 0 <=> able to continue
      isbuf(cpl_fields_ibuf_cdate)   = idate   ! model date, coded: yyyymmdd
      isbuf(cpl_fields_ibuf_sec)     = sec     ! elapsed seconds on model date
      isbuf(cpl_fields_ibuf_lisize)  = nx_block-2*nghost  ! local size wrt i-index
      isbuf(cpl_fields_ibuf_ljsize)  = ny_block-2*nghost  ! local size wrt j-index
      isbuf(cpl_fields_ibuf_ncpl)    = nadv_i  ! number of msg-pairs per day
      isbuf(cpl_fields_ibuf_lsize) &
         = (nx_block-2*nghost)*(ny_block-2*nghost)*max_blocks 
      isbuf(cpl_fields_ibuf_dead)    = 0       ! not a dead model

!     call ice_timer_start(18)      ! Time spent packing

      !--- pack & send msg buffer ---

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
               flag = .true.
            endif
         end do
         end do
      end do
      if (flag) then
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                 write(nu_diag,*) &
                  ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j,iblk)
                 call shr_sys_flush(nu_diag)
               endif
            end do
            end do
         end do
      endif

      buffs(:,:)=spval_dbl
      n=0
      do iblk = 1, nblocks

       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo,jhi
       do i = ilo,ihi

         n=n+1

         !--- ice states
         buffs(n,index_i2c_Si_ifrac) = ailohi(i,j,iblk)     ! frac 

         !---  compression index
         buffs(n,index_i2c_index) = rndex_global(i,j,iblk)  ! global index

         if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then
            buffs(n,index_i2c_Si_t      ) = Tsrf(i,j,iblk)  ! temperature

            buffs(n,index_i2c_Si_avsdr  ) = alvdr(i,j,iblk) ! alb, visible/dir
            buffs(n,index_i2c_Si_anidr  ) = alidr(i,j,iblk) ! alb, near-ir/dir
            buffs(n,index_i2c_Si_avsdf  ) = alvdf(i,j,iblk) ! alb, visible/dif
            buffs(n,index_i2c_Si_anidf  ) = alidf(i,j,iblk) ! alb, near-ir/dif
            buffs(n,index_i2c_Si_tref   ) = Tref(i,j,iblk)  ! diagnostic: 2m ref temp
            buffs(n,index_i2c_Si_qref   ) = Qref(i,j,iblk)  ! diagnostic: 2m ref sp hum
            !--- a/i fluxes computed by ice
            buffs(n,index_i2c_Faii_taux ) = tauxa(i,j,iblk) ! stress: a/i zonal
            buffs(n,index_i2c_Faii_tauy ) = tauya(i,j,iblk) ! stress: a/i meridional
            buffs(n,index_i2c_Faii_lat  ) = flat(i,j,iblk)  ! latent heat flux
            buffs(n,index_i2c_Faii_sen  ) = fsens(i,j,iblk) ! sensible heat flux
            buffs(n,index_i2c_Faii_lwup ) = flwout(i,j,iblk)! upward longwave flux
            buffs(n,index_i2c_Faii_evap ) = evap(i,j,iblk)  ! evaporation h2o flux
            buffs(n,index_i2c_Faii_swnet) = fswabs(i,j,iblk)! sw net absorbed hf
            !--- i/o fluxes computed by ice
            buffs(n,index_i2c_Fioi_swpen) = fswthru(i,j,iblk) ! solar thru ice to ocn hf
            buffs(n,index_i2c_Fioi_melth) = fhocn(i,j,iblk) ! hf from melting
            buffs(n,index_i2c_Fioi_meltw) = fresh(i,j,iblk) ! h2o flux from melting
            buffs(n,index_i2c_Fioi_salt ) = fsalt(i,j,iblk) ! salt flux from melting
            buffs(n,index_i2c_Fioi_taux ) = tauxo(i,j,iblk) ! stress : i/o zonal
            buffs(n,index_i2c_Fioi_tauy ) = tauyo(i,j,iblk) ! stress : i/o meridional
         endif  ! tmask and ailohi > c0
       end do
       end do
      end do

!     call ice_timer_stop(18)      ! Time spent packing

      call ice_timer_start(timer_cplsend)     ! Time spent sending
      call cpl_interface_contractSend &
           (cpl_fields_cplname, contractS, isbuf, buffs)
      call ice_timer_stop(timer_cplsend)      ! Time spent sending

      !-----------------------------------------------------------------
      ! diagnostics
      !-----------------------------------------------------------------

      if (info_dbug == 1 .AND. stop_now /= 1) then

         do n=1,nsend
            work1 = c0
            n2 = 0
            do iblk = 1, nblocks

               this_block = get_block(blocks_ice(iblk),iblk)         
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

               do j = jlo,jhi
               do i = ilo,ihi
                  n2 = n2 + 1
                  if (ailohi(i,j,iblk) > c0 .and. &
                      ailohi(i,j,iblk) <= c1) then
                     work1(i,j,iblk) = buffs(n2,n)
                  endif
               enddo
               enddo
            enddo

            call update_ghost_cells (work1,            bndy_info, &
                                     field_loc_center, field_type_scalar)

            gsum = global_sum(work1, distrb_info, field_loc_center, &
                              tarea)

            if (my_task == master_task) then
               write (nu_diag,100) 'ice','send',n,gsum
            endif
         enddo
      endif
 100  format('comm_diag',1x,a3,1x,a4,1x,i3,es26.19)
      call shr_sys_flush(nu_diag)

      call ice_timer_stop(timer_couple)    ! time spent coupling

      end subroutine to_coupler

!=======================================================================
!BOP
!
! !IROUTINE: exit_coupler - exit from coupled/mpi environment
!
! !INTERFACE:
!
      subroutine exit_coupler
!
! !DESCRIPTION:
!
! Exit from coupled/MPI environment
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      include "mpif.h"         ! MPI library definitions

      integer (kind=int_kind) :: &
        ierr  ! error flag

      if (my_task == master_task) then
         if (irbuf(cpl_fields_ibuf_stopnow) == 1) then
            write (nu_diag,*) '(ice) received final coupler msg'
         else
            write (nu_diag,*) '(ice) terminating before coupler'
            call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
         endif
      endif

      call cpl_interface_finalize (cpl_fields_icename)

      if (my_task == master_task) then
         write(nu_diag,*) '(ice) exit_coupler finished',my_task
      endif

      end subroutine exit_coupler

!=======================================================================

#endif
! ifdef CCSM (previous five subroutines)
!=======================================================================

      end module ice_coupling

!=======================================================================
