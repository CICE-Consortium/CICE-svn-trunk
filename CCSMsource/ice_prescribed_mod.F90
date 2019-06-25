!=======================================================================
! CVS: $Id: ice_prescribed_mod.F90 782 2006-04-25 22:42:24Z dbailey $
! CVS: $Source$
! CVS: $Name$
!=======================================================================
!BOP ===================================================================
!
! !MODULE: ice_prescribed_mod - Prescribed Ice Model
!
! !DESCRIPTION:
!
! The prescribed ice model reads in ice concentration data from a netCDF
! file.  Ice thickness, temperature, the ice temperature profile are
! prescribed.  Air/ice fluxes are computed to get surface temperature,
! Ice/ocean fluxes are set to zero, and ice dynamics are not calculated.
! Regridding and data cycling capabilities are included.
!
! !REVISION HISTORY:
!
! 2005-May-19 - J. Schramm - first version
! 2005-Apr-19 - B. Kauffman, J. Schramm, M. Vertenstein, NCAR - design
!
! !INTERFACE: ----------------------------------------------------------
 
module ice_prescribed_mod
 
! !USES:

   use shr_stream_mod
   use shr_map_mod
   use shr_ncread_mod
   use shr_sys_mod

   use ice_kinds_mod
   use ice_fileunits,  only : nu_diag
   use ice_exit,       only : abort_ice
   use ice_model_size, only : imt_global, jmt_global, ncat, nilyr, ntilay
   use ice_constants,  only : c0, c1, c2, eps04, c360
   use ice_domain,     only : imt_local, jmt_local, ilo, ihi, jlo, jhi, my_task, master_task
   use ice_grid,       only : TLAT,TLON,hm,tmask
   use ice_calendar,   only : idate, sec
   use ice_itd,        only : ilyr1, hin_max
   use ice_work,       only : work_g1, work_g2

   implicit none
   save

   private ! except

! !PUBLIC TYPES:

! !PUBLIC MEMBER FUNCTIONS:

   public :: ice_prescribed_init      ! initialize input data stream
   public :: ice_prescribed_run       ! get time slices and time interp
   public :: ice_prescribed_readField ! reads field from input data file
   public :: ice_prescribed_phys      ! set prescribed ice state and fluxes

! !PUBLIC DATA MEMBERS:

   character(len=char_len_long), public :: stream_info_file  ! file with stream info
   logical(kind=log_kind), public :: prescribed_ice      ! true if prescribed ice
   integer(kind=int_kind), public :: stream_year_first   ! first year in stream to use
   integer(kind=int_kind), public :: stream_year_last    ! last year in stream to use
   integer(kind=int_kind), public :: model_year_align    ! align stream_year_first
                                                         ! with this model year
   logical(kind=log_kind), public :: prescribed_ice_fill ! true if data fill required

!EOP

   real(kind=dbl_kind),    allocatable :: dataXCoord(:,:)  ! input data longitudes 
   real(kind=dbl_kind),    allocatable :: dataYCoord(:,:)  ! input data latitudes
   integer(kind=int_kind), allocatable :: dataMask(:,:)    ! input data mask
   real(kind=dbl_kind),    allocatable :: dataArea(:,:)    ! input data area
   real(kind=dbl_kind),    allocatable :: dataUB(:,:)      ! upper bound on model domain
   real(kind=dbl_kind),    allocatable :: dataLB(:,:)      ! lower bound on model domain

   type(shr_stream_streamType), save  :: csim_stream        ! csim data stream
   type(shr_map_mapType)        :: csim_map           ! used by shr_map_mapSet 
   type(shr_map_mapType)        :: csim_fill          ! used by shr_map_mapSet

   real(kind=dbl_kind), allocatable :: dataInLB(:,:,:)  ! input data LB
   real(kind=dbl_kind), allocatable :: dataInUB(:,:,:)  ! input data UB
   real(kind=dbl_kind), allocatable :: dataSrcLB(:,:)   ! reformed data for mapping
   real(kind=dbl_kind), allocatable :: dataSrcUB(:,:)   ! reformed data for mapping
   real(kind=dbl_kind), allocatable :: dataDstLB(:,:)   ! output from mapping
   real(kind=dbl_kind), allocatable :: dataDstUB(:,:)   ! output from mapping
   real(kind=dbl_kind), allocatable :: dataOutLB(:,:,:) ! output for model use
   real(kind=dbl_kind), allocatable :: dataOutUB(:,:,:) ! output for model use

   real(kind=dbl_kind), allocatable :: ice_cov_global(:,:)      ! ice cover for model
   real(kind=dbl_kind)              :: ice_cov(ilo:ihi,jlo:jhi) ! scattered ice cover 

   logical(kind=log_kind) :: regrid                      ! true if remapping required

   integer(kind=int_kind) :: dateUB, dateLB, secUB, secLB
   integer(kind=int_kind) :: nlon           ! longitudes in netCDF data file
   integer(kind=int_kind) :: nlat           ! latitudes  in netCDF data file
   integer(kind=int_kind) :: nflds          ! number of fields in list

   character(len=char_len_long) :: fldList      ! list of fields in data stream
   character(len=char_len)      :: fldName      ! name of field in stream

!=======================================================================
contains
!=======================================================================
!BOP ===================================================================
!
! !IROUTINE: ice_prescribed_init --  initialize data stream information
!
! !DESCRIPTION:
! (1) Initialize data stream information.  
! (2) Check input data domain - currently supports a regular lat-lon grid
!     or cpl6 history output.   If input data domain does not match ice
!     model domain, initialize mapping weights and set 'regrid' to true.
!
! !REVISION HISTORY:
!     2005-May-19 - J. Schramm - first version
!
! !INTERFACE: -----------------------------------------------------------

subroutine ice_prescribed_init(stream_info_file, stream_year_first, &
&                              stream_year_last, model_year_align,  &
&                              prescribed_ice_fill)

! !USES:

   use shr_string_mod
   use ice_mpi_internal, only : global_gather

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(len=char_len_long), intent(in) :: stream_info_file   ! file with stream info
   integer(kind=int_kind), intent(in) :: stream_year_first  ! first year in stream to use
   integer(kind=int_kind), intent(in) :: stream_year_last   ! last year in stream to use
   integer(kind=int_kind), intent(in) :: model_year_align   ! align stream_year_first
                                                            ! with this model year
   logical(kind=log_kind), intent(in) :: prescribed_ice_fill ! true if data fill required

!EOP
   !----- Local ------
   real(kind=dbl_kind), allocatable :: work4d(:,:,:,:) ! 4D work array
   integer(kind=int_kind), allocatable :: csim_mask_g(:,:) ! global csim land mask
   real(kind=dbl_kind)     :: aice_max                 ! maximun ice concentration
   character(len=char_len_long) :: first_data_file     ! first data file in stream
   character(len=char_len_long) :: domain_info_fn      ! file with domain info
   integer(kind=int_kind)  :: ndims                    ! # dimensions in domain info
   character(len=char_len) :: timeName                 ! domain time var name
   character(len=char_len) :: latName                  ! domain latitude var name
   character(len=char_len) :: lonName                  ! domain longitude var name
   character(len=char_len) :: maskName                 ! domain mask var name
   character(len=char_len) :: areaName                 ! domain area var name
   logical(kind=log_kind)  :: check                    ! true if field is found

   !----- formats -----
   character(*),parameter :: subName = "('ice_prescribed_init')"
   character(*),parameter :: F00 = "('(ice_prescribed_init) ',4a)"
   character(*),parameter :: F01 = "('(ice_prescribed_init) ',a,2i6)"
   character(*),parameter :: F02 = "('(ice_prescribed_init) ',a,2g20.13)"

   !------------------------------------------------------------------
   ! Create integer CSIM mask with global dimensions
   !------------------------------------------------------------------
      allocate (work_g1(imt_global,jmt_global),work_g2(imt_global,jmt_global),&
      &  csim_mask_g(imt_global,jmt_global))

      call global_gather(work_g1,hm(ilo:ihi,jlo:jhi))
      csim_mask_g = work_g1     ! Convert to integer array

      call global_gather(work_g1,TLAT(ilo:ihi,jlo:jhi))
      call global_gather(work_g2,TLON(ilo:ihi,jlo:jhi))

   if (my_task == master_task) then


      !---------------------------------------------------------------------
      ! Parse info file, initialize csim_stream datatype and load with info
      !---------------------------------------------------------------------
      call shr_stream_init(csim_stream, stream_info_file, stream_year_first, &
   &                       stream_year_last, model_year_align)

      !---------------------------------------------------------------------
      ! Check that ice cover data exists
      ! Assumes that stream has one field, ice fraction
      !---------------------------------------------------------------------
      call shr_stream_getFileFieldList(csim_stream,fldList)
      nflds = shr_string_listGetNum(fldList)                ! How many fields?

      if (nflds == 1) then
         call shr_string_listGetName(fldList,nflds,fldName) ! Get name of 1st field
      else
         write(nu_diag,F00) "ERROR: Only one field can exist in stream"
         write(nu_diag,F01)  trim(fldList)//' has ',nflds
         call abort_ice(subName)
      end if

      call shr_stream_getFirstFileName(csim_stream,first_data_file)
      check = shr_ncread_varExists(first_data_file,fldName) 

      if (.not.check) then
         write(nu_diag,F00) "ERROR: ice concentration field does not exist"
         call abort_ice(subName)
      end if

      !---------------------------------------------------------------------
      ! Get size of the input data domain and allocate arrays
      !---------------------------------------------------------------------

      call shr_stream_getDomainInfo(csim_stream, domain_info_fn, timeName, &
      &                             lonName, latName, maskName, areaName)

      call shr_ncread_varDimSizes(domain_info_fn,areaName,nlon,nlat)
      write (nu_diag,F01) 'dimsizes for areaName', nlon,nlat ; call shr_sys_flush(6)

      allocate(dataXCoord(nlon,nlat)) 
      allocate(dataYCoord(nlon,nlat))
      allocate(dataMask(nlon,nlat))
      allocate(dataArea(nlon,nlat))

      !------------------------------------------------------------------
      ! Read in lat-lon grid info from input data file, ALL output is 2D
      !------------------------------------------------------------------
      call shr_ncread_domain(domain_info_fn, lonName, dataXCoord, latName, &
      &                dataYCoord, maskName, dataMask, areaName, dataArea)

      write (nu_diag,F02) 'min/max dataXCoord = ',minval(dataXCoord), maxval(dataXCoord)
      write (nu_diag,F02) 'min/max dataYCoord = ',minval(dataYCoord), maxval(dataYCoord)
      write (nu_diag,F01) 'min/max dataMask = ',minval(dataMask), maxval(dataMask)
      write (nu_diag,F02) 'min/max dataArea = ',minval(dataArea), maxval(dataArea)
      write (nu_diag,F02) 'min/max TLON = ',minval(work_g2), maxval(work_g2)
      write (nu_diag,F02) 'min/max TLAT = ',minval(work_g1), maxval(work_g1)
      write (nu_diag,F01) 'min/max csim_mask_g = ',minval(csim_mask_g), maxval(csim_mask_g)
      !------------------------------------------------------------------
      ! Check that data domain matches csim domain
      !------------------------------------------------------------------
      call ice_prescribed_checkDomain(work_g2, work_g1, csim_mask_g, dataXCoord, &
      &                               dataYCoord, dataMask, regrid)

      if (regrid) then

         if (prescribed_ice_fill) then
            call shr_map_mapSet(csim_fill, work_g2, work_g1, csim_mask_g, &
            &                              work_g2, work_g1, csim_mask_g, &
            &                   name='csim_map',type='fill',algo='nnoni', &
            &                   mask='nomask')
         end if

         write (nu_diag,F00) 'Computing mapping weights'

         call shr_map_mapSet(csim_map, dataXCoord, dataYCoord, dataMask, &
         &                             work_g2,   work_g1,  csim_mask_g, &
         &                   name='csim_map',type='remap',algo='bilinear', &
         &                   mask='dstmask',vect='scalar')
      end if

      allocate(dataInLB(nlon,nlat,nflds))      ! netCDF input data size
      allocate(dataInUB(nlon,nlat,nflds))
      allocate(dataSrcLB(nflds,nlon*nlat))     ! reformed data for mapping
      allocate(dataSrcUB(nflds,nlon*nlat))

      allocate(dataDstLB(nflds,imt_global*jmt_global))  ! output from mapping
      allocate(dataDstUB(nflds,imt_global*jmt_global))
      allocate(dataOutLB(imt_global,jmt_global,nflds))  ! output for model use
      allocate(dataOutUB(imt_global,jmt_global,nflds))
      allocate(ice_cov_global(imt_global,jmt_global))   ! Assumes 1 field

      !--------------------------------------------------------------------
      ! Check to see that ice concentration is in fraction, not percent
      !--------------------------------------------------------------------
      allocate(work4d(nlon,nlat,1,1))
      call shr_ncread_field4dG(first_data_file, fldName, rfld=work4d, dim3i=1)
      aice_max = maxval(work4d)
      deallocate(work4d)

      if (aice_max > c2) then
         write(nu_diag,F02) "ERROR: Ice conc data must be in fraction, aice_max= ",aice_max
         call abort_ice(subName)
      end if

   end if           ! master_task

   deallocate(dataXCoord,dataYCoord)
   deallocate(dataMask,dataArea)
   deallocate(work_g1,work_g2)
   deallocate(csim_mask_g)

   !-----------------------------------------------------------------
   ! For one ice category, set hin_max(1) to something big
   !-----------------------------------------------------------------
   if (ncat == 1) then
      hin_max(1) = 999._dbl_kind
   end if
    
end subroutine ice_prescribed_init
  
!=======================================================================
!BOP ===================================================================
!
! !IROUTINE: ice_prescribed_run -- Obtain two time slices of data for 
!                                  current time
!
! !DESCRIPTION:
! 
!  Finds two time slices bounding current model time, remaps if necessary
!
! !REVISION HISTORY:
!     2005-May-19 - J. Schramm - first version
!
! !INTERFACE: -----------------------------------------------------------

subroutine ice_prescribed_run(mDateIn, secIn, prescribed_ice_fill)

! !USES:

   use shr_tInterp_mod
   use ice_mpi_internal, only : global_scatter

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(kind=int_kind), intent(in) :: mDateIn  ! Current model date (yyyymmdd)
   integer(kind=int_kind), intent(in) :: secIn    ! Elapsed seconds on model date
   logical(kind=log_kind), intent(in) ::  prescribed_ice_fill  ! true if data fill required

!EOP

   real(kind=dbl_kind)    :: fLB         ! weight for lower bound
   real(kind=dbl_kind)    :: fUB         ! weight for upper bound
   integer(kind=int_kind) :: mDateLB     ! model date    of lower bound
   integer(kind=int_kind) :: mDateUB     ! model date    of upper bound
   integer(kind=int_kind) :: dDateLB     ! data  date    of lower bound
   integer(kind=int_kind) :: dDateUB     ! data  date    of upper bound
   integer(kind=int_kind) ::   secUB     ! elap sec      of upper bound
   integer(kind=int_kind) ::   secLB     ! elap sec      of lower bound
   integer(kind=int_kind) ::    n_lb     ! t-coord index of lower bound
   integer(kind=int_kind) ::    n_ub     ! t-coord index of upper bound
   character(len=char_len_long) ::  fileLB     ! file containing  lower bound
   character(len=char_len_long) ::  fileUB     ! file containing  upper bound

   integer(kind=int_kind) :: mDateLB_old = -999
   integer(kind=int_kind) :: secLB_old = -999
   integer(kind=int_kind) :: i,j,n,icnt  ! loop indices and counter

   !------------------------------------------------------------------------
   ! get two time slices of monthly ice coverage data
   ! check that upper and lower time bounds have not changed
   !------------------------------------------------------------------------

   if (my_task ==  master_task) then

      call shr_stream_findBounds(csim_stream, mDateIn,secIn,             &
      &                          mDateLB, dDateLB, secLB, n_lb, fileLB, &
      &                          mDateUB, dDateUB, secUB, n_ub, fileUB)

      if (mDateLB_old /= mDateLB .or. secLB_old /= secLB) then

         call ice_prescribed_readField(fileLB, fldName, n_lb, dataInLB)
         call ice_prescribed_readField(fileUB, fldName, n_ub, dataInUB)

         if (regrid) then
            !-------------------------------------------------
            ! copy input data to arrays ordered for mapping
            !-------------------------------------------------
            do n=1,nflds
              icnt = 0
              do j=1,nlat
              do i=1,nlon
                 icnt = icnt + 1
                 dataSrcLB(n,icnt) = dataInLB(i,j,n)
                 dataSrcUB(n,icnt) = dataInUB(i,j,n)
              enddo
              enddo
            enddo

            !-------------------------------------------------
            ! map the ordered arrays
            !-------------------------------------------------
            if (prescribed_ice_fill) then
               call shr_map_mapData(dataSrcLB, dataDstLB, csim_fill)
               call shr_map_mapData(dataSrcUB, dataDstUB, csim_fill)
            end if
            call shr_map_mapData(dataSrcLB, dataDstLB, csim_map)
            call shr_map_mapData(dataSrcUB, dataDstUB, csim_map)

            !-------------------------------------------------
            ! copy mapped fields back to general 3d arrays
            !-------------------------------------------------
            do n=1,nflds
               icnt = 0
               do j=1,jmt_global
               do i=1,imt_global
                  icnt = icnt + 1
                  dataOutLB(i,j,n) = dataDstLB(n,icnt)
                  dataOutUB(i,j,n) = dataDstUB(n,icnt)
               enddo
               enddo
            enddo
         else       ! no regrid
            dataOutLB = dataInLB
            dataOutUB = dataInUB

         end if    ! regrid

         mDateLB_old = mDateLB
         secLB_old = secLB

      end if       ! mDate
    
      call shr_tInterp_getFactors(mDateLB, secLB, mDateUB, secUB, &
      &                           mDateIn, secIN, fLB, fUB)

      ice_cov_global(:,:) = fLB*dataOutLB(:,:,1) + fUB*dataOutUB(:,:,1)

   end if    ! master_task

   deallocate(dataInLB,dataInUB)
   deallocate(dataOutLB,dataOutUB)
   deallocate(dataSrcLB,dataSrcUB)
   deallocate(dataDstLB,dataDstUB)

  !-----------------------------------------------------------------
  ! Scatter ice concentration to all processors
  !-----------------------------------------------------------------
   call global_scatter(ice_cov_global,ice_cov(ilo:ihi,jlo:jhi))
   deallocate(ice_cov_global)

  !-----------------------------------------------------------------
  ! Set prescribed ice state and fluxes
  !-----------------------------------------------------------------

   call ice_prescribed_phys

end subroutine ice_prescribed_run

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ice_prescribed_readField -- Read a field from input data
!
! !DESCRIPTION:
!     Read a field from input data, assumes prescribed ice uses only one field
!
! !REVISION HISTORY:
!     2005-May-23 - J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine ice_prescribed_readField(fileName, fldName, nTime, fldRead)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*), intent(in)  :: fileName     ! netCDF input file
   character(*), intent(in)  :: fldName      ! name of field in stream
   integer(kind=int_kind), intent(in)  :: nTime          ! index of time slice to read
   real   (kind=dbl_kind), intent(out) :: fldRead(:,:,:) ! field read in

!EOP

   real(kind=dbl_kind), allocatable :: A4d(:,:,:,:)    ! 4D field array
   integer(kind=int_kind)           :: ni, nj          ! size of field to read

   ni = size(fldRead,1)
   nj = size(fldRead,2)

   allocate(A4d(ni,nj,1,1))

   call shr_ncread_field4dG(fileName, fldName, rfld=A4d, dim3i=nTime)
   fldRead(:,:,1) = A4d(:,:,1,1)

   deallocate(A4d)

end subroutine ice_prescribed_readField

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ice_prescribed_checkDomain -- Compare input data domain and ice domain 
!
! !DESCRIPTION:
!     Check that input data domain matches that of ice model domain
!
! !REVISION HISTORY:
!     2005-May-16 - J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine ice_prescribed_checkDomain(csimXCoord, csimYCoord, csimMask, &
   &                                  dataXCoord, dataYCoord, dataMask, regrid)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   real(kind=dbl_kind),    intent(in) :: csimXCoord(:,:) ! csim longitudes (degrees)
   real(kind=dbl_kind),    intent(in) :: csimYCoord(:,:) ! csim latitudes  (degrees)
   integer(kind=int_kind), intent(in) :: csimMask(:,:)   ! csim land mask
   real(kind=dbl_kind),    intent(in) :: dataXCoord(:,:) ! data longitudes (degrees)
   real(kind=dbl_kind),    intent(in) :: dataYCoord(:,:) ! data latitudes  (degrees)
   integer(kind=int_kind), intent(in) :: dataMask(:,:)   ! data land mask
   logical(kind=log_kind), intent(out):: regrid          ! true if remapping required

!EOP

   integer(kind=int_kind)           :: ni_csim, nj_csim     ! size of csim domain
   integer(kind=int_kind)           :: ni_data, nj_data     ! size of data domain
   real(kind=dbl_kind)              :: max_XDiff, max_YDiff ! max x,y grid diffs
   real(kind=dbl_kind), allocatable :: csimXTemp(:,:)       ! csim grid 0 to 360

   !----- formats -----
   character(*),parameter :: F00 = "('(ice_prescribed_checkDomain) ',a)"
   character(*),parameter :: F01 = "('(ice_prescribed_checkDomain) ',a,2i6)"
   character(*),parameter :: F02 = "('(ice_prescribed_checkDomain) ',a,g20.13)"

   regrid = .false.
   !------------------------------------------------------------------------
   ! Check domain size
   !------------------------------------------------------------------------
   ni_csim = size(csimXCoord,1)
   nj_csim = size(csimXCoord,2)
   ni_data = size(dataXCoord,1)
   nj_data = size(dataXCoord,2)

   !------------------------------------------------------------------------
   ! Convert T grid longitude from -180 -> 180 to 0 to 360
   !------------------------------------------------------------------------
   allocate(csimXTemp(ni_csim,nj_csim)) 
   csimXTemp = csimXCoord
   where (csimXTemp >= c360) csimXTemp = csimXTemp - c360
   where (csimXTemp < c0 )   csimXTemp = csimXTemp + c360

   if (ni_data == ni_csim .and. nj_data == nj_csim) then
      max_XDiff = maxval(abs(csimXTemp  - dataXCoord))
      max_YDiff = maxval(abs(csimYCoord - dataYCoord))
      if (max_XDiff > eps04) regrid = .true.
      if (max_YDiff > eps04) regrid = .true.

      write(nu_diag,F02) 'Maximum X difference = ', max_XDiff
      write(nu_diag,F02) 'Maximum Y difference = ', max_YDiff
      write(nu_diag,F00) 'CSIM and input data domain sizes match'
      write(nu_diag,F01) 'CSIM and data grids are', ni_csim,nj_csim
   else
      regrid = .true.
   end if
   deallocate(csimXTemp) 

   if (regrid) then
      write(nu_diag,F00) 'Ice concentration will be interpolated to CSIM grid'
   else
      write(nu_diag,F00) 'Ice concentration grid and CSIM grid the same'
   end if

end subroutine ice_prescribed_checkDomain

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ice_prescribed_phys -- set prescribed ice state and fluxes
!
! !DESCRIPTION:
!
! Set prescribed ice state using input ice concentration;
! set surface ice temperature to atmospheric value; use
! linear temperature gradient in ice to ocean temperature.
!
! !REVISION HISTORY:
!     2005-May-23 - J. Schramm - Updated with data models
!     2004-July   - J. Schramm - Modified to allow variable snow cover
!     2001-May    - B. P. Briegleb - Original version
!
! !INTERFACE: ------------------------------------------------------------------
 
subroutine ice_prescribed_phys

! !USES:
 
   use ice_flux
   use ice_grid, only : bound
   use ice_state
   use ice_itd, only  : aggregate
   use ice_dyn_evp

   implicit none
 
! !INPUT/OUTPUT PARAMETERS:
 
!EOP

   !----- Local ------
   integer(kind=int_kind) :: layer    ! level index
   integer(kind=int_kind) :: nc       ! ice category index
   integer(kind=int_kind) :: i,j,k    ! longitude, latitude and level indices

   real(kind=dbl_kind) :: slope     ! diff in underlying ocean tmp and ice surface tmp
   real(kind=dbl_kind) :: Ti        ! ice level temperature
   real(kind=dbl_kind) :: hi        ! ice prescribed (hemispheric) ice thickness 
   real(kind=dbl_kind) :: hs        ! snow cover
   real(kind=dbl_kind) :: dz        ! distance freeboard below SL (m)
   real(kind=dbl_kind) :: dhs       ! snow to remove 
   real(kind=dbl_kind) :: zn        ! normalized ice thickness
   real(kind=dbl_kind) :: salin(nilyr)  ! salinity (ppt) 

   real(kind=dbl_kind), parameter :: nsal    = 0.407_dbl_kind
   real(kind=dbl_kind), parameter :: msal    = 0.573_dbl_kind
   real(kind=dbl_kind), parameter :: saltmax = 3.2_dbl_kind   ! max salinity at ice base (ppm)

   !-----------------------------------------------------------------
   ! Initialize ice state
   !-----------------------------------------------------------------
   do nc=1,ncat
      do j=1,jmt_local
      do i=1,imt_local
         aicen(i,j,nc) = c0
         vicen(i,j,nc) = c0
      enddo
      enddo
   enddo

   do nc=1,ncat
      Tsfcn(ilo:ihi,jlo:jhi,nc) = Tf(ilo:ihi,jlo:jhi)
      call bound(Tsfcn(:,:,nc))
   enddo

   do layer=1,ntilay
      do j=1,jmt_local
      do i=1,imt_local
         eicen(i,j,layer) = c0
      enddo
      enddo
   enddo

   !-----------------------------------------------------------------
   ! Set ice cover over land to zero, not sure if this should be
   ! be done earier, before time/spatial interp??????
   !-----------------------------------------------------------------
   do j = jlo,jhi
   do i = ilo,ihi
      if (tmask(i,j)) then
         if (ice_cov(i,j) .lt. eps04) ice_cov(i,j) = c0
         if (ice_cov(i,j) .gt. c1)    ice_cov(i,j) = c1
      else
         ice_cov(i,j) = c0
      end if
   enddo
   enddo

   do j = jlo,jhi
   do i = ilo,ihi

      if (tmask(i,j)) then   ! Over ocean points

         !--------------------------------------------------------------
         ! Place ice where ice concentration > .0001
         !--------------------------------------------------------------
         if (ice_cov(i,j) >= eps04) then

            hi = 0.0_dbl_kind
            !----------------------------------------------------------
            ! Set ice thickness in each hemisphere
            !----------------------------------------------------------
            if(TLAT(i,j)*rad_to_deg > 40.0_dbl_kind) then
              hi  = 2.0_dbl_kind
            else if(TLAT(i,j)*rad_to_deg < -40.0_dbl_kind) then
              hi  = 1.0_dbl_kind
            end if

            !----------------------------------------------------------
            ! All ice in appropriate thickness category
            !----------------------------------------------------------
            do nc = 1,ncat
               if(hin_max(nc-1) < hi .and. hi < hin_max(nc)) then
                  aicen(i,j,nc) = ice_cov(i,j)
                  vicen(i,j,nc) = hi*aicen(i,j,nc) 

                  !---------------------------------------------------------
                  ! keep snow/ice boundary above sea level by reducing snow
                  !---------------------------------------------------------
                  hs = vsnon(i,j,nc)/aicen(i,j,nc)
                  dz = hs - hi*(rhow-rhoi)/rhos

                  if (dz > puny .and. hs > puny) then ! snow below freeboard
                     dhs = min(dz*rhoi/rhow, hs) ! snow to remove
                     hs = hs - dhs
                     vsnon(i,j,nc) = hs * aicen(i,j,nc)
                  end if

                  !---------------------------------------------------------
                  ! make linear temp profile and compute enthalpy
                  !---------------------------------------------------------
                  Tsfcn(i,j,nc) = min(Tair(i,j)-Tffresh,-p2)   ! deg C       
                  slope = Tf(i,j) - Tsfcn(i,j,nc)
                  do k = 1, nilyr
                     zn = (real(k,kind=dbl_kind)-p5) / real(nilyr,kind=dbl_kind)
                     Ti = Tsfcn(i,j,nc) + slope*zn
                     salin(k) = (saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
                     eicen(i,j,ilyr1(nc)+k-1) =                             &
                     &    (-rLfi - rcpi*(-depressT*salin(k)-Ti)             &
                     &     -rLfidepressT*salin(k)/Ti) *vicen(i,j,nc)/nilyr
                  enddo
               end if    ! hin_max
            enddo        ! ncat
         end if          ! ice_cov >= eps04
      end if             ! tmask
   enddo                 ! i
   enddo                 ! j

   !--------------------------------------------------------------------
   ! compute aggregate ice state and open water area
   !--------------------------------------------------------------------
   call aggregate 

   do j = 1, jmt_local
     do i = 1, imt_local
       aice_init(i,j) = aice(i,j)
     enddo
   enddo

   !--------------------------------------------------------------------
   ! set non-computed fluxes, ice velocities, ice-ocn stresses to zero
   !--------------------------------------------------------------------

   do j=jlo,jhi
   do i=ilo,ihi
      frzmlt(i,j)   = c0     ! freezing/melting potential (W/m^2)
      uvel(i,j)     = c0     ! ice velocity u
      vvel(i,j)     = c0     ! ice velocity v
      strocnxT(i,j) = c0     ! ice-ocean stress, x-direction (T-cell)
      strocnyT(i,j) = c0     ! ice-ocean stress, y-direction (T-cell)
   enddo
   enddo

   !-----------------------------------------------------------------
   ! other atm and ocn fluxes
   !-----------------------------------------------------------------
   call init_flux_atm
   call init_flux_ocn

end subroutine ice_prescribed_phys

!=====================================================================================

end module ice_prescribed_mod

!=====================================================================================
