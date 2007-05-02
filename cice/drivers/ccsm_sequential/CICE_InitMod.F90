!=======================================================================
!
!BOP
!
! !MODULE: CICE_InitMod - performs CICE initialization
!
! !DESCRIPTION:
!
!  This module contains the CICE initialization routine that sets model
!  parameters and initializes the grid and CICE state variables.
!
! !REVISION HISTORY:
!  SVN:$Id: CICE_InitMod.F90 52 2007-01-30 18:04:24Z eclare $
!
!  authors Elizabeth C. Hunke, LANL
!          William H. Lipscomb, LANL
!          Philip W. Jones, LANL
!
! 2006: Converted to free form source (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module CICE_InitMod
!
! !USES:
!
      use ice_calendar
      use ice_communicate
      use ice_diagnostics
      use ice_domain
      use ice_dyn_evp
      use ice_exit
      use ice_fileunits
      use ice_flux
      use ice_grid
      use ice_history
      use ice_restart
      use ice_init
      use ice_itd
      use ice_kinds_mod
      use ice_mechred
      use ice_ocean
      use ice_orbital
      use ice_shortwave
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap
      use ice_work
      use ice_prescribed_mod

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: cice_init

!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: cice_init - initialize CICE model
!
! !DESCRIPTION:
!
!  Initialize CICE model.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine cice_init( mpicom_ice )
!
! !USES:
!
      use ice_domain_size
      use ice_blocks
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         mpicom_ice ! communicator for sequential ccsm
!EOP
!
!     local temporary variables

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         alvdrni      , & 
         alidrni      , &
         alvdfni      , &
         alidfni      , &
         alvdrns      , & 
         alidrns      , &
         alvdfns      , &
         alidfns

      integer (kind=int_kind) :: i, j, ij, n, iblk, ilo, ihi, jlo, jhi

      call init_communicate( mpicom_ice ) ! initial setup for message passing
      call input_data           ! namelist variables
      call init_work            ! work arrays

      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution
      call init_ice_timers      ! initialize all timers
      call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! grid variables

      if (advection == 'remap') &
         call init_remap        ! grid variables for remapping transport
      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file
      call init_evp (dt)        ! define evp dynamics parameters, variables
      call init_coupler_flux    ! initialize fluxes exchanged with coupler
      call init_thermo_vertical ! initialize vertical thermodynamics
      if (shortwave == 'dEdd') then
         call init_orbit        ! initialize orbital parameters
         call init_dEdd         ! initialize delta-Eddington scheme
      endif
      call init_itd             ! initialize ice thickness distribution
      call calendar(time)       ! determine the initial date
      call init_state           ! initialize the ice state
      call ice_prescribed_init
      if (restart) call restartfile      ! start from restart file

      ! Need to compute albedos before init_cpl in CCSM

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      alvdr   (:,:,:) = c0
      alidr   (:,:,:) = c0
      alvdf   (:,:,:) = c0
      alidf   (:,:,:) = c0

      do iblk=1,nblocks
      do n=1,ncat

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n,iblk) > puny) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

         call compute_albedos (nx_block,   ny_block, &
                               icells,               &
                               indxi,      indxj,    &
                               aicen(:,:,n,iblk), vicen(:,:,n,iblk),    &
                               vsnon(:,:,n,iblk), trcrn(:,:,1,n,iblk),  &
                               alvdrni,           alidrni,  &
                               alvdfni,           alidfni,  &
                               alvdrns,           alidrns,  &
                               alvdfns,           alidfns,  &
                               alvdrn(:,:,n,iblk),alidrn(:,:,n,iblk),   &
                               alvdfn(:,:,n,iblk),alidfn(:,:,n,iblk),   &
                               apondn(:,:,n,iblk),hpondn(:,:,n,iblk))

         ! Aggregate albedos for coupler

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            alvdf(i,j,iblk) = alvdf(i,j,iblk) &
               + alvdfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidf(i,j,iblk) = alidf(i,j,iblk) &
               + alidfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alvdr(i,j,iblk) = alvdr(i,j,iblk) &
               + alvdrn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidr(i,j,iblk) = alidr(i,j,iblk) &
               + alidrn(i,j,n,iblk)*aicen(i,j,n,iblk)
         enddo

      enddo
      enddo

      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables

      write_ic = .true.        ! write initial conditions
      if(.not.prescribed_ice) call ice_write_hist(dt)
      write_ic = .false.

      end subroutine cice_init

!=======================================================================

      end module CICE_InitMod

!=======================================================================
