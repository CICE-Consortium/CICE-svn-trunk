c $Id: ice_init.F,v 1.10 2004/02/25 17:35:13 eclare Exp $
!=======================================================================
!BOP
!
! !MODULE:   ice_init - parameter and variable initializations
!
! !DESCRIPTION:
! 
! parameter and variable initializations
!
! !REVISION HISTORY:
! 
! authors Elizabeth C. Hunke, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      module ice_init
!
! !USES:
!
      use ice_domain
!
!EOP
!
      implicit none

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: input_data - namelist variables
!
! !INTERFACE:
!
      subroutine input_data
!
! !DESCRIPTION:
!
! Namelist variables, set to default values; may be altered
! at run time
!
! !REVISION HISTORY:
! 
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_albedo
      use ice_mechred_cice
      use ice_diagnostics
      use ice_history
      use ice_calendar
      use ice_dyn_evp
      use ice_itd
      use ice_transport_remap, only: advection
      use ice_ocean, only: oceanmixed_ice
      use ice_flux_in, only:  
     &            ycycle, fyear_init, atm_data_dir, ocn_data_dir

      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag

      character (len=6) :: chartmp

      !-----------------------------------------------------------------
      ! namelist variables used in CSIM and not defined elsewhere in CICE
      !-----------------------------------------------------------------
      integer (kind=int_kind) ::
     &   kcolumn              ! 1 for column model

      logical (kind=log_kind) ::
     &   snow_into_ocn  ! if true, throw snow into ocean
     &,  prescribed_ice ! if true, use prescribed ice
     &,  prescribed_ice_climatology ! if true, use climate-prescribed ice
     &,  oceanmixed_ice_sst_init    ! if true, use ocn frc sst for init
     &,  prntdiag_oceanmixed        ! if true, print out some diagnostics

      character(len=char_len_long) :: 
     &   prescribed_ice_file   ! netCDF file name for prescribed ice
     &,  oceanmixed_ice_file   ! netCDF file name for ocean forcing data
 
      !-----------------------------------------------------------------
      ! Namelist variables.
      ! NOTE: This list includes all the namelist variables in CSIM.
      !       Not all of these are used by standalone CICE.
      !-----------------------------------------------------------------

      namelist /ice_nml/          runid,       runtype, 
     & istep0,      dt,           ndte,        ndyn_dt,        npt, 
     & diagfreq,    histfreq,     dumpfreq,    dumpfreq_n,
     & restart,     print_points, kcolumn,     kitd, 
     & kdyn,        kstrength,    evp_damping, snow_into_ocn,
     & grid_file,   kmt_file,     
     & incond_dir,  incond_file,  restart_dir, 
     & dump_file,   pointer_file, history_dir, history_file, 
     & grid_type,   advection,    hist_avg, 
     & prescribed_ice, prescribed_ice_file, prescribed_ice_climatology,
     & oceanmixed_ice, oceanmixed_ice_file, oceanmixed_ice_sst_init,
     & prntdiag_oceanmixed
     &,albicev, albicei, albsnowv, albsnowi 
     &,year_init, ksmooth, print_global, diag_file     ! CICE only 
     &,ycycle, fyear_init, atm_data_dir, ocn_data_dir, dump_dir

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------
      runid   = 'unknown_runid'     ! default run ID
      runtype = 'unknown_runtype'   ! default runtype
      istep0 = 0         ! number of steps taken in previous integrations,
                         ! real (dumped) or imagined (use to set calendar)
      dt = 3600.0_dbl_kind  ! time step, s 
      ndte = 120         ! subcycles per dynamics timestep:  ndte=dyn_dt/dte
      ndyn_dt = 1        ! timestep for dynamics = dt/ndyn_dt
      npt = 99999        ! total number of time steps (dt) 
      diagfreq = 24      ! how often diag output is written
      histfreq='m'       ! output frequency
      dumpfreq='y'       ! restart frequency option
      dumpfreq_n= 1      ! restart frequency
      hist_avg = .true.  ! write snapshots rather than time-averages
      restart = .false.  ! do not read restart files for initialization
      print_points = .false. ! if true, print point data
      kcolumn = 0        ! 1 = column model
      kitd = 1           ! type of itd conversions (0 = delta f, 1 = linear)
      kdyn = 1           ! type of dynamics (1 = evp)
      kstrength = 1      ! 1 = use Rothrock 1975 pressure formulation
      evp_damping = .false.  ! if true, use damping procedure in evp dynamics
      advection  = 'remap'   ! incremental remapping transport scheme
      snow_into_ocn = .false.! if true, throw snow into ocn during ridging
      grid_type    = 'displaced_pole'  ! define rectangular grid internally

      grid_file    = 'data.domain.grid'
      kmt_file     = 'data.domain.kmt'
      incond_dir   = ' '             ! Write to executable dir for default
      incond_file  = 'incond'
      restart_dir  = ' '             ! Write to executable dir for default
      dump_file    = 'iced'          ! This is the restart file name prefix
      pointer_file = 'ice.restart_file'
      history_dir  = ' '             ! Write to executable dir for default
      history_file = 'iceh'
      prescribed_ice          = .false.             ! if true, prescribe ice
      prescribed_ice_file     = 'prescribed_ice.nc' ! local prescribed ice file name
      prescribed_ice_climatology = .true.           ! if true, use climatology
      oceanmixed_ice          = .false.             ! if true, use internal ocean mixed layer
      oceanmixed_ice_file     = 'oceanmixed_ice.nc' ! local ocean forcing file name
      oceanmixed_ice_sst_init = .false.             ! if true, initialize sst from ocn frc data
      prntdiag_oceanmixed     = .false.             ! if true, print out some diagnostics

      albicev   = 0.78_dbl_kind   ! visible ice albedo for h > ahmax
      albicei   = 0.36_dbl_kind   ! near-ir ice albedo for h > ahmax
      albsnowv  = 0.98_dbl_kind   ! cold snow albedo, visible
      albsnowi  = 0.70_dbl_kind   ! cold snow albedo, near IR

      ! The following are in CICE but not CSIM.
      year_init = 1997   ! initial year
      print_global = .false. ! if true, print global diagnostic data
      ksmooth = 0        ! 1 = smooth the ice strength; 0 = no smoothing
      diag_file    = 'ice_diag.d'
      ycycle = 1                   ! number of years in forcing cycle
      fyear_init = 1900            ! first year of forcing cycle
      dump_dir = './'              ! directory in which restart files are written

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------
      if (my_task == master_task) then
        open (nu_nml, file='ice_in', status='old')
   10   continue  !*** keep reading until right namelist is found
        read(nu_nml, nml=ice_nml, iostat=nml_error)
        if (nml_error > 0) goto 10    ! An error occurred
        if (nml_error < 0) goto 20    ! End of file condition
        close(nu_nml)
   20   continue
      endif
      call ice_bcast_iscalar(nml_error)

      if (nml_error /= 0) then
         call abort_ice ('ice: Namelist read error in ice_init.F')
      endif

      call ice_bcast_char(runid)
      call ice_bcast_char(runtype)
      call ice_bcast_iscalar(istep0)
      call ice_bcast_rscalar(dt)
      call ice_bcast_iscalar(ndte)
      call ice_bcast_iscalar(ndyn_dt)
      call ice_bcast_iscalar(npt)
      call ice_bcast_iscalar(diagfreq)
      call ice_bcast_char(histfreq)
      call ice_bcast_char(dumpfreq)
      call ice_bcast_iscalar(dumpfreq_n)
      call ice_bcast_logical(hist_avg)
      call ice_bcast_logical(print_points)
      call ice_bcast_iscalar(kcolumn)
      call ice_bcast_iscalar(kitd)
      call ice_bcast_iscalar(kdyn)
      call ice_bcast_logical(restart)
      call ice_bcast_iscalar(kstrength)
      call ice_bcast_logical(evp_damping)
      call ice_bcast_logical(snow_into_ocn)
      call ice_bcast_char(advection)
      call ice_bcast_char(grid_type)
      call ice_bcast_char(grid_file)
      call ice_bcast_char(kmt_file)
      call ice_bcast_char(restart_dir)
      call ice_bcast_char(dump_file)
      call ice_bcast_char(pointer_file)
      call ice_bcast_char(history_dir)
      call ice_bcast_char(history_file)
      call ice_bcast_logical(prescribed_ice)
      call ice_bcast_char(prescribed_ice_file)
      call ice_bcast_logical(prescribed_ice_climatology)
      call ice_bcast_logical(oceanmixed_ice)
      call ice_bcast_char(oceanmixed_ice_file)
      call ice_bcast_logical(oceanmixed_ice_sst_init)
      call ice_bcast_logical(prntdiag_oceanmixed)
      call ice_bcast_rscalar(albicev)
      call ice_bcast_rscalar(albicei)
      call ice_bcast_rscalar(albsnowv)
      call ice_bcast_rscalar(albsnowi)

      call ice_bcast_iscalar(year_init)
      call ice_bcast_logical(print_global)
      call ice_bcast_iscalar(ksmooth)
      call ice_bcast_char   (diag_file)

      if (histfreq.eq.'1') hist_avg = .false. ! potential conflict
      if (kcolumn.eq.1) grid_type    = 'column' 

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------
      if (nu_diag /= 6) then
         if (my_task == master_task) 
     &        write(6,*) 'Diagnostic output will be in file ', diag_file
         open (nu_diag, file=diag_file, status='unknown')
      endif

      if (my_task == master_task) then
        write (nu_diag,*) '--------------------------------'
        write (nu_diag,*) '  CICE model diagnostic output  '
        write (nu_diag,*) '--------------------------------'

        write (nu_diag,900) dt, ndte, npt, year_init, istep0
        write (nu_diag,901) diagfreq, histfreq, dumpfreq
        write (nu_diag,902) kitd, kdyn, kstrength, kcolumn
        write (nu_diag,*) ' '
        if (evp_damping) then
          write (nu_diag,*) 'Damping procedure used in evp dynamics'
        else
          write (nu_diag,*) 'Damping procedure NOT used in evp dynamics'
        endif
        write (nu_diag,*) ' '
#ifdef coupled
       if( oceanmixed_ice ) then
         write (nu_diag,*) 'WARNING:  coupled and oceanmixed flags are '
         write (nu_diag,*) 'BOTH ON.  Ocean data received from coupler '
         write (nu_diag,*) 'will be altered by mixed layer routine!    '
         write (nu_diag,*) 'WARNING WARNING WARNING '
         write (nu_diag,*) ' '
       endif
#endif

        chartmp = advection(1:6)
        if (chartmp /= 'upwind' .and. chartmp /= 'mpdata')
     &       advection = 'remap'
        write (nu_diag,*) 'Advection algorithm is ', advection
        write (nu_diag,*) ' '
        write (nu_diag,*) 'Diagnostic data for grid points specified'
        if (print_points) then
          write (nu_diag,*) 'in ice_diagnostics.F will be printed'
        else
          write (nu_diag,*) 'in ice_diagnostics.F will NOT be printed'
        endif
        write (nu_diag,*) ' '
        if (grid_type .eq. 'displaced_pole') then
          write (nu_diag,*) 'Displaced pole grid'
          write (nu_diag,*) 'The grid file is ',grid_file
          write (nu_diag,*) 'The land mask file is ', kmt_file
        elseif (grid_type .eq. 'column') then
          write (nu_diag,*) 'Column model grid'
        elseif (grid_type .eq. 'rectangular') then
          write (nu_diag,*) 'Rectangular grid'
        else
          call abort_ice('ice_init: unknown grid_type')
        endif

        write (nu_diag,*) 'The dump file prefix is ', dump_file
        if (restart) then
          write (nu_diag,*) 'The restart file pointer is ', pointer_file
        else
          write (nu_diag,*) 'Starting from default initialization '
          write (nu_diag,*) 'rather than restart file'
        endif
        write (nu_diag,*) 'The history file prefix is ', history_file
        if (hist_avg) then
          write (nu_diag,*) 'History data will be averaged over 1 ',
     &                      histfreq
        else
          write (nu_diag,*) 'History data will be snapshots'
        endif
        write (nu_diag,*) 'Albedo, albicev: ',albicev
        write (nu_diag,*) 'Albedo, albicei: ',albicei
        write (nu_diag,*) 'Albedo, albsnowv: ',albsnowv
        write (nu_diag,*) 'Albedo, albsnowi: ',albsnowi
      endif

  900 format ('  dt=',f7.1,'  ndte=',i3,'  npt=',i6,
     &    '  year_init=',i6,'  istep0=',i8)
  901 format ('  diagfreq=',i3,'  histfreq=',a2,'  dumpfreq=',a2)
  902 format ('  kitd=',i2,'  kdyn=',i2,'  kstrength=',i2,
     &    '  kcolumn=',i2)

      end subroutine input_data

!=======================================================================
!BOP
!
! !IROUTINE: init_state - initialize state for itd 
!
! !INTERFACE:
!
      subroutine init_state
!
! !DESCRIPTION:
!
! Initialize state for the itd model
!
! !REVISION HISTORY:
! 
! author C. M. Bitz
!
! !USES:
!
      use ice_model_size
      use ice_constants
      use ice_flux
      use ice_therm_vertical
      use ice_grid
      use ice_state
      use ice_itd
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  ij              ! horizontal index, combines i and j loops
     &,  k               ! ice layer index
     &,  n               ! thickness category index
     &,  icells          ! number of cells initialized with ice

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) :: 
     &   indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind) ::
     &   slope, Ti
     &,  ainit(ncat)

      !-----------------------------------------------------------------
      ! Initialize grid with no ice.
      !-----------------------------------------------------------------
      do n = 1,ncat
         do j = 1,jmt_local
         do i = 1,imt_local
            aicen(i,j,n) = c0
            vicen(i,j,n) = c0
            vsnon(i,j,n) = c0
            esnon(i,j,n) = c0
         enddo
         enddo
         do j = jlo,jhi
         do i = ilo,ihi
            Tsfcn(i,j,n) = Tf(i,j)  ! Tf not defined for ghost cells
         enddo
         enddo
         call bound(Tsfcn(:,:,n))
      enddo

      do k = 1,ntilay
         do j = 1,jmt_local
         do i = 1,imt_local
            eicen(i,j,k) = c0
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Place ice where ocean surface is cold.
      !-----------------------------------------------------------------

      ! initial category areas in cells with ice
      if (ncat == 1) then
         ainit(1) = .8_dbl_kind
      else
         do n = 1, ncat
            ainit(n) = c1 / real(ncat)
         enddo
      endif

      icells = 0
      do j = jlo,jhi
      do i = ilo,ihi
         if (tmask(i,j)) then
          ! place ice where ocean sfc is cold
            if ( (sst (i,j) <= Tf(i,j)+p2) .and.
     &           (ULAT(i,j) < -40.0_dbl_kind/rad_to_deg .or.
     &            ULAT(i,j) >  40.0_dbl_kind/rad_to_deg)    ) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif               ! cold surface
         endif                  ! tmask
      enddo                     ! i
      enddo                     ! j

      do n = 1,ncat
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            aicen(i,j,n) = ainit(n)
            if (n < ncat) then
               vicen(i,j,n) = p5*(hin_max(n-1) + hin_max(n))
     &                      * aicen(i,j,n) ! m
            else                ! n=ncat
               vicen(i,j,n) = (hin_max(n-1) + c1)*aicen(i,j,n) ! m
            endif
!!!            vsnon(i,j,n) = p2*vicen(i,j,n)   ! m
            
            Tsfcn(i,j,n) = min(Tsmelt, Tair(i,j) - Tffresh) ! deg C

            ! snow
            Ti = min(c0, Tsfcn(i,j,n))
            esnon(i,j,n) = -rhos*(Lfresh - cp_ice*Ti)*vsnon(i,j,n)
         enddo                  ! ij
       
         do k = 1, nilyr
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               ! assume linear temp profile and compute enthalpy
               slope = Tf(i,j) - Tsfcn(i,j,n)
               Ti = Tsfcn(i,j,n) + slope*(real(k)-p5)/real(nilyr)

               eicen(i,j,ilyr1(n)+k-1) = 
     &              -(rhoi * (cp_ice*(Tmlt(k)-Ti) 
     &              + Lfresh*(c1-Tmlt(k)/Ti) - cp_ocn*Tmlt(k)))
     &              * vicen(i,j,n)/real(nilyr)

            enddo               ! ij
         enddo                  ! k
      enddo                     ! n

      ! compute aggregate ice state and open water area
      call aggregate

      call bound_aggregate

      end subroutine init_state

!=======================================================================
!BOP
!
! !IROUTINE: init_flux - initialize fluxes exchanged with coupler
!
! !INTERFACE:
!
      subroutine init_flux
!
! !DESCRIPTION:
!
! Initialize all fluxes exchanged with flux coupler
! and some data derived fields
!
! !REVISION HISTORY:
! 
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_constants
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer i,j

      do j=jlo,jhi
      do i=ilo,ihi
      !-----------------------------------------------------------------
      ! fluxes received
      !-----------------------------------------------------------------
        zlvl  (i,j) = c10             ! atm level height (m)
        uatm  (i,j) = c5              ! wind velocity    (m/s)
        vatm  (i,j) = c5
        potT  (i,j) = 273.0_dbl_kind  ! air potential temperature  (K)
        rhoa  (i,j) = 1.3_dbl_kind    ! air density (kg/m^3)
        fsnow (i,j) = 3.3e-6_dbl_kind ! snowfall rate (kg/m2/s)
        frain (i,j) = c0              ! rainfall rate (kg/m2/s)
        fsw   (i,j) = c0              ! shortwave radiation (W/m^2)
        swvdr (i,j) = 0._dbl_kind     ! shortwave radiation (W/m^2)
        swvdf (i,j) = 0._dbl_kind     ! shortwave radiation (W/m^2)
        swidr (i,j) = 0._dbl_kind     ! shortwave radiation (W/m^2)
        swidf (i,j) = 0._dbl_kind     ! shortwave radiation (W/m^2)
        Tair  (i,j) = 273._dbl_kind   ! air temperature  (K)
        Qa    (i,j) = 0.014_dbl_kind  ! specific humidity (kg/kg)
        flw   (i,j) = 280.0_dbl_kind  ! incoming longwave radiation (W/m^2)
        sss   (i,j) = 34.0_dbl_kind   ! sea surface salinity (o/oo)
        uocn  (i,j) = c0              ! surface ocean currents (m/s)
        vocn  (i,j) = c0
        frzmlt(i,j) = -c2             ! freezing/melting potential (W/m^2)

      !-----------------------------------------------------------------
      ! derived or computed fields
      !-----------------------------------------------------------------

        Tf      (i,j) = -depressT*sss(i,j) ! freezing temp (C)
        sst     (i,j) = Tf(i,j)            ! sea surface temp (C)

        wind    (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, (m/s)

        strocnx (i,j) = c0        ! ice-ocean stress, x-direction (U-cell)
        strocny (i,j) = c0        ! ice-ocean stress, y-direction (U-cell)
        strocnxT(i,j) = c0        ! ice-ocean stress, x-direction (T-cell)
        strocnyT(i,j) = c0        ! ice-ocean stress, y-direction (T-cell)

      enddo
      enddo

      call init_flux_atm
      call init_flux_ocn

      end subroutine init_flux

!=======================================================================
!BOP
!
! !IROUTINE: setup_mpi - initialize mpi
!
! !INTERFACE:
!
      subroutine setup_mpi
!
! !DESCRIPTION:
!
! This routine initializes mpi for either internal parallel
! processing or for message passing with the coupler
!
! !REVISION HISTORY:
! 
! author Elizabeth C. Hunke, LANL
! code originally based on POP routine
!
! !USES:
!
      use ice_mpi_internal
      use ice_coupling
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &  coords1, coords2, n, ilen, jlen
     &, interior_i, interior_j  ! dummies for interior blocks

      master_task = 0

#ifdef coupled
      !  if running in coupled mode
      call ice_coupling_setup('ice',MPI_COMM_ICE)
#else
      !  if running in stand-alone MPI mode
#ifdef _MPI
#if fcd_coupled
#else
      call MPI_INIT(ierr)
#endif
      call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_ICE, ierr)
#endif
#endif

#ifdef _MPI
      call MPI_COMM_SIZE (MPI_COMM_ICE, nb_tasks, ierr)
      call MPI_COMM_RANK (MPI_COMM_ICE, my_task, ierr)

      if (nb_tasks /= nproc_s) then
         write (nu_diag,*) 'nb_tasks, nproc_s =', nb_tasks, nproc_s
         call abort_ice ('nb_tasks must equal nproc_s')
      endif

      if (real(imt_global)/real(nproc_x) /= int(imt_global/nproc_x))
     &  then
        write (nu_diag,*) 'nproc_x, imt_global =', nproc_x, imt_global
        call abort_ice
     &       ('number of pes in x must evenly divide imt_global')
      endif

      if (real(jmt_global)/real(nproc_y) /= int(jmt_global/nproc_y))
     &  then
        write (nu_diag,*) 'nproc_y, jmt_global =', nproc_y, jmt_global
        call abort_ice
     &       ('number of pes in y must evenly divide jmt_global')
      endif

      if ( ierr /= MPI_SUCCESS ) then 
         call abort_ice('(setup_mpi) ERROR after MPI_COMM_xxx')
      endif

      coords1 = mod(my_task,nproc_x)
      coords2 = my_task/nproc_x
      nbr_east = coords2*nproc_x+mod(my_task+1,nproc_x)
      nbr_west = coords2*nproc_x+mod(my_task-1+nproc_x,nproc_x)
      nbr_north = my_task+nproc_x
      nbr_south = my_task-nproc_x
      if (nbr_south < 0) nbr_south = -1
      if (nbr_north > nproc_s-1) nbr_north=-1

      ilen = ihi-ilo+1
      jlen = jhi-jlo+1

      do n=1,nproc_s

      local_start(1,n)=((imt_global-1)/nproc_x+1)*mod((n-1),nproc_x)+1
      local_start(2,n)=((jmt_global-1)/nproc_y+1)*((n-1)/nproc_x)+1

      call MPI_TYPE_VECTOR(jlen, ilen, ilen, 
     &     mpi_integer, mpi_interior_int(n), ierr)
      call MPI_TYPE_COMMIT(mpi_interior_int(n), ierr)

      call MPI_TYPE_VECTOR(jlen, ilen, ilen, 
     &     mpi_real8, mpi_interior_real(n), ierr)
      call MPI_TYPE_COMMIT(mpi_interior_real(n), ierr)

      call MPI_TYPE_VECTOR(jlen, ilen, imt_global, 
     &     mpi_integer, mpi_interior_int_global(n), ierr)
      call MPI_TYPE_COMMIT(mpi_interior_int_global(n), ierr)

      call MPI_TYPE_VECTOR(jlen, ilen, imt_global, 
     &     mpi_real8, mpi_interior_real_global(n), ierr)
      call MPI_TYPE_COMMIT(mpi_interior_real_global(n), ierr)

      enddo

      do n=1,nproc_s
         if (my_task == n-1) then
            write (nu_diag,*) ' my_task,e,w,n,s ', my_task,
     &           nbr_east, nbr_west, nbr_north, nbr_south
         endif
      enddo
      write (nu_diag,*) ''

#else
      ! not MPI
      local_start(1,1)= 1
      local_start(2,1)= 1
      my_task = master_task
      nb_tasks = 1 
#endif

      end subroutine setup_mpi

!=======================================================================

      end module ice_init

!=======================================================================
