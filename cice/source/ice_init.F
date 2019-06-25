c $Id: $
c=======================================================================
!---!   parameter and variable initializations
!---!
!---! authors Elizabeth C. Hunke, LANL
!---!         C. M. Bitz
c=======================================================================

      module ice_init

      use ice_domain
      use ice_therm_driver

      implicit none

c=======================================================================

      contains

c=======================================================================

      subroutine input_data

!---!-------------------------------------------------------------------
!---! Namelist variables, set to default values; may be altered
!---! at run time
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_mechred_cice
      use ice_diagnostics
      use ice_history
      use ice_calendar
      use ice_dyn_evp
      use ice_itd
      use ice_transport

      integer (kind=int_kind) :: nml_error ! namelist i/o error flag

      namelist /ice_nml/ year, istep0, dt, ndte, npt, 
     & diagfreq, histfreq, dumpfreq, restart, print_points,
     & kitd, kdyn, kstrength, evp_damping, 
     & grid_file, kmt_file, dump_file, restrt_file, history_file,
     & grid_type, advection, hist_avg

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------
      year = 1997       ! initial year
      istep0 = 0        ! number of steps taken in previous integrations,
                        ! real (dumped) or imagined (use to set calendar)
      dt = 7200.0_dbl_kind  ! time step, s 
      ndte = 120        ! number of subcycles:  ndte=dt/dte
      npt = 20          ! total number of time steps (dt) 
      diagfreq = 1      ! how often diag output is written
      histfreq='0'       ! output frequency
      dumpfreq='0'      ! restart frequency
      hist_avg = .false. ! write snapshots rather than time-averages
      print_points = .false.     ! if true, print point data
      kitd = 0          ! type of itd conversions (0 = delta f, 1 = linear)
      kdyn = 1          ! type of dynamics (1 = evp)
      restart = .false. ! do not read restart files for initialization
      kstrength = 0     ! 1 = use Rothrock 1975 pressure formulation
      evp_damping = .false. ! if true, use damping procedure in evp dynamics
      advection  = 'mpdata2'  ! second order advection using mpdata
      grid_type    = 'rectangular'  ! define rectangular grid internally

      grid_file    = 'grid.d'
      kmt_file     = 'kmt.d'
      dump_file    = 'ice_dump.d'
      restrt_file  = 'ice_restart.d'
      history_file = 'history_file.d'

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------
        nml_error = -1
        open (21, file='ice_in', status='old')
   10   continue  !*** keep reading until right namelist is found
        read(21, nml=ice_nml,err=10,end=20)
        close(21)
        nml_error = 0
   20   continue

        if (histfreq.eq.'1') hist_avg = .false. ! potential conflict

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) then
        write (6,900) dt, ndte, npt, year, istep0
        write (6,901) diagfreq, histfreq, dumpfreq
        write (6,902) kitd, kdyn, kstrength
        write (6,*) ' '
        if (evp_damping) then
          write (6,*) 'Damping procedure used in evp dynamics'
        else
          write (6,*) 'Damping procedure NOT used in evp dynamics'
        endif
        write (6,*) ' '
        write (6,*) ' '
#ifdef coupled
#ifdef oceanmixed
        write (6,*) 'WARNING:  coupled and oceanmixed flags are BOTH ON'
        write (6,*) 'Ocean data received from coupler will be altered '
        write (6,*) 'by mixed layer routine! WARNING WARNING WARNING '
        write (6,*) ' '
#endif
#endif
        write (6,*) 'Advection algorithm is ',advection
        write (6,*) ' '
        write (6,*) 'Diagnostic data for the grid points specified'
        if (print_points) then
          write (6,*) 'in ice_diagnostics.F will be printed'
        else
          write (6,*) 'in ice_diagnostics.F will NOT be printed'
        endif
        write (6,*) ' '
        if (grid_type .eq. 'displaced_pole') then
          write (6,*) 'Displaced pole grid'
          write (6,*) 'The grid file is ',grid_file
          write (6,*) 'The land mask file is ', kmt_file
        else
          write (6,*) 'Rectangular grid'
        endif

        write (6,*) 'The dump file prefix is ', dump_file
        if (restart) then
          write (6,*) 'The restart file pointer is ', restrt_file
        else
          write (6,*) 'Starting from default initialization '
          write (6,*) 'rather than restart file'
        endif
        write (6,*) 'The history file prefix is ', history_file
        if (hist_avg) then
          write (6,*) 'History data will be averaged over 1 ',histfreq
        else
          write (6,*) 'History data will be snapshots'
        endif
      endif

 900  format ('  dt=',f7.1,'  ndte=',i3,'  npt=',i6,
     &    '  year=',i8,'  istep0=',i8)
 901  format ('  diagfreq=',i3,'  histfreq=',a2,'  dumpfreq=',a2)
 902  format ('  kitd=',i2,'  kdyn=',i2,'  kstrength=',i2)

      end subroutine input_data

c=======================================================================

      subroutine init_state

!---!-------------------------------------------------------------------
!---! Initialize state for the itd model
!---!
!---! author C. M. Bitz
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_constants
      use ice_flux
      use ice_therm_cice
      use ice_grid
      use ice_state
      use ice_itd

      integer (kind=int_kind) :: layer,nc,i,j,k
      real (kind=dbl_kind) ::
     &   slope, Ti
     &,  ainit(ncat)

      do nc=1,ncat
        do j=1,jmt_local
        do i=1,imt_local
          aicen(i,j,nc) = c0
          vicen(i,j,nc) = c0
          vsnon(i,j,nc) = c0
          Tsfcn(i,j,nc) = c0
        enddo
        enddo
        do j=jlo,jhi
        do i=ilo,ihi
          Tsfcn(i,j,nc) = Tf(i,j)
        enddo
        enddo
      enddo

      do layer=1,ntilay
        do j=1,jmt_local
        do i=1,imt_local
            eicen(i,j,layer) = c0
        enddo
        enddo
      enddo

      if (ncat.eq.1) then
      ainit(1)=.8_dbl_kind
      else
      ainit(1)=.2_dbl_kind
      ainit(2)=.3_dbl_kind
      ainit(3)=.3_dbl_kind
      do nc=4,ncat
         ainit(nc)=0.1_dbl_kind
      enddo
      endif

      do j = jlo,jhi
      do i = ilo,ihi
        if (tmask(i,j)) then
        ! place ice where ocean sfc is cold
        if ( (sst (i,j).le.Tf(i,j)+p2) .and.
     &       (ULAT(i,j).lt.-40.0_dbl_kind/rad_to_deg .or.
     &        ULAT(i,j).gt. 40.0_dbl_kind/rad_to_deg)    )    then

          do nc = 1,ncat
            aicen(i,j,nc) = ainit(nc)
            if (nc.lt.ncat) then
              vicen(i,j,nc) = p5*(hin_max(nc-1)
     $           + hin_max(nc))*aicen(i,j,nc)   ! m
            else ! nc=ncat
              vicen(i,j,nc) = (hin_max(nc-1)+0.9_dbl_kind)*aicen(i,j,nc) ! m
            endif
!!!            vsnon(i,j,nc) = p2*vicen(i,j,nc)   ! m
            
            Tsfcn(i,j,nc) = min(Tsmelt,Tair(i,j) - Tffresh)    ! deg C       
            ! make linear temp profile and compute enthalpy
            slope = Tf(i,j) - Tsfcn(i,j,nc)
            k=1
            do layer=layer1(nc),layern(nc)
              Ti=Tsfcn(i,j,nc)+slope*(real(k)-p5)/nlyr

c         eicen(i,j,layer)=(-rLfi - rcpi*(-depressT*salin(layer)-Ti)
c     &   -rLfidepressT*salin(layer)/Ti) *vicen(i,j,nc)/nlyr

              eicen(i,j,layer)=-(rhoi * (cp_ice*(Tmlt(layer)-Ti) 
     &           + Lfresh*(c1-Tmlt(layer)/Ti) - cp_ocn*Tmlt(layer)))
     &   * vicen(i,j,nc)/nlyr

              k=k+1
            enddo ! layer
          ! snow
          Ti = min(c0, Tsfcn(i,j,nc))
          esnon(i,j,nc) = -rhos*(Lfresh - cp_ice*Ti)*vsnon(i,j,nc)
          enddo  ! nc
        endif
        endif
      enddo
      enddo

      ! compute aggregate ice state and open water area
      call aggregate 

      end subroutine init_state

c=======================================================================

      subroutine init_flux

!---!-------------------------------------------------------------------
!---! initialize all fluxes exchanged with flux coupler
!---! and some data derived fields
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_constants
      use ice_flux

      integer i,j

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
      do i=ilo,ihi
      !-----------------------------------------------------------------
      ! fluxes received
      !-----------------------------------------------------------------
        zlvl(i,j) = c10    ! atm level height (m)
        uatm(i,j) = c5     ! wind velocity    (m/s)
        vatm(i,j) = c5
        potT(i,j) = 273.0_dbl_kind   ! air potential temperature  (K)
        rhoa(i,j) = 1.3_dbl_kind    ! air density (kg/m^3)
        Fsnow(i,j) = 3.3e-6_dbl_kind ! snowfall rate (kg/m2/s)
        Frain(i,j) = c0     ! rainfall rate (kg/m2/s)
        Fsw(i,j) = c0      ! shortwave radiation (W/m^2)
        swvdr(i,j) = 0._dbl_kind      ! shortwave radiation (W/m^2)
        swvdf(i,j) = 0._dbl_kind      ! shortwave radiation (W/m^2)
        swidr(i,j) = 0._dbl_kind      ! shortwave radiation (W/m^2)
        swidf(i,j) = 0._dbl_kind      ! shortwave radiation (W/m^2)
        Tair(i,j) = 273._dbl_kind   ! air temperature  (K)
        Qa(i,j) = 0.014_dbl_kind   ! specific humidity (kg/kg)
        Flw(i,j) = 280.0_dbl_kind    ! incoming longwave radiation (W/m^2)
!!!        sst(i,j) = -1.9_dbl_kind    ! sea surface temperature (C) 
        sss(i,j) = 34.0_dbl_kind     ! sea surface salinity (o/oo)
        uocn(i,j) = c0    ! surface ocean currents (m/s)
        vocn(i,j) = c0
        frzmlt(i,j)=-c2    ! freezing/melting potential (W/m^2)

      !-----------------------------------------------------------------
      ! derived or computed fields
      !-----------------------------------------------------------------

        Tf(i,j) = -depressT*sss(i,j)                   ! freezing temp (C)
        sst(i,j) = Tf(i,j)

        wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, (m/s)

        strocnx(i,j) = c0        ! ice-ocean stress, x-direction (U-cell)
        strocny(i,j) = c0        ! ice-ocean stress, y-direction (U-cell)
        strocnxT(i,j) = c0        ! ice-ocean stress, x-direction (T-cell)
        strocnyT(i,j) = c0        ! ice-ocean stress, y-direction (T-cell)

        Fhocn(i,j) = c0  

      call init_flux_atm(i,j)
      call init_flux_ocn(i,j)

      enddo
      enddo

      end subroutine init_flux

c=======================================================================

      subroutine setup_mpi

!---!-------------------------------------------------------------------
!---! this routine initializes mpi for either internal parallel
!---! processing or for message passing with the coupler
!---!
!---! author Elizabeth C. Hunke, LANL
!---! code originally based on POP routine
!---!-------------------------------------------------------------------

      use ice_mpi_internal
      use ice_coupling

      integer (kind=int_kind) ::
     &  ice_task           ! master task for ice
     &, coords1, coords2, n, ilen, jlen
     &, interior_i, interior_j  ! dummies for interior blocks

      master_task = 0

#ifdef coupled
      !  if running in coupled mode
      call MPI_INIT(ierr)
      call mpi_coupled ('ice', cpl_task, ice_task, MPI_COMM_ICE)
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
      if (nb_tasks.ne.nproc_s) 
     &  write (6,*) ' error in nb_tasks ',nb_tasks,nproc_s
      if ( ierr /= MPI_SUCCESS ) then 
         write (6,*) '(setup_mpi) ERROR after MPI_COMM_xxx'
         call MPI_FINALIZE(ierr)
         stop
      endif

      coords1=mod(my_task,nproc_x)
      coords2=my_task/nproc_x
      nbr_east = coords2*nproc_x+mod(my_task+1,nproc_x)
      nbr_west = coords2*nproc_x+mod(my_task-1+nproc_x,nproc_x)
      nbr_north = my_task+nproc_x
      nbr_south = my_task-nproc_x
      if (nbr_south.lt.0) nbr_south = -1
      if (nbr_north.gt.nproc_s-1) nbr_north=-1

      ilen=ihi-ilo+1
      jlen=jhi-jlo+1

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
      if (my_task.eq.n-1) then
      write (6,*) ' my_task,e,w,n,s ',my_task,nbr_east,nbr_west,
     &  nbr_north,nbr_south
      endif
      enddo
#else
      my_task = master_task
      nb_tasks = 1 
#endif

      end subroutine setup_mpi

c=======================================================================

      end module ice_init

c=======================================================================

