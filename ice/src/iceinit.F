c $Id: iceinit.F,v 2.0 1999/12/09 22:25:35 eclare Exp $

c-----------------------------------------------------------------------
c   parameter and variable initializations
c-----------------------------------------------------------------------
      subroutine input_data
c.. Namelist variables, set to default values; may be altered
c.. at run time
      implicit none
      include "ice.h"    ! ice code common blocks
      include "history.h"
      integer n

      namelist /inp/ hice0, hsnow0, hzero, 
     & u0, v0, c0, pst, cst, eyc, dt, 
     & ndte, npt, diagfreq, ktherm, kdata, kread,
     & year, istep0, outfreq, dumpfreq,
     & grid_file, kmt_file, dump_file, restrt_file, history_file

c.. default values
      dt = 7200.        ! time step dt_zeta, s 
      ndte = 120        ! number of subcycles:  ndte=dt/dte
      npt = 20          ! total number of time steps (dt) 
      diagfreq = 1      ! how often diag output is written
      outfreq='0'       ! output frequency
      dumpfreq='0'      ! restart frequency
      u0 = 0.           ! initial u velocity component, m/s
      v0 = 0.           ! initial v velocity component, m/s
      pst = 2.75e4      ! P* pressure coefficient, N/m^2
      cst = 20.         ! c* pressure parameter
      eyc = 0.36        ! E_o coefficient in formula for E
      kread = 0         ! type of read for restart
      ktherm = 0        ! type of thermodynamics 
      kdata = 0         ! type of dynamics data - stand-alone only
      year = 1997       ! initial year
      istep0 = 0        ! number of steps taken in previous integrations,
                        ! real (dumped) or imagined (use to set calendar)

      do n=1,ncat-1
      hzero(n) = 0.10*n ! category discriminating thickness, m
      enddo
      do n=1,ncat
      c0(n) = 0.0       ! initial ice compactness
      hice0(n) = 0.     ! initial ice thicknesses, m
      hsnow0(n) = 0.    ! initial snow thicknesses, m
      enddo
      c0(1) = 0.2       ! initial thin ice compactness
      c0(2) = 0.8       ! initial thick ice compactness
      hice0(1) = 0.0    ! initial thin ice thickness, m
      hice0(2) = 1.0    ! initial thick ice thickness, m
      hsnow0(1) = 0.    ! initial snow thickness, m
      hsnow0(2) = 0.01  ! initial snow thickness, m

      grid_file    = 'grid.d'
      kmt_file     = 'kmt.d'
      dump_file    = 'ice_dump.d'
      restrt_file  = 'ice_restart.d'
      history_file = 'history_file.d'

c.. read from input files
      open  (20,file='ice.inp',status='unknown')       ! parameters
      read  (20,inp)
      close (20)

      iout=.false.
      open  (21,file='icefields.inp',status='unknown')
      read  (21,icefields)
      close (21)

      if (my_task.eq.master_task) then
        write (6,900) dt, ndte, npt, year
        write (6,901) diagfreq, outfreq, dumpfreq
        write (6,910) ktherm, kdata, kread
        write (6,905) u0, v0, eyc, cst, pst
        write (6,914)
        do n=1,ncat-1
          write (6,915) n ,c0(n), hice0(n),hzero(n),hsnow0(n)
        enddo
        write (6,916) ncat,c0(ncat),hice0(ncat),          hsnow0(ncat)

        write (6,*) 'The grid file is ',grid_file
        write (6,*) 'The land mask file is ', kmt_file
        write (6,*) 'The following variables will be',
     &    ' written to the history tape: '
        write (6,*) ' Description                             Units',
     &    '     netCDF variable'
        do n=1,avgsiz
           if (iout(n))   write (6,1000) vdesc(n), vunit(n), vname(n)
        enddo
        write (6,*) ' '
      endif

 900  format ('  dt=',f7.1,'  ndte=',i3,'  npt=',i6,
     1    '  year=',i8)
 901  format ('  diagfreq=',i3,'  outfreq=',a2,'  dumpfreq=',a2)
 910  format ('  ktherm=',i2,'  kdata=',i2,'  kread=',i2)
 905  format ('  u0=',f5.1,'  v0=',f5.1,'  Eyc=',f8.2,
     1    '  cst=',f4.1,'  pst=',1pe8.2)
 914  format (2x,'n',7x,'  c0 ',7x,'hice0',7x,'hzero',7x,'hsnow0')
 915  format (2x,i2,4(5x,f7.2))
 916  format (2x,i2,2(5x,f7.2),17x,f7.2)
1000  format (1x,A40,2x,A10,2x,A10)

      return
      end
c-----------------------------------------------------------------------
      subroutine init_par
c.. Initialize variables
c.. The dynamics component is initialized in init_dyn
      implicit none
      include "ice.h"            ! ice code common blocks
      include "thermw.h"         ! thermodynamics common blocks
      include "history.h"
      integer i, j, k, n
      real S_winton

      do n=1,ncat
       iceruf(n) = 0.0005        ! ice surface roughness, m
      enddo
      
      S_winton = 1.              ! ppt
      Tmlt(1) = -mu_Tf*S_winton  ! C
      Tmlt(2) = 0.               ! C (Winton thermodynamics)

      ! initialize some calendar stuff
      istep1 = 0        ! number of steps at current timestep
                        ! real (dumped) or imagined (use to set calendar)
      istep = 1         ! local timestep number
      time=0.
      yday=1.           ! absolute day number
      mday=1            ! day of the month
      month=1           ! month
      nyr=0             ! year
      idate=00000101    ! date
      sec=0             ! seconds into date

      ! initialize the history array
      aa = 0.0
      avgct = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine init_var
c.. Initialize variables
c.. The dynamics component is initialized in init_dyn
      implicit none
      include "ice.h"            ! ice code common blocks
      include "thermw.h"         ! thermodynamics common blocks
      integer i, j, k, n, nk
      real slope
      real Tij (0:nlyr)          ! temperature of ice layers, C

c$OMP PARALLEL DO PRIVATE(i,j,n)
      do j = jlo,jhi
      do i = ilo,ihi
        strairx(i,j) = 0.        ! wind stress
        strairy(i,j) = 0.
        strocnx(i,j) = 0.        ! ice-ocean stress, x-direction (U-cell)
        strocny(i,j) = 0.        ! ice-ocean stress, y-direction (U-cell)
        do n = 1,ncat
          albedo(i,j,n) = 0.     ! surface albedo
          flatent(i,j,n) = 0.    ! latent heat flux, W/m^2
          fsensible(i,j,n) = 0.  ! sensible heat flux, W/m^2
          flwout(i,j,n) = 0.     ! outgoing longwave radiation, W/m^2
          fswthru(i,j,n) = 0.    ! shortwave through ice to ocean, W/m^2
          fresh(i,j,n) = 0.      ! fresh water flux, kg/m^2 s
          fhnet(i,j,n) = 0.      ! net heat flux to ocean, W/m^2
          fsalt(i,j,n) = 0.      ! salt flux to ocean, kg/m^2 s
        enddo
      enddo
      enddo

c$OMP PARALLEL DO PRIVATE(n,i,j)
      do n=0,ncat
      do j=1,jmt_local
      do i=1,imt_local
       hice (i,j,n) = 0.         ! ice thickness, m
       hsnow(i,j,n) = 0.         ! snow depth, m
       compact(i,j,n) = 0.       ! ice fractional area
       Tsfc(i,j,n) = 0.          ! surface temperature, C
       do k=0,nlyr
         qice(i,j,n,k) = 0.      ! ice enthalpy, J/kg
       enddo
      enddo
      enddo
      enddo

c$OMP PARALLEL DO PRIVATE(i,j,n,k)
      do j = jlo,jhi
      do i = ilo,ihi
        do n = 1,ncat
          if (ULAT(i,j).lt.-65.*3.14/180.or.ULAT(i,j).gt.70.*3.14/180.)
     &    then
             compact(i,j,n) = c0  (n)*hm(i,j)             
             hice (i,j,n) = hice0 (n)*hm(i,j)   ! m
             hsnow(i,j,n) = hsnow0(n)*hm(i,j)   ! m 
          endif
        enddo
      enddo
      enddo

      do n = 0,ncat
        call bound(hice(:,:,n))
        call bound(hsnow(:,:,n))
        call bound(compact(:,:,n))
        call bound(Tsfc(:,:,n))
        call bound(qice(:,:,n,1))
        call bound(qice(:,:,n,2))
      enddo

c..  update for consistency ............
c$OMP PARALLEL DO PRIVATE(i,j,slope,n,Tij)
      do j=jlo,jhi
      do i=ilo,ihi
      do n=1,ncat
      ! thick/thin conversions if necessary
       if (compact(i,j,n).lt.tiny) then
        hice(i,j,n) = 0.
        hsnow(i,j,n) = 0.
       endif
      enddo
      do n=1,ncat-1
       if (hice(i,j,n).gt.hzero(n)) then
        if (compact(i,j,n+1).gt.tiny) then
         hice(i,j,n+1) = hice(i,j,n+1) +
     &    compact(i,j,n)/compact(i,j,n+1)*(hice(i,j,n)-hzero(n))
        else
         hice(i,j,n+1) = hzero(n)
         compact(i,j,n+1) = compact(i,j,n)
        endif
        hice(i,j,n) = hzero(n)
       endif
      enddo

      do n=1,ncat
      ! initialize layers and other thermo info
       if (hice(i,j,n).gt.hzero(1)) then 
        Tsfc(i,j,n) = Tair(i,j) - Tf0kel 
        ! linear temperature profile, deg C
        slope = Tf(i,j) - Tsfc(i,j,n)
        do k=0,nlyr
         Tij(k)=Tsfc(i,j,n)+slope*(k-0.5)/nlyr
        enddo
       elseif (hice(i,j,n).gt.tiny) then
        Tsfc(i,j,n) = Tf(i,j)            
        do k=0,nlyr
         Tij(k) = Tf(i,j)
        enddo
       endif
       if (hice(i,j,n).gt.tiny) then
       ! enthalpy, J/kg     !!! 2-layer ice, Winton model !!!
       qice(i,j,n,1) = -cwater*Tmlt(1) + cfresh*(Tmlt(1)-Tij(1))
     $    + Lfresh*(1.-Tmlt(1)/Tij(1))
       qice(i,j,n,2) = -cfresh*Tij(2) 
       endif
      enddo
      enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine flux_init
c.. initialize all fluxes received from flux coupler
c.. and some data derived fields
      implicit none
      include "ice.h"      ! ice code common blocks
      integer i,j

c$OMP PARALLEL DO PRIVATE(i,j)
         do j=jlo,jhi
         do i=ilo,ihi
        zlvl(i,j) = 10.    ! atm level height (m)
        uatm(i,j) = 5.     ! wind velocity    (m/s)
        vatm(i,j) = 5.
        potT(i,j) = 260.   ! air potential temperature  (K)
        rhoa(i,j) = 1.3    ! air density (kg/m^3)
        snow(i,j) = 1.e-8  ! snowfall rate (m/s)
        fsw(i,j) = 0.      ! shortwave radiation (W/m^2)
        Tair(i,j) = 260.   ! air temperature  (K)
        Qa(i,j) = 0.0014   ! specific humidity (kg/kg)
        flw(i,j) = 140.    ! incoming longwave radiation (W/m^2)
        sst(i,j) = -1.9    ! sea surface temperature (C)
        sss(i,j) = 34.     ! sea surface salinity (o/oo)
        gwatx(i,j) = 0.    ! surface ocean currents (m/s)
        gwaty(i,j) = 0.
        frzmlt(i,j)=-2.    ! freezing/melting potential (W/m^2)

        fw(i,j)=-2.        ! oceanic heat flux (W/m^2)
        Tf(i,j) = -mu_Tf*sss(i,j)                      ! freezing temp (C)
        wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, (m/s)
         enddo
         enddo

       return
       end
c----------------------------------------------------------------------
      subroutine setup_mpi

      implicit none
      include "ice.h"               ! ice code common blocks
      include "ice_mpi.h"           ! common blocks for mpi
#ifdef coupled
      include "mpif.h"              ! MPI library definitions
#else
#ifdef _MPI
      include "mpif.h"              ! MPI library definitions
#endif
#endif
      integer ::
     &  cpl_task                    ! master task for coupler
     &, ice_task                    ! master task for ice
      integer ::
     &  coords1, coords2, n, ilen, jlen

      master_task = 0

#ifdef coupled
      !  if running in coupled mode
      call MPI_INIT(ierr)
      call mpi_coupled ('ice', cpl_task, ice_task, MPI_COMM_ICE)
#else
      !  if running in stand-alone MPI mode
#ifdef _MPI
      call MPI_INIT(ierr)
      call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_ICE, ierr)
#endif
#endif

#ifdef _MPI
      call MPI_COMM_SIZE (MPI_COMM_ICE, nb_tasks, ierr)
      call MPI_COMM_RANK (MPI_COMM_ICE, my_task, ierr)
      if (nb_tasks.ne.nproc_s) 
     &  write (6,*) ' error in nb_tasks ',nb_tasks,nproc_s
      if ( ierr /= MPI_SUCCESS ) then 
         write (6,*) '(ice) ERROR after mpi_init'
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

      return
      end
c----------------------------------------------------------------------
