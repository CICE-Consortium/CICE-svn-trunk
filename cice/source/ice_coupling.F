c $Id: $
c=======================================================================
!---! message passing to/from the flux coupler
!---!
!---! author Elizabeth C. Hunke
!---!        others at NCAR
!---!
!---! Copyright, 2001.  The Regents of the University of California.
!---! This software was produced under a U.S. Government contract 
!---! (W-7405-ENG-36) by Los Alamos National Laboratory, which is 
!---! operated by the University of California for the U.S. Department 
!---! of Energy.  The U.S. Government is licensed to use, reproduce, and 
!---! distribute this software.  Permission is granted to the public to 
!---! copy and use this software without charge, provided that this 
!---! Notice and any statement of authorship are reproduced on all 
!---! copies.  Neither the Government nor the University makes any 
!---! warranty, express or implied, or assumes any liability or 
!---! responsibility for the use of this software.
c=======================================================================

      module ice_coupling

      use ice_model_size
      use ice_fileunits

      implicit none

      integer (kind=int_kind), parameter ::
     &   ncbuffi      = 100 ! size of integer control buffer
     &,  nsnd         = 21  ! number of fields sent to coupler
     &,  nrcv         = 21  ! number of fields sent from coupler
     &,  msgtype_c2ii = 31  ! message type for cpl->ice init
     &,  msgtype_i2ci = 40  ! message type for ice->cpl init
     &,  msgtype_c2i  = 30  ! message type for cpl->ice
     &,  msgtype_i2c  = 41  ! message type for ice->cpl
     &,  nrbuff = nrcv*imt_global*jmt_global ! size of receive buffer
     &,  nsbuff = nsnd*imt_global*jmt_global ! size of send    buffer

      integer (kind=int_kind) ::
     &   cbuffi(ncbuffi)    ! control buffer from cpl
     &,  cpl_task           ! master task id for coupler
     &,  nadv_i             ! number of coupler comms per day

      real (kind=dbl_kind) ::
     &   sbuff(nsnd,imt_global*jmt_global)  ! buffer sent to   coupler

      ! buffer location indexes
      integer (kind=int_kind), parameter ::
        !  ocean states
     &   kbufr_sst     =  1 ! sea surface temp         (K)
     &,  kbufr_sss     =  2 ! sea surface salinity     (o/oo)
     &,  kbufr_uocn    =  3 ! ocean current            (m/s)
     &,  kbufr_vocn    =  4 ! 
        !  atm states
     &,  kbufr_zlvl    =  5 ! atm level height         (m)
     &,  kbufr_uatm    =  6 ! wind                     (m/s)
     &,  kbufr_vatm    =  7 ! 
     &,  kbufr_potT    =  8 ! atm potential T          (K)
     &,  kbufr_Tair    =  9 ! atm temperature          (K) 
     &,  kbufr_Qa      = 10 ! atm specific humidity    (kg/kg)
     &,  kbufr_rhoa    = 11 ! atm air density          (kg/m^3)
        !  more ocean states
     &,  kbufr_tiltx   = 12 ! sea surface slope  
     &,  kbufr_tilty   = 13 !
        !  ocean -> ice flux
     &,  kbufr_fw      = 14 ! ptntl to form/melt ice   (W/m^2)
        !  atm -> ice fluxes
     &,  kbufr_swdidr  = 15 !          near IR,   drct
     &,  kbufr_swdvdr  = 16 ! sw down, vsbl,      drct (W/m^2)
     &,  kbufr_swdidf  = 17 !          near IR,   dffs
     &,  kbufr_swdvdf  = 18 !          vsbl,      dffs
     &,  kbufr_flw     = 19 ! longwave down            (W/m^2)
     &,  kbufr_rain    = 20 ! precip, rain             (kg/m^2 s)
     &,  kbufr_snow    = 21 ! precip, snow             (kg/m^2 s)
        !  ice states
     &,  kbufs_ifrc    =  1 ! ice fraction
     &,  kbufs_tsfc    =  2 ! surface temperature      (K)
     &,  kbufs_alb1    =  3 ! albedo, visible, direct
     &,  kbufs_alb2    =  4 !         near IR, direct
     &,  kbufs_alb3    =  5 !         visible, diffuse
     &,  kbufs_alb4    =  6 !         near IR, diffuse
        !  ice -> atm fluxes
     &,  kbufs_tauxa   =  7 ! wind stress              (N/m^2)
     &,  kbufs_tauya   =  8
     &,  kbufs_lat     =  9 ! latent heat flux         (W/m^2)
     &,  kbufs_sens    = 10 ! sensible heat flux       (W/m^2)
     &,  kbufs_lwup    = 11 ! outgoing longwave radiation (W/m^2)
     &,  kbufs_evap    = 12 ! evaporated water         (kg/m^2 s)
        !  2m atm reference temperature
     &,  kbufs_Tref    = 13 ! (K)
        !  ice -> ocean fluxes
     &,  kbufs_netsw   = 14 ! penetrating shortwave -> ocean (W/m^2)
     &,  kbufs_melth   = 15 ! net ocean heat used      (W/m^2)
     &,  kbufs_meltw   = 16 ! water flux -> ocean      (kg/m^2 s)
     &,  kbufs_salt    = 17 ! salt flux -> ocean       (kg/m^2 s)
     &,  kbufs_tauxo   = 18 ! ice/ocean stress         (N/m^2)
     &,  kbufs_tauyo   = 19
     &,  kbufs_swabs   = 20 ! absorbed shortwave
     &,  kbufs_index   = 21 ! a counter

c=======================================================================

#if coupled

      contains

c=======================================================================

      subroutine init_cpl

!---!-------------------------------------------------------------------
!---! initializes message passing between ice and coupler
!---!-------------------------------------------------------------------

      use ice_domain
      use ice_constants
      use ice_mpi_internal
      use ice_grid
      use ice_calendar

      integer (kind=int_kind), dimension(MPI_STATUS_SIZE,2) :: 
     &   status             ! status array for communications

      integer (kind=int_kind) ::
     &   i, j                              ! generic indices
     &,  imask(imt_global,jmt_global)      ! temporary boundary mask

      real (kind=dbl_kind) ::
     &   TLAT_V(4,imt_global,jmt_global)   ! latitude  of cell vertices
     &,  TLON_V(4,imt_global,jmt_global)   ! longitude of cell vertices
     &,  work(ilo:ihi,jlo:jhi)             ! tmp array for gathers
     &,  workg(imt_global,jmt_global)      ! temporary global work array

      !-----------------------------------------------------------------
      ! receive initial buffer from coupler and ignore content
      !-----------------------------------------------------------------

      if (my_task == master_task) then
      call MPI_RECV(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                       msgtype_c2ii, MPI_COMM_WORLD, status, ierr)

      if (ierr /= MPI_SUCCESS ) then
        write (nu_diag,*) '(init_cpl) ERROR after initial recv'
        call exit_coupler
        stop
      endif

      write(nu_diag,*) '(init_cpl) Received control buffer from coupler'
      endif

      !-----------------------------------------------------------------
      ! initialize and send buffer
      !-----------------------------------------------------------------

      sec = 0                  ! elapsed seconds into date
      nadv_i = nint(secday/dt) ! number coupling steps per day

      cbuffi = 0
      cbuffi( 3) = stop_now    ! stop now flag
      cbuffi( 4) = idate       ! initial date (coded: yyyymmdd)
      cbuffi( 5) = sec         ! elapsed seconds into date
      cbuffi( 7) = imt_global  ! grid size in x-direction
      cbuffi( 8) = jmt_global  ! grid size in y-direction
      cbuffi( 9) = nadv_i      ! ice steps per day
      cbuffi(35) = 0           ! use model restart data for initial state

      if (my_task == master_task) 
     & call MPI_SEND(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                     msgtype_i2ci, MPI_COMM_WORLD, ierr)

      if (ierr /= MPI_SUCCESS ) then
         write (nu_diag,*)'(init_cpl) ERROR after inital send_init'
         call exit_coupler
         stop
      endif

      TLON_V = 1.e30  !!! ech:  WHY ARE WE SENDING THESE?
      TLAT_V = 1.e30

      !-----------------------------------------------------------------
      ! send grid center and vertex coordinates in degrees 
      !-----------------------------------------------------------------

      if (my_task == master_task) then
      call MPI_SEND(TLON_G, imt_global*jmt_global, MPI_DOUBLE_PRECISION,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ierr)

      call MPI_SEND(TLAT_G, imt_global*jmt_global, MPI_DOUBLE_PRECISION,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ierr)

      call MPI_SEND(TLON_V,4*imt_global*jmt_global,MPI_DOUBLE_PRECISION,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ierr)

      call MPI_SEND(TLAT_V,4*imt_global*jmt_global,MPI_DOUBLE_PRECISION,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ierr)
      endif

      !-----------------------------------------------------------------
      ! create and send grid cell area in square radians
      !-----------------------------------------------------------------

      do j=jlo,jhi
        do i=ilo,ihi
          work(i,j) = tarea(i,j)/(radius*radius)
        end do
      end do
      
      call global_gather(workg,work)

      call MPI_SEND(workg, imt_global*jmt_global, 
     & MPI_DOUBLE_PRECISION,cpl_task,msgtype_i2ci,MPI_COMM_WORLD,ierr)

      !-----------------------------------------------------------------
      ! create and send integer mask
      !-----------------------------------------------------------------

      do j=jlo,jhi
        do i=ilo,ihi
          work(i,j) = hm(i,j)
        end do
      end do

      call global_gather(workg,work)

      if (my_task == master_task) then

      do j=1,jmt_global
        do i=1,imt_global
          imask(i,j) = nint(workg(i,j))
        end do
      end do

      call MPI_SEND(imask, imt_global*jmt_global, MPI_INTEGER,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ierr)

      if (ierr /= MPI_SUCCESS ) then
        write (nu_diag,*) '(init_cpl) ERROR after grid info send'
        call exit_coupler
        stop
      endif

      write(nu_diag,*) '(init_cpl) Sent grid info to coupler'
      endif

      !-----------------------------------------------------------------
      ! send initial state info to coupler
      !-----------------------------------------------------------------

      call to_coupler

      end subroutine init_cpl

c=======================================================================

      subroutine from_coupler

!---!-------------------------------------------------------------------
!---! flux coupler -> ice data
!---!-------------------------------------------------------------------

      use ice_domain
      use ice_constants
      use ice_flux
      use ice_timers
      use ice_mpi_internal
      use ice_grid
      use ice_calendar
      use ice_state

      integer (kind=int_kind) ::
     &   i, j, k, n, i1, i2                  ! generic indices
     &,  xlen, cnt, nrecv, trecv, nrs, nrn   ! for packed receives
     &,  request1(nproc_s),status1(MPI_STATUS_SIZE,nproc_s)
     &,  request1a(nproc_s),status2(MPI_STATUS_SIZE)
     &,  status(MPI_STATUS_SIZE)

      real (kind=dbl_kind) ::  
     &   rbuffy (nrcv+1,(ihi-ilo+1)*(jhi-jlo+1),nproc_s) 
     &,  rbuffy1(nrcv+1,(ihi-ilo+1)*(jhi-jlo+1))
     &,  workx, worky
      real (kind=dbl_kind), allocatable :: rbuffp(:,:)  ! packed buffer 

         call ice_timer_start(8)  ! time spent coupling

      !-----------------------------------------------------------------
      ! receive message
      !-----------------------------------------------------------------

      if (my_task == master_task) then
      call MPI_RECV(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                       msgtype_c2i, MPI_COMM_WORLD, status, ierr)

        nrecv = cbuffi(13)
        nrs   = cbuffi(14)
        nrn   = cbuffi(15)
        trecv = nrecv*nrcv

        allocate(rbuffp(nrecv,nrcv))

        call MPI_RECV(rbuffp, trecv, MPI_DOUBLE_PRECISION, cpl_task,
     &                       msgtype_c2i, MPI_COMM_WORLD, status, ierr)

      if (ierr /= MPI_SUCCESS ) then
        write (nu_diag,*) '(from_coupler) ERROR after receive ',ierr
        call exit_coupler
        stop
      endif
c      write (nu_diag,*) '(from_coupler) received message'

      endif

      !-----------------------------------------------------------------
      ! broadcast write_restart flag
      !-----------------------------------------------------------------
      if (cbuffi(21) == 1 .and. new_day) then
        if (my_task == master_task) 
     &  write (nu_diag,*) '(from_coupler) received write restart signal'
        write_restart = 1
      endif
      call ice_bcast_iscalar(write_restart)

      !-----------------------------------------------------------------
      ! broadcast cpl_write_history flag
      !-----------------------------------------------------------------
      if (cbuffi(23) == 1 .and. new_day) then
        if (my_task == master_task) 
     &  write (nu_diag,*) '(from_coupler) received write history signal'
        cpl_write_history = 1
      endif
      call ice_bcast_iscalar(cpl_write_history)

      !-----------------------------------------------------------------
      ! broadcast stop_now flag
      !-----------------------------------------------------------------
      if (cbuffi(3) == 1) then
        if (my_task == master_task) 
     &     write (nu_diag,*) '(from_coupler) received terminate signal'
        stop_now = 1
      endif
      call ice_bcast_iscalar(stop_now)

      !-----------------------------------------------------------------
      ! unpack rbuff into full size
      ! manually distribute data to subdomains for performance
      !-----------------------------------------------------------------

      xlen = ihi-ilo+1

      if (my_task == master_task) then
        do n = 1,nproc_s
        cnt = 0
        do j = jlo,jhi
        do i = ilo,ihi
          i2 = (local_start(2,n)+j-jlo-1)*imt_global + 
     &         (local_start(1,n)+i-ilo)
          i1 = (j-jlo)*xlen + i-ilo+1
          if (i2 <= nrs) then
            cnt = cnt + 1
            do k = 1,nrcv
              rbuffy(k,cnt,n) = rbuffp(i2,k)
            enddo
            rbuffy(nrcv+1,cnt,n) = float(i1)
          elseif (i2 >= nrn) then
            cnt = cnt + 1
            do k = 1,nrcv
              rbuffy(k,cnt,n) = rbuffp(nrs+i2-nrn+1,k)
            enddo
            rbuffy(nrcv+1,cnt,n) = float(i1)
          endif
        enddo
        enddo
        call MPI_ISEND(cnt,1,MPI_INTEGER,n-1,n-1+1000,
     &    MPI_COMM_ICE,request1a(n),ierr)
        call MPI_ISEND(rbuffy(1,1,n),cnt*(nrcv+1),MPI_REAL8,
     &    n-1,n-1,MPI_COMM_ICE,request1(n),ierr)
        enddo
      endif

      ! zero stuff while waiting, only filling in active cells
      do j=jlo,jhi
       do i=ilo:ihi
        zlvl   (i,j) = c0
        uatm   (i,j) = c0
        vatm   (i,j) = c0
        potT   (i,j) = c0
        Tair   (i,j) = c0
        Qa     (i,j) = c0
        rhoa   (i,j) = c0
        swvdr  (i,j) = c0
        swvdf  (i,j) = c0
        swidr  (i,j) = c0
        swidf  (i,j) = c0
        Flw    (i,j) = c0
        Frain  (i,j) = c0
        Fsnow  (i,j) = c0
        sst    (i,j) = c0
        sss    (i,j) = c0
        uocn   (i,j) = c0
        vocn   (i,j) = c0
        strtltx(i,j) = c0
        strtlty(i,j) = c0
        frzmlt (i,j) = c0
       enddo
      enddo

      call MPI_RECV(cnt,1,MPI_INTEGER,master_task,my_task+1000,
     &  MPI_COMM_ICE,status2,ierr)
      call MPI_RECV(rbuffy1(1,1),cnt*(nrcv+1),MPI_REAL8,
     &    master_task,my_task,MPI_COMM_ICE,status2,ierr)

      ! distribute data to subdomains
      do i1 = 1,cnt
        i = mod(nint(rbuffy1(nrcv+1,i1))-1,xlen)+ilo
        j = (nint(rbuffy1(nrcv+1,i1))-1)/xlen+jlo
        !  atm states                                     ! arrival units
        zlvl   (i,j) = rbuffy1(kbufr_zlvl,i1)    ! m
        uatm   (i,j) = rbuffy1(kbufr_uatm,i1)    ! m/s
        vatm   (i,j) = rbuffy1(kbufr_vatm,i1)    ! m/s
        potT   (i,j) = rbuffy1(kbufr_potT,i1)    ! K
        Tair   (i,j) = rbuffy1(kbufr_Tair,i1)    ! K
        Qa     (i,j) = rbuffy1(kbufr_Qa,  i1)    ! kg/kg
        rhoa   (i,j) = rbuffy1(kbufr_rhoa,i1)    ! kg/m^3
        !  atm -> ice fluxes
        swvdr  (i,j) = rbuffy1(kbufr_swdvdr,i1)  ! W/m^2
        swvdf  (i,j) = rbuffy1(kbufr_swdvdf,i1)  ! W/m^2
        swidr  (i,j) = rbuffy1(kbufr_swdidr,i1)  ! W/m^2
        swidf  (i,j) = rbuffy1(kbufr_swdidf,i1)  ! W/m^2
        Flw    (i,j) = rbuffy1(kbufr_flw,   i1)  ! W/m^2
        Frain  (i,j) = rbuffy1(kbufr_rain,  i1)  ! kg/m^2 s
        Fsnow  (i,j) = rbuffy1(kbufr_snow,  i1)  ! kg/m^2 s liquid
        !  ocean states
        sst    (i,j) = rbuffy1(kbufr_sst,  i1)   ! K
        sss    (i,j) = rbuffy1(kbufr_sss,  i1)   ! o/oo
        uocn   (i,j) = rbuffy1(kbufr_uocn, i1)   ! m/s
        vocn   (i,j) = rbuffy1(kbufr_vocn, i1)   ! m/s
        strtltx(i,j) = rbuffy1(kbufr_tiltx,i1)   ! m/m
        strtlty(i,j) = rbuffy1(kbufr_tilty,i1)   ! m/m
        !  ocean -> ice flux
        frzmlt (i,j) = rbuffy1(kbufr_fw,   i1)   ! W/m^2
      enddo

      if (my_task == master_task) then
        deallocate(rbuffp)
        call MPI_WAITALL(nproc_s,request1a,status1,ierr)
        call MPI_WAITALL(nproc_s,request1,status1,ierr)
      endif

      !-----------------------------------------------------------------
      ! rotate zonal/meridional vectors to local coordinates
      ! compute data derived quantities
      !-----------------------------------------------------------------

      do j=jlo,jhi
       do i=ilo,ihi
        ! ocean
!popgrid:   ocn sending fields on u-grid - don't need to rotate
!popgrid        workx      = uocn  (i,j)              ! currents, m/s 
!popgrid        worky      = vocn  (i,j)
!popgrid        uocn(i,j) = workx*cos(ANGLET(i,j))    ! rotate to POP grid 
!popgrid     1            + worky*sin(ANGLET(i,j))    ! T-cell centers
!popgrid        vocn(i,j) = worky*cos(ANGLET(i,j))
!popgrid     1            - workx*sin(ANGLET(i,j))

!popgrid        workx      = strtltx  (i,j)           ! sea sfc tilt, m/m
!popgrid        worky      = strtlty  (i,j)      
!popgrid        strtltx(i,j) = workx*cos(ANGLET(i,j)) ! rotate to POP grid 
!popgrid     1               + worky*sin(ANGLET(i,j)) ! T-cell centers
!popgrid        strtlty(i,j) = worky*cos(ANGLET(i,j))
!popgrid     1               - workx*sin(ANGLET(i,j))

        sst   (i,j) = sst(i,j) - Tffresh    ! sea src temperature (C)
        Tf    (i,j) = -1.8_dbl_kind ! hardwired to ignore SSS (mismatches SST)
c        Tf    (i,j) = -depressT*sss(i,j)      ! freezing temperature (C)
c        Tf    (i,j) = -depressT*max(sss(i,j),ice_ref_salinity) ! freezing T(C)
        enddo
      enddo
      ! interpolate from T-cell centers to U-cell centers
!popgrid:   ocn sending fields on u-grid - don't need to interpolate
!popgrid      call t2ugrid(uocn)
!popgrid      call t2ugrid(vocn)
!popgrid      call t2ugrid(strtltx)
!popgrid      call t2ugrid(strtlty)


      ! atmo variables are needed in T cell centers in subroutine stability,
      ! and are interpolated to the U grid later as necessary
      do j=jlo,jhi
       do i=ilo,ihi
        ! atmosphere
        workx      = uatm(i,j)                ! wind velocity, m/s
        worky      = vatm(i,j) 
        uatm (i,j) = workx*cos(ANGLET(i,j))   ! convert to POP grid
     1             + worky*sin(ANGLET(i,j))   ! note uatm, vatm, wind
        vatm (i,j) = worky*cos(ANGLET(i,j))   !  are on the T-grid here
     1             - workx*sin(ANGLET(i,j))

        wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s

        Fsw  (i,j) = swvdr(i,j) + swvdf(i,j)
     1             + swidr(i,j) + swidf(i,j)
       enddo
      enddo

      time_forc=time

      call ice_timer_stop(8)   ! time spent coupling

      end subroutine from_coupler

c=======================================================================

      subroutine to_coupler

!---!-------------------------------------------------------------------
!---! ice -> flux coupler data
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_domain
      use ice_constants
      use ice_flux
      use ice_timers
      use ice_mpi_internal
      use ice_albedo
      use ice_grid
      use ice_calendar
      use ice_state
      use ice_history

      integer (kind=int_kind) ::
     &   i, j, n, n2, n3           ! generic indices
     &,  nsend, tsend, nsize(nproc_s)
     &,  recv_req1(nproc_s),recv_status1(MPI_STATUS_SIZE,nproc_s)
     &,  recv_req2(nproc_s),recv_status2(MPI_STATUS_SIZE,nproc_s)
     &,  request1,status1(MPI_STATUS_SIZE)
     &,  request2,status2(MPI_STATUS_SIZE)

      real (kind=dbl_kind) ::  
     &   workx, worky              ! tmps for converting grid
     &,  Tsrf  (ilo:ihi,jlo:jhi)   ! surface temperature
     &,  tauxa (ilo:ihi,jlo:jhi)   ! atmo/ice stress
     &,  tauya (ilo:ihi,jlo:jhi)               
     &,  tauxo (ilo:ihi,jlo:jhi)   ! ice/ocean stress
     &,  tauyo (ilo:ihi,jlo:jhi)               
     &,  ailohi(ilo:ihi,jlo:jhi)               

         call ice_timer_start(8)   ! time spent coupling

      do j=jlo,jhi
       do i=ilo,ihi
        ! surface temperature
        Tsrf(i,j)  = Tffresh + Tsfc(i,j)                      !K

        ! wind stress  (on POP T-grid:  convert to lat-lon)
        workx = strairxT(i,j)                               ! N/m^2
        worky = strairyT(i,j)                               ! N/m^2
        tauxa(i,j) = workx*cos(ANGLET(i,j)) - worky*sin(ANGLET(i,j))
        tauya(i,j) = worky*cos(ANGLET(i,j)) + workx*sin(ANGLET(i,j))
        ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
        workx = -strocnxT(i,j)                               ! N/m^2
        worky = -strocnyT(i,j)                               ! N/m^2
        tauxo(i,j) = workx*cos(ANGLET(i,j)) - worky*sin(ANGLET(i,j))
        tauyo(i,j) = worky*cos(ANGLET(i,j)) + workx*sin(ANGLET(i,j))
       enddo
      enddo

      !-----------------------------------------------------------------
      ! gather coupling variables from subdomains
      !-----------------------------------------------------------------

      do j=jlo,jhi
        do i=ilo,ihi
          ailohi(i,j) = aice(i,j)
        enddo
      enddo
      n2 = 0
      do j=jlo,jhi
      do i=ilo,ihi
      if (tmask(i,j) .and. ailohi(i,j) /= c0) then  ! departure units 
         n2 = n2 + 1
         sbuff(kbufs_ifrc, n2) = ailohi      (i,j)  ! none
         sbuff(kbufs_tsfc, n2) = Tsrf        (i,j)  ! K
         sbuff(kbufs_alb1, n2) = alvdr       (i,j)  ! none
         sbuff(kbufs_alb2, n2) = alidr       (i,j)  ! none
         sbuff(kbufs_alb3, n2) = alvdf       (i,j)  ! none
         sbuff(kbufs_alb4, n2) = alidf       (i,j)  ! none
         sbuff(kbufs_lat , n2) = Flatent     (i,j)  ! W/m^2
         sbuff(kbufs_sens, n2) = Fsensible   (i,j)  ! W/m^2
         sbuff(kbufs_lwup, n2) = Flwout      (i,j)  ! W/m^2
         sbuff(kbufs_evap, n2) = evap        (i,j)  ! kg/m^2 s
         sbuff(kbufs_tauxa,n2) = tauxa       (i,j)  ! N/m^2
         sbuff(kbufs_tauya,n2) = tauya       (i,j)  ! N/m^2
         sbuff(kbufs_netsw,n2) = Fswthru     (i,j)  ! W/m^2
         sbuff(kbufs_melth,n2) = Fhnet       (i,j)  ! W/m^2
         sbuff(kbufs_meltw,n2) = Fresh       (i,j)  ! kg/m^2 s
         sbuff(kbufs_salt ,n2) = Fsalt       (i,j)  ! kg/m^2 s
         sbuff(kbufs_tauxo,n2) = tauxo       (i,j)  ! N/m^2
         sbuff(kbufs_tauyo,n2) = tauyo       (i,j)  ! N/m^2
         sbuff(kbufs_Tref, n2) = Tref        (i,j)  ! K
         sbuff(kbufs_swabs,n2) = sabs        (i,j)  ! W/m^2
         sbuff(kbufs_index,n2) = rndex_global(i,j)  ! none
      endif
      enddo
      enddo

      if (my_task .eq. master_task) then
        do n=1,nproc_s
           call MPI_IRECV(nsize(n),1,MPI_INTEGER,n-1,n-1+1000,
     &                    MPI_COMM_ICE,recv_req1(n),ierr)
        end do
      endif

      call MPI_ISEND(n2,1,MPI_INTEGER,master_task,my_task+1000,
     &               MPI_COMM_ICE,request1,ierr)
      call MPI_ISEND(sbuff(1,1),n2*nsnd,
     &     MPI_REAL8,master_task,my_task,MPI_COMM_ICE,request2,ierr)

      if (my_task .eq. master_task) then
        nsend = 0
        do n=1,nproc_s
         call MPI_WAITANY(nproc_s,recv_req1,n3,recv_status1(1,n),ierr)
         call MPI_IRECV(sbuff(1,nsend+1),nsize(n3)*nsnd,
     &        MPI_REAL8,n3-1,n3-1,MPI_COMM_ICE,recv_req2(n),ierr)
         nsend = nsend + nsize(n3)
        end do
        tsend = nsend*nsnd
      endif

      if (my_task .eq. master_task) then
        call MPI_WAITALL(nproc_s,recv_req2,recv_status2,ierr)
      endif
      call MPI_WAIT(request1,status1,ierr)
      call MPI_WAIT(request2,status2,ierr)

      cbuffi( 4) = idate          ! date (coded: yyyymmdd)
      cbuffi( 5) = sec            ! elapsed seconds into date
      cbuffi( 9) = nadv_i         ! ice comm pairs per day
      cbuffi(10) = nsend          ! size of real buffer to send

      !-----------------------------------------------------------------
      ! send message to coupler
      !-----------------------------------------------------------------
   
      if (my_task == master_task) then
      call MPI_SEND(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                     msgtype_i2c, MPI_COMM_WORLD, ierr)

      call MPI_SEND(sbuff, tsend, MPI_DOUBLE_PRECISION, cpl_task,
     &                     msgtype_i2c, MPI_COMM_WORLD, ierr)

      if (ierr .ne. MPI_SUCCESS ) then
        write (nu_diag,*) '(to_coupler) ERROR after coupler send ',ierr
      endif
c        write (nu_diag,*) '(to_coupler) sent message'
      endif

      call ice_timer_stop(8)   ! time spent coupling

      ! for debugging
c      call ice_global_real_minmax((jhi-jlo+1)*(ihi-ilo+1),Fresh,'Fresh')

      end subroutine to_coupler

c=======================================================================

      subroutine exit_coupler

!---!-------------------------------------------------------------------
!---! exit from coupled/MPI environment
!---!-------------------------------------------------------------------

      use ice_kinds_mod
      use ice_model_size
      use ice_domain

      include "mpif.h"         ! MPI library definitions

      integer (kind=int_kind) ::
     &   i, j                  ! generic indices
     &,  ierr                  ! error flag
     &,  status(MPI_STATUS_SIZE)

      if (my_task == master_task) then
      if (cbuffi(3) .eq. 1) then
         write (6,*) '(ice) received final coupler msg'
      else
         write (6,*) '(ice) terminating before coupler'
         call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
      endif
      endif

      call MPI_FINALIZE(ierr)

      write(6,*) '(ice) exit_coupler finished',my_task

      end subroutine exit_coupler

c=======================================================================

      subroutine mpi_coupled (in_model_name, 
     &                        cpl_task, model_task, model_comm)

!---!-------------------------------------------------------------------
!---! this routine queries all the components of the full coupled
!---! system and sets up proper communicators and task ids for each
!---! component of the coupled system
!---!
!---! this routine should be called after mpi_init, but before 
!---! setting up any internal mpi setups (since these will require
!---! the internal communicators returned by this routine)
!---!
!---! code originally based on POP routine
!---!-------------------------------------------------------------------

      include "mpif.h"                  ! MPI library definitions

      character (3), intent(in) :: in_model_name   
                        ! 3-letter identifier (atm, lnd, ocn, ice, cpl)
                        ! for the model calling this routine

      integer, intent(out) ::
     &  cpl_task        ! master task of coupler
     &, model_task      ! master task of model (in MPI_COMM_WORLD)
     &, model_comm      ! communicator for internal model comms

      character (3) :: cmodel   ! model name temporary

      integer, dimension(3) :: range  ! array for creating groups for
                                      !  each coupled model

      integer :: 
     &  MPI_GROUP_WORLD ! group id for MPI_COMM_WORLD
     &, MPI_GROUP_CPL   ! group of processors assigned to cpl
     &, MPI_GROUP_ATM   ! group of processors assigned to atm
     &, MPI_GROUP_OCN   ! group of processors assigned to ocn
     &, MPI_GROUP_ICE   ! group of processors assigned to ice
     &, MPI_GROUP_LND   ! group of processors assigned to lnd
     &, MPI_COMM_CPL    ! communicator for processors assigned to cpl
     &, MPI_COMM_ATM    ! communicator for processors assigned to atm
     &, MPI_COMM_OCN    ! communicator for processors assigned to ocn
     &, MPI_COMM_ICE    ! communicator for processors assigned to ice
     &, MPI_COMM_LND    ! communicator for processors assigned to lnd
     &, n               ! dummy loop counter
     &, ierr            ! error flag for MPI comms
     &, nprocs_cpl      ! total processor count
     &, my_task_coupled ! rank of process in coupled domain
     &, cpl_rank_min, cpl_rank_max   ! processor range for each
     &, ocn_rank_min, ocn_rank_max   !  component model
     &, atm_rank_min, atm_rank_max
     &, ice_rank_min, ice_rank_max
     &, lnd_rank_min, lnd_rank_max

      !-----------------------------------------------------------------
      !     determine processor rank in coupled domain
      !-----------------------------------------------------------------

      call MPI_COMM_RANK  (MPI_COMM_WORLD, my_task_coupled, ierr)

      !-----------------------------------------------------------------
      !     determine which group of processes assigned to each model
      !-----------------------------------------------------------------

      call MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs_cpl, ierr)

      atm_rank_min = nprocs_cpl
      atm_rank_max = 0
      ocn_rank_min = nprocs_cpl
      ocn_rank_max = 0
      ice_rank_min = nprocs_cpl
      ice_rank_max = 0
      lnd_rank_min = nprocs_cpl
      lnd_rank_max = 0
      cpl_rank_min = nprocs_cpl
      cpl_rank_max = 0

      !***
      !*** each processor broadcasts its model to all the processors
      !*** in the coupled domain
      !***

      do n=0,nprocs_cpl-1
        if (n == my_task_coupled) then
          cmodel = in_model_name
        else
          cmodel = 'unk'
        endif

        call MPI_BCAST(cmodel, 3, MPI_CHARACTER, n,
     &                            MPI_COMM_WORLD, ierr)

        select case(cmodel)
        case ('atm')
          atm_rank_min = min(atm_rank_min, n)
          atm_rank_max = max(atm_rank_max, n)
        case ('ocn')
          ocn_rank_min = min(ocn_rank_min, n)
          ocn_rank_max = max(ocn_rank_max, n)
        case ('ice')
          ice_rank_min = min(ice_rank_min, n)
          ice_rank_max = max(ice_rank_max, n)
        case ('lnd')
          lnd_rank_min = min(lnd_rank_min, n)
          lnd_rank_max = max(lnd_rank_max, n)
        case ('cpl')
          cpl_rank_min = min(cpl_rank_min, n)
          cpl_rank_max = max(cpl_rank_max, n)
        case ('dud')
        case default
          write (nu_diag,*)
     &      'Unknown model ',cmodel,' in coupled model domain'
          write (nu_diag,*)
     &      'Model must be atm, cpl, ice, lnd, ocn'
          stop
        end select

      end do

      !-----------------------------------------------------------------
      !     create subroup and communicators for each models internal 
      !     communciations, note that MPI_COMM_CREATE must be called by
      !     all processes in MPI_COMM_WORLD so this must be done by all
      !     models consistently and in the same order.
      !-----------------------------------------------------------------

      call MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr)

      range(3) = 1

      range(1) = atm_rank_min
      range(2) = atm_rank_max
      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_ATM, ierr)

      range(1) = ocn_rank_min
      range(2) = ocn_rank_max
      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_OCN, ierr)

      range(1) = ice_rank_min
      range(2) = ice_rank_max
      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_ICE, ierr)

      range(1) = lnd_rank_min
      range(2) = lnd_rank_max
      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_LND, ierr)

      range(1) = cpl_rank_min
      range(2) = cpl_rank_max
      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_CPL, ierr)

      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_ATM,
     &                      MPI_COMM_ATM, ierr)
 
      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_OCN,
     &                      MPI_COMM_OCN, ierr)
 
      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_ICE,
     &                      MPI_COMM_ICE, ierr)
 
      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_LND,
     &                      MPI_COMM_LND, ierr)
 
      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_CPL,
     &                      MPI_COMM_CPL, ierr)
 
      !-----------------------------------------------------------------
      !     determine coupler process and model processes
      !     assume the first processor in each domain is the task that 
      !     will communicate coupled model messages
      !-----------------------------------------------------------------

      cpl_task = cpl_rank_min

      select case(in_model_name)
      case ('atm')
        model_task = atm_rank_min
        model_comm = MPI_COMM_ATM
      case ('ocn')
        model_task = ocn_rank_min
        model_comm = MPI_COMM_OCN
      case ('ice')
        model_task = ice_rank_min
        model_comm = MPI_COMM_ICE
      case ('lnd')
        model_task = lnd_rank_min
        model_comm = MPI_COMM_LND
      case ('cpl')
        model_task = cpl_rank_min
        model_comm = MPI_COMM_CPL
      case default
        write (nu_diag,*)
     &    'Unknown model ',in_model_name,' in coupled model'
        write (nu_diag,*)
     &    'Model must be atm, cpl, ice, lnd, ocn'
        stop
      end select

      end subroutine mpi_coupled

#endif

c=======================================================================

      end module ice_coupling

c=======================================================================
