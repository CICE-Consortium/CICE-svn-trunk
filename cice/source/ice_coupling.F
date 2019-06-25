c $Id: $
c=======================================================================
!---! message passing to/from the flux coupler
!---!
!---! author Elizabeth C. Hunke
!---!
!---! Copyright, 1999.  The Regents of the University of California.
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

      implicit none

      integer (kind=int_kind), parameter ::
     &   ncbuffi = 100       ! size of integer control buffer
     &,  nsnd = 19          ! number of fields sent to coupler
     &,  nrcv = 21          ! number of fields sent from coupler
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
     &   rbuff(imt_global,jmt_global,nrcv)  ! buffer rcvd from coupler
     &,  sbuff(imt_global,jmt_global,nsnd)  ! buffer sent to   coupler

      ! buffer location indexes
      integer (kind=int_kind), parameter ::
        !  ocean states
     &    kbufr_sst     =  1     ! sea surface temp         (K)
     &,   kbufr_sss     =  2     ! sea surface salinity     (o/oo)
     &,   kbufr_uocn    =  3     ! ocean current            (m/s)
     &,   kbufr_vocn    =  4     ! 
        !  atm states
     &,   kbufr_zlvl    =  5     ! atm level height         (m)
     &,   kbufr_uatm    =  6     ! wind                     (m/s)
     &,   kbufr_vatm    =  7     ! 
     &,   kbufr_potT    =  8     ! atm potential T          (K)
     &,   kbufr_Tair    =  9     ! atm temperature          (K) 
     &,   kbufr_Qa      = 10     ! atm specific humidity    (kg/kg)
     &,   kbufr_rhoa    = 11     ! atm air density          (kg/m^3)
        !  more ocean states
     &,   kbufr_tiltx   = 12     ! sea surface slope  
     &,   kbufr_tilty   = 13     !
        !  ocean -> ice flux
     &,   kbufr_fw      = 14     ! ptntl to form/melt ice   (W/m^2)
        !  atm -> ice fluxes
     &,   kbufr_swdidr  = 15     !          near IR,   drct
     &,   kbufr_swdvdr  = 16     ! sw down, vsbl,      drct (W/m^2)
     &,   kbufr_swdidf  = 17     !          near IR,   dffs
     &,   kbufr_swdvdf  = 18     !          vsbl,      dffs
     &,   kbufr_flw     = 19     ! longwave down            (W/m^2)
     &,   kbufr_rain    = 20     ! precip, rain             (kg/m^2 s)
     &,   kbufr_snow    = 21     ! precip, snow             (kg/m^2 s)
        !  ice states
     &,   kbufs_ifrc    =  1     ! ice fraction
     &,   kbufs_tsfc    =  2     ! surface temperature      (K)
     &,   kbufs_alb1    =  3     ! albedo, visible, direct
     &,   kbufs_alb2    =  4     !         near IR, direct
     &,   kbufs_alb3    =  5     !         visible, diffuse
     &,   kbufs_alb4    =  6     !         near IR, diffuse
        !  ice -> atm fluxes
     &,   kbufs_tauxa   =  7     ! wind stress              (N/m^2)
     &,   kbufs_tauya   =  8
     &,   kbufs_lat     =  9     ! latent heat flux         (W/m^2)
     &,   kbufs_sens    = 10     ! sensible heat flux       (W/m^2)
     &,   kbufs_lwup    = 11     ! outgoing longwave radiation (W/m^2)
     &,   kbufs_evap    = 12     ! evaporated water         (kg/m^2 s)
        !  2m atm reference temperature
     &,   kbufs_Tref    = 13     ! (K)
        !  ice -> ocean fluxes
     &,   kbufs_netsw   = 14     ! penetrating shortwave -> ocean (W/m^2)
     &,   kbufs_melth   = 15     ! net ocean heat used      (W/m^2)
     &,   kbufs_meltw   = 16     ! water flux -> ocean      (kg/m^2 s)
     &,   kbufs_salt    = 17     ! salt flux -> ocean       (kg/m^2 s)
     &,   kbufs_tauxo   = 18     ! ice/ocean stress         (N/m^2)
     &,   kbufs_tauyo   = 19

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
     &   status                 ! status array for communications

      integer (kind=int_kind) ::
     &   i, j                          ! generic indices
     &,  imask(imt_global,jmt_global)  ! temporary boundary mask

      real (kind=dbl_kind) ::
     &   TLAT_V(4,imt_global,jmt_global) ! latitude  of cell vertices
     &,  TLON_V(4,imt_global,jmt_global) ! longitude of cell vertices
     &,  hmask(ilo:ihi,jlo:jhi)          ! tmp array for tmask gather
     &,  tareag_rad(imt_global,jmt_global) ! Grid cell area in rad^2 
     &,  work(ilo:ihi,jlo:jhi)           ! tmp array for gathers
     &,  workg(imt_global,jmt_global)    ! temporary global work array

      !-----------------------------------------------------------------
      ! receive initial buffer from coupler and ignore content
      !-----------------------------------------------------------------

      if (my_task == master_task) then
      call MPI_RECV(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                       msgtype_c2ii, MPI_COMM_WORLD, status, ierr)

      if (ierr /= MPI_SUCCESS ) then
        write (6,*) '(init_cpl) ERROR after initial recv'
        call exit_coupler
        stop
      endif

      write(*,*) '(init_cpl) Received control buffer from coupler'
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
         write (6,*)'(init_cpl) ERROR after inital send_init'
         call exit_coupler
         stop
      endif

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
      ! create tareag_rad - grid cell area in square radians
      !-----------------------------------------------------------------

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
        do i=ilo,ihi
          work(i,j) = tarea(i,j)/(radius*radius)
        end do
      end do
      
      call global_gather(tareag_rad,work)

      !-----------------------------------------------------------------
      ! create imask - use workg as temp space
      !-----------------------------------------------------------------

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
        do i=ilo,ihi
          hmask(i,j) = hm(i,j)
        end do
      end do

      call global_gather(workg,hmask)

      if (my_task == master_task) then

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_global
        do i=1,imt_global
          imask(i,j) = nint(workg(i,j))
        end do
      end do

      !-----------------------------------------------------------------
      ! send tareag_rad and imask
      !-----------------------------------------------------------------

      call MPI_SEND(tareag_rad, imt_global*jmt_global, 
     & MPI_DOUBLE_PRECISION,cpl_task,msgtype_i2ci,MPI_COMM_WORLD,ierr)

      call MPI_SEND(imask, imt_global*jmt_global, MPI_INTEGER,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ierr)

      if (ierr /= MPI_SUCCESS ) then
        write (6,*) '(init_cpl) ERROR after initial send (grid info)'
        call exit_coupler
        stop
      endif

      write(*,*) '(init_cpl) Sent grid info to coupler'
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
     &   i, j                      ! generic indices
     &,  status(MPI_STATUS_SIZE)

      real (kind=dbl_kind) ::  
     &   workx, worky

         call ice_timer_start(8)  ! time spent coupling

      !-----------------------------------------------------------------
      ! receive message
      !-----------------------------------------------------------------

      if (my_task == master_task) then
      call MPI_RECV(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                       msgtype_c2i, MPI_COMM_WORLD, status, ierr)

      call MPI_RECV(rbuff, nrbuff, MPI_DOUBLE_PRECISION, cpl_task,
     &                       msgtype_c2i, MPI_COMM_WORLD, status, ierr)

      if (ierr /= MPI_SUCCESS ) then
        write (6,*) '(from_coupler) ERROR after receive ',ierr
        call exit_coupler
        stop
      endif
c      write (6,*) '(from_coupler) received message'

      endif

      !-----------------------------------------------------------------
      ! broadcast write_restart flag
      !-----------------------------------------------------------------
      if (cbuffi(21) == 1 .and. new_day) then
        if (my_task == master_task) 
     &     write (*,*) '(from_coupler) received write restart signal'
        write_restart = 1
      endif
      call ice_bcast_iscalar(write_restart)


      !-----------------------------------------------------------------
      ! broadcast cpl_write_history flag
      !-----------------------------------------------------------------
      if (cbuffi(23) == 1 .and. new_day) then
        if (my_task == master_task) 
     &     write (*,*) '(from_coupler) received write history signal'
        cpl_write_history = 1
      endif
      call ice_bcast_iscalar(cpl_write_history)

      !-----------------------------------------------------------------
      ! broadcast stop_now flag
      !-----------------------------------------------------------------
      if (cbuffi(3) == 1) then
        if (my_task == master_task) 
     &     write (*,*) '(from_coupler) received terminate signal'
        stop_now = 1
      endif
      call ice_bcast_iscalar(stop_now)

      !-----------------------------------------------------------------
      ! distribute data to subdomains
      !-----------------------------------------------------------------

        !  atm states                                     ! arrival units
      call global_scatter(rbuff(1,1,kbufr_zlvl)  ,zlvl )  ! m
      call global_scatter(rbuff(1,1,kbufr_uatm)  ,uatm )  ! m/s
      call global_scatter(rbuff(1,1,kbufr_vatm)  ,vatm )  ! m/s
      call global_scatter(rbuff(1,1,kbufr_potT)  ,potT )  ! K
      call global_scatter(rbuff(1,1,kbufr_Tair)  ,Tair )  ! K
      call global_scatter(rbuff(1,1,kbufr_Qa)    ,Qa   )  ! kg/kg
      call global_scatter(rbuff(1,1,kbufr_rhoa)  ,rhoa )  ! kg/m^3
        !  atm -> ice fluxes
      call global_scatter(rbuff(1,1,kbufr_swdvdr),swvdr)  ! W/m^2
      call global_scatter(rbuff(1,1,kbufr_swdvdf),swvdf)  ! W/m^2
      call global_scatter(rbuff(1,1,kbufr_swdidr),swidr)  ! W/m^2
      call global_scatter(rbuff(1,1,kbufr_swdidf),swidf)  ! W/m^2
      call global_scatter(rbuff(1,1,kbufr_flw)   ,Flw  )  ! W/m^2
      call global_scatter(rbuff(1,1,kbufr_rain)  ,Frain)  ! kg/m^2 s
      call global_scatter(rbuff(1,1,kbufr_snow)  ,Fsnow)  ! kg/m^2 s liquid
        !  ocean states
      call global_scatter(rbuff(1,1,kbufr_sst)   ,sst  )  ! K
      call global_scatter(rbuff(1,1,kbufr_sss)   ,sss  )  ! o/oo
      call global_scatter(rbuff(1,1,kbufr_uocn)  ,uocn)   ! m/s
      call global_scatter(rbuff(1,1,kbufr_vocn)  ,vocn)   ! m/s
      call global_scatter(rbuff(1,1,kbufr_tiltx) ,strtltx)! m/m
      call global_scatter(rbuff(1,1,kbufr_tilty) ,strtlty)! m/m
        !  ocean -> ice flux
      call global_scatter(rbuff(1,1,kbufr_fw)    ,frzmlt) ! W/m^2

      !-----------------------------------------------------------------
      ! rotate zonal/meridional vectors to local coordinates
      ! compute data derived quantities
      !-----------------------------------------------------------------

      ! interpolate ocean dynamics variables from T-cell centers to 
      ! U-cell centers    ! NOTE ANGLE is on the U grid !
      call t2ugrid(uocn)
      call t2ugrid(vocn)
      call t2ugrid(strtltx)
      call t2ugrid(strtlty)

c$OMP PARALLEL DO PRIVATE(i,j,workx,worky)
      do j=jlo,jhi
       do i=ilo,ihi
        ! ocean
        workx      = uocn  (i,j)            ! currents, m/s 
        worky      = vocn  (i,j)
        uocn(i,j) = workx*cos(ANGLE(i,j))   ! convert to POP grid 
     1             + worky*sin(ANGLE(i,j))
        vocn(i,j) = worky*cos(ANGLE(i,j))
     1             - workx*sin(ANGLE(i,j))

        workx      = strtltx  (i,j)            ! sea sfc tilt, m/m
        worky      = strtlty  (i,j)      
        strtltx(i,j) = workx*cos(ANGLE(i,j))   ! convert to POP grid 
     1             + worky*sin(ANGLE(i,j))
        strtlty(i,j) = worky*cos(ANGLE(i,j))
     1             - workx*sin(ANGLE(i,j))

        sst   (i,j) = sst(i,j) - Tffresh    ! sea src temperature (C)
        Tf    (i,j) = -1.8_dbl_kind ! hardwired to ignore SSS (mismatches SST)
c        Tf    (i,j) = -depressT*sss(i,j)      ! freezing temperature (C)
c        Tf    (i,j) = -depressT*max(sss(i,j),ice_ref_salinity) ! freezing T(C)
        enddo
      enddo

      ! atmo variables are needed in T cell centers in subroutine stability,
      ! and are interpolated to the U grid later as necessary
c$OMP PARALLEL DO PRIVATE(i,j,workx,worky)
      do j=jlo,jhi
       do i=ilo,ihi
        ! atmosphere
        workx      = uatm(i,j)               ! wind velocity, m/s
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

      integer (kind=int_kind) ::
     &   i, j                      ! generic indices

      real (kind=dbl_kind) ::  
     &   workx, worky              ! tmps for converting grid
     &,  Tsrf (ilo:ihi,jlo:jhi)       ! surface temperature
     &,  tauxa(ilo:ihi,jlo:jhi)       ! atmo/ice stress
     &,  tauya(ilo:ihi,jlo:jhi)               
     &,  tauxo(ilo:ihi,jlo:jhi)       ! ice/ocean stress
     &,  tauyo(ilo:ihi,jlo:jhi)               
     &,  ailohi (ilo:ihi,jlo:jhi)               

         call ice_timer_start(8)  ! time spent coupling

c$OMP PARALLEL DO PRIVATE(i,j,workx,worky)
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
                                                        ! departure units 
      do j=jlo,jhi
        do i=ilo,ihi
          ailohi(i,j) = aice(i,j)
        enddo
      enddo
      call global_gather(sbuff(1,1,kbufs_ifrc), ailohi   )  ! none
      call global_gather(sbuff(1,1,kbufs_tsfc), Tsrf     )  ! K
      call global_gather(sbuff(1,1,kbufs_alb1), alvdr    )  ! none
      call global_gather(sbuff(1,1,kbufs_alb2), alidr    )  ! none
      call global_gather(sbuff(1,1,kbufs_alb3), alvdf    )  ! none
      call global_gather(sbuff(1,1,kbufs_alb4), alidf    )  ! none
      call global_gather(sbuff(1,1,kbufs_lat ), Flatent  )  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_sens), Fsensible)  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_lwup), Flwout   )  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_evap), evap )  ! kg/m^2 s
      call global_gather(sbuff(1,1,kbufs_tauxa),tauxa)  ! N/m^2
      call global_gather(sbuff(1,1,kbufs_tauya),tauya)  ! N/m^2
      call global_gather(sbuff(1,1,kbufs_netsw),Fswthru)  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_melth),Fhnet  )  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_meltw),Fresh  )  ! kg/m^2 s
      call global_gather(sbuff(1,1,kbufs_salt ),Fsalt  )  ! kg/m^2 s
      call global_gather(sbuff(1,1,kbufs_tauxo),tauxo)  ! N/m^2
      call global_gather(sbuff(1,1,kbufs_tauyo),tauyo)  ! N/m^2
      call global_gather(sbuff(1,1,kbufs_Tref), Tref )  ! K

      cbuffi( 4) = idate          ! date (coded: yyyymmdd)
      cbuffi( 5) = sec            ! elapsed seconds into date
      cbuffi( 7) = imt_local      ! full size in x-dir of ice grid
      cbuffi( 8) = jmt_local      ! full size in y-dir of ice grid
      cbuffi( 9) = nadv_i         ! ice comm pairs per day

      !-----------------------------------------------------------------
      ! send message to coupler
      !-----------------------------------------------------------------
   
      if (my_task == master_task) then
      call MPI_SEND(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                     msgtype_i2c, MPI_COMM_WORLD, ierr)

      call MPI_SEND(sbuff, nsbuff, MPI_DOUBLE_PRECISION, cpl_task,
     &                     msgtype_i2c, MPI_COMM_WORLD, ierr)

      if (ierr .ne. MPI_SUCCESS ) then
         write (6,*) '(to_coupler) ERROR after coupler send ',ierr
      endif
c      write (6,*) '(to_coupler) sent message'
      endif

      call ice_timer_stop(8)   ! time spent coupling

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
     &   i, j                      ! generic indices
     &,  ierr  ! error flag
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
     &  cpl_task           ! master task of coupler
     &, model_task         ! master task of model (in MPI_COMM_WORLD)
     &, model_comm         ! communicator for internal model comms

      character (3) :: cmodel   ! model name temporary

      integer, dimension(3) :: range  ! array for creating groups for
                                      !  each coupled model

      integer :: 
     &  MPI_GROUP_WORLD  ! group id for MPI_COMM_WORLD
     &, MPI_GROUP_CPL    ! group of processors assigned to cpl
     &, MPI_GROUP_ATM    ! group of processors assigned to atm
     &, MPI_GROUP_OCN    ! group of processors assigned to ocn
     &, MPI_GROUP_ICE    ! group of processors assigned to ice
     &, MPI_GROUP_LND    ! group of processors assigned to lnd
     &, MPI_COMM_CPL     ! communicator for processors assigned to cpl
     &, MPI_COMM_ATM     ! communicator for processors assigned to atm
     &, MPI_COMM_OCN     ! communicator for processors assigned to ocn
     &, MPI_COMM_ICE     ! communicator for processors assigned to ice
     &, MPI_COMM_LND     ! communicator for processors assigned to lnd
     &, n                ! dummy loop counter
     &, ierr             ! error flag for MPI comms
     &, nprocs_cpl       ! total processor count
     &, my_task_coupled  ! rank of process in coupled domain
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
          write (6,*)'Unknown model ',cmodel,' in coupled model domain'
          write (6,*)'Model must be atm, cpl, ice, lnd, ocn'
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
        write (6,*)'Unknown model ',in_model_name,' in coupled model'
        write (6,*)'Model must be atm, cpl, ice, lnd, ocn'
        stop
      end select

      end subroutine mpi_coupled

#endif

c=======================================================================

      end module ice_coupling

c=======================================================================
