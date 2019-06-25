c $Id: coupling.F,v 2.0 1999/12/09 22:25:34 eclare Exp $

c-----------------------------------------------------------------------
c   message passing to/from the flux coupler
c-----------------------------------------------------------------------
      subroutine init_cpl
c.. initializes message passing between ice and coupler
      implicit none
      include "ice.h"          ! ice code common blocks
      include "mpif.h"         ! MPI library definitions
      include "coupling.h"     ! coupling common blocks

      integer, dimension(MPI_STATUS_SIZE,2) :: status 
                                            ! status array for communications
      integer stop_now
      integer ier                           ! error flag for messages
      integer i, j                          ! generic indices
      integer imask(imt_global,jmt_global)  ! temporary boundary mask

      real pi2,                             ! 2*pi
     &     rad2deg                          ! radian to degrees factor

      real LONSW(imt_local,jmt_local),      ! longitude of SW corner
     &     LONSE(imt_local,jmt_local),      ! longitude of SE corner
     &     LONNW(imt_local,jmt_local),      ! longitude of NW corner
     &     LONNE(imt_local,jmt_local)       ! longitude of NE corner

      real TLAT_G(imt_global,jmt_global),   ! latitude  of cell center
     &     TLON_G(imt_global,jmt_global),   ! longitude of cell center
     &     TLAT_V(4,imt_global,jmt_global), ! latitude  of cell vertices
     &     TLON_V(4,imt_global,jmt_global)  ! longitude of cell vertices

c.. exchange initial messages ..........................................

c.. receive coupler message and ignore content

      if (my_task == master_task) then
      call MPI_RECV(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                       msgtype_c2ii, MPI_COMM_WORLD, status, ier)

      if (ier /= MPI_SUCCESS ) then
        write (6,*) '(ice) ERROR after initial recv'
        call MPI_FINALIZE(ierr)
        stop
      endif
      endif

      write(*,*) '(ice) Received control buffer from coupler'

c.. fill buffer and send initial message to coupler

      sec = 0                  ! elapsed seconds into date
      nadv_i = nint(86400./dt) ! number coupling steps per day

      cbuffi = 0
      cbuffi( 3) = stop_now    ! stop now flag
      cbuffi( 4) = idate       ! initial date (coded: yyyymmdd)
      cbuffi( 5) = sec         ! elapsed seconds into date
      cbuffi( 7) = imt_global  ! grid size in x-direction
      cbuffi( 8) = jmt_global  ! grid size in y-direction
      cbuffi( 9) = nadv_i      ! ice steps per day
      cbuffi(35) = 0           ! use model restart data for initial state

c.. send buffer to coupler

      if (my_task == master_task) 
     & call MPI_SEND(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                     msgtype_i2ci, MPI_COMM_WORLD, ier)

      if (ierr /= MPI_SUCCESS ) then
         write (6,*)'(ice) ERROR after inital send_init'
         call MPI_FINALIZE(ierr)
         stop
      endif

c.. compute grid vertex coordinates
c.. for longitude, must be careful of longitude branch cut (eg 0,2pi)

      pi2 = 8.0*atan(1.0)
      rad2deg = 360./pi2

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
        do i=ilo,ihi
          LONSE(i,j) = ULONG(i  ,j-1)
          LONNW(i,j) = ULONG(i-1,j  )
          LONSW(i,j) = ULONG(i-1,j-1)

          if (LONSE(i,j) - ULONG(i,j) .gt.  4.0) LONSE(i,j) = 
     &                                           LONSE(i,j) - pi2
          if (LONSE(i,j) - ULONG(i,j) .lt. -4.0) LONSE(i,j) = 
     &                                           LONSE(i,j) + pi2
          if (LONNW(i,j) - ULONG(i,j) .gt.  4.0) LONNW(i,j) = 
     &                                           LONNW(i,j) - pi2
          if (LONNW(i,j) - ULONG(i,j) .lt. -4.0) LONNW(i,j) = 
     &                                           LONNW(i,j) + pi2
          if (LONSW(i,j) - ULONG(i,j) .gt.  4.0) LONSW(i,j) = 
     &                                           LONSW(i,j) - pi2
          if (LONSW(i,j) - ULONG(i,j) .lt. -4.0) LONSW(i,j) = 
     &                                           LONSW(i,j) + pi2
 
          LONNE(i,j) = rad2deg*ULONG(i,j)
          LONNW(i,j) = rad2deg*LONNW(i,j)
          LONSW(i,j) = rad2deg*LONSW(i,j)
          LONSE(i,j) = rad2deg*LONSE(i,j)
        end do
      end do

      call global_gather(TLON_G,LONSW)

      if (my_task == master_task) then
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_global
        do i=1,imt_global
          TLON_V(1,i,j) = TLON_G(i,j)
        end do
      end do
      endif

      call global_gather(TLON_G,LONSE)

      if (my_task == master_task) then
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_global
        do i=1,imt_global
          TLON_V(2,i,j) = TLON_G(i,j)
        end do
      end do
      endif

      call global_gather(TLON_G,LONNE)

      if (my_task == master_task) then
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_global
        do i=1,imt_global
          TLON_V(3,i,j) = TLON_G(i,j)
        end do
      end do
      endif

      call global_gather(TLON_G,LONNW)

      if (my_task == master_task) then
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_global
        do i=1,imt_global
          TLON_V(4,i,j) = TLON_G(i,j)
          TLON_G(i,j)   = 0.25*(TLON_V(1,i,j) + TLON_V(2,i,j) +
     &                          TLON_V(3,i,j) + TLON_V(4,i,j))
        end do
      end do
  
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_local
        do i=1,imt_local
          LONNE(i,j) = rad2deg*ULAT(i  ,j  )
          LONSE(i,j) = rad2deg*ULAT(i  ,j-1)
          LONNW(i,j) = rad2deg*ULAT(i-1,j  )
          LONSW(i,j) = rad2deg*ULAT(i-1,j-1)
        end do
      end do
      endif

      call global_gather(TLAT_G,LONSW)

      if (my_task == master_task) then
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_global
        do i=1,imt_global
          TLAT_V(1,i,j) = TLAT_G(i,j)
        end do
      end do
      endif

      call global_gather(TLAT_G,LONSE)

      if (my_task == master_task) then
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_global
        do i=1,imt_global
          TLAT_V(2,i,j) = TLAT_G(i,j)
        end do
      end do
      endif

      call global_gather(TLAT_G,LONNE)

      if (my_task == master_task) then
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_global
        do i=1,imt_global
          TLAT_V(3,i,j) = TLAT_G(i,j)
        end do
      end do
      endif

      call global_gather(TLAT_G,LONNW)

      if (my_task == master_task) then
c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_global
        do i=1,imt_global
          TLAT_V(4,i,j) = TLAT_G(i,j)
          TLAT_G(i,j)   = 0.25*(TLAT_V(1,i,j) + TLAT_V(2,i,j) +
     &                          TLAT_V(3,i,j) + TLAT_V(4,i,j))
        end do
      end do
      endif

c.. now send grid information to coupler

      if (my_task == master_task) then
      call MPI_SEND(TLON_G, imt_global*jmt_global, MPI_DOUBLE_PRECISION,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ier)

      call MPI_SEND(TLAT_G, imt_global*jmt_global, MPI_DOUBLE_PRECISION,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ier)

      call MPI_SEND(TLON_V,4*imt_global*jmt_global,MPI_DOUBLE_PRECISION,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ier)

      call MPI_SEND(TLAT_V,4*imt_global*jmt_global,MPI_DOUBLE_PRECISION,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ier)
      endif

c.. create and send imask - use TLON_G as temp space

      call global_gather(TLON_G,hm)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jmt_global
        do i=1,imt_global
          imask(i,j) = nint(TLON_G(i,j))
        end do
      end do

      if (my_task == master_task) then
      call MPI_SEND(imask, imt_global*jmt_global, MPI_INTEGER,
     &              cpl_task, msgtype_i2ci, MPI_COMM_WORLD, ier)

      write(*,*) '(ice) Sent grid info to coupler'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine from_coupler
c.. flux coupler -> ice data
      implicit none
      include "ice.h"          ! ice code common blocks
      include "mpif.h"         ! MPI library definitions
      include "coupling.h"     ! coupling common blocks
      include "thermw.h"       ! thermodynamics common blocks

      integer i, j                      ! generic indices
      integer ier                       ! error flag
      integer status(MPI_STATUS_SIZE)
      real    workx, worky
      real    cdepth, chtd(ilo:ihi,jlo:jhi), cht, tscale

c.. receive message

      if (my_task == master_task) then
      call MPI_RECV(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                       msgtype_c2i, MPI_COMM_WORLD, status, ier)

      call MPI_RECV(rbuff, nrbuff, MPI_DOUBLE_PRECISION, cpl_task,
     &                       msgtype_c2i, MPI_COMM_WORLD, status, ier)

      if (ier.ne.MPI_SUCCESS) then
        write (6,*) '(ice) ERROR after receive ',ier
        STOP
      end if

      if (cbuffi(3) .eq. 1) then
         write (*,*) '(ice) received terminate signal'
         stop
      endif
      endif

c.. distribute data to subdomains

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
      call global_scatter(rbuff(1,1,kbufr_flw)   ,flw  )  ! W/m^2
      call global_scatter(rbuff(1,1,kbufr_rain)  ,rain )  ! kg/m^2 s
      call global_scatter(rbuff(1,1,kbufr_snow)  ,snow )  ! kg/m^2 s liquid
        !  ocean states
      call global_scatter(rbuff(1,1,kbufr_sst)   ,sst  )  ! K
      call global_scatter(rbuff(1,1,kbufr_sss)   ,sss  )  ! o/oo
      call global_scatter(rbuff(1,1,kbufr_gwatx) ,gwatx)  ! m/s
      call global_scatter(rbuff(1,1,kbufr_gwaty) ,gwaty)  ! m/s
      call global_scatter(rbuff(1,1,kbufr_tiltx) ,tiltx)  ! m/m
      call global_scatter(rbuff(1,1,kbufr_tilty) ,tilty)  ! m/m
        !  ocean -> ice flux
      call global_scatter(rbuff(1,1,kbufr_fw)    ,frzmlt)  ! W/m^2

c.. rotate zonal/meridional vectors to local coordinates
c.. compute data derived quantities

      cht = 4.19e4            ! heat capacity of water (J kg-1 deg-1)
      tscale = dt             ! time scale for fw, s
      tscale = max(tscale,dt)

c$OMP PARALLEL DO PRIVATE(i,j)
      do j=jlo,jhi
       do i=ilo,ihi
        sst  (i,j) = sst(i,j) - Tf0kel      ! sea sfc temperature, deg C
        Tf   (i,j) = -mu_Tf*sss(i,j)        ! freezing temp, deg C
c.. determine the maximum ocean layer thickness used to compute fw
        if (hm(i,j).gt.tiny.and.frzmlt(i,j).ne.0.) then
         chtd(i,j) = abs(frzmlt(i,j)*dt/(Tf(i,j) - sst(i,j))) ! depth * cht
        else
         chtd(i,j) = 0.
        endif
       enddo
      enddo
      cdepth = ice_global_real_maxval((jhi-jlo+1)*(ihi-ilo+1),chtd)
      if (my_task.eq.master_task) write (6,*) 'depth = ',cdepth/cht,'m'

c$OMP PARALLEL DO PRIVATE(i,j,workx,worky)
      do j=jlo,jhi
       do i=ilo,ihi
c..  atmosphere
        workx      = uatm(i,j)               ! wind velocity, m/s
        worky      = vatm(i,j) 
        uatm (i,j) = workx*cos(ANGLE(i,j))   ! convert to POP grid
     1             + worky*sin(ANGLE(i,j))
        vatm (i,j) = worky*cos(ANGLE(i,j))
     1             - workx*sin(ANGLE(i,j))

        wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s
        fsw  (i,j) = swvdr(i,j) + swvdf(i,j)
     1             + swidr(i,j) + swidf(i,j)
        snow (i,j) = snow(i,j)/rhos ! snowfall rate, m/s 
                ! rain drains directly into the ocean with meltwater

c..  ocean
        workx      = gwatx  (i,j)            ! currents, m/s 
        worky      = gwaty  (i,j)
        gwatx(i,j) = workx*cos(ANGLE(i,j))   ! convert to POP grid 
     1             + worky*sin(ANGLE(i,j))
        gwaty(i,j) = worky*cos(ANGLE(i,j))
     1             - workx*sin(ANGLE(i,j))

        workx      = tiltx  (i,j)            ! sea sfc tilt, m/m
        worky      = tilty  (i,j)      
        tiltx(i,j) = workx*cos(ANGLE(i,j))   ! convert to POP grid 
     1             + worky*sin(ANGLE(i,j))
        tilty(i,j) = worky*cos(ANGLE(i,j))
     1             - workx*sin(ANGLE(i,j))

        fw  (i,j) = -(sst(i,j)-Tf(i,j))*cdepth/tscale 
        fw  (i,j) = min(fw(i,j),0.)         ! oceanic heat flux, W/m^2
        if (fw(i,j).lt.frzmlt(i,j)) fw(i,j) = frzmlt(i,j)
       enddo
      enddo

c.. interpolate ocean dynamics variables from T-cell centers to U-cell 
c.. centers 
c.. NOTE:  atmo variables are interpolated in subroutine roughness as
c.. necessary
      call t2ugrid(gwatx)
      call t2ugrid(gwaty)
      call t2ugrid(tiltx)
      call t2ugrid(tilty)

      time_forc=time

      return
      end
c     ------------------------------------------------------      
      subroutine to_coupler
c.. ice -> flux coupler data
      implicit none
      include "ice.h"          ! ice code common blocks
      include "mpif.h"         ! MPI library definitions
      include "coupling.h"     ! coupling common blocks
      include "thermw.h"       ! thermodynamics common blocks

      integer i, j                      ! generic indices
      integer n                         ! category index
      integer ier                       ! error flag
      real    workx, worky              ! tmps for converting grid

      real Tsrf (ilo:ihi,jlo:jhi)       ! surface temperature
      real alb1 (ilo:ihi,jlo:jhi)       ! albedo, visible, direct
      real alb2 (ilo:ihi,jlo:jhi)       ! albedo, near IR, direct
      real alb3 (ilo:ihi,jlo:jhi)       ! albedo, visible, diffuse
      real alb4 (ilo:ihi,jlo:jhi)       ! albedo, near IR, diffuse
      real lat  (ilo:ihi,jlo:jhi)       ! latent heat flux
      real sens (ilo:ihi,jlo:jhi)       ! sensible heat flux
      real lwup (ilo:ihi,jlo:jhi)       ! outgoing longwave radiation
      real evap (ilo:ihi,jlo:jhi)       ! evaporated water 
      real meltw(ilo:ihi,jlo:jhi)       ! water flux -> ocean
      real melth(ilo:ihi,jlo:jhi)       ! net heat flux -> ocean
      real netsw(ilo:ihi,jlo:jhi)       ! penetrating shortwave -> ocean
      real salt (ilo:ihi,jlo:jhi)       ! salt flux -> ocean
      real tauxa(ilo:ihi,jlo:jhi)       ! atmo/ice stress
      real tauya(ilo:ihi,jlo:jhi)               
      real tauxo(ilo:ihi,jlo:jhi)       ! ice/ocean stress
      real tauyo(ilo:ihi,jlo:jhi)               
      real work (ilo:ihi,jlo:jhi)               

c.. merge across categories and scale fluxes by final ice fraction

      call merge_state(Tsfc,   Tsrf)
      call merge_state(albedo, alb1)
c      call merge_state(albd2, alb2)    ! for sophisticated albedo paramet.
c      call merge_state(albd3, alb3)
c      call merge_state(albd4, alb4)
      call merge_flux(fhnet,    ifrc,melth)
      call merge_flux(fresh,    ifrc,meltw)
      call merge_flux(fsalt,    ifrc,salt)
      call merge_flux(flatent,  ifrc,lat)
      call merge_flux(fsensible,ifrc,sens)
      call merge_flux(flwout,   ifrc,lwup)
      call merge_flux(fswthru,  ifrc,netsw)

c$OMP PARALLEL DO PRIVATE(i,j,workx,worky)
      do j=jlo,jhi
       do i=ilo,ihi
c.. surface temperature
        Tsrf(i,j)  = Tf0kel + Tsrf(i,j)                    !K
c.. albedo                     ! simple version for now
        alb2(i,j) = alb1(i,j)
        alb3(i,j) = alb1(i,j)
        alb4(i,j) = alb2(i,j)
c.. evaporated water -> atmosphere
        evap(i,j) = lat(i,j)/(Lfresh + Lvap)               ! kg/m^2 s
c.. water -> ocean, due to melt/freezing at top and bottom surfaces
c..                 and rainfall on ice surface
        meltw(i,j) = evap(i,j) + meltw(i,j) + rain(i,j)    ! kg/m^2 s
c.. wind stress  (on POP grid:  convert to lat-lon)
        workx = strairx(i,j)                               ! N/m^2
        worky = strairy(i,j)                               ! N/m^2
        tauxa(i,j) = workx*cos(ANGLE(i,j)) - worky*sin(ANGLE(i,j))
        tauya(i,j) = worky*cos(ANGLE(i,j)) + workx*sin(ANGLE(i,j))
c.. ice/ocean stress (on POP grid)
        workx = -strocnx(i,j)                               ! N/m^2
        worky = -strocny(i,j)                               ! N/m^2
        tauxo(i,j) = workx*cos(ANGLE(i,j)) - worky*sin(ANGLE(i,j)) ! convert
        tauyo(i,j) = worky*cos(ANGLE(i,j)) + workx*sin(ANGLE(i,j)) ! to lat-lon
       enddo
      enddo

c.. interpolate from U-cell centers to T-cell centers
      call u2tgrid(tauxa)
      call u2tgrid(tauya)
      call u2tgrid(tauxo)
      call u2tgrid(tauyo)

c.. gather coupling variables from subdomains
                                                        ! departure units 
      call global_gather(sbuff(1,1,kbufs_ifrc), ifrc )  ! none
      call global_gather(sbuff(1,1,kbufs_tsfc), Tsrf )  ! K
      call global_gather(sbuff(1,1,kbufs_alb1), alb1 )  ! none
      call global_gather(sbuff(1,1,kbufs_alb2), alb2 )  ! none
      call global_gather(sbuff(1,1,kbufs_alb3), alb3 )  ! none
      call global_gather(sbuff(1,1,kbufs_alb4), alb4 )  ! none
      call global_gather(sbuff(1,1,kbufs_lat),  lat  )  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_sens), sens )  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_lwup), lwup )  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_evap), evap )  ! kg/m^2 s
      call global_gather(sbuff(1,1,kbufs_tauxa),tauxa)  ! N/m^2
      call global_gather(sbuff(1,1,kbufs_tauya),tauya)  ! N/m^2

      call global_gather(sbuff(1,1,kbufs_netsw),netsw)  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_melth),melth)  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_meltw),meltw)  ! kg/m^2 s
      call global_gather(sbuff(1,1,kbufs_salt), salt )  ! kg/m^2 s
      call global_gather(sbuff(1,1,kbufs_tauxo),tauxo)  ! N/m^2
      call global_gather(sbuff(1,1,kbufs_tauyo),tauyo)  ! N/m^2

      call global_gather(sbuff(1,1,kbufs_Tref), Tref )  ! K

      cbuffi( 4) = idate          ! date (coded: yyyymmdd)
      cbuffi( 5) = sec            ! elapsed seconds into date
      cbuffi( 7) = imt_local      ! full size in x-dir of ice grid
      cbuffi( 8) = jmt_local      ! full size in y-dir of ice grid
      cbuffi( 9) = nadv_i         ! ice comm pairs per day

c.. send message to coupler
   
      if (my_task == master_task) then
      call MPI_SEND(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                     msgtype_i2c, MPI_COMM_WORLD, ier)

      call MPI_SEND(sbuff, nsbuff, MPI_DOUBLE_PRECISION, cpl_task,
     &                     msgtype_i2c, MPI_COMM_WORLD, ier)

      if (ier .ne. MPI_SUCCESS ) then
         write (6,*) '(ice) ERROR after coupler send ',ier
      endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine merge_flux(work1,ifr,work2)
c.. sums total flux through ice categories
c.. work1 has already been multiplied by its compactness in thermo routine
      implicit none
      include "ice.h"          ! ice code common blocks
      integer i, j, n
      real work1(ilo:ihi,jlo:jhi,ncat),work2(ilo:ihi,jlo:jhi)
      real ifr  (ilo:ihi,jlo:jhi)       ! ice fraction

c$OMP PARALLEL DO PRIVATE(i,j,n)
      do j=jlo,jhi
       do i=ilo,ihi
        work2(i,j) = 0.
        if (ifr(i,j).ne.0.) then
         do n = 1,ncat
          work2(i,j) = work2(i,j) + work1(i,j,n)
         enddo
         work2(i,j) = work2(i,j)/ifr(i,j)
        endif
       enddo
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine exit_coupler
c.. final message from flux coupler 
      implicit none
      include "ice.h"          ! ice code common blocks
      include "mpif.h"         ! MPI library definitions
      include "coupling.h"     ! coupling common blocks

      integer i, j                      ! generic indices
      integer ier                       ! error flag
      integer status(MPI_STATUS_SIZE)

      if (my_task == master_task) then
      call MPI_RECV(cbuffi, ncbuffi, MPI_INTEGER, cpl_task,
     &                       msgtype_c2i, MPI_COMM_WORLD, status, ier)

      call MPI_RECV(rbuff, nrbuff, MPI_DOUBLE_PRECISION, cpl_task,
     &                       msgtype_c2i, MPI_COMM_WORLD, status, ier)

      if (ier.ne.MPI_SUCCESS) then
        write (6,*) '(ice) ERROR after receive ',ier
        STOP
      end if

      if (cbuffi(3) .eq. 1) then
         write (*,*) '(ice) received final coupler msg'
      else
         write (*,*) '(ice) terminating before coupler'
         call MPI_ABORT(MPI_COMM_WORLD,-1,ier)
      endif
      endif

      call MPI_FINALIZE(ierr)

      return
      end
c-----------------------------------------------------------------------

      subroutine mpi_coupled (in_model_name, 
     &                        cpl_task, model_task, model_comm)

!-----------------------------------------------------------------------
!
!     this routine queries all the components of the full coupled
!     system and sets up proper communicators and task ids for each
!     component of the coupled system
!
!     this routine should be called after mpi_init, but before 
!     setting up any internal mpi setups (since these will require
!     the internal communicators returned by this routine)
!
!-----------------------------------------------------------------------

      implicit none
      include "mpif.h"                  ! MPI library definitions

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character (3), intent(in) :: in_model_name   
                        ! 3-letter identifier (atm, lnd, ocn, ice, cpl)
                        ! for the model calling this routine

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer, intent(out) ::
     &  cpl_task           ! master task of coupler
     &, model_task         ! master task of model (in MPI_COMM_WORLD)
     &, model_comm         ! communicator for internal model comms

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

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

!-----------------------------------------------------------------------
!
!     determine processor rank in coupled domain
!
!-----------------------------------------------------------------------

      call MPI_COMM_RANK  (MPI_COMM_WORLD, my_task_coupled, ierr)

!-----------------------------------------------------------------------
!
!     determine which group of processes assigned to each model
!
!-----------------------------------------------------------------------

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
        case default
          write (6,*)'Unknown model ',cmodel,' in coupled model domain'
          write (6,*)'Model must be atm, cpl, ice, lnd, ocn'
          stop
        end select

      end do

!-----------------------------------------------------------------------
!
!     create subroup and communicators for each models internal 
!     communciations, note that MPI_COMM_CREATE must be called by
!     all processes in MPI_COMM_WORLD so this must be done by all
!     models consistently and in the same order.
!
!-----------------------------------------------------------------------

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
 
!-----------------------------------------------------------------------
!
!     determine coupler process and model processes
!     assume the first processor in each domain is the task that will
!     communicate coupled model messages
!
!-----------------------------------------------------------------------

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

      end

!***********************************************************************
