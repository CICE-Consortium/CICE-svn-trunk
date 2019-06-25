c=======================================================================
c SCCS module : msg_ice.MCL1.0.f   r3.2
c SCCS get    : 04 Jun 1996 / 11:26:54
c SCCS edit   : 04 Jun 1996 / 11:26:49
c SCCS path   : /net/crestone.u1/nieman/csm/drv-3.n/SCCS/s.msg_ice.MCL1.0.f
c=======================================================================
c  edited by E. C. Hunke, 8 April 1998

      SUBROUTINE msg_icei ()

      implicit none

      !--- global ---
      include "driver.h"           ! control flags
      include "msg.MCL1.0.h"       ! msg passing stuff
      include "dims.h"             ! defines domain dimensions
      include "grids.h"            ! domain grids & masks

      !--- local ---
      integer info                 ! msg passing error code
      real, dimension(nx_i,ny_i) :: rmask_i
c-----------------------------------------------------------------------
c PURPOSE:
c   send/recv initial msg to/from ice via NCAR MCL 1.0 msg passing
c INPUT:
c   o flags contained in driver.h
c OUTPUT:
c   o ice info-buffer data         : ibuff_i
c   o ice spatial coordinate arrays: xe_i,xc_i,ye_i,yc_i
c   o ice model domain mask        : mask_i
c MESSAGES:
c   o data  sent : drv info-buffer flags
c   o data  recvd: ice info-buffer flags, grids, domain mask
c ASSUMPTIONS: 
c   o ice domain extends to north pole: ye_i(ny_i+1) = 90.0
c   o ice domain is periodic in longitude: xe_i(1)+360 = xe_i(nx_i+1)
c-----------------------------------------------------------------------
 
      !-----------------------------------------------------------------
      ! send initial msg
      !-----------------------------------------------------------------
      ibuff_d( 1) = ierr      ! error code
      ibuff_d( 2) = stop_now  ! stop now flag
      ibuff_d( 3) = stop_endd ! stop at end of day flag
      ibuff_d( 4) = 1         ! drv steps per day
      ibuff_d( 5) = nstep     ! drv time step
      ibuff_d( 6) = cdate     ! drv date (coded: yyyymmdd)
      ibuff_d( 7) = sec       ! drv elapsed sec into date
      ibuff_d( 8) = hist_tavg ! data flag: time  avg  history data
      ibuff_d( 9) = hist_endd ! data flag: fixed time history data
      ibuff_d(10) = rest_endd ! data flag: restart data
      ibuff_d(15) = info_time ! mp timing flag
      ibuff_d(16) = info_dbug ! dbug diagnostic level

      write(*,9001) '(msg_icei) send d->i msg, id = ',msgid_d2i_init
      call mcl_send_init(                   tid_ice,dd_icec,info)
      call mcl_pack (ibuff_d,nibuff,MCL_INT,tid_ice,dd_null,info)
      call mcl_send (msgid_d2i_init        ,tid_ice,dd_icec,info)
      write(*,9001) '(drvicei) send d->i msg, id = ',msgid_d2i_init,info
      if (info.lt.MCL_SUCCESS) then
         write (*,*) '(msg_icei) ERROR: mcl error code: ',info
         call mcl_perror(info,'(msg_icei) after mcl_send')
         call mcl_finalize(info)
         STOP
      end if

      write(*,*) '(msg_icei) Sent control buffer'

      !-----------------------------------------------------------------
      ! receive initial msg
      !-----------------------------------------------------------------
      write(*,9001) '(msg_icei) recv i->d msg, id = ',msgid_i2d_init
      call mcl_recv(msgid_i2d_init            ,tid_ice,dd_icei,info)
      call mcl_unpk(ibuff_i,nibuff   ,MCL_INT ,tid_ice,dd_null,info)
      call mcl_unpk(  lat_i,nx_i*ny_i,MCL_REAL,tid_ice,dd_icei,info)
      call mcl_unpk(  lon_i,nx_i*ny_i,MCL_REAL,tid_ice,dd_icei,info)
      call mcl_unpk(rmask_i,nx_i*ny_i,MCL_REAL,tid_ice,dd_icei,info)
      if (info.lt.MCL_SUCCESS) then
         write (*,*) '(msg_icei) ERROR: mcl error code: ',info
         call mcl_perror(info,'(msg_icei) after mcl_recv')
         call mcl_finalize(info)
         STOP
      end if

      mask_i = rmask_i

      write(*,*) '(msg_icei) Received grid info'

      !-----------------------------------------------------------------
      ! extract number of ice timesteps per day
      !-----------------------------------------------------------------
      nadv_i = ibuff_i(4)
      write(*,*) '(drv) Ice advanced ',nadv_i,' times/day'

      RETURN
 9001 format(a,i4)
 9010 format(a,' ice domain: [',f7.3,',',f7.3,']x[',f7.3,',',f7.3,']')
      END

c=======================================================================
c  edited by E. C. Hunke, 8 April 1998

      SUBROUTINE msg_icer (Si,Fioi,Faii)

      implicit none

      !--- global variables ---
      include "driver.h"              ! control flags
      include "msg.MCL1.0.h"          ! msg passing stuff
      include "dims.h"                ! defines domain dimensions
      include "index.h"               ! decodes buffer contents
      include "grids.h"               ! domain grids & masks

      !--- output ---
      real    Si  (nx_i,ny_i,nSi  )   ! ice state
      real    Fioi(nx_i,ny_i,nFioi)   ! ice/ocn fluxes from ice
      real    Faii(nx_i,ny_i,nFaii)   ! ice/atm fluxes from ice

      !--- local variables ---
      integer i,j                     ! generic loop indicies
      integer info                    ! msg passing return error code

c-----------------------------------------------------------------------
c PURPOSE:
c   receive msg from ice model via NCAR MCL 1.0 msg passing
c INPUT:
c   o none
c OUTPUT:
c   o ibuff_i, Si_i, Fioi_i
c MESSAGES:
c   o data received: ice info-buffer, state variables, output fluxes
c-----------------------------------------------------------------------
      
      if ((info_dbug.gt.2 .and. sec.eq.0) .or. info_dbug.gt.3)
     &  write(*,9001) '(msg_icer) recv i->d msg, id = ',msgid_i2d

      !--- recv ice msg and unpack control buffer ----------
      call mcl_recv(msgid_i2d                    ,tid_ice,dd_icer,info)
      call mcl_unpk(ibuff_i,nibuff      ,MCL_INT ,tid_ice,dd_null,info)
      if (info.lt.MCL_SUCCESS) then
         write (*,*) '(msg_icer) ERROR: mcl error code: ',info
         call mcl_perror(info,'(msg_icer) after mcl_recv')
         call mcl_finalize(info)
         STOP
      end if
      print *,'(drv) Received ice message ',ibuff_i(6),ibuff_i(7)

      !--- unpack data ----------
      !--- ice states ----------
      call mcl_unpk(Si(:,:,kSi_ifrac), nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_icer, info)
      call mcl_unpk(Si(:,:,kSi_u    ), nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_icer, info)
      call mcl_unpk(Si(:,:,kSi_v    ), nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_icer, info)
      call mcl_unpk(Si(:,:,kSi_t    ), nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_icer, info)
      call mcl_unpk(Si(:,:,kSi_avsdr), nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_icer, info)
      call mcl_unpk(Si(:,:,kSi_anidr), nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_icer, info)
      call mcl_unpk(Si(:,:,kSi_avsdf), nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_icer, info)
      call mcl_unpk(Si(:,:,kSi_anidf), nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_icer, info)
      !--- ice/atm output fluxes ---
      call mcl_unpk(Faii(:,:,kFaii_lat ), nx_i*ny_i, MCL_REAL, 
     &                                    tid_ice, dd_icer, info)
      call mcl_unpk(Faii(:,:,kFaii_sen ), nx_i*ny_i, MCL_REAL, 
     &                                    tid_ice, dd_icer, info)
      call mcl_unpk(Faii(:,:,kFaii_lwup), nx_i*ny_i, MCL_REAL, 
     &                                    tid_ice, dd_icer, info)
      call mcl_unpk(Faii(:,:,kFaii_evap), nx_i*ny_i, MCL_REAL, 
     &                                    tid_ice, dd_icer, info)
      !--- ice/ocn output fluxes ---
      call mcl_unpk(Fioi(:,:,kFioi_meltw), nx_i*ny_i, MCL_REAL, 
     &                                     tid_ice, dd_icer, info)
      call mcl_unpk(Fioi(:,:,kFioi_melth), nx_i*ny_i, MCL_REAL, 
     &                                     tid_ice, dd_icer, info)
      call mcl_unpk(Fioi(:,:,kFioi_netsw), nx_i*ny_i, MCL_REAL, 
     &                                     tid_ice, dd_icer, info)
      call mcl_unpk(Fioi(:,:,kFioi_tauxo), nx_i*ny_i, MCL_REAL, 
     &                                     tid_ice, dd_icer, info)
      call mcl_unpk(Fioi(:,:,kFioi_tauyo), nx_i*ny_i, MCL_REAL, 
     &                                     tid_ice, dd_icer, info)

      !--- set out-of-domain ice T to allow a/i flux calcs there ---
c$    doacross local(i,j)
      do j=1,ny_i
      do i=1,nx_i
      if ( mask_i(i,j).eq.0 .or. Si(i,j,kSi_ifrac).le.0.0 )
     &    Si(i,j,kSi_t)=275.0  
      end do
      end do

      RETURN
 9001 format(a,i4)
      END

c=======================================================================
c  edited by E. C. Hunke, 8 April 1998

      SUBROUTINE msg_ices (Sa, So, Fixx)

      implicit none

      !--- global ---
      include "driver.h"              ! control flags
      include "msg.MCL1.0.h"          ! msg passing stuff
      include "dims.h"                ! defines domain dimensions
      include "index.h"               ! decodes buffer contents

      !--- input ---
      real    So  (nx_i,ny_i,nSo  )   ! ocn state variables on ice grid
      real    Sa  (nx_i,ny_i,nSa  )   ! atm state variables on ice grid
      real    Fixx(nx_i,ny_i,nFixx)   ! ice input fluxes

      !--- local ---
      integer i,j                     ! generic loop indicies
      integer info                    ! msg passing error code

      real    wind(nx_i,ny_i)         ! wind speed
c-----------------------------------------------------------------------
c PURPOSE:
c   send msg to ice model via NCAR MCL 1.0 msg passing
c INPUT:
c   o Foxx_i & various flags contained in driver.h 
c OUTPUT:
c   o none
c MESSAGES:
c   o data sent: driver info-buffer, ice input fluxes
c-----------------------------------------------------------------------
      
      if ((info_dbug.gt.2 .and. sec.eq.0) .or. info_dbug.gt.3)
     &  write(*,9001) '(msg_ices) send d->i msg, id = ', msgid_d2i

      !--- fill info-buffer array ---
      ibuff_d( 1) = ierr      ! error code
      ibuff_d( 2) = stop_now  ! stop now flag
      ibuff_d( 3) = stop_endd ! stop at end of day flag
      ibuff_d( 4) = 1         ! drv steps per day
      ibuff_d( 5) = nstep     ! drv time step
      ibuff_d( 6) = cdate     ! drv date (coded: yyyymmdd)
      ibuff_d( 7) = sec       ! drv elapsed sec into date
      ibuff_d( 8) = hist_tavg ! data flag: time  avg  history data
      ibuff_d( 9) = hist_endd ! data flag: fixed time history data
      ibuff_d(10) = rest_endd ! data flag: restart data
      ibuff_d(15) = info_time ! mp timing flag
      ibuff_d(16) = info_dbug ! dbug diagnostic level

      !--- initialize msg ---
      call mcl_send_init(                          tid_ice,dd_ices,info)
      call mcl_pack (ibuff_d,nibuff      ,MCL_INT ,tid_ice,dd_null,info)

      !--- pack msg buffer ---
      !--- ocn states ---
      call mcl_pack (So(:,:,kSo_t), nx_i*ny_i, MCL_REAL, 
     &                              tid_ice, dd_ices, info)
      call mcl_pack (So(:,:,kSo_s), nx_i*ny_i, MCL_REAL, 
     &                              tid_ice, dd_ices, info)
      call mcl_pack (So(:,:,kSo_u), nx_i*ny_i, MCL_REAL, 
     &                              tid_ice, dd_ices, info)
      call mcl_pack (So(:,:,kSo_v), nx_i*ny_i, MCL_REAL, 
     &                              tid_ice, dd_ices, info)
      !--- atm states ---
c$    doacross local(i,j)
      do j=1,ny_i
        do i=1,nx_i
          wind(i,j) = sqrt(Sa(i,j,kSa_u)**2 + Sa(i,j,kSa_v)**2)
        end do
      end do
      call mcl_pack (wind,             nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_ices, info)
      call mcl_pack (Sa(:,:,kSa_ptem), nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_ices, info)
      call mcl_pack (Sa(:,:,kSa_shum), nx_i*ny_i, MCL_REAL, 
     &                                 tid_ice, dd_ices, info)
      !--- ocn fluxes ---
      call mcl_pack (Fixx(:,:,kFixx_q   ), nx_i*ny_i, MCL_REAL, 
     &                                     tid_ice, dd_ices, info)
      call mcl_pack (Fixx(:,:,kFixx_dhdx), nx_i*ny_i, MCL_REAL, 
     &                                     tid_ice, dd_ices, info)
      call mcl_pack (Fixx(:,:,kFixx_dhdy), nx_i*ny_i, MCL_REAL, 
     &                                     tid_ice, dd_ices, info)
      !--- atm fluxes ---
      call mcl_pack (Fixx(:,:,kFixx_netsw), nx_i*ny_i, MCL_REAL, 
     &                                      tid_ice, dd_ices, info)
      call mcl_pack (Fixx(:,:,kFixx_lwd  ), nx_i*ny_i, MCL_REAL, 
     &                                      tid_ice, dd_ices, info)
      call mcl_pack (Fixx(:,:,kFixx_water), nx_i*ny_i, MCL_REAL, 
     &                                      tid_ice, dd_ices, info)
      call mcl_pack (Fixx(:,:,kFixx_taux ), nx_i*ny_i, MCL_REAL, 
     &                                      tid_ice, dd_ices, info)
      call mcl_pack (Fixx(:,:,kFixx_tauy ), nx_i*ny_i, MCL_REAL, 
     &                                      tid_ice, dd_ices, info)

      !--- send msg ---
      call mcl_send (     msgid_d2i              ,tid_ice,dd_ices,info)
      if (info.lt.MCL_SUCCESS) then
         write (*,*) '(msg_ices) ERROR: mcl error code: ',info
         call mcl_perror(info,'(msg_ices) after mcl_send')
         call mcl_finalize(info)
         STOP
      end if
      
      RETURN
 9001 format(a,i4)
      END

c=======================================================================
c edited by E. C. Hunke 9 April 1998

      SUBROUTINE msg_icef ( )

      implicit none

      !--- global ---
      include "driver.h"              ! control flags
      include "msg.MCL1.0.h"          ! msg passing stuff
      include "dims.h"                ! defines domain dimensions
      include "index.h"               ! decodes buffer contents

      !--- local ---
      integer i,j,k                   ! generic loop indicies
      integer info                    ! msg passing error code
      real    buff(nx_i,ny_i,11)      ! outgoing msg buffer

c-----------------------------------------------------------------------
c PURPOSE:
c   send final msg to ice model via NCAR MCL 1.0 msg passing
c INPUT:
c   o various flags contained in driver.h 
c OUTPUT:
c   o none
c MESSAGES:
c   o data sent: driver info-buffer
c-----------------------------------------------------------------------
      
      if ((info_dbug.gt.2 .and. sec.eq.0) .or. info_dbug.gt.3)
     &  write(*,9001) '(msg_icef) send d->i msg, id = ', msgid_d2i

      !--- fill info-buffer array ---
      ibuff_d( 1) = ierr      ! error code
      ibuff_d( 2) = stop_now  ! stop now flag
      ibuff_d( 3) = stop_endd ! stop at end of day flag
      ibuff_d( 4) = 1         ! drv steps per day
      ibuff_d( 5) = nstep     ! drv time step
      ibuff_d( 6) = cdate     ! drv date (coded: yyyymmdd)
      ibuff_d( 7) = sec       ! drv elapsed sec into date
      ibuff_d( 8) = hist_tavg ! data flag: time  avg  history data
      ibuff_d( 9) = hist_endd ! data flag: fixed time history data
      ibuff_d(10) = rest_endd ! data flag: restart data
      ibuff_d(15) = info_time ! mp timing flag
      ibuff_d(16) = info_dbug ! dbug diagnostic level

      !--- msg buffer contains no useful data ---
      do k=1,11
c$    doacross local(i,j)
      do j=1,ny_i
      do i=1,nx_i
         buff(i,j,k) = 1.0e36
      end do
      end do
      end do

      !--- send msg ---
      call mcl_send_init(                          tid_ice,dd_ices,info)
      call mcl_pack (ibuff_d,nibuff      ,MCL_INT ,tid_ice,dd_null,info)
      call mcl_pack ( buff  ,nx_i*ny_i*11,MCL_REAL,tid_ice,dd_ices,info)
      call mcl_send (     msgid_d2i               ,tid_ice,dd_ices,info)
      if (info.lt.MCL_SUCCESS) then
         write (*,*) '(msg_icef) ERROR: mcl error code: ',info
         call mcl_perror(info,'(msg_icef) after mcl_send')
      end if
      
      RETURN
 9001 format(a,i4)
      END

c=======================================================================
c  created by E. C. Hunke, 10 April 1998

      SUBROUTINE flux_ices (Sa, So, Fixx)

      implicit none

      !--- global ---
      include "driver.h"              ! control flags
      include "dims.h"                ! defines domain dimensions
      include "index.h"               ! decodes buffer contents
      include "grids.h"               ! ice grid lat/lon and mask

      !--- input ---
      real    So  (nx_i,ny_i,nSo  )   ! ocn state variables on ice grid
      real    Sa  (nx_i,ny_i,nSa  )   ! atm state variables on ice grid
      real    Fixx(nx_i,ny_i,nFixx)   ! ice input fluxes

      !--- local ---
      integer i,j                     ! generic loop indicies

c-----------------------------------------------------------------------
c PURPOSE:
c   create bogus ocn/atm fluxes to send to ice model 
c INPUT:
c   o latitude, longitude
c OUTPUT:
c   o fluxes to be sent to the ice model
c Sign convention:
c    positive value <=> downward flux
c Unit convention:
c    T             ~ Kelvin
c    u,v           ~ m/s
c    heat flux     ~ W/m^2 
c    momentum flux ~ N/m^2
c    water flux    ~ (kg/s)/m^2
c-----------------------------------------------------------------------
      
c$    doacross local(i,j)
      do j=1,ny_i
        do i=1,nx_i
      !--- ocn states ---
          So(i,j,kSo_s) = 34.8         ! Greenland Sea salinity 
          So(i,j,kSo_t) = (25.*cos(9.*lat_i(i,j)/7.) ! lat-dependent
     &                  - 0.0544*34.8) ! equals freezing temperature at 70N,S
     &                  * mask_i(i,j)  ! land mask
     &                  + 273.15       ! Kelvin
          So(i,j,kSo_u) = 0.1          ! ocean currents
     &                  * mask_i(i,j)  ! land mask
          So(i,j,kSo_v) = 0.               ! this is not necessarily zonal
      !--- atm states ---
          Sa(i,j,kSa_u) = 1.           ! wind
          Sa(i,j,kSa_v) = 0.
          Sa(i,j,kSa_ptem) = -10+35.*cos(9.*lat_i(i,j)/7.) ! lat-dependent
     &                  + 273.         ! Kelvin
          Sa(i,j,kSa_shum) = 0.        ! nice and dry
      !--- ocn fluxes ---
          Fixx(i,j,kFixx_q   ) = 0.    ! no heat flux 
          Fixx(i,j,kFixx_dhdx) = 0.    ! no swell
          Fixx(i,j,kFixx_dhdy) = 0.
      !--- atm fluxes ---
          Fixx(i,j,kFixx_netsw) = 0.   ! dark
          Fixx(i,j,kFixx_lwd  ) = 0.
          Fixx(i,j,kFixx_water) = 0.   ! and dry
          Fixx(i,j,kFixx_taux ) = 1.56e-3  ! zonal wind stress
          Fixx(i,j,kFixx_tauy ) = 0.
        end do
      end do

      RETURN
      END
