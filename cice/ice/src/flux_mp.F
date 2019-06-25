c $Id: flux_mp.F,v 1.3.2.2 1998/05/04 16:43:24 eclare Exp $

c     program flux_mp.F
c   message passing to/from the flux coupler
c-----------------------------------------------------------------------
      subroutine init_mp
c.. initializes message passing
      implicit none
#include "ice.h"
#include "mcl.h"                         ! MCL specific definitions
#include "ice_mp.h"                    

      integer ier                        ! error flag for message calls
      integer geom                       ! unused argument to mcl_alloc() 
      integer nbytes                     ! number of bytes for allocation calls
      integer maxsec                     ! max time to wait for message 
      parameter (maxsec=300)
      integer model_num                  ! model number of ice in coupler
      integer msgtype_d2ii               ! message type for drv->ice init
      integer msgtype_i2di               ! message type for ice->drv init
      character*3 group_name 

      integer i, j                       ! generic indices
      real    rmask(imt,jmt)             ! temporary boundary mask

      nsnd = 17                          ! number of fields sent to coupler
      nrcv = 15                          ! number of fields sent from coupler
      nadv_i = int(86400./dt)            ! number of timesteps per day

      group_name='LCM'

      model_num=8                        ! model number of ice in coupler
      msgtype_d2ii = 31                  ! message type for drv->ice init
      msgtype_i2di = 40                  ! message type for ice->drv init
      msgtype_d2i  = 30                  ! message type for drv->ice
      msgtype_i2d  = 41                  ! message type for ice->drv

c.. connect this model to the flux coupler ...............................

      call mcl_con2coupler(model_num,group_name,maxsec,msg_id,ier)
      if ( ier .lt. MCL_SUCCESS ) then 
         call mcl_perror(ier,'(ice) ERROR after mcl_con2coupler')
         call mcl_finalize(ier)
         stop
      endif
      write(*,*) '(ice) Connected to coupler '

c.. allocate message buffer space ........................................
c   Note: MCL_SIZEOF_INT and MCL_SIZEOF_REAL definitions provided in mcl.h 

c   control message contains nibuff int flags
c      nbytes = nibuff*MCL_SIZEOF_INT
c      call mcl_alloc(bufidc, nbytes, msg_id, geom, ier)
      bufidc(1) = 0 ! control message not buffered

c   initial message contains grid coordinates and mask
        nbytes = 3*imt*jmt*MCL_SIZEOF_REAL
        call mcl_alloc(bufidi, nbytes, msg_id, geom, ier)

c   send nsnd fields to coupler
	nbytes=imt*jmt*nsnd*MCL_SIZEOF_REAL
        call mcl_alloc(bufids, nbytes, msg_id, geom, ier)

c   receive nrcv fields from coupler
	nbytes=imt*jmt*nrcv*MCL_SIZEOF_REAL
        call mcl_alloc(bufidr, nbytes, msg_id, geom, ier)

      if ( ier .lt. MCL_SUCCESS ) then 
         call mcl_perror(ier,'(ice) ERROR allocating MCL buffers')
         call mcl_finalize(ier)
         stop
      endif
      write(*,*) '(ice) Allocated buffer space'

c.. exchange initial messages ............................................

c   receive coupler message and ignore content
      call mcl_recv(msgtype_d2ii,msg_id,bufidc,ier)
      if (ier .lt. MCL_SUCCESS ) then
        call mcl_perror(ier,'(ice) ERROR after initial recv')
        call mcl_finalize(ier)
        stop
      endif

      call mcl_unpk(ibuff_d,nibuff,MCL_INT,msg_id,bufidc,ier)
      if (ier .lt. MCL_SUCCESS ) then
        call mcl_perror(ier,'(ice) ERROR after initial unpk')
        call mcl_finalize(ier)
        stop
      endif
      write(*,*) '(ice) Received message from coupler '

c.. fill buffer and send initial message to coupler

      idate = year*10000 + 101 ! initial date (yyyymmdd) - January 1
      sec = 0                  ! elapsed seconds into date

      ibuff( 1) = 0        ! error code
      ibuff( 2) = 0        ! stop now flag
      ibuff( 3) = 0        ! stop at end of day flag
      ibuff( 4) = nadv_i   ! ice steps per day
      ibuff( 5) = 1        ! initial time step
      ibuff( 6) = idate    ! initial date (coded: yyyymmdd)
      ibuff( 7) = sec      ! elapsed seconds into date
      ibuff( 8) = 0        ! data flag: time  avg  history data
      ibuff( 9) = 0        ! data flag: fixed time history data
      ibuff(10) = 0        ! data flag: restart data
      ibuff(15) = 0        ! mp timing flag
      ibuff(16) = 0        ! dbug diagnostic level

c.. send grid to coupler
c$DOACROSS LOCAL(i,j)
      do j=1,jmt
       do i=1,imt
        rmask(i,j)  = hm(i,j)        ! temporary mask array: REAL(imt,jmt)
       enddo
      enddo

      call mcl_send_init(msg_id,bufidi,ier) 

      print*,' call mcl_pack(ibuff,nibuff,MCL_INT,msg_id,bufidi,ier)'
      call mcl_pack(ibuff,nibuff,MCL_INT,msg_id,bufidc,ier)
      call mcl_pack(ULAT, imt*jmt,MCL_REAL,msg_id,bufidi,ier)
      call mcl_pack(ULONG,imt*jmt,MCL_REAL,msg_id,bufidi,ier)
      call mcl_pack(rmask,imt*jmt,MCL_REAL,msg_id,bufidi,ier)

         if (ier .lt. MCL_SUCCESS ) then
           call mcl_perror(ier,'(ice) ERROR after initial pack')
           call mcl_finalize(ier)
           stop
         endif

      print*,'call mcl_send(msgtype_i2di,msg_id,bufidi,ier)'
      call mcl_send(msgtype_i2di,msg_id,bufidi,ier)

      if (ier .lt. MCL_SUCCESS ) then
        call mcl_perror(ier,'(ice) ERROR after inital send')
        call mcl_finalize(ier)
        stop
      endif

      write(*,*) '(ice) Sent grid info'

      return
      end
c-----------------------------------------------------------------------
      subroutine from_coupler
c.. flux coupler -> ice data
      implicit none
#include "ice.h"
#include "mcl.h"                        ! MCL specific definitions
#include "ice_mp.h"                    

      integer i, j
      integer ier                       ! error flag
      real    workx, worky

      call mcl_recv(msgtype_d2i,msg_id,bufidr,ier)
      if (ier .lt. MCL_SUCCESS ) then
        call mcl_perror(ier,'(ice) ERROR after recv')
        call mcl_finalize(ier)
        stop
      endif

      call mcl_unpk(ibuff_d, nibuff, MCL_INT,  msg_id, bufidc, ier)
      call mcl_unpk(sst,    imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(sss,    imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(gwatx,  imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(gwaty,  imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(wind,   imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(tair,   imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(Qa,     imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(fw,     imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(tiltx,  imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(tilty,  imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(fsw,    imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(flw,    imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(snow,   imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(strairx,imt*jmt, MCL_REAL, msg_id, bufidr, ier)
      call mcl_unpk(strairy,imt*jmt, MCL_REAL, msg_id, bufidr, ier)

      if (ier.ne.MCL_SUCCESS) then
        call mcl_perror(ier,'(ice) ERROR after unpack')
        call mcl_finalize(ier)
        STOP
      end if

      if (ibuff_d(2) .eq. 1) then
         write (*,*) '(ice) received terminate signal'
         call mcl_finalize(ier)
         stop
      endif

c$DOACROSS LOCAL(i,j,workx,worky)
c....... Convert units MKS -> CGS and
c....... rotate zonal/meridional vectors to local coordinates
      do j=1,jmt
       do i=1,imt
c..  state variables from ocean
        sst    (i,j) = sst    (i,j) - 273.15   ! sea sfc temperature, deg C
        sss    (i,j) = sss    (i,j)            ! salinity, o/oo         
        workx        = gwatx  (i,j) * 100.     ! currents, cm/s 
        worky        = gwaty  (i,j) * 100.     
        gwatx  (i,j) = workx*cos(ANGLE(i,j))   ! convert to POP grid 
     1               + worky*sin(ANGLE(i,j))
        gwaty  (i,j) = worky*cos(ANGLE(i,j))
     1               - workx*sin(ANGLE(i,j))
c..  state variables from atmosphere
        wind   (i,j) = wind   (i,j) * 100.     ! wind speed, cm/s
        tair   (i,j) = tair   (i,j)            ! surface air temperature, K
        Qa     (i,j) = Qa     (i,j)            ! specific humidity (10m)
c..  fluxes from ocean
        fw     (i,j) = fw     (i,j) * 2.388e-5 ! heat flux, cal/cm^2 s
        workx        = tiltx  (i,j) * 0.1      ! sea sfc tilt, g/cm s^2 
        worky        = tilty  (i,j) * 0.1      
        tiltx  (i,j) = workx*cos(ANGLE(i,j))   ! convert to POP grid 
     1               + worky*sin(ANGLE(i,j))
        tilty  (i,j) = worky*cos(ANGLE(i,j))
     1               - workx*sin(ANGLE(i,j))
c..  fluxes from atmosphere
        fsw    (i,j) = fsw    (i,j) * 2.388e-5 ! incoming shortwave, cal/cm^2 s
        flw    (i,j) = flw    (i,j) * 2.388e-5 ! incoming longwave, cal/cm^2 s
        snow   (i,j) = snow   (i,j) /(10.*rhos)! snowfall rate, cm/s 

        workx        = strairx(i,j) * 10.0     ! wind stress, g/cm s^2
        worky        = strairy(i,j) * 10.0    
        strairx(i,j) = workx*cos(ANGLE(i,j))   ! convert to POP grid
     1               + worky*sin(ANGLE(i,j))
        strairy(i,j) = worky*cos(ANGLE(i,j))
     1               - workx*sin(ANGLE(i,j))

c....... Initialize freezing temperature
        tf     (i,j) = -0.0544*sss(i,j)        ! freezing temp, deg C
c....... Set maximum ocean heat flux to -2 W/m^2
        fw     (i,j) = min(fw(i,j),-4.776e-5)  ! heat flux, cal/cm^2 s
       enddo
      enddo

c.. albedo of the snow
      data salbedo /0.,0.,.83,.81,.82,.78,.64,.69,.84,.85,0.,0./         
      albsnow = salbedo(month)        ! Semtner's
c**      albsnow = 0.                 ! if already included in shortwave

      return
      end
c     ------------------------------------------------------      
      subroutine to_coupler
c.. ice -> flux coupler data
      implicit none
#include "ice.h"
#include "mcl.h"                        ! MCL specific definitions
#include "ice_mp.h"                    

      integer i, j                      ! generic indices
      real ttmerge, a0, a1, b           ! merging function
      integer ier                       ! error flag
      real    tsfc                      ! melt/freeze flag
      real    sfrac                     ! temporary snow fraction
      real    aicevs                    ! ice albedo, visible
      real    aiceir                    ! ice albedo, near IR
      real    asnovs                    ! snow albedo, visible
      real    asnoir                    ! snow albedo, near IR
      real    avsdr1                    ! surface albedo, visible, thick
      real    anidr1                    ! surface albedo, near IR, thick
      real    avsdr0                    ! surface albedo, visible, thin
      real    anidr0                    ! surface albedo, near IR, thin
      real    workx, worky              ! tmps for converting grid

      real ifrc (imt,jmt)               ! ice fraction
      real iceu (imt,jmt)               ! ice velocity, zonal
      real icev (imt,jmt)               ! ice velocity, meridional
      real tsrf (imt,jmt)               ! surface temperature
      real alb1 (imt,jmt)               ! albedo, visible, direct
      real alb2 (imt,jmt)               ! albedo, near IR, direct
      real alb3 (imt,jmt)               ! albedo, visible, diffuse
      real alb4 (imt,jmt)               ! albedo, near IR, diffuse
      real lat  (imt,jmt)               ! latent heat flux
      real sens (imt,jmt)               ! sensible heat flux
      real lwup (imt,jmt)               ! outgoing longwave radiation
      real evap (imt,jmt)               ! evaporated water 
      real meltw(imt,jmt)               ! water flux -> ocean
      real melth(imt,jmt)               ! net heat flux -> ocean
      real netsw(imt,jmt)               ! penetrating shortwave -> ocean
      real tauxo(imt,jmt)               ! ice/ocean stress
      real tauyo(imt,jmt)               

      ttmerge(a0,a1,b) = b*a1 + (1.-b)*a0 ! merge thick/thin characteristics

c$DOACROSS LOCAL(i,j,workx,worky)
      do j=1,jmt
       do i=1,imt
c.. ice fraction
        ifrc(i,j) = compact(i,j)
        if (hithin(i,j).gt.0.) ifrc(i,j) = 1.
c.. ice velocity (on POP grid)
        workx = u(i,j)*0.01                                ! m/s
        worky = v(i,j)*0.01                                ! m/s
        iceu(i,j) = workx*cos(ANGLE(i,j)) - worky*sin(ANGLE(i,j)) ! convert
        icev(i,j) = worky*cos(ANGLE(i,j)) + workx*sin(ANGLE(i,j)) ! to lat-lon
c.. surface temperature
        tsrf(i,j)  = 273.15 +
     1   ttmerge(tsfcij(i,j,0),tsfcij(i,j,1),compact(i,j))     ! K
c.. albedo                     ! NCAR CSM 
cc..   thick ice
c        sfrac = hsthick(i,j)/(hsthick(i,j) + 1.)
c        tsfc = 0.                                   ! melting
c        if (tsfcij(i,j,1).lt.0) tsfc = 1.           ! freezing
c        aicevs = tsfc*0.7  + (1.-tsfc)*0.5          ! ice, visible
c        aiceir = 0.5                                ! ice, near IR
c        asnovs = tsfc*0.95 + (1.-tsfc)*0.85         ! snow, visible
c        asnoir = tsfc*0.7  + (1.-tsfc)*0.55         ! snow, near IR
cc      direct
c        avsdr1 = ttmerge(aicevs,asnovs,sfrac)
c        anidr1 = ttmerge(aiceir,asnoir,sfrac)
cc..   thin ice
c        sfrac = hsthin(i,j)/(hsthin(i,j) + 1.)
c        tsfc = 0.                                   ! melting
c        if (tsfcij(i,j,0).lt.0) tsfc = 1.           ! freezing
c        aicevs = tsfc*0.7  + (1.-tsfc)*0.5          ! ice, visible
c        aiceir = 0.5                                ! ice, near IR
c        asnovs = tsfc*0.95 + (1.-tsfc)*0.85         ! snow, visible
c        asnoir = tsfc*0.7  + (1.-tsfc)*0.55         ! snow, near IR
cc      direct
c        avsdr0 = ttmerge(aicevs,asnovs,sfrac)
c        anidr0 = ttmerge(aiceir,asnoir,sfrac)
cc      effective albedo
c        alb1(i,j) = ttmerge(avsdr0,avsdr1,compact(i,j))
c        alb2(i,j) = ttmerge(anidr0,anidr1,compact(i,j))
c.. albedo                     ! simple version for now
        alb1(i,j) = 
     1   ttmerge(albij(i,j,0),albij(i,j,1),compact(i,j))  
        alb2(i,j) = alb1(i,j)
        alb3(i,j) = alb1(i,j)
        alb4(i,j) = alb2(i,j)
c.. latent heat flux
        lat(i,j) = 4.1876e4*
     1   ttmerge(flatent(i,j,0),flatent(i,j,1),compact(i,j))     ! W/m^2   
c.. sensible heat flux
        sens(i,j) = 4.1876e4*
     1   ttmerge(fsensible(i,j,0),fsensible(i,j,1),compact(i,j)) ! W/m^2
c.. outgoing longwave radiation
        lwup(i,j) = 4.1876e4*
     1   ttmerge(flwout(i,j,0),flwout(i,j,1),compact(i,j))       ! W/m^2
c.. evaporated water -> atmosphere
c   2.834e6 = 2.5e6 + 0.334e6 = latent heats of vaporization + fusion 
        evap(i,j) = lat(i,j)/2.834e6                             ! kg/m^2 s
c.. water -> ocean, due to melt/freezing at top and bottom surfaces
        meltw(i,j) = evap(i,j) + 10.*
     1   ttmerge(fresh(i,j,0),fresh(i,j,1),compact(i,j))         ! kg/m^2 s
c.. net heat flux -> ocean
        melth(i,j) = 4.1876e4*
     1   ttmerge(fhnet(i,j,0),fhnet(i,j,1),compact(i,j))         ! W/m^2
c.. penetrating shortwave radiation (through ice into ocean)
        netsw(i,j) = 0.                                          ! W/m^2 
c.. ice/ocean stress (on POP grid)
        workx = -strocnx(i,j)*0.1                                ! m/s
        worky = -strocny(i,j)*0.1                                ! m/s
        tauxo(i,j) = workx*cos(ANGLE(i,j)) - worky*sin(ANGLE(i,j)) ! convert
        tauyo(i,j) = worky*cos(ANGLE(i,j)) + workx*sin(ANGLE(i,j)) ! to lat-lon
       enddo
      enddo

      ibuff( 4) = nadv_i   ! ice steps per day
      ibuff( 5) = istep    ! time step
      ibuff( 6) = idate    ! date (coded: yyyymmdd)
      ibuff( 7) = sec      ! elapsed seconds into date
   
      call mcl_send_init(msg_id,bufids,ier) 

      call mcl_pack(ibuff, nibuff,MCL_INT ,msg_id,bufidc,ier)
      call mcl_pack( ifrc,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack( iceu,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack( icev,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack( tsrf,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack( alb1,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack( alb2,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack( alb3,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack( alb4,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack(  lat,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack( sens,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack( lwup,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack( evap,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack(meltw,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack(melth,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack(netsw,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack(tauxo,imt*jmt,MCL_REAL,msg_id,bufids,ier)
      call mcl_pack(tauyo,imt*jmt,MCL_REAL,msg_id,bufids,ier)

      if (ier .lt. MCL_SUCCESS ) then
         call mcl_perror(ier,'(ice) ERROR after pack')
         call mcl_finalize(ier)
         stop
      endif

      call mcl_send(msgtype_i2d,msg_id,bufids,ier)

      if (ier .lt. MCL_SUCCESS ) then
         call mcl_perror(ier,'(ice) ERROR after send')
         call mcl_finalize(ier)
         stop
      endif

      return
      end
c-----------------------------------------------------------------------
