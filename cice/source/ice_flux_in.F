c $Id: $
c=======================================================================
!---! reads and interpolates input forcing data
!---!
!---! author Elizabeth C. Hunke, LANL
c=======================================================================

      module ice_flux_in

      use ice_kinds_mod
      use ice_domain
      use ice_constants
      use ice_flux
      use ice_calendar
      use ice_read_write

      implicit none

      character (char_len) ::   ! input data file names
     & height_file
     &, uwind_file
     &, vwind_file
     &,  potT_file
     &,  tair_file
     &, humid_file
     &,  rhoa_file
     &,   fsw_file
     &,   flw_file
     &,  rain_file
     &,   sst_file
     &,   sss_file

      real (kind=dbl_kind) :: c0intp, c1intp, c2intp, c3intp

c=======================================================================

      contains

c=======================================================================

      subroutine getflux

!---!-------------------------------------------------------------------
!---! gets forcing data and interpolates as necessary
!---!-------------------------------------------------------------------

      call sss_sst_restore
c      if (istep.eq.1) call sss_clim_sst_ic

      if (istep.eq.1.or.mod(time,21600.).eq.0) call NCAR_bulk_dat

      end subroutine getflux

c=======================================================================

      subroutine monthly_dat

!---!-------------------------------------------------------------------
!---! interpolate monthly data to 1x daily
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j, im,ip,ipp

      !-----------------------------------------------------------------
      ! Compute index in table
      !-----------------------------------------------------------------
      im  = mod(month+10,12) + 1
      ip  = mod(month,12) + 1
      ipp = mod(month+1,12) + 1

      call interpolate         ! interpolation coefficients

      !-----------------------------------------------------------------
      ! Read and interpolate
      !-----------------------------------------------------------------
      ! Atmo data
      call intp_data(0,im,ip,ipp,height_file,zlvl)
      call intp_data(0,im,ip,ipp, uwind_file,uatm)
      call intp_data(0,im,ip,ipp, vwind_file,vatm)
      call intp_data(0,im,ip,ipp,  potT_file,potT)
      call intp_data(0,im,ip,ipp,  tair_file,Tair)
      call intp_data(0,im,ip,ipp, humid_file,Qa)
      call intp_data(0,im,ip,ipp,  rhoa_file,rhoa)
      call intp_data(0,im,ip,ipp,   fsw_file,Fsw)
      call intp_data(0,im,ip,ipp,   flw_file,Flw)
      call intp_data(0,im,ip,ipp,  rain_file,Fsnow)
      ! Ocean data
      call intp_data(0,im,ip,ipp,sst_file,sst)
      call intp_data(0,im,ip,ipp,sss_file,sss)

      do j=jlo,jhi
       do i=ilo,ihi
        wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s
       enddo
      enddo

      call complete_getflux_ocn

      end subroutine monthly_dat

c=======================================================================

      subroutine interpolate

!---!-------------------------------------------------------------------
!---!  NOTE this scales time to 12 months of equal length
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j
      real (kind=dbl_kind) ::
     & tt, t0, t1, t2, t3, secmo, secyr, s3, time360

      ! make time cyclic
      secyr = c360 * secday
      time360 = time*c360/c365        ! scale to 360-day year
      tt = mod(time360,secyr)
 
      ! Find neighboring times
      secmo = secyr/c12
      t1 = (real(month)-p5)*secmo        ! midpoint, current month
      t0 = t1 - secmo                     !  - 1 month
      t2 = t1 + secmo                     !  + 1 month
      t3 = t1 + c2*secmo                  !  + 2 months
 
      ! Compute coefficients
      s3 = c2*secmo**3
c     c0intp = - (tt - t1)*(tt - t2)*(tt - t3)/(c3*s3)
      c1intp =   (tt - t0)*(tt - t2)*(tt - t3)/s3
      c2intp = - (tt - t0)*(tt - t1)*(tt - t3)/s3
      c3intp =   (tt - t0)*(tt - t1)*(tt - t2)/(c3*s3)
      c0intp =   c1 - c1intp - c2intp - c3intp

      end subroutine interpolate

c=======================================================================

      subroutine intp_data(recd,imx,ipx,ipp,data_file,field)

      use ice_diagnostics  !! debugging

!---!-------------------------------------------------------------------
!---! performs the interpolation
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j, nbits, recd,imx,ipx,ipp, nrec
      character (72) ::  data_file
      logical (kind=log_kind) :: scatter, diag
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::
     &  fieldm(ilo:ihi,jlo:jhi)
     &, field (ilo:ihi,jlo:jhi)
     &, fieldp(ilo:ihi,jlo:jhi)
     &, fieldpp(ilo:ihi,jlo:jhi)

      nbits = 64
      scatter = .true.
      diag = .false.

      if (istep1.gt.check_step) diag = .true.  !! debugging

      if (my_task.eq.master_task.and.(diag)) write (6,*) '  ',data_file

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------
      call ice_open (49,data_file,nbits)
      nrec = recd+imx
      call ice_read(49,nrec,fieldm,'rda8',scatter,diag)
      nrec = recd+month
      call ice_read(49,nrec,field,'rda8',scatter,diag)
      nrec = recd+ipx
      call ice_read(49,nrec,fieldp,'rda8',scatter,diag)
      nrec = recd+ipp
      call ice_read(49,nrec,fieldpp,'rda8',scatter,diag)
      ! close file
      if (my_task.eq.master_task) close (49)

      !-----------------------------------------------------------------
      ! interpolate
      !-----------------------------------------------------------------
      do j=jlo,jhi
       do i=ilo,ihi
         field(i,j) = c0intp*fieldm(i,j) + c1intp*field(i,j) 
     &              + c2intp*fieldp(i,j) + c3intp*fieldpp(i,j)
       enddo
      enddo

      end subroutine intp_data

c=======================================================================

      subroutine complete_getflux_ocn

!---!-------------------------------------------------------------------
!---! compute remaining ocean forcing fields
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j

      do j=jlo,jhi
       do i=ilo,ihi
         Tf(i,j) = -depressT*sss(i,j)   ! freezing temp, deg C
       enddo
      enddo

      ! interpolate ocean dynamics variables from T-cell centers to 
      ! U-cell centers   
c      call t2ugrid(uocn)
c      call t2ugrid(vocn)
c      call t2ugrid(strtltx)
c      call t2ugrid(strtlty)

      end subroutine complete_getflux_ocn

c=======================================================================

      subroutine sss_sst_restore

      use ice_ocean

!---!-------------------------------------------------------------------
!---! interpolate monthly sss, sst data to timestep
!---! restore sst computed by ice model to data
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j, n, imx,ipx,ipp
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: sstdat
      real (kind=dbl_kind) :: Trestore
      logical (kind=log_kind) :: scatter, diag
      scatter = .true.
      diag = .true.

c      Trestore = dt       ! STRONG restoring 
c      Trestore = 1.296e6  ! 15 days
c      Trestore = 2.592e6  ! 30 days
      Trestore = 7.776e6  ! 90 days

      sss_file = './data_in/sss.mm.100x116.da'
      sst_file = './data_in/sst.mm.100x116.da'

      ! Compute index in table
      imx  = mod(month+10,12) + 1
      ipx  = mod(month,12) + 1
      ipp = mod(month+1,12) + 1

      call interpolate         ! interpolation coefficients

      ! Read and interpolate
      call intp_data(0,imx,ipx,ipp,   sss_file,sss)
      call intp_data(0,imx,ipx,ipp,   sst_file,sstdat)

      ! restore sst to data
      do j = jlo,jhi
       do i = ilo,ihi
        if (sss(i,j).lt.c0) sss(i,j) = c0
        sst(i,j) = sst(i,j) + (sstdat(i,j)-sst(i,j))*dt/Trestore 
#ifdef oceanmixed
        sst_mixedlayer(i,j) = sst(i,j)
#endif
       enddo
      enddo

      call complete_getflux_ocn

      end subroutine sss_sst_restore

c=======================================================================

      subroutine sss_clim_sst_ic

!---!-------------------------------------------------------------------
!---! annual mean climatology for Levitus sss
!---! sst initial condition
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j, nbits, k
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: ssstmp
      logical (kind=log_kind) :: scatter, diag
      scatter = .true.
      diag = .true.

      nbits = 64

      sss_file = './data_in/sss.mm.100x116.da'
      sst_file = './data_in/sst.mm.100x116.da'

      !-----------------------------------------------------------------
      ! surface salinity
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) call ice_open(49,sss_file,nbits)

      do j = jlo,jhi
       do i = ilo,ihi
        ssstmp(i,j) = c0
       enddo
      enddo

      do k = 1,12

      call ice_read(49,k,ssstmp,'rda8',scatter,diag)
      do j = jlo,jhi
       do i = ilo,ihi
        sss(i,j) = ssstmp(i,j) + sss(i,j)
       enddo
      enddo

      enddo  ! k

      do j = jlo,jhi
       do i = ilo,ihi
        sss(i,j) = sss(i,j)/c12       ! annual average salinity
        if (sss(i,j).lt.c0) sss(i,j) = c0
       enddo
      enddo

      !-----------------------------------------------------------------
      ! Shea, Trenberth and Reynolds SST
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) call ice_open(50,sst_file,nbits)
      call ice_read(50,month,sst,'rda8',scatter,diag)

      ! close files
      if (my_task.eq.master_task) close (49)
      if (my_task.eq.master_task) close (50)

      do j=jlo,jhi
       do i=ilo,ihi
         ! Initialize freezing temperature
         Tf     (i,j) = -depressT*sss(i,j)   ! freezing temp, deg C
         ! Make sure sst is not less than Tf
         sst(i,j) = max(sst(i,j),Tf(i,j))
       enddo
      enddo

      end subroutine sss_clim_sst_ic

c=======================================================================

      subroutine NCAR_bulk_dat

!---!-------------------------------------------------------------------
!---! reads monthly NCAR_bulk atmospheric data
!---!-------------------------------------------------------------------

      integer (kind=int_kind) :: i, j, n, imx,ipx,ipp, nbits, recnum
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) :: cldf
      logical (kind=log_kind) :: scatter, diag
      scatter = .true.
      diag = .false.

c      if (istep1.gt.13000) diag = .true.  !! debugging

       if (my_task.eq.master_task) then
         flw_file   = './data_in/cldf_gx3_85-88.dat'
         rain_file  = './data_in/prec_gx3_85-88.dat'
         fsw_file   = './data_in/swdn_gx3_85-88.dat'
         uwind_file = './data_in/u_10_gx3_85-88.dat'
         vwind_file = './data_in/v_10_gx3_85-88.dat'
         tair_file  = './data_in/t_10_gx3_85-88.dat'
         humid_file = './data_in/q_10_gx3_85-88.dat'
         rhoa_file  = './data_in/dn10_gx3_85-88.dat'

         if (istep.eq.1) then
          write (6,*) flw_file
          write (6,*) rain_file
          write (6,*) fsw_file
          write (6,*) uwind_file
          write (6,*) vwind_file
          write (6,*) tair_file
          write (6,*) humid_file
          write (6,*) rhoa_file
         endif
       endif   ! master_task

10     continue

      nbits = 64

      !-----------------------------------------------------------------
      ! read monthly data 
      !-----------------------------------------------------------------
      ! Compute index in table
      imx = mod(month+10,12) + 1
      ipx = mod(month,12) + 1
      ipp = mod(month+1,12) + 1

      call interpolate         ! interpolation coefficients

      ! Read and interpolate
      call intp_data(0,imx,ipx,ipp,   fsw_file,Fsw)
      call intp_data(0,imx,ipx,ipp,   flw_file,cldf)
      call intp_data(0,imx,ipx,ipp,  rain_file,Fsnow)

      !-----------------------------------------------------------------
      ! read 6-hourly forcing
      !-----------------------------------------------------------------
      ! open files
      call ice_open (34,tair_file, nbits)     !sfc air temperature
      call ice_open (35,uwind_file,nbits)     !wind velocity, x
      call ice_open (36,vwind_file,nbits)     !wind velocity, y
      call ice_open (37,rhoa_file, nbits)     !wind speed
      call ice_open (38,humid_file,nbits)     !humidity
      ! read data
      recnum = 4*int(yday) - 3 + int(sec/21600.)
c      if (my_task.eq.master_task) print*,yday,recnum

      call ice_read (34,recnum,Tair,'rda8',scatter,diag)
      call ice_read (35,recnum,uatm,'rda8',scatter,diag)
      call ice_read (36,recnum,vatm,'rda8',scatter,diag)
      call ice_read (37,recnum,rhoa,'rda8',scatter,diag)
      call ice_read (38,recnum,Qa,  'rda8',scatter,diag)
      ! close files
      if (my_task.eq.master_task) then
       close (34)
       close (35)
       close (36)
       close (37)
       close (38)
      endif

      !-----------------------------------------------------------------
      ! clean up
      !-----------------------------------------------------------------
      do j=jlo,jhi
       do i=ilo,ihi
        ! fix interpolation errors
        Fsw (i,j) = max(Fsw(i,j),c0)
        cldf(i,j) = max(min(cldf(i,j),c1),c0)
        Fsnow(i,j) = max(Fsnow(i,j),c0)   
        rhoa(i,j) = max(rhoa(i,j),c0)
        Qa  (i,j) = max(Qa(i,j),c0)

        ! set other necessary fields
        zlvl(i,j) = c10
        swvdr(i,j) = Fsw(i,j)*(.875)*(.28)   ! as in the dummy atm (latm)
        swvdf(i,j) = Fsw(i,j)*(.875)*(.24)   ! as in the dummy atm (latm)
        swidr(i,j) = Fsw(i,j)*(.875)*(.31)   ! as in the dummy atm (latm)
        swidf(i,j) = Fsw(i,j)*(.875)*(.17)   ! as in the dummy atm (latm)

        Qa(i,j) = Qa(i,j) * 0.94_dbl_kind    ! as in the dummy atm (latm)
        Frain(i,j) = c0
        Fsnow(i,j) = Fsnow(i,j)/2.592e+06_dbl_kind  ! mm/month -> kg/m^2 s
        if (Tair(i,j).ge.Tffresh) then
            Frain(i,j) = Fsnow(i,j)
            Fsnow(i,j) = c0
        endif

        wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s
        potT(i,j) = Tair(i,j)
           ! flw as in Parkinson and Washington (1979)
        Flw(i,j) = stefan_boltzmann*Tair(i,j)**4   !!! downward longwave !!!
     &   * (c1-
     &   0.261_dbl_kind*exp(-7.77e-4_dbl_kind*(Tffresh - Tair(i,j))**2))
     &   * (c1 + 0.275_dbl_kind*cldf(i,j))
       enddo
      enddo

      end subroutine NCAR_bulk_dat

c=======================================================================

      end module ice_flux_in

c=======================================================================
