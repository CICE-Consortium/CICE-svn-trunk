c $Id: $
c=======================================================================
!---! diagnostic information output during run
!---!
!---! authors  (?)
!---!         Elizabeth C. Hunke, LANL
c=======================================================================

      module ice_diagnostics

      use ice_domain
      use ice_constants
      use ice_calendar

      implicit none

      ! point print data
      logical (kind=log_kind) ::
     &    print_points             ! if true, print point data
      integer (kind=int_kind), parameter ::
     &          npnt = 2  ! total number of points to be printed

      real (kind=dbl_kind), dimension(npnt) :: 
     &     latpnt,  !  latitude of desired points
     $     lonpnt   ! longitude of desired points
      integer (kind=int_kind) ::
     &      iindx(npnt)  ! i index for points
     &,     jindx(npnt)  ! j index for points

      ! for water and heat budgets
      real (kind=dbl_kind), dimension(npnt), save ::
     &   pFs_tot  ! total snowfall (m)
     &,  hefi   ! ice thickness, initial (m)
     &,  hefs   ! snow thickness, initial (m)
     &,  tote   ! ice and snow energy, initial
     &,  plat, plon           ! latitude, longitude of points
     &,  piloc, pjloc, pmloc  ! location of points on processors

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::
     &   dhefi  ! thermo diagnostic

      data latpnt /   90.,  -70. /
      data lonpnt /    0.,   40. /

      ! printing info for routine print_state
      character (len=20) :: plabel
      integer (kind=int_kind), parameter ::
     &  check_step = 87700000
     &, ip = 25
     &, jp = 111
     &, mtask = 1

c=======================================================================

      contains

c=======================================================================

      subroutine runtime_diags

!---!-------------------------------------------------------------------
!---! Writes diagnostic info (max, min, global sums, etc) to standard out
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_flux
      use ice_albedo
      use ice_mpi_internal
      use ice_history
      use ice_grid
      use ice_state
      use ice_dyn_evp
      use ice_coupling
      use ice_itd

      real (kind=dbl_kind) ::
     &   umaxn, hmaxn, shmaxn, arean, snwmxn, extentn, atotn
     &,  umaxs, hmaxs, shmaxs, areas, snwmxs, extents, atots
     &,  workn(imt_local,jmt_local), works(imt_local,jmt_local)
     &,  ptmp, etot

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::  
     &   Tsrf          ! surface temperature
     &,  ailohi        ! ice area
     &,  work1, work   ! temporary arrays

      real (kind=dbl_kind), dimension(imt_global,jmt_global) ::
     &   Fsw_g   ! shortwave
     &,  Tf_g    ! freezing temperature
     &,  vice_g  ! ice volume
     &,  vsno_g  ! snow volume
     &,  meltb_g
     &,  meltt_g
     &,  snoice_g
     &,  frazil_g
     &,  congel_g
     &,  eice_g
     &,  dhefi_g

      real (kind=dbl_kind), dimension(2) ::
     &   par, pTair, pQa, pFs, pFr, pFsw, pFlw, pTsf, psub, pFsb
     &,  pdhi, pdhs, pde, psst, pTf, hiavg, hsavg, pFw
     &,  pFlb, pmeltb, pmeltt, psnoi, pfraz, pcong, pdhit

      integer (kind=int_kind) :: n, layer, i,j, ii,jj

      !-----------------------------------------------------------------
      ! state of the ice
      !-----------------------------------------------------------------
      ! hemispheric quantities
      ! maximum effective thickness
      do j=1,jmt_local
       do i=1,imt_local
       workn(i,j) = mask_n(i,j)*vice(i,j)
       works(i,j) = mask_s(i,j)*vice(i,j)
       enddo
      enddo
      hmaxn = ice_global_real_maxval(imt_local*jmt_local,WORKN)
      hmaxs = ice_global_real_maxval(imt_local*jmt_local,WORKS)
      ! total ice area
      call get_sum(0,tarean,one,aice,arean)
      call get_sum(0,tareas,one,aice,areas)
      arean = arean * m2_to_km2
      areas = areas * m2_to_km2
      ! ice extent
      do j=1,jmt_local
       do i=1,imt_local
       workn(i,j) = c0
       if (aice(i,j).gt.puny) workn(i,j) = c1
       enddo
      enddo
      call get_sum(0,tarean,one,workn,extentn)
      call get_sum(0,tareas,one,workn,extents)
      extentn = extentn * m2_to_km2
      extents = extents * m2_to_km2
      ! total ice volume
      call get_sum(0,tarean,one,vice,shmaxn)
      call get_sum(0,tareas,one,vice,shmaxs)
      ! total snow volume
      call get_sum(0,tarean,one,vsno,snwmxn)
      call get_sum(0,tareas,one,vsno,snwmxs)
      ! average ice albedo
       do j=jlo,jhi
        do i=ilo,ihi
         work1(i,j) = alvdr(i,j)*awtvdr + alidr(i,j)*awtidr 
     $              + alvdf(i,j)*awtvdf + alidf(i,j)*awtidf
        enddo
       enddo
      call get_sum(1,tarean,work1,aice,atotn)
      call get_sum(1,tareas,work1,aice,atots)
      atotn = atotn / arean * m2_to_km2  ! convert arean back to m2
      atots = atots / areas * m2_to_km2
      ! maximum velocity  
      do j=1,jmt_local
       do i=1,imt_local
       workn(i,j) = max(abs(mask_n(i,j)*u(i,j)),abs(mask_n(i,j)*v(i,j)))
       works(i,j) = max(abs(mask_s(i,j)*u(i,j)),abs(mask_s(i,j)*v(i,j)))
       enddo
      enddo
      umaxn = ice_global_real_maxval(imt_local*jmt_local,WORKN)
      umaxs = ice_global_real_maxval(imt_local*jmt_local,WORKS)

      if (print_points) then
      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------
      do j=jlo,jhi
       do i=ilo,ihi
        ! ice/snow enthalpy change
        work(i,j) = esno(i,j)
        do layer = 1,ntilay
         work(i,j) = work(i,j) + eicen(i,j,layer)
        enddo
       enddo
      enddo
      call global_gather(eice_g,work)

      call global_gather(Fsw_g,Fsw)
      call global_gather(Tf_g,Tf)
      call global_gather(vice_g,vice(ilo:ihi,jlo:jhi))
      call global_gather(vsno_g,vsno(ilo:ihi,jlo:jhi))
      call global_gather(meltb_g,meltb)
      call global_gather(meltt_g,meltt)
      call global_gather(snoice_g,snoice)
      call global_gather(congel_g,congel)
      call global_gather(frazil_g,frazil)
      call global_gather(dhefi_g,dhefi)

#ifndef coupled
      call global_gather(rbuff(1,1,kbufr_Tair)  ,Tair )  ! K
      call global_gather(rbuff(1,1,kbufr_Qa)    ,Qa   )  ! kg/kg
      call global_gather(rbuff(1,1,kbufr_rain)  ,Frain)  ! kg/m^2 s
      call global_gather(rbuff(1,1,kbufr_snow)  ,Fsnow)  ! kg/m^2 s liquid
      call global_gather(rbuff(1,1,kbufr_flw)   ,Flw  )  ! W/m^2
      call global_gather(rbuff(1,1,kbufr_sst)   ,sst  )  ! K
      do j=jlo,jhi
        do i=ilo,ihi
          ailohi(i,j) = aice(i,j)
          Tsrf(i,j) = Tsfc(i,j)
        enddo
      enddo
      call global_gather(sbuff(1,1,kbufs_ifrc), ailohi   )  ! none
      call global_gather(sbuff(1,1,kbufs_tsfc), Tsrf     )  ! K
      call global_gather(sbuff(1,1,kbufs_evap), evap )      ! kg/m^2 s
      call global_gather(sbuff(1,1,kbufs_sens), Fsensible)  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_lat ), Flatent  )  ! W/m^2
      call global_gather(sbuff(1,1,kbufs_melth),Foht_hist)  ! W/m^2
      do j=1,jmt_global
        do i=1,imt_global
          rbuff(i,j,kbufr_sst) = rbuff(i,j,kbufr_sst) + Tffresh
          sbuff(i,j,kbufs_Tsfc) = sbuff(i,j,kbufs_Tsfc) + Tffresh
        enddo
      enddo
#endif

      if (my_task.eq.master_task) then
      do n = 1,npnt
        ii = iindx(n)
        jj = jindx(n)

        ! air temperature
        pTair(n) = rbuff(ii,jj,kbufr_Tair) - Tffresh
        ! specific humidity
        pQa(n) = rbuff(ii,jj,kbufr_Qa)
        ! snowfall
        pFs(n) = rbuff(ii,jj,kbufr_snow)*dt/rhos
        pFs_tot(n) = pFs_tot(n) + pFs(n)*diagfreq ! approx total to date
        ! rainfall
        pFr(n) = rbuff(ii,jj,kbufr_rain)*dt/rhow
        ! shortwave radiation
        pFsw(n) = Fsw_g(ii,jj)
        ! longwave radiation
        pFlw(n) = rbuff(ii,jj,kbufr_flw)
        ! ice area
        par(n) = sbuff(ii,jj,kbufs_ifrc)
        ! average snow/ice thicknesses
        hiavg(n) = c0
        hsavg(n) = c0
        if (par(n).ne.c0) then 
          hiavg(n) = vice_g(ii,jj)/par(n)
          hsavg(n) = vsno_g(ii,jj)/par(n)
        endif
        ! ice/snow surface temperature
        pTsf(n) = sbuff(ii,jj,kbufs_Tsfc) - Tffresh
        ! sublimation/condensation
        psub(n) = sbuff(ii,jj,kbufs_evap)*dt/rhoi
        ! latent heat flux
        pFlb(n) = sbuff(ii,jj,kbufs_lat)
        ! sensible heat flux
        pFsb(n) = sbuff(ii,jj,kbufs_sens)
        ! top melt
        pmeltt(n) = meltt_g(ii,jj)
        ! bottom melt
        pmeltb(n) = meltb_g(ii,jj)
        ! snoice
        psnoi(n) = snoice_g(ii,jj)
        ! frazil ice
        pfraz(n) = frazil_g(ii,jj)
        ! congelation ice
        pcong(n) = congel_g(ii,jj)
        ! ice thickness change
        pdhi(n) = vice_g(ii,jj)-hefi(n)
        ! misc thermo diagnostic
        pdhit(n) = dhefi_g(ii,jj)
        ! snow thickness change
        pdhs(n) = (vsno_g(ii,jj)-hefs(n))
        ! ice/snow enthalpy change
        pde(n) = -(eice_g(ii,jj)-tote(n))/dt
        ! sea surface temperature
        psst(n) = rbuff(ii,jj,kbufr_sst) - Tffresh
        ! freezing temperature
        pTf(n) = Tf_g(ii,jj)
        ! ocean heat used by ice
        pFw(n) = -sbuff(ii,jj,kbufs_melth)
      enddo ! n
      endif ! master_task
      endif ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------
      if (my_task.eq.master_task) then

        write (6,899) 'Arctic','Antarctic'

        write (6,900) 'max cH (m)             = ',hmaxn,hmaxs
        write (6,901) 'total area (km^2)      = ',arean,areas
        write (6,901) 'total ice volume (m^3) = ',shmaxn,shmaxs
        write (6,901) 'total snw volume (m^3) = ',snwmxn,snwmxs
        write (6,900) 'average albedo         = ',atotn,atots
        write (6,900) 'max u, v (m/s)         = ',umaxn,umaxs

        if (print_points) then  

          write (6,*) '                         '
        write (6,902) '       Lat, Long         ',plat(1),plon(1),
     &                                            plat(2),plon(2)
        write (6,903) '     my_task, i, j       ',
     &                                      pmloc(1),piloc(1),pjloc(1),
     &                                      pmloc(2),piloc(2),pjloc(2)
          write (6,*) '----------atm----------'
        write (6,900) 'air temperature (C)    = ',pTair(1),pTair(2)
        write (6,900) 'specific humidity      = ',pQa(1),pQa(2)
        write (6,900) 'snowfall (m)           = ',pFs(1),pFs(2)
        write (6,900) 'total snow to date (m) = ',pFs_tot(1),pFs_tot(2)
        write (6,900) 'rainfall (m)           = ',pFr(1),pFr(2)
        write (6,900) 'shortwave radiation sum= ',pFsw(1),pFsw(2)
        write (6,900) 'longwave radiation     = ',pFlw(1),pFlw(2)
          write (6,*) '----------ice----------'
        write (6,900) 'area fraction          = ',par(1),par(2)
        write (6,900) 'avg ice thickness (m)  = ',hiavg(1),hiavg(2)
        write (6,900) 'avg snow depth (m)     = ',hsavg(1),hsavg(2)
        write (6,900) 'surface temperature(C) = ',pTsf(1),pTsf(2)
        write (6,900) 'sensible heat flx      = ',pFsb(1),pFsb(2)
        write (6,900) 'latent heat flx        = ',pFlb(1),pFlb(2)
        write (6,900) 'subl/cond (m ice)      = ',psub(1),psub(2)
        write (6,900) 'top melt (m)           = ',pmeltt(1),pmeltt(2)
        write (6,900) 'bottom melt (m)        = ',pmeltb(1),pmeltb(2)
        write (6,900) 'new ice (m)            = ',pfraz(1),pfraz(2)
        write (6,900) 'congelation (m)        = ',pcong(1),pcong(2)
        write (6,900) 'snow-ice (m)           = ',psnoi(1),psnoi(2)
        write (6,900) 'effective dhi (m)      = ',pdhi(1),pdhi(2)
        write (6,900) 'effective dhs (m)      = ',pdhs(1),pdhs(2)
        write (6,900) 'intnl enrgy chng(W/m^2)= ',pde(1),pde(2)
        write (6,900) 'thermo diag (m)        = ',pdhi(1),pdhi(2)
          write (6,*) '----------ocn----------'
        write (6,900) 'sst (C)                = ',psst(1),psst(2)
        write (6,900) 'freezing temp (C)      = ',pTf(1),pTf(2)
        write (6,900) 'heat used (W/m^2)      = ',pFw(1),pFw(2)

        endif
      endif

 899  format (27x,a24,2x,a24)
 900  format (a25,2x,f24.17,2x,f24.17)
 901  format (a25,2x,1pe24.17,2x,1pe24.17)
 902  format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
 903  format (a25,10x,f4.0,1x,f4.0,1x,f4.0,9x,f4.0,1x,f4.0,1x,f4.0)

      end subroutine runtime_diags

c=======================================================================

      subroutine init_mass_diags

!---!-------------------------------------------------------------------
!---! computes global combined ice and snow mass sum
!---!-------------------------------------------------------------------

      use ice_mpi_internal
      use ice_grid
      use ice_state

      integer (kind=int_kind) :: n, layer, ii,jj, i,j
      real (kind=dbl_kind) ::
     &   vice_g(imt_global,jmt_global)
     &,  vsno_g(imt_global,jmt_global)
     &,  work_g(imt_global,jmt_global)
     &,  work(ilo:ihi,jlo:jhi)

      if (print_points) then
      call global_gather(vice_g,vice(ilo:ihi,jlo:jhi))
      call global_gather(vsno_g,vsno(ilo:ihi,jlo:jhi))
      do n = 1,npnt
       if (my_task.eq.master_task) then
        ii = iindx(n)
        jj = jindx(n)

        ! save north/south ice effective thickness (m)
        hefi(n) = vice_g(ii,jj)
        ! save north/south snow effective thickness (m)
        hefs(n) = vsno_g(ii,jj)
       endif
        ! save north/south effective energy (J/m^2)
        do j=jlo,jhi
         do i=ilo,ihi
          dhefi(i,j) = c0       ! convenient place to initialize dhi
          work(i,j) = esno(i,j)
          do layer = 1,ntilay
           work(i,j) = work(i,j) + eicen(i,j,layer)
          enddo
         enddo
        enddo
       call global_gather(work_g,work)
       if (my_task.eq.master_task) tote(n) = work_g(ii,jj)        
       if (istep.eq.1) pFs_tot(n) = c0
      enddo
      endif

      end subroutine init_mass_diags

c=======================================================================

      subroutine init_diags

!---!-------------------------------------------------------------------
!---! initialize diagnostic output
!---!-------------------------------------------------------------------

      use ice_grid
      use ice_mpi_internal

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::  
     &   ispace,jspace  ! work arrays for finding points on processors
      real (kind=dbl_kind), dimension(imt_global,jmt_global) ::
     &   hm_g        ! ocean mask
     &,  ispace_g,jspace_g ! work arrays for finding points on processors
      real (kind=dbl_kind) ::
     &   minbth(1)   ! sum of distances in lat and lon
     &,  mindis      ! min(minbth) over process tasks
     &,  latdis      ! latitude distance
     &,  londis      ! longitude distance
     &,  totdis      ! total distance
      integer (kind=int_kind) ::
     &   n   ! index for search
     &,  i,j ! grid indices

      character (char_len) :: label(npnt)
      logical (kind=log_kind) :: prnt     ! for prints
      data prnt / .true. /

      if (print_points) then

      ! initialize labels
      if (my_task.eq.master_task) 
     &   write(6,*) ' initialize diagnostic prints '

      label(1)(1:40)  = 'Near North Pole pack ice                '
      label(2)(1:40)  = 'Weddell Sea                             '

      ! gather mask on global grid
      call global_gather(hm_g,hm(ilo:ihi,jlo:jhi))

      if (my_task.eq.master_task) then

      write(6,*) ' find indices of points '
      do n=1,npnt
        ! compute indices
        minbth(n) = 540.0_dbl_kind    !  360.0 + 180.
        iindx(n)   = 0
        jindx(n)   = 0
        do j=1,jmt_global
          do i=1,imt_global
            if (hm_g(i,j).gt.p5) then
             latdis = abs(latpnt(n)-TLAT_G(i,j)*rad_to_deg)
             londis = abs(lonpnt(n)-TLON_G(i,j)*rad_to_deg)
             totdis = sqrt(latdis**2 + londis**2)
             if( totdis .lt. minbth(n) ) then
              minbth(n) = totdis
              jindx(n)   = j
              iindx(n)   = i
             endif
            endif
          enddo
        enddo

        write(6,2234) n,latpnt(n),lonpnt(n),
     &    TLAT_G(iindx(n),jindx(n))*rad_to_deg,
     &    TLON_G(iindx(n),jindx(n))*rad_to_deg,
     &    iindx(n),jindx(n)
 2234   format(' found point',i4/
     &  '   lat     lon     TLAT    TLON    iindx jindx '/
     &  4(f7.1,1x),4x,2(i3,3x))

      enddo ! n
      endif ! master_task

      ! find points on local domain processors
      do j=jlo,jhi
       do i=ilo,ihi
        ispace(i,j) = real(i + my_task*10000)
        jspace(i,j) = real(j + my_task*10000)
       enddo
      enddo
      call global_gather(ispace_g,ispace)
      call global_gather(jspace_g,jspace)
      if (my_task.eq.master_task) then
       do n=1,npnt
        ! Lat, Long
        plat(n) = TLAT_G(iindx(n),jindx(n))*rad_to_deg
        plon(n) = TLON_G(iindx(n),jindx(n))*rad_to_deg
        ! location on processor
        piloc(n) = mod(ispace_g(iindx(n),jindx(n)),10000.)       ! i
        pjloc(n) = mod(jspace_g(iindx(n),jindx(n)),10000.)       ! j
        pmloc(n) = (ispace_g(iindx(n),jindx(n))-piloc(n))/10000 ! my_task
       enddo
      endif ! master_task

      endif ! print_points

      end subroutine init_diags

c=======================================================================

      subroutine print_state(plabel,i,j)

!---!-------------------------------------------------------------------
!---! this routine is useful for debugging
!---! call to it should be inserted in the form (after thermo, for example)
!---!      do j=jlo,jhi
!---!       do i=ilo,ihi
!---!         call to_column(i,j)
!---!         plabel = 'post thermo'
!---!         if (istep1.ge.check_step.and.i.eq.ip.and.j.eq.jp
!---!     &       .and.my_task.eq.mtask) 
!---!     &   call print_state(plabel,i,j)
!---!       enddo
!---!      enddo
!---! 'use ice_diagnostics' may need to be inserted also
!---!
!---! author Elizabeth C. Hunke, LANL
!---!-------------------------------------------------------------------

      use ice_model_size
      use ice_kinds_mod
      use ice_state
      use ice_itd
      use ice_flux

      character (len=20), intent(in) :: plabel
      integer (kind=int_kind), intent(in) :: i,j
      real (kind=dbl_kind) :: edebug, qs,Tsnow
      integer (kind=int_kind) :: nc, layer

      write (6,*) plabel
      write (6,*) 'my_task',my_task,' i, j=',i,j,' istep1 = ',istep1

      write (6,*) 'ai0', ai0
      do nc=1,ncat
        write (6,*) ' '
        write (6,*) 'nc=',nc
        write (6,*) 'ain', ain(nc)
        write (6,*) 'vin',vin(nc)
        write (6,*) 'vsn',vsn(nc)
        if (ain(nc).gt.puny) write (6,*) 'hin',vin(nc)/ain(nc)
        if (ain(nc).gt.puny) write (6,*) 'hsn',vsn(nc)/ain(nc)
        write (6,*) 'Tsfn',Tsfn(nc)

        write (6,*) 'esn',esn(nc)
        if (vsn(nc).gt.puny) then
          qs = -esn(nc)/vsn(nc)
          Tsnow = (Lfresh - qs/rhos) / cp_ice
          write (6,*) 'qs',qs
          write (6,*) 'Tsnow',Tsnow
          write (6,*) ' '
        endif
      enddo

c      print*,'hsnow avg',
c     & (vsn(1)+vsn(2)+vsn(3)+vsn(4)+vsn(5))/
c     & (ain(1)+ain(2)+ain(3)+ain(4)+ain(5))

      edebug = c0
      do nc = 1,ncat
        do layer = 1,nlyr
          write (6,*) 'ein, cat ',nc,' layer ',layer,ein(layer,nc)
          edebug = edebug+ein(layer,nc)
        enddo
      enddo
      write (6,*) 'eice(i,j)',edebug 
      write (6,*) ' '
c      write (6,*) eicen(i-1,j+1,1),eicen(i,j+1,1),eicen(i+1,j+1,1)
c      write (6,*) eicen(i-1,j  ,1),eicen(i,j  ,1),eicen(i+1,j  ,1)
c      write (6,*) eicen(i-1,j-1,1),eicen(i,j-1,1),eicen(i+1,j-1,1)
c      write (6,*) ' '


      write (6,*) 'Forcing'
      write (6,*) 'uatm',uatm(i,j)
      write (6,*) 'vatm',vatm(i,j)
      write (6,*) 'potT',potT(i,j)
      write (6,*) 'Tair',Tair(i,j)
      write (6,*) 'Qa',Qa(i,j)
      write (6,*) 'rhoa',rhoa(i,j)
      write (6,*) 'swvdr',swvdr(i,j)
      write (6,*) 'swvdf',swvdf(i,j)
      write (6,*) 'swidr',swidr(i,j)
      write (6,*) 'swidf',swidf(i,j)
      write (6,*) 'Flw',Flw(i,j)
      write (6,*) 'Frain',Frain(i,j)
      write (6,*) 'Fsnow',Fsnow(i,j)
      write (6,*) 'sst',sst(i,j)
      write (6,*) 'sss',sss(i,j)
      write (6,*) 'uocn',uocn(i,j)
      write (6,*) 'vocn',vocn(i,j)
      write (6,*) 'strtltx',strtltx(i,j)
      write (6,*) 'strtlty',strtlty(i,j)
      write (6,*) 'frzmlt',frzmlt(i,j)

      end subroutine print_state

c=======================================================================

      end module ice_diagnostics

c=======================================================================
