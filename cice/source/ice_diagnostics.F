! $Id: ice_diagnostics.F,v 1.16 2004/03/01 16:56:23 eclare Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_diagnostics - diagnostic information output during run
!
! !DESCRIPTION:
!
! Diagnostic information output during run
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL 
!          Bruce P. Briegleb, NCAR
!
! !INTERFACE:
!
      module ice_diagnostics
!
! !USES:
!
      use ice_domain
      use ice_constants
      use ice_calendar
      use ice_fileunits
      use ice_work, only:  work_g2, work_l1, work_l2, worka, workb
!
!EOP
!
      implicit none

      ! dignostic output file
      character (len=char_len) :: diag_file

      ! point print data

      logical (kind=log_kind) ::
     &   print_points         ! if true, print point data
     &,  print_global         ! if true, print global data

      integer (kind=int_kind), parameter ::
     &   npnt = 2             ! total number of points to be printed

      real (kind=dbl_kind), dimension(npnt) ::
     &   latpnt               !  latitude of desired points
     &,  lonpnt               ! longitude of desired points

      integer (kind=int_kind) ::
     &   iindx(npnt)          ! i index for points
     &,  jindx(npnt)          ! j index for points

      ! for water and heat budgets
      real (kind=dbl_kind), dimension(npnt), save ::
     &   pFs_tot              ! total snowfall (m)
     &,  hefi                 ! ice thickness, initial (m)
     &,  hefs                 ! snow thickness, initial (m)
     &,  tote                 ! ice and snow energy, initial
     &,  plat, plon           ! latitude, longitude of points
     &,  piloc, pjloc, pmloc  ! location of points on processors

      ! for hemispheric water and heat budgets
      real (kind=dbl_kind) ::
     &   totmn                ! total ice/snow water mass (nh)
     &,  totms                ! total ice/snow water mass (sh)
     &,  totmin               ! total ice water mass (nh)
     &,  totmis               ! total ice water mass (sh)
     &,  toten                ! total ice/snow energy (J)
     &,  totes                ! total ice/snow energy (J)

      data latpnt /   90.,  -65. /
      data lonpnt /    0.,  -45. /

      ! printing info for routine print_state
      character (char_len) :: plabel
      integer (kind=int_kind), parameter ::
     &   check_step = 99999999
     &,  ip = 15
     &,  jp = 14
     &,  mtask = 1

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: runtime_diags - writes max,min,global sums to standard out
!
! !INTERFACE:
!
      subroutine runtime_diags
!
! !DESCRIPTION:
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!          Bruce P. Briegleb, NCAR
!          Cecilia M. Bitz, UW
!
! !USES:
!
      use ice_model_size
      use ice_flux
      use ice_albedo
      use ice_mpi_internal
      use ice_grid
      use ice_state
      use ice_coupling
      use ice_itd
#ifdef coupled
      use shr_sys_mod, only : shr_sys_flush
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind), dimension(:,:), allocatable ::
     &   fsw_g  ,  snoice_g,  rain_g,  sens_g
     &,  Tf_g   ,  frazil_g,  snow_g,  lat_g 
     &,  vice_g ,  congel_g,  flw_g          
     &,  vsno_g ,  eice_g  ,  sst_g          
     &,  meltt_g,  meltb_g,   meltl_g        
     &,  aice_g,   Tair_g  ,  tsfc_g         
     &,  fhnet_g,  Qa_g    ,  evap_g
     &,  swabs_g,  lwout_g

      real (kind=dbl_kind), dimension(2) ::
     &   par, pTair, pQa, pFs, pFr, pFsw, pFlw, pTsf, psub, pFsb
     &,  pdhi, pdhs, pde, psst, pTf, hiavg, hsavg, pFw
     &,  pFlb, pmeltt, pmeltb, pmeltl, psnoi, pfraz, pcong
     &,  pswabs, plwout

      integer (kind=int_kind) :: i, j, k, n, ii,jj

      ! Hemispheric state quantities
      real (kind=dbl_kind) ::
     &   umaxn, hmaxn, shmaxn, arean, snwmxn, extentn, atotn
     &,  umaxs, hmaxs, shmaxs, areas, snwmxs, extents, atots
     &,  etotn, mtotn, micen,  msnwn
     &,  etots, mtots, mices,  msnws

      ! Hemispheric flux quantities
      real (kind=dbl_kind) ::
     &   rnn, snn, frzn,  hnetn, fhocnn, fhatmn,  fhfrzn
     &,  rns, sns, frzs,  hnets, fhocns, fhatms,  fhfrzs
     &,  sfsaltn, sfreshn, evpn, fluxn , delmxn,  delmin
     &,  sfsalts, sfreshs, evps, fluxs , delmxs,  delmis
     &,  delein, werrn, herrn, msltn, delmsltn, serrn
     &,  deleis, werrs, herrs, mslts, delmslts, serrs
     &,  ftmp

      !-----------------------------------------------------------------
      ! state of the ice
      !-----------------------------------------------------------------
      ! hemispheric quantities

      ! maximum effective thickness
      do j=1,jmt_local
      do i=1,imt_local
         work_l1(i,j) = mask_n(i,j)*vice(i,j)
         work_l2(i,j) = mask_s(i,j)*vice(i,j)
      enddo
      enddo
      hmaxn = ice_global_real_maxval(imt_local*jmt_local,work_l1)
      hmaxs = ice_global_real_maxval(imt_local*jmt_local,work_l2)

      ! total ice area
      call get_sum(0,tarean,one,aice,arean)
      call get_sum(0,tareas,one,aice,areas)
      arean = arean * m2_to_km2
      areas = areas * m2_to_km2

      ! ice extent
      do j=1,jmt_local
      do i=1,imt_local
         work_l1(i,j) = c0
         if (aice(i,j) > puny) work_l1(i,j) = c1
      enddo
      enddo
      call get_sum(0,tarean,one,work_l1,extentn)
      call get_sum(0,tareas,one,work_l1,extents)
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
         worka(i,j) = alvdr(i,j)*awtvdr + alidr(i,j)*awtidr
     &              + alvdf(i,j)*awtvdf + alidf(i,j)*awtidf
      enddo
      enddo
      call get_sum(1,tarean,worka,aice,atotn)
      call get_sum(1,tareas,worka,aice,atots)

      if (arean > c0) then 
         atotn = atotn / arean * m2_to_km2  ! convert arean back to m2
      else
         atotn = c0
      endif

      if (areas > c0) then 
         atots = atots / areas * m2_to_km2
      else
         atots = c0
      endif

      ! maximum velocity  
      do j=1,jmt_local
      do i=1,imt_local
         work_l1(i,j) = max(abs(mask_n(i,j)*uvel(i,j)),
     &                      abs(mask_n(i,j)*vvel(i,j)))
         work_l2(i,j) = max(abs(mask_s(i,j)*uvel(i,j)),
     &                      abs(mask_s(i,j)*vvel(i,j)))
      enddo
      enddo
      umaxn = ice_global_real_maxval(imt_local*jmt_local,work_l1)
      umaxs = ice_global_real_maxval(imt_local*jmt_local,work_l2)

      ! total ice/snow internal energy
      do j=1,jmt_local
      do i=1,imt_local
         work_l2(i,j) = esno(i,j) + eice(i,j)
      enddo
      enddo
      call get_sum(0,tarean,one,work_l2,etotn)
      call get_sum(0,tareas,one,work_l2,etots)

      !-----------------------------------------------------------------
      ! various fluxes
      !-----------------------------------------------------------------
      ! evap, fsalt, fresh, fhnet, fswthru, fsens, and flwout 
      ! need to be multiplied by aice because 
      ! regretfully they have been divided by aice for the coupler
      !-----------------------------------------------------------------

      ! evaporation
      call get_sum(1,tarean,evap,aice_init,evpn)	
      call get_sum(1,tareas,evap,aice_init,evps)	
      evpn = evpn*dt
      evps = evps*dt

      ! salt flux
      call get_sum(1,tarean,fsalt_hist,aice_init,sfsaltn) 
      call get_sum(1,tareas,fsalt_hist,aice_init,sfsalts) 
      sfsaltn = sfsaltn*dt
      sfsalts = sfsalts*dt

      ! fresh water flux
      call get_sum(1,tarean,fresh_hist,aice_init,sfreshn) 
      call get_sum(1,tareas,fresh_hist,aice_init,sfreshs) 
      sfreshn = sfreshn*dt
      sfreshs = sfreshs*dt

      ! ocean heat
      ! Note: fswthru not included because it does not heat ice
      call get_sum(1,tarean,fhnet_hist,aice_init,fhocnn)
      call get_sum(1,tareas,fhnet_hist,aice_init,fhocns)

      ! latent heat
      ! You may be wondering, where is the latent heat flux? 
      ! It is not included here because it cancels with 
      ! the evaporative flux times the enthalpy of the  
      ! ice/snow that evaporated.

      ! atmo heat flux
      ! Note: flwout includes the reflected longwave down, needed by the
      !  atmosphere as an upwards radiative boundary condition.
      ! Also note: fswabs includes solar radiation absorbed in ocean,
      !  which must be subtracted here.
      do j=jlo,jhi
      do i=ilo,ihi
         ftmp = (fswabs(i,j) - fswthru_hist(i,j) + flw(i,j) + 
     &           flwout(i,j) + fsens(i,j)) * aice_init(i,j)
         worka(i,j) = tarean(i,j)*ftmp
         workb(i,j) = tareas(i,j)*ftmp
      enddo
      enddo
      fhatmn = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),worka)
      fhatms = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),workb)

      ! freezing potential
      do j=jlo,jhi
      do i=ilo,ihi
         ftmp = max(c0,frzmlt(i,j))
         worka(i,j) = tarean(i,j)*ftmp
         workb(i,j) = tareas(i,j)*ftmp
      enddo
      enddo
      fhfrzn = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),worka)
      fhfrzs = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),workb)

      ! rain
      call get_sum(1,tarean,frain,aice_init,rnn)	
      call get_sum(1,tareas,frain,aice_init,rns)	
      rnn = rnn*dt
      rns = rns*dt

      ! snow
      call get_sum(1,tarean,fsnow,aice_init,snn)	
      call get_sum(1,tareas,fsnow,aice_init,sns)	
      snn = snn*dt
      sns = sns*dt

      ! frazil ice growth !! should not be multiplied by aice
      do j=jlo,jhi
      do i=ilo,ihi
         worka(i,j) = tarean(i,j)*frazil(i,j)*rhoi/dt ! m/step->kg/m^2/s
         workb(i,j) = tareas(i,j)*frazil(i,j)*rhoi/dt
      enddo
      enddo
      frzn = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),worka)
      frzs = ice_global_real_sum((jhi-jlo+1)*(ihi-ilo+1),workb)
      frzn = frzn*dt
      frzs = frzs*dt

      ! ice and snow mass
      micen = rhoi*shmaxn
      msnwn = rhos*snwmxn
      mices = rhoi*shmaxs
      msnws = rhos*snwmxs

      mtotn = micen + msnwn
      mtots = mices + msnws

      ! mass change since beginning of time step
      delmin = mtotn - totmn
      delmis = mtots - totms

      ! ice mass change excluding frazil ice formation
      delmxn = micen - totmin - frzn
      delmxs = mices - totmis - frzs

      ! total water flux
      fluxn  = c0
      fluxs  = c0
      if( arean > c0) then
         fluxn = rnn + snn + evpn - sfreshn + frzn
      endif
      if( areas > c0) then
         fluxs = rns + sns + evps - sfreshs + frzs
      endif
      
      werrn = (fluxn-delmin)/(mtotn+c1)
      werrs = (fluxs-delmis)/(mtots+c1)

      ! energy change
      delein = etotn - toten
      deleis = etots - totes

      fhatmn = fhatmn + ( - snn * Lfresh + evpn * Lvap ) / dt 
      fhatms = fhatms + ( - sns * Lfresh + evps * Lvap ) / dt 

      hnetn = (fhatmn - fhocnn - fhfrzn) * dt
      hnets = (fhatms - fhocns - fhfrzs) * dt

      herrn = (hnetn - delein) / (etotn - c1)
      herrs = (hnets - deleis) / (etots - c1)

      ! salt mass
      msltn = micen*ice_ref_salinity*p001
      mslts = mices*ice_ref_salinity*p001

      ! change in salt mass, exclusing frazil
      delmsltn = delmxn*ice_ref_salinity*p001
      delmslts = delmxs*ice_ref_salinity*p001

      ! salt error
      serrn = (sfsaltn + delmsltn) / (msltn + c1)
      serrs = (sfsalts + delmslts) / (mslts + c1)


      if (print_points) then

         allocate(
     &   fsw_g   (imt_global,jmt_global),
     &   snoice_g(imt_global,jmt_global),
     &   rain_g  (imt_global,jmt_global),
     &   sens_g  (imt_global,jmt_global),
     &   Tf_g    (imt_global,jmt_global),
     &   frazil_g(imt_global,jmt_global),
     &   snow_g  (imt_global,jmt_global),
     &   lat_g   (imt_global,jmt_global),
     &   vice_g  (imt_global,jmt_global),
     &   congel_g(imt_global,jmt_global),
     &   flw_g   (imt_global,jmt_global),
     &   vsno_g  (imt_global,jmt_global),
     &   eice_g  (imt_global,jmt_global),
     &   sst_g   (imt_global,jmt_global),
     &   meltt_g (imt_global,jmt_global),
     &   meltb_g (imt_global,jmt_global),
     &   meltl_g (imt_global,jmt_global),
     &   aice_g  (imt_global,jmt_global),
     &   Tair_g  (imt_global,jmt_global),
     &   tsfc_g  (imt_global,jmt_global),
     &   fhnet_g (imt_global,jmt_global),
     &   Qa_g    (imt_global,jmt_global),
     &   evap_g  (imt_global,jmt_global),
     &   swabs_g (imt_global,jmt_global),
     &   lwout_g (imt_global,jmt_global) )

      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------
         ! gather all information from all processors to master_task
         call global_gather(  Tair_g,Tair )     ! K
         call global_gather(    Qa_g,Qa   )     ! kg/kg
         call global_gather(  snow_g,fsnow)     ! kg/m^2 s liquid
         call global_gather(  rain_g,frain)     ! kg/m^2 s
         call global_gather(   fsw_g,fsw)       ! W/m^2
         call global_gather(   flw_g,flw  )     ! W/m^2
         call global_gather(  aice_g,aice(ilo:ihi,jlo:jhi)) ! none
         call global_gather(  vice_g,vice(ilo:ihi,jlo:jhi)) ! m
         call global_gather(  vsno_g,vsno(ilo:ihi,jlo:jhi)) ! m
         call global_gather(  Tsfc_g,Tsfc(ilo:ihi,jlo:jhi)) ! C
         call global_gather( swabs_g,fswabs)    ! W/m^2
         call global_gather( lwout_g,flwout)    ! W/m^2
         call global_gather(  sens_g,fsens)     ! W/m^2
         call global_gather(   lat_g,flat)      ! W/m^2
         call global_gather(  evap_g,evap )     ! kg/m^2 s
         call global_gather( meltt_g,meltt)     ! m
         call global_gather( meltb_g,meltb)     ! m
         call global_gather( meltl_g,meltl)     ! m
         call global_gather(frazil_g,frazil)    ! m
         call global_gather(congel_g,congel)    ! m
         call global_gather(snoice_g,snoice)    ! m

         do j=jlo,jhi
         do i=ilo,ihi
            worka(i,j) = esno(i,j) + eice(i,j)
         enddo
         enddo

         call global_gather( eice_g,worka)      ! J/m^2
         call global_gather(  sst_g,sst  )      ! C
         call global_gather(   Tf_g,Tf)         ! C
         call global_gather(fhnet_g,fhnet_hist) ! W/m^2

         if (my_task == master_task) then
         ! determine diagnostic for desired grid points from global array
            do n = 1,npnt
               ii = iindx(n)
               jj = jindx(n)

               ! air temperature
               pTair(n) = Tair_g(ii,jj) - Tffresh

               ! specific humidity
               pQa(n) = Qa_g(ii,jj)

               ! snowfall 
               pFs(n) = snow_g(ii,jj)*dt/rhos
               pFs_tot(n) = pFs_tot(n) + pFs(n)*diagfreq ! approx total to date

               ! rainfall
               pFr(n) = rain_g(ii,jj)*dt/rhow

               ! shortwave radiation
               pFsw(n) = fsw_g(ii,jj)

               ! longwave radiation
               pFlw(n) = flw_g(ii,jj)

               ! ice area
               par(n) = aice_g(ii,jj)

               ! average snow/ice thicknesses
               hiavg(n) = c0
               hsavg(n) = c0
               if (par(n) /= c0) then
                  hiavg(n) = vice_g(ii,jj)/par(n)
                  hsavg(n) = vsno_g(ii,jj)/par(n)
               endif

               ! ice/snow surface temperature
               pTsf(n) = Tsfc_g(ii,jj)

               ! sublimation/condensation
               psub(n) = evap_g(ii,jj)*dt/rhoi

               ! absorbed solar flux
               pswabs(n) = swabs_g(ii,jj)

               ! outward longwave flux
               plwout(n) = lwout_g(ii,jj)

               ! latent heat flux
               pFlb(n) = lat_g(ii,jj)

               ! sensible heat flux
               pFsb(n) = sens_g(ii,jj)

               ! top melt
               pmeltt(n) = meltt_g(ii,jj)

               ! bottom melt
               pmeltb(n) = meltb_g(ii,jj)

               ! lateral melt
               pmeltl(n) = meltl_g(ii,jj)

               ! snoice
               psnoi(n) = snoice_g(ii,jj)

               ! frazil ice
               pfraz(n) = frazil_g(ii,jj)

               ! congelation ice
               pcong(n) = congel_g(ii,jj)

               ! ice thickness change
               pdhi(n) = vice_g(ii,jj)-hefi(n)

               ! snow thickness change
               pdhs(n) = (vsno_g(ii,jj)-hefs(n))

               ! ice/snow energy change
               pde(n) = -(eice_g(ii,jj)-tote(n))/dt

               ! sea surface temperature
               psst(n) = sst_g(ii,jj)

               ! freezing temperature
               pTf(n) = Tf_g(ii,jj)

               ! ocean heat used by ice
               pFw(n) = -fhnet_g(ii,jj)

            enddo               ! n
         endif                  ! master_task

         deallocate (
     &   fsw_g   ,  snoice_g,  rain_g,  sens_g
     &,  Tf_g   ,  frazil_g,  snow_g,  lat_g
     &,  vice_g ,  congel_g,  flw_g         
     &,  vsno_g ,  eice_g  ,  sst_g         
     &,  meltt_g,  meltb_g,   meltl_g       
     &,  aice_g,   Tair_g  ,  tsfc_g        
     &,  fhnet_g,  Qa_g    ,  evap_g
     &,  swabs_g,  lwout_g )

      endif                     ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------
      if (my_task == master_task) then

       write (nu_diag,899) 'Arctic','Antarctic'

       write (nu_diag,900) 'max cH (m)             = ',hmaxn,hmaxs
       write (nu_diag,901) 'total area (km^2)      = ',arean,areas
       write (nu_diag,901) 'total ice volume (m^3) = ',shmaxn,shmaxs
       write (nu_diag,901) 'total snw volume (m^3) = ',snwmxn,snwmxs
       write (nu_diag,900) 'average albedo         = ',atotn,atots
       write (nu_diag,900) 'max u, v (m/s)         = ',umaxn,umaxs

       if (print_global) then

       write (nu_diag,*) '----------------------------'
       write (nu_diag,901) 'arwt rain h2o kg in dt = ',rnn,rns
       write (nu_diag,901) 'arwt snow h2o kg in dt = ',snn,sns
       write (nu_diag,901) 'arwt evap h2o kg in dt = ',evpn,evps
       write (nu_diag,901) 'arwt frzl h2o kg in dt = ',frzn,frzs
       write (nu_diag,901) 'arwt frsh h2o kg in dt = ',sfreshn,sfreshs
         
       write (nu_diag,901) 'arwt ice mass (kg)     = ',micen,mices
       write (nu_diag,901) 'arwt snw mass (kg)     = ',msnwn,msnws

       write (nu_diag,901) 'arwt tot mass (kg)     = ',mtotn,mtots
       write (nu_diag,901) 'arwt tot mass chng(kg) = ',delmin,delmis
       write (nu_diag,901) 'arwt water flux        = ',fluxn,fluxs
       write (nu_diag,*) '(=rain+snow+evap+frzl-fresh)  '
       write (nu_diag,901) 'water flux error       = ',werrn,werrs
       
       write (nu_diag,*) '----------------------------'
       write (nu_diag,901) 'arwt atm heat flux (W) = ',fhatmn,fhatms
       write (nu_diag,901) 'arwt ocn heat flux (W) = ',fhocnn,fhocns
       write (nu_diag,901) 'arwt frzl heat flux(W) = ',fhfrzn,fhfrzs
       write (nu_diag,901) 'arwt tot energy    (J) = ',etotn,etots
       write (nu_diag,901) 'arwt net heat      (J) = ',hnetn,hnets
       write (nu_diag,901) 'arwt tot energy chng(J)= ',delein,deleis
       write (nu_diag,901) 'arwt heat error        = ',herrn,herrs
       
       write (nu_diag,*) '----------------------------'
       write (nu_diag,901) 'arwt salt mass (kg)    = ',msltn,mslts
       write (nu_diag,901) 'arwt salt mass chng(kg)= ',delmsltn,
     &                                                 delmslts
       write (nu_diag,901) 'arwt salt flx in dt(kg)= ',sfsaltn,
     &                                                 sfsalts
       write (nu_diag,901) 'arwt salt flx error    = ',serrn,serrs
       write (nu_diag,*) '----------------------------'

       endif                     ! print_global

#ifdef coupled
       call shr_sys_flush(nu_diag)
#endif

       if (print_points) then

       write(nu_diag,*) '                         '
       write(nu_diag,902) '       Lat, Long         ',plat(1),plon(1),
     &                                                plat(2),plon(2)
       write(nu_diag,903) '     my_task, i, j       ',
     &                                     pmloc(1),piloc(1),pjloc(1),
     &                                     pmloc(2),piloc(2),pjloc(2)
       write(nu_diag,*) '----------atm----------'
       write(nu_diag,900) 'air temperature (C)    = ',pTair(1),pTair(2)
       write(nu_diag,900) 'specific humidity      = ',pQa(1),pQa(2)
       write(nu_diag,900) 'snowfall (m)           = ',pFs(1),pFs(2)
       write(nu_diag,900) 'total snow to date (m) = ',pFs_tot(1)
     &                                               ,pFs_tot(2)
       write(nu_diag,900) 'rainfall (m)           = ',pFr(1),pFr(2)
       write(nu_diag,900) 'shortwave radiation sum= ',pFsw(1),pFsw(2)
       write(nu_diag,900) 'longwave radiation     = ',pFlw(1),pFlw(2)
       write(nu_diag,*) '----------ice----------'
       write(nu_diag,900) 'area fraction          = ',par(1),par(2)
       write(nu_diag,900) 'avg ice thickness (m)  = ',hiavg(1),hiavg(2)
       write(nu_diag,900) 'avg snow depth (m)     = ',hsavg(1),hsavg(2)
       write(nu_diag,900) 'surface temperature(C) = ',pTsf(1),pTsf(2)
       write(nu_diag,900) 'absorbed shortwave flx = ',pswabs(1),
     &                                                pswabs(2)
       write(nu_diag,900) 'outward longwave flx   = ',plwout(1),
     &                                                plwout(2)
       write(nu_diag,900) 'sensible heat flx      = ',pFsb(1),pFsb(2)
       write(nu_diag,900) 'latent heat flx        = ',pFlb(1),pFlb(2)
       write(nu_diag,900) 'subl/cond (m ice)      = ',psub(1),psub(2)
       write(nu_diag,900) 'top melt (m)           = ',pmeltt(1)
     &                                               ,pmeltt(2)
       write(nu_diag,900) 'bottom melt (m)        = ',pmeltb(1)
     &                                               ,pmeltb(2)
       write(nu_diag,900) 'lateral melt (m)       = ',pmeltl(1)
     &                                               ,pmeltl(2)
       write(nu_diag,900) 'new ice (m)            = ',pfraz(1),pfraz(2)
       write(nu_diag,900) 'congelation (m)        = ',pcong(1),pcong(2)
       write(nu_diag,900) 'snow-ice (m)           = ',psnoi(1),psnoi(2)
       write(nu_diag,900) 'effective dhi (m)      = ',pdhi(1),pdhi(2)
       write(nu_diag,900) 'effective dhs (m)      = ',pdhs(1),pdhs(2)
       write(nu_diag,900) 'intnl enrgy chng(W/m^2)= ',pde(1),pde(2)
       write(nu_diag,*) '----------ocn----------'
       write(nu_diag,900) 'sst (C)                = ',psst(1),psst(2)
       write(nu_diag,900) 'freezing temp (C)      = ',pTf(1),pTf(2)
       write(nu_diag,900) 'heat used (W/m^2)      = ',pFw(1),pFw(2)

       endif                    ! print_points
      endif                     ! my_task = master_task

  899 format (27x,a24,2x,a24)
  900 format (a25,2x,f24.17,2x,f24.17)
  901 format (a25,2x,1pe24.17,2x,1pe24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
  903 format (a25,10x,f4.0,1x,f4.0,1x,f4.0,9x,f4.0,1x,f4.0,1x,f4.0)

      end subroutine runtime_diags

!=======================================================================
!BOP
!
! !IROUTINE: init_mass_diags - computes global combined ice and snow mass sum
!
! !INTERFACE:
!
      subroutine init_mass_diags
!
! !DESCRIPTION:
!
! Computes global combined ice and snow mass sum
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_mpi_internal
      use ice_grid
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: n, k, ii, jj, i, j

      real (kind=dbl_kind) ::
     &   shmaxn, snwmxn,  shmaxs, snwmxs

      ! total ice volume
      call get_sum(0,tarean,one,vice,shmaxn)
      call get_sum(0,tareas,one,vice,shmaxs)

      ! total snow volume
      call get_sum(0,tarean,one,vsno,snwmxn)
      call get_sum(0,tareas,one,vsno,snwmxs)

      ! north/south ice mass
      totmin = rhoi*shmaxn
      totmis = rhoi*shmaxs

      ! north/south ice+snow mass
      totmn = totmin + rhos*snwmxn
      totms = totmis + rhos*snwmxs

      ! north/south ice+snow energy
      ! total ice/snow energy
      do j=1,jmt_local
      do i=1,imt_local
         work_l2(i,j) = esno(i,j) + eice(i,j)
      enddo
      enddo

      call get_sum(0,tarean,one,work_l2,toten)
      call get_sum(0,tareas,one,work_l2,totes)

      if (print_points) then

         allocate (work_g2(imt_global,jmt_global))

         call global_gather(work_g2,vice(ilo:ihi,jlo:jhi))
         do n = 1,npnt
            if (my_task == master_task) then
               ii = iindx(n)
               jj = jindx(n)
               ! save north/south ice effective thickness (m)
               hefi(n) = work_g2(ii,jj)
            endif
         enddo                  ! n

         call global_gather(work_g2,vsno(ilo:ihi,jlo:jhi))
         do n = 1,npnt
            if (my_task == master_task) then
               ii = iindx(n)
               jj = jindx(n)
               ! save north/south snow effective thickness (m)
               hefs(n) = work_g2(ii,jj)
            endif
         enddo                  ! n

         do n = 1,npnt
            ! save north/south effective energy (J/m^2)
            do j=jlo,jhi
            do i=ilo,ihi
               worka(i,j) = esno(i,j) + eice(i,j)
            enddo
            enddo
            call global_gather(work_g2,worka)
            if (my_task == master_task) then
               ii = iindx(n)
               jj = jindx(n)
               tote(n) = work_g2(ii,jj)
            endif
            if (istep==1) pFs_tot(n) = c0
         enddo                  ! n

         deallocate (work_g2)

      endif                     ! print_points

      end subroutine init_mass_diags

!=======================================================================
!BOP
!
! !IROUTINE: init_diags - find tasks for requested points
!
! !INTERFACE:
!
      subroutine init_diags
!
! !DESCRIPTION:
!
!  Find tasks for requested points.
! 
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
      use ice_grid
      use ice_mpi_internal
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind) ::
     &   minbth      ! sum of distances in lat and lon
     &,  mindis      ! min(minbth) over process tasks
     &,  latdis      ! latitude distance
     &,  londis      ! longitude distance
     &,  totdis      ! total distance

      integer (kind=int_kind) ::
     &   n           ! index for search
     &,  i,j         ! grid indices

      character (char_len) :: label(npnt)

      if (print_points) then

         allocate (work_g2(imt_global,jmt_global))

         ! initialize labels
         if (my_task == master_task)
     &        write(nu_diag,*) ' initialize diagnostic prints '

         label(1)(1:40)  = 'Near North Pole pack ice                '
         label(2)(1:40)  = 'Weddell Sea                             '

         ! gather mask on global grid
         call global_gather(work_g2,hm(ilo:ihi,jlo:jhi))

         if (my_task == master_task) then

            write(nu_diag,*) ' find indices of points '
            do n=1,npnt
               ! compute indices
               minbth   = 540.0_dbl_kind !  360.0 + 180.
               iindx(n) = 0
               jindx(n) = 0
               do j=1,jmt_global
               do i=1,imt_global
                  if (work_g2(i,j) > p5) then
                     latdis = abs(latpnt(n)-TLAT_G(i,j)*rad_to_deg)
                     londis = abs(lonpnt(n)-TLON_G(i,j)*rad_to_deg)
                     totdis = sqrt(latdis**2 + londis**2)
                     if( totdis < minbth ) then
                        minbth   = totdis
                        jindx(n) = j
                        iindx(n) = i
                     endif      ! totdis < minbth
                  endif         ! work_g2 > p5
               enddo            ! i
               enddo            ! j

               write(nu_diag,2234) n,latpnt(n),lonpnt(n),
     &              TLAT_G(iindx(n),jindx(n))*rad_to_deg,
     &              TLON_G(iindx(n),jindx(n))*rad_to_deg,
     &              iindx(n),jindx(n)
 2234          format(' found point',i4/
     &              '   lat     lon     TLAT    TLON    iindx jindx '/
     &              4(f7.1,1x),4x,2(i3,3x))

            enddo               ! n
         endif                  ! master_task

         ! find points on local domain processors
         do j=jlo,jhi
         do i=ilo,ihi
            worka(i,j) = real(i + my_task*10000)
            workb(i,j) = real(j + my_task*10000)
         enddo
         enddo
         call global_gather(work_g2,workb)
         if (my_task == master_task) then
           do n=1,npnt
             ! Lat, Long
             plat(n) = TLAT_G(iindx(n),jindx(n))*rad_to_deg
             plon(n) = TLON_G(iindx(n),jindx(n))*rad_to_deg

             ! location on processor
             pjloc(n) = mod(work_g2(iindx(n),jindx(n)),10000._dbl_kind) ! j
           enddo
         endif                  ! master_task
         call global_gather(work_g2,worka)
         if (my_task == master_task) then
           do n=1,npnt
             ! location on processor
             piloc(n) = mod(work_g2(iindx(n),jindx(n)),10000._dbl_kind) ! i
             pmloc(n) = (work_g2(iindx(n),jindx(n))-piloc(n))/10000 ! my_task
           enddo
         endif                  ! master_task

         deallocate (work_g2)

      endif                     ! print_points

      end subroutine init_diags

!=======================================================================
!BOP
!
! !IROUTINE: print_state - print ice state for specified grid point
!
! !INTERFACE:
!
      subroutine print_state(plabel,i,j)
!
! !DESCRIPTION:
!
! This routine is useful for debugging \\
! Calls to it should be inserted in the form (after thermo, for example) \\
!      do j=jlo,jhi                                         \\
!       do i=ilo,ihi                                        \\
!         plabel = 'post thermo'                            \\
!         if (istep1 >= check\_step .and. i==ip .and. j==jp  \\
!     &       .and. my\_task == mtask)                       \\
!     &   call print\_state(plabel,i,j)                      \\
!       enddo                                               \\
!      enddo                                                \\
!
! 'use ice\_diagnostics' may need to be inserted also \\
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_model_size
      use ice_kinds_mod
      use ice_state
      use ice_itd
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (len=20), intent(in) :: plabel
      integer (kind=int_kind), intent(in) :: i,j
!
!EOP
!
      real (kind=dbl_kind) :: 
     &     edebug
     &,    qs, Tsnow

      integer (kind=int_kind) :: n, k

      write(nu_diag,*) plabel,' istep1 = ',istep1
      write(nu_diag,*) 'my_task',my_task,' i, j=',i,j 

      write(nu_diag,*) 'aice0', aice0(i,j)
      do n = 1,ncat
         write(nu_diag,*) ' '
         write(nu_diag,*) 'n =',n
         write(nu_diag,*) 'aicen', aicen(i,j,n)
         write(nu_diag,*) 'vicen', vicen(i,j,n)
         write(nu_diag,*) 'vsnon', vsnon(i,j,n)
         if (aicen(i,j,n) > puny) then
            write(nu_diag,*) 'hin', vicen(i,j,n)/aicen(i,j,n)
            write(nu_diag,*) 'hsn', vsnon(i,j,n)/aicen(i,j,n)
         endif
         write(nu_diag,*) 'Tsfcn',Tsfcn(i,j,n)
         write(nu_diag,*) ' '
         write(nu_diag,*) 'esnon',esnon(i,j,n)
         if (vsnon(i,j,n) > puny) then
            qs = -esnon(i,j,n)/vsnon(i,j,n)
            Tsnow = (Lfresh - qs/rhos) / cp_ice
            write(nu_diag,*) 'qs', qs
            write(nu_diag,*) 'Tsnow', Tsnow
            write(nu_diag,*) ' '
         endif
      enddo                     ! n

      edebug = c0
      do n = 1,ncat
         do k = 1,nilyr
            write(nu_diag,*) 'eicen, cat ',n,' layer ',k,
     &           eicen(i,j,ilyr1(n)+k-1)
            edebug = edebug + eicen(i,j,ilyr1(n)+k-1)
         enddo
      enddo
      write(nu_diag,*) 'eice(i,j)',edebug 
      write(nu_diag,*) ' '
      write(nu_diag,*) 'atm states and fluxes'
      write(nu_diag,*) '            uatm  = ',uatm(i,j)
      write(nu_diag,*) '            vatm  = ',vatm(i,j)
      write(nu_diag,*) '            potT  = ',potT(i,j)
      write(nu_diag,*) '            Tair  = ',Tair(i,j)
      write(nu_diag,*) '            Qa    = ',Qa(i,j)
      write(nu_diag,*) '            rhoa  = ',rhoa(i,j)
      write(nu_diag,*) '            swvdr = ',swvdr(i,j)
      write(nu_diag,*) '            swvdf = ',swvdf(i,j)
      write(nu_diag,*) '            swidr = ',swidr(i,j)
      write(nu_diag,*) '            swidf = ',swidf(i,j)
      write(nu_diag,*) '            flw   = ',flw(i,j)
      write(nu_diag,*) '            frain = ',frain(i,j)
      write(nu_diag,*) '            fsnow = ',fsnow(i,j)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'ocn states and fluxes'
      write(nu_diag,*) '            frzmlt  = ',frzmlt(i,j)
      write(nu_diag,*) '            sst     = ',sst(i,j)
      write(nu_diag,*) '            sss     = ',sss(i,j)
      write(nu_diag,*) '            Tf      = ',Tf(i,j)
      write(nu_diag,*) '            uocn    = ',uocn(i,j)
      write(nu_diag,*) '            vocn    = ',vocn(i,j)
      write(nu_diag,*) '            strtltx = ',strtltx(i,j)
      write(nu_diag,*) '            strtlty = ',strtlty(i,j)
      write(nu_diag,*) ' '
      write(nu_diag,*) 'srf states and fluxes'
      write(nu_diag,*) '            Tref    = ',Tref(i,j)
      write(nu_diag,*) '            Qref    = ',Qref(i,j)
      write(nu_diag,*) '            fsens   = ',fsens(i,j)
      write(nu_diag,*) '            flat    = ',flat(i,j)
      write(nu_diag,*) '            evap    = ',evap(i,j)
      write(nu_diag,*) '            flwout  = ',flwout(i,j)
      write(nu_diag,*) ' '
      
      end subroutine print_state

!=======================================================================

      end module ice_diagnostics

!=======================================================================
