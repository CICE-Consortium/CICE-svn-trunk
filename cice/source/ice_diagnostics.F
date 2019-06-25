! $Id: ice_diagnostics.F,v 1.27 2006/08/02 19:45:56 eclare Exp $
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
      use ice_work, only: work_l1, work_l2, worka, workb
!
!EOP
!
      implicit none
      save

      ! diagnostic output file
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
     &   pdhi                 ! change in mean ice thickness (m)
     &,  pdhs                 ! change in mean snow thickness (m)
     &,  pde                  ! change in ice and snow energy (J m-2)
     &,  plat, plon           ! latitude, longitude of points

      integer (kind=int_kind), dimension(npnt), save ::
     &   piloc, pjloc, pmloc  ! location of points on processors

      ! for hemispheric water and heat budgets
      real (kind=dbl_kind) ::
     &   totmn                ! total ice/snow water mass (nh)
     &,  totms                ! total ice/snow water mass (sh)
     &,  totmin               ! total ice water mass (nh)
     &,  totmis               ! total ice water mass (sh)
     &,  toten                ! total ice/snow energy (J)
     &,  totes                ! total ice/snow energy (J)

      data latpnt /   90._dbl_kind,  -65._dbl_kind /
      data lonpnt /    0._dbl_kind,  -45._dbl_kind /

      ! printing info for routine print_state
      character (char_len) :: plabel
      integer (kind=int_kind), parameter ::
     &   check_step = 999999999
     &,  ip = 91
     &,  jp = 382
     &,  mtask = 0

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
      use ice_itd
#ifdef CCSMcoupled
      use shr_sys_mod, only : shr_sys_flush
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
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

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) ::
     &   paice, pTair, pQa, pfsnow, pfrain, pfsw, pflw 
     &,  pTsfc, pevap, pfswabs, pflwout, pflat, pfsens 
     &,  psst,  pTf, hiavg, hsavg, pfhocn 
     &,  pmeltt, pmeltb, pmeltl, psnoice, pfrazil, pcongel

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
!      do j=1,jmt_local
!      do i=1,imt_local
!         work_l1(i,j) = c0
!         if (aice(i,j) >= p15) work_l1(i,j) = c1
!      enddo
!      enddo
!      call get_sum(0,tarean,one,work_l1,extentn)
!      call get_sum(0,tareas,one,work_l1,extents)
!      extentn = extentn * m2_to_km2
!      extents = extents * m2_to_km2

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
         atotn = atotn / arean * m2_to_km2 ! convert arean back to m2
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

      if (print_global) then

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

      endif ! print_global

      if (print_points) then

      !-----------------------------------------------------------------
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)

               pTair(n) = Tair(i,j) - Tffresh    ! air temperature
               pQa(n) = Qa(i,j)                  ! specific humidity
               pfsnow(n) = fsnow(i,j)*dt/rhos    ! snowfall
               pfrain(n) = frain(i,j)*dt/rhow    ! rainfall
               pfsw(n) = fsw(i,j)                ! shortwave radiation
               pflw(n) = flw(i,j)                ! longwave radiation
               paice(n) = aice(i,j)              ! ice area
               
               hiavg(n) = c0                     ! avg snow/ice thickness
               hsavg(n) = c0
               if (paice(n) /= c0) then
                  hiavg(n) = vice(i,j)/paice(n)
                  hsavg(n) = vsno(i,j)/paice(n)
               endif
               pTsfc(n) = Tsfc(i,j)              ! ice/snow sfc temperature
               pevap(n) = evap(i,j)*dt/rhoi      ! sublimation/condensation
               pfswabs(n) = fswabs(i,j)          ! absorbed solar flux
               pflwout(n) = flwout(i,j)          ! outward longwave flux
               pflat(n) = flat(i,j)              ! latent heat flux
               pfsens(n) = fsens(i,j)            ! sensible heat flux
               pmeltt(n) = meltt(i,j)            ! top melt
               pmeltb(n) = meltb(i,j)            ! bottom melt
               pmeltl(n) = meltl(i,j)            ! lateral melt
               psnoice(n) = snoice(i,j)          ! snow ice
               pfrazil(n) = frazil(i,j)          ! frazil ice
               pcongel(n) = congel(i,j)          ! congelation ice
               pdhi(n) = vice(i,j) - pdhi(n)     ! ice thickness change
               pdhs(n) = vsno(i,j) - pdhs(n)     ! snow thickness change
               pde(n) = -(eice(i,j)              ! ice/snow energy change 
     &                  + esno(i,j) - pde(n)) / dt
               psst(n) = sst(i,j)                ! sea surface temperature
               pTf(n) = Tf(i,j)                  ! freezing temperature
               pfhocn(n) = -fhnet_hist(i,j)      ! ocean heat used by ice

            endif  ! my_task = pmloc

            call ice_bcast_rscalar(pTair   (n), pmloc(n))
            call ice_bcast_rscalar(pQa     (n), pmloc(n))
            call ice_bcast_rscalar(pfsnow  (n), pmloc(n))
            call ice_bcast_rscalar(pfrain  (n), pmloc(n))
            call ice_bcast_rscalar(pfsw    (n), pmloc(n))
            call ice_bcast_rscalar(pflw    (n), pmloc(n))
            call ice_bcast_rscalar(paice   (n), pmloc(n))
            call ice_bcast_rscalar(hsavg   (n), pmloc(n))
            call ice_bcast_rscalar(hiavg   (n), pmloc(n))
            call ice_bcast_rscalar(pTsfc   (n), pmloc(n))
            call ice_bcast_rscalar(pevap   (n), pmloc(n))
            call ice_bcast_rscalar(pfswabs (n), pmloc(n))
            call ice_bcast_rscalar(pflwout (n), pmloc(n))
            call ice_bcast_rscalar(pflat   (n), pmloc(n))
            call ice_bcast_rscalar(pfsens  (n), pmloc(n))
            call ice_bcast_rscalar(pmeltt  (n), pmloc(n))
            call ice_bcast_rscalar(pmeltb  (n), pmloc(n))
            call ice_bcast_rscalar(pmeltl  (n), pmloc(n))
            call ice_bcast_rscalar(psnoice (n), pmloc(n))
            call ice_bcast_rscalar(pfrazil (n), pmloc(n))
            call ice_bcast_rscalar(pcongel (n), pmloc(n))
            call ice_bcast_rscalar(pdhi    (n), pmloc(n))
            call ice_bcast_rscalar(pdhs    (n), pmloc(n))
            call ice_bcast_rscalar(pde     (n), pmloc(n))
            call ice_bcast_rscalar(psst    (n), pmloc(n))
            call ice_bcast_rscalar(pTf     (n), pmloc(n))
            call ice_bcast_rscalar(pfhocn  (n), pmloc(n))
            
         enddo                  ! npnt
      endif                     ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------
      if (my_task == master_task) then

       write (nu_diag,899) 'Arctic','Antarctic'

       write (nu_diag,900) 'max cH (m)             = ',hmaxn,hmaxs
!      write (nu_diag,901) 'ice extent (km^2)      = ',extentn, extents
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

       if (print_points) then

        write(nu_diag,*) '                         '
        write(nu_diag,902) '       Lat, Long         ',plat(1),plon(1),
     &                                                 plat(2),plon(2)
        write(nu_diag,903) '  my_task, i, j     ',
     &                        pmloc(1),piloc(1),pjloc(1),
     &                        pmloc(2),piloc(2),pjloc(2)
        write(nu_diag,*) '----------atm----------'
        write(nu_diag,900) 'air temperature (C)    = ',pTair(1),pTair(2)
        write(nu_diag,900) 'specific humidity      = ',pQa(1),pQa(2)
        write(nu_diag,900) 'snowfall (m)           = ',pfsnow(1),
     &                                                 pfsnow(2)
        write(nu_diag,900) 'rainfall (m)           = ',pfrain(1),
     &                                                 pfrain(2)
        write(nu_diag,900) 'shortwave radiation sum= ',pfsw(1),pfsw(2)
        write(nu_diag,900) 'longwave radiation     = ',pflw(1),pflw(2)
        write(nu_diag,*) '----------ice----------'
        write(nu_diag,900) 'area fraction          = ',paice(1),paice(2)
        write(nu_diag,900) 'avg ice thickness (m)  = ',hiavg(1),hiavg(2)
        write(nu_diag,900) 'avg snow depth (m)     = ',hsavg(1),hsavg(2)
        write(nu_diag,900) 'surface temperature(C) = ',pTsfc(1),pTsfc(2)
        write(nu_diag,900) 'absorbed shortwave flx = ',pfswabs(1),
     &                                                 pfswabs(2)
        write(nu_diag,900) 'outward longwave flx   = ',pflwout(1),
     &                                                 pflwout(2)
        write(nu_diag,900) 'sensible heat flx      = ',pfsens(1),
     &                                                 pfsens(2)
        write(nu_diag,900) 'latent heat flx        = ',pflat(1),pflat(2)
        write(nu_diag,900) 'subl/cond (m ice)      = ',pevap(1),pevap(2)
        write(nu_diag,900) 'top melt (m)           = ',pmeltt(1)
     &                                                ,pmeltt(2)
        write(nu_diag,900) 'bottom melt (m)        = ',pmeltb(1)
     &                                                ,pmeltb(2)
        write(nu_diag,900) 'lateral melt (m)       = ',pmeltl(1)
     &                                                ,pmeltl(2)
        write(nu_diag,900) 'new ice (m)            = ',pfrazil(1),
     &                                                 pfrazil(2)
        write(nu_diag,900) 'congelation (m)        = ',pcongel(1),
     &                                                 pcongel(2)
        write(nu_diag,900) 'snow-ice (m)           = ',psnoice(1),
     &                                                 psnoice(2)
        write(nu_diag,900) 'effective dhi (m)      = ',pdhi(1),pdhi(2)
        write(nu_diag,900) 'effective dhs (m)      = ',pdhs(1),pdhs(2)
        write(nu_diag,900) 'intnl enrgy chng(W/m^2)= ',pde (1),pde (2)
        write(nu_diag,*) '----------ocn----------'
        write(nu_diag,900) 'sst (C)                = ',psst(1),psst(2)
        write(nu_diag,900) 'freezing temp (C)      = ',pTf(1),pTf(2)
        write(nu_diag,900) 'heat used (W/m^2)      = ',pfhocn(1),
     &                                                 pfhocn(2)

       endif                    ! print_points
      endif                     ! my_task = master_task

#ifdef CCSMcoupled
       call shr_sys_flush(nu_diag)
#endif

  899 format (27x,a24,2x,a24)
  900 format (a25,2x,f24.17,2x,f24.17)
  901 format (a25,2x,1pe24.17,2x,1pe24.17)
  902 format (a25,10x,f6.1,1x,f6.1,9x,f6.1,1x,f6.1)
  903 format (a25,10x,i4,1x,i4,1x,i4,9x,i4,1x,i4,1x,i4)

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

         do n = 1, npnt

            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)

               pdhi(n) = vice(i,j)
               pdhs(n) = vsno(i,j)
               pde(n)  = esno(i,j) + eice(i,j)
            endif

         enddo  ! npnt

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
     &,  latdis      ! latitude distance
     &,  londis      ! longitude distance
     &,  totdis      ! total distance
     &,  mindis      ! minimum distance from desired location
     &,  mindis_g    ! global minimum distance from desired location

      real (kind=dbl_kind), dimension(1) ::
     &   wrk         ! for maxval and minval calls that require an array

      integer (kind=int_kind) ::
     &   n           ! index for point search
     &,  i,j         ! grid indices

      character (char_len) :: label(npnt)

      if (print_points) then

         ! initialize labels
         if (my_task == master_task) then
            write(nu_diag,*) ' '
            write(nu_diag,*) 'Find indices of diagnostic points '
         endif

         label(1)(1:40)  = 'Near North Pole pack ice                '
         label(2)(1:40)  = 'Weddell Sea                             '

         piloc(:) = 0
         pjloc(:) = 0
         pmloc(:) = -999
         plat(:)  = -999._dbl_kind
         plon(:)  = -999._dbl_kind

         ! find minimum distance to diagnostic points on this processor 
         do n = 1, npnt
            iindx(n) = 0
            jindx(n) = 0
            mindis = 540.0_dbl_kind !  360. + 180.

            do j = jlo, jhi
               do i = ilo, ihi
                  if (hm(i,j) > p5) then
                     latdis = abs(latpnt(n)-TLAT(i,j)*rad_to_deg)
                     londis = abs(lonpnt(n)-TLON(i,j)*rad_to_deg)
     &                      * cos(TLON(i,j))
                     totdis = sqrt(latdis**2 + londis**2)
                     if (totdis < mindis) then
                        mindis = totdis
                        jindx(n) = j
                        iindx(n) = i
                     endif      ! totdis < mindis
                  endif         ! hm > p5
               enddo            ! i
            enddo               ! j

            ! find global minimum distance to diagnostic points 
            wrk(1) = mindis
            mindis_g = ice_global_real_minval(1,wrk)

            ! save indices of minimum-distance grid cell
            if (abs(mindis_g - mindis) < puny) then
               piloc(n) = iindx(n)
               pjloc(n) = jindx(n)
               pmloc(n) = my_task
               plat(n) = TLAT(iindx(n),jindx(n))*rad_to_deg
               plon(n) = TLON(iindx(n),jindx(n))*rad_to_deg
            endif

            ! communicate to all processors
            ! need integer/real conversions because there is no
            !  maxval routine for integers in version 3.14
            wrk(1) = real(piloc(n),kind=dbl_kind)
            piloc(n) = nint(ice_global_real_maxval(1,wrk))
            wrk(1) = real(pjloc(n),kind=dbl_kind)
            pjloc(n) = nint(ice_global_real_maxval(1,wrk))
            wrk(1) = real(pmloc(n),kind=dbl_kind)
            pmloc(n) = nint(ice_global_real_maxval(1,wrk))
            wrk(1) = plat(n)
            plat(n)  = ice_global_real_maxval(1,wrk)
            wrk(1) = plon(n)
            plon(n)  = ice_global_real_maxval(1,wrk)

            ! write to log file
            if (my_task==master_task) then
               write(nu_diag,*) ''
               write(nu_diag,100) n,latpnt(n),lonpnt(n),plat(n),plon(n),
     &              piloc(n), pjloc(n), pmloc(n)
            endif
 100        format(' found point',i4/
     &         '   lat    lon   TLAT   TLON     i     j  my_task'/
     &          4(f6.1,1x),1x,3(i4,2x) )

         enddo                  ! npnt
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
