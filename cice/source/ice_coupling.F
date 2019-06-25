! CVS: $Id: ice_coupling.F,v 1.13 2004/02/25 17:35:13 eclare Exp $
! CVS: $Source: /home/climate/CVS-COSIM/cice/source/ice_coupling.F,v $
! CVS: $Name:  $
!=======================================================================
!
!BOP
!
! !MODULE: ice_coupling - message passing to and from the coupler
!
! !DESCRIPTION:
!
! Message passing to and from the coupler
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!         Tony Craig, NCAR, Dec-30-2002, modified for cpl6
!
! !INTERFACE:
!
      module ice_coupling
!
! !USES:
!
      use ice_kinds_mod
      use ice_model_size
      use ice_constants
      use ice_calendar
      use ice_grid
      use ice_state
      use ice_flux
      use ice_albedo
      use ice_mpi_internal
      use ice_timers
      use ice_fileunits
      use ice_work, only: worka, work_l1
#ifdef coupled
      use shr_sys_mod, only : shr_sys_flush
      use ice_history, only : runtype
      use cpl_contract_mod
      use cpl_interface_mod
      use cpl_fields_mod
#endif
!
!EOP
!

#ifdef coupled  
! #endif at end of module

      implicit none

      integer (kind=int_kind), dimension (cpl_fields_ibuf_total) ::
     &   isbuf
     &,  irbuf

      real (kind=dbl_kind), allocatable ::
     &   sbuf(:,:)

      real (kind=dbl_kind) ::
     &   buffs((jhi-jlo+1)*(ihi-ilo+1),cpl_fields_i2c_total)
     &,  buffr((jhi-jlo+1)*(ihi-ilo+1),cpl_fields_c2i_total)

      type(cpl_contract) ::
     &   contractS
     &,  contractR

      integer(kind=int_kind), save ::
     &   nadv_i
     &,  info_dbug

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: ice_coupling_setup - sets mpi communicators and task ids
!
! !INTERFACE:
!
      subroutine ice_coupling_setup(in_model_name,model_comm)
!
! !DESCRIPTION:
!
! This routine uses get the model communicator from ccsm share code
!
! !REVISION HISTORY:
!
! author: T. Craig, NCAR, Dec 30, 2002: for cpl6
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (3), intent(in) :: in_model_name   

      integer, intent(out) ::  model_comm     ! communicator for model
!
!EOP
!
      write(nu_diag,*) 'calling cpl_interface_init for model: ',
     &     in_model_name,' ', trim(cpl_fields_icename)

      call cpl_interface_init(cpl_fields_icename,model_comm)

      call shr_sys_flush(nu_diag)

      end subroutine ice_coupling_setup

!=======================================================================
!BOP
!
! !IROUTINE: init_cpl - initializes message passing between ice and coupler
!
! !INTERFACE:
!
      subroutine init_cpl
!
! !DESCRIPTION:
!
! Initializes message passing between ice and coupler
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer(kind=int_kind) :: i,j,n     ! local loop indices

      write(nu_diag,*) 
     &     '(ice_coupling,init_cpl) send initial msg. set contract'
      call shr_sys_flush(nu_diag)

      nadv_i = nint(secday/dt)

      isbuf                          = 0         ! default info-buffer value
      isbuf(cpl_fields_ibuf_cdate  ) = idate     ! initial date (coded: yyyymmdd)
      isbuf(cpl_fields_ibuf_sec    ) = sec       ! elapsed seconds into date
      isbuf(cpl_fields_ibuf_stopnow) = stop_now  ! stop now flag
      isbuf(cpl_fields_ibuf_userest) = 0         ! use model restart data initally
      isbuf(cpl_fields_ibuf_ncpl   ) = nadv_i    ! number of comms per day
      isbuf(cpl_fields_ibuf_lsize  ) = (ihi-ilo+1)*(jhi-jlo+1) ! size of local grid
      isbuf(cpl_fields_ibuf_lisize ) = (ihi-ilo+1) ! local size wrt i-index
      isbuf(cpl_fields_ibuf_ljsize ) = (jhi-jlo+1) ! local size wrt i-index
      isbuf(cpl_fields_ibuf_gsize  ) = imt_global*jmt_global ! size of global grid
      isbuf(cpl_fields_ibuf_gisize ) = imt_global  ! global size wrt i-index
      isbuf(cpl_fields_ibuf_gjsize ) = jmt_global  ! global size wrt j-index
      isbuf(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
      isbuf(cpl_fields_ibuf_dead   ) = 0           ! not a dead model

      allocate(sbuf((ihi-ilo+1)*(jhi-jlo+1),cpl_fields_grid_total))
      sbuf = -888.0
      n=0
      do j=jlo,jhi
      do i=ilo,ihi
         n=n+1
         sbuf(n,cpl_fields_grid_lon  ) = TLON(i,j)*rad_to_deg
         sbuf(n,cpl_fields_grid_lat  ) = TLAT(i,j)*rad_to_deg
         sbuf(n,cpl_fields_grid_area ) = tarea(i,j)/(radius*radius)
         sbuf(n,cpl_fields_grid_mask ) = float(nint(hm(i,j)))
         sbuf(n,cpl_fields_grid_index) = rndex_global(i,j)
      enddo
      enddo

      call cpl_interface_contractInit
     &     (contractS, cpl_fields_icename, cpl_fields_cplname,
     &      cpl_fields_i2c_fields, isbuf, sbuf)

      call cpl_interface_contractInit
     &     (contractR, cpl_fields_icename, cpl_fields_cplname,
     &      cpl_fields_c2i_fields, isbuf, sbuf)

      write(nu_diag,*) '(init_cpl) Initialized contracts with coupler'
      call shr_sys_flush(nu_diag)

      !-----------------------------------------------------------------
      ! Receive initial message from coupler.
      !-----------------------------------------------------------------

      call cpl_interface_ibufRecv(cpl_fields_cplname,irbuf)

      if (my_task==master_task) then
         write(nu_diag,*)
     &        '(init_cpl) Received control buffer from coupler'
         call shr_sys_flush(nu_diag)

         if (trim(runtype)=='startup' .or.
     &       trim(runtype)== 'hybrid') then
            idate = irbuf(cpl_fields_ibuf_cdate)
            write(nu_diag,*) '(init_cpl) idate from coupler = ',idate
            nyr   = (idate/10000)               ! integer year of basedate
            month = (idate-nyr*10000)/100       ! integer month of basedate
            mday  = idate-nyr*10000-month*100-1 ! day of year of basedate
            time  = ((nyr-1)*daycal(13)+daycal(month)+mday)*secday
            call calendar(time)                 ! recompute calendar info
            time_forc = time
            call shr_sys_flush(nu_diag)
         endif

      endif                     ! my_task==master_task

      call ice_bcast_iscalar(idate)
      call ice_bcast_rscalar(time)
      call ice_bcast_rscalar(time_forc)

      deallocate(sbuf)

      write(nu_diag,*) '(ice_coupling,init_cpl) done setting contract'

      !-----------------------------------------------------------------
      ! Send initial state info to coupler.
      !-----------------------------------------------------------------

      call to_coupler

      end subroutine init_cpl

!=======================================================================
!BOP
!
! !IROUTINE: from_coupler - input from coupler to sea ice model
!
! !INTERFACE:
!
      subroutine from_coupler
!
! !DESCRIPTION:
!
! Reads input data from coupler to sea ice model
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i,j,n,n2   ! local loop indices

      real (kind=dbl_kind) ::  
     &   gsum, workx, worky

      call ice_timer_start(8)  ! time spent coupling

      !-----------------------------------------------------------------
      ! Zero stuff while waiting, only filling in active cells.
      !-----------------------------------------------------------------

      zlvl(:,:)    = c0
      uatm(:,:)    = c0
      vatm(:,:)    = c0
      potT(:,:)    = c0
      Tair(:,:)    = c0
      Qa(:,:)      = c0
      rhoa(:,:)    = c0
      swvdr(:,:)   = c0
      swvdf(:,:)   = c0
      swidr(:,:)   = c0
      swidf(:,:)   = c0
      flw(:,:)     = c0
      frain(:,:)   = c0
      fsnow(:,:)   = c0
      sst(:,:)     = c0
      sss(:,:)     = c0
      uocn(:,:)    = c0
      vocn(:,:)    = c0
      ss_tltx(:,:) = c0
      ss_tlty(:,:) = c0
      frzmlt(:,:)  = c0

      !-----------------------------------------------------------------
      ! recv input field msg
      !-----------------------------------------------------------------
     
      call ice_timer_start(16)  ! time spent receiving

      call cpl_interface_contractRecv
     &     (cpl_fields_cplname, contractR, irbuf, buffr)

      call ice_timer_stop(16)
      call ice_timer_start(17)  ! time spent cr-unpacking

      !--- unpack message
      n=0
      do j=jlo,jhi
      do i=ilo,ihi
         n=n+1

         !--- ocn states--
         sst  (i,j) = buffr(n,cpl_fields_c2i_ot)
         sss  (i,j) = buffr(n,cpl_fields_c2i_os)
         uocn (i,j) = buffr(n,cpl_fields_c2i_ou)
         vocn (i,j) = buffr(n,cpl_fields_c2i_ov)

         !--- atm states-
         zlvl (i,j) = buffr(n,cpl_fields_c2i_z)
         uatm (i,j) = buffr(n,cpl_fields_c2i_u)
         vatm (i,j) = buffr(n,cpl_fields_c2i_v)
         potT (i,j) = buffr(n,cpl_fields_c2i_ptem)
         Tair (i,j) = buffr(n,cpl_fields_c2i_tbot)
         Qa   (i,j) = buffr(n,cpl_fields_c2i_shum)
         rhoa (i,j) = buffr(n,cpl_fields_c2i_dens)

         !--- ocn states--
         ss_tltx(i,j) = buffr(n,cpl_fields_c2i_dhdx)
         ss_tlty(i,j) = buffr(n,cpl_fields_c2i_dhdy)
         frzmlt (i,j) = buffr(n,cpl_fields_c2i_q)

         !--- atm fluxes--
         swvdr(i,j) = buffr(n,cpl_fields_c2i_swvdr)
         swidr(i,j) = buffr(n,cpl_fields_c2i_swndr)
         swvdf(i,j) = buffr(n,cpl_fields_c2i_swvdf)
         swidf(i,j) = buffr(n,cpl_fields_c2i_swndf)
         flw  (i,j) = buffr(n,cpl_fields_c2i_lwdn)
         frain(i,j) = buffr(n,cpl_fields_c2i_rain)
         fsnow(i,j) = buffr(n,cpl_fields_c2i_snow)

      end do
      end do

      call ice_timer_stop(17)  ! time spent cr-unpacking

      !-----------------------------------------------------------------
      ! broadcast dbug diagnostic level
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_infobug) >= 2 ) then
         if (my_task == master_task) write (nu_diag,*)
     &        '(from_coupler) dbug level >= 2'
         info_dbug = 1
      endif

      !-----------------------------------------------------------------
      ! broadcast write_restart flag
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_resteod) == 1 .AND. new_day) then
         if (my_task == master_task) write (nu_diag,*)
     &        '(from_coupler) received write restart signal'
         write_restart = 1
      endif

      !-----------------------------------------------------------------
      ! broadcast cpl_write_history flag
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_histeod) == 1 .AND. new_day) then
         if (my_task == master_task) write (nu_diag,*)
     &        '(from_coupler) received write history signal'
         cpl_write_history = 1
      endif

      !-----------------------------------------------------------------
      ! broadcast stop_now flag
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_stopnow) == 1) then
         if (my_task==master_task) write (nu_diag,*)
     &        '(from_coupler) received terminate signal'
         stop_now = 1
      endif

      if (info_dbug == 1 .AND. stop_now /= 1) then

        do j=jlo,jhi
        do i=ilo,ihi
           worka(i,j) = tarea(i,j)
        enddo
        enddo

        do n=1,cpl_fields_c2i_total
           work_l1 = c0
           n2   = 0
           do j=jlo,jhi
           do i=ilo,ihi
              n2 = n2 + 1 
              if (hm(i,j) > p5) work_l1(i,j) = buffr(n2,n)
           enddo
           enddo
           call bound(work_l1)
           call get_sum(0, worka, one, work_l1, gsum)
           if (my_task == master_task) then
              write (nu_diag,100) 'ice', 'recv', n, gsum
           endif
        enddo                   ! cpl_fields_c2i_total
      endif
 100  format ('comm_diag',1x,a3,1x,a4,1x,i3,es26.19)

      !-----------------------------------------------------------------
      ! rotate zonal/meridional vectors to local coordinates
      ! compute data derived quantities
      !-----------------------------------------------------------------

      ! Vector fields come in on T grid, but are oriented geographically
      ! need to rotate to pop-grid FIRST using ANGLET
      ! then interpolate to the U-cell centers  (otherwise we
      ! interpolate across the pole)
      ! use ANGLET which is on the T grid !

      do j=jlo,jhi
      do i=ilo,ihi
         ! ocean
         workx      = uocn  (i,j) ! currents, m/s 
         worky      = vocn  (i,j)
         uocn(i,j) = workx*cos(ANGLET(i,j))    ! convert to POP grid 
     &             + worky*sin(ANGLET(i,j))
         vocn(i,j) = worky*cos(ANGLET(i,j))
     &             - workx*sin(ANGLET(i,j))

         workx      = ss_tltx  (i,j)           ! sea sfc tilt, m/m
         worky      = ss_tlty  (i,j)
         ss_tltx(i,j) = workx*cos(ANGLET(i,j)) ! convert to POP grid 
     &                + worky*sin(ANGLET(i,j))
         ss_tlty(i,j) = worky*cos(ANGLET(i,j))
     &                - workx*sin(ANGLET(i,j))

         sst(i,j) = sst(i,j) - Tffresh         ! sea sfc temp (C)
         Tf (i,j) = -1.8_dbl_kind              ! hardwired for NCOM
c        Tf (i,j) = -depressT*sss(i,j)         ! freezing temp (C)
c        Tf (i,j) = -depressT*max(sss(i,j),ice_ref_salinity)

      enddo
      enddo

      ! Interpolate ocean dynamics variables from T-cell centers to 
      ! U-cell centers.

      call t2ugrid(uocn)
      call t2ugrid(vocn)
      call t2ugrid(ss_tltx)
      call t2ugrid(ss_tlty)

      ! Atmosphere variables are needed in T cell centers in
      ! subroutine stability and are interpolated to the U grid
      ! later as necessary.

      do j=jlo,jhi
      do i=ilo,ihi
         ! atmosphere
         workx      = uatm(i,j) ! wind velocity, m/s
         worky      = vatm(i,j) 
         uatm (i,j) = workx*cos(ANGLET(i,j)) ! convert to POP grid
     &              + worky*sin(ANGLET(i,j)) ! note uatm, vatm, wind
         vatm (i,j) = worky*cos(ANGLET(i,j)) !  are on the T-grid here
     &              - workx*sin(ANGLET(i,j))

         wind (i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2) ! wind speed, m/s
         fsw  (i,j) = swvdr(i,j) + swvdf(i,j)
     &              + swidr(i,j) + swidf(i,j)
      enddo
      enddo

      time_forc=time

      call ice_timer_stop(8)   ! time spent coupling

      end subroutine from_coupler

!=======================================================================
!BOP
!
! !IROUTINE: to_coupler - send data from sea ice model to coupler
!
! !INTERFACE:
!
      subroutine to_coupler
!
! !DESCRIPTION:
!
! Sea ice model to coupler
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP

      integer(kind=int_kind) :: i,j,n,n2     ! local loop indices

      real (kind=dbl_kind) ::  
     &   gsum, workx, worky           ! tmps for converting grid
     &,  Tsrf (ilo:ihi,jlo:jhi)       ! surface temperature
     &,  tauxa(ilo:ihi,jlo:jhi)       ! atmo/ice stress
     &,  tauya(ilo:ihi,jlo:jhi)               
     &,  tauxo(ilo:ihi,jlo:jhi)       ! ice/ocean stress
     &,  tauyo(ilo:ihi,jlo:jhi)               
     &,  ailohi(ilo:ihi,jlo:jhi)      ! fractional ice area
!YY+
      logical :: flag
      flag=.false.
!YY-

      call ice_timer_start(8)  ! time spent coupling

      do j=jlo,jhi
      do i=ilo,ihi

        ! ice fraction
        ailohi(i,j) = aice(i,j)

        ! surface temperature
        Tsrf(i,j)  = Tffresh + Tsfc(i,j)                    !K

        ! wind stress  (on POP T-grid:  convert to lat-lon)
        workx = strairxT(i,j)                               ! N/m^2
        worky = strairyT(i,j)                               ! N/m^2
        tauxa(i,j) = workx*cos(ANGLET(i,j)) - worky*sin(ANGLET(i,j))
        tauya(i,j) = worky*cos(ANGLET(i,j)) + workx*sin(ANGLET(i,j))

        ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
        workx = -strocnxT(i,j)                              ! N/m^2
        worky = -strocnyT(i,j)                              ! N/m^2
        tauxo(i,j) = workx*cos(ANGLET(i,j)) - worky*sin(ANGLET(i,j))
        tauyo(i,j) = worky*cos(ANGLET(i,j)) + workx*sin(ANGLET(i,j))

      enddo
      enddo

      !--- set info buffer flags ---
      isbuf                          = 0       ! unused
      isbuf(cpl_fields_ibuf_stopnow) = 0       ! stop flag: 0 <=> able to continue
      isbuf(cpl_fields_ibuf_cdate)   = idate   ! model date, coded: yyyymmdd
      isbuf(cpl_fields_ibuf_sec)     = sec     ! elapsed seconds on model date
      isbuf(cpl_fields_ibuf_lisize)  = (ihi-ilo+1) 
      isbuf(cpl_fields_ibuf_ljsize)  = (jhi-jlo+1) 
      isbuf(cpl_fields_ibuf_ncpl)    = nadv_i  ! number of msg-pairs per day
      isbuf(cpl_fields_ibuf_lsize)   = (jhi-jlo+1)*(ihi-ilo+1)
      isbuf(cpl_fields_ibuf_dead)    = 0       ! not a dead model

      call ice_timer_start(18)      ! Time spent packing

      !--- pack & send msg buffer ---
!YY+
      do j=jlo,jhi
      do i=ilo,ihi
         if (tmask(i,j) .and. ailohi(i,j) < c0 ) then
            flag = .true.
         endif
      end do
      end do
      if (flag) then
        do j=jlo,jhi
        do i=ilo,ihi
          if (tmask(i,j) .and. ailohi(i,j) < c0 ) then
            write(nu_diag,*)
     &           ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j)
            call shr_sys_flush(nu_diag)
          endif
        end do
        end do
      endif
!YY-
      buffs(:,:)=spval
      n=0
      do j=jlo,jhi
      do i=ilo,ihi
         n=n+1
!YY         if (tmask(i,j) .and. ailohi(i,j) < c0 ) then
!YY            write(nu_diag,*)
!YY     &           ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j)
!YY            call shr_sys_flush(nu_diag)
!YY         endif

         !--- ice states
         buffs(n,cpl_fields_i2c_ifrac) = ailohi(i,j)    ! frac 

         !---  compression index
         buffs(n,cpl_fields_i2c_index) = rndex_global(i,j)  ! global index
         if (tmask(i,j) .and. ailohi(i,j) > c0 ) then
            buffs(n,cpl_fields_i2c_t    ) = Tsrf(i,j)  ! temperature

            buffs(n,cpl_fields_i2c_avsdr) = alvdr(i,j) ! alb, visible/dir
            buffs(n,cpl_fields_i2c_anidr) = alidr(i,j) ! alb, near-ir/dir
            buffs(n,cpl_fields_i2c_avsdf) = alvdf(i,j) ! alb, visible/dif
            buffs(n,cpl_fields_i2c_anidf) = alidf(i,j) ! alb, near-ir/dif
           !--- a/i fluxes computed by ice
            buffs(n,cpl_fields_i2c_taux ) = tauxa(i,j) ! stress: a/i zonal
            buffs(n,cpl_fields_i2c_tauy ) = tauya(i,j) ! stress: a/i meridional
            buffs(n,cpl_fields_i2c_lat  ) = flat(i,j)  ! latent heat flux
            buffs(n,cpl_fields_i2c_sen  ) = fsens(i,j) ! sensible heat flux
            buffs(n,cpl_fields_i2c_lwup ) = flwout(i,j)! upward longwave flux
            buffs(n,cpl_fields_i2c_evap ) = evap(i,j)  ! evaporation h2o flux
            buffs(n,cpl_fields_i2c_tref ) = Tref(i,j)  ! diagnostic: 2m ref temp
            buffs(n,cpl_fields_i2c_qref ) = Qref(i,j)  ! diagnostic: 2m ref sp hum
            buffs(n,cpl_fields_i2c_swnet) = fswabs(i,j)! sw net absorbed hf
           !--- i/o fluxes computed by ice
            buffs(n,cpl_fields_i2c_swpen) = fswthru(i,j) ! solar thru ice to ocn hf
            buffs(n,cpl_fields_i2c_melth) = fhnet(i,j) ! hf from melting
            buffs(n,cpl_fields_i2c_meltw) = fresh(i,j) ! h2o flux from melting
            buffs(n,cpl_fields_i2c_salt ) = fsalt(i,j) ! salt flux from melting
            buffs(n,cpl_fields_i2c_otaux) = tauxo(i,j) ! stress : i/o zonal
            buffs(n,cpl_fields_i2c_otauy) = tauyo(i,j) ! stress : i/o meridional
         endif  ! tmask and ailohi > c0
      end do
      end do

      call ice_timer_stop(18)      ! Time spent packing

      call ice_timer_start(19)     ! Time spent sending
      call cpl_interface_contractSend
     &     (cpl_fields_cplname, contractS, isbuf, buffs)
      call ice_timer_stop(19)      ! Time spent sending

      !-----------------------------------------------------------------
      ! diagnostics
      !-----------------------------------------------------------------

      if (info_dbug == 1 .AND. stop_now /= 1) then

         do j=jlo,jhi
         do i=ilo,ihi
            worka(i,j) = tarea(i,j)
         enddo
         enddo

         do n=1,cpl_fields_i2c_total
            work_l1 = c0
            n2 = 0
            do j=jlo,jhi
            do i=ilo,ihi
               n2 = n2 + 1
               if(ailohi(i,j) > c0 .and. ailohi(i,j) <= c1) then
                  work_l1(i,j) = buffs(n2,n)
               endif
            enddo
            enddo
            call bound(work_l1)
            call get_sum(0, worka, one, work_l1, gsum)
            if (my_task == master_task) then
               write (nu_diag,100) 'ice','send',n,gsum
            endif
         enddo
      endif
 100  format('comm_diag',1x,a3,1x,a4,1x,i3,es26.19)
      call shr_sys_flush(nu_diag)

      call ice_timer_stop(8)    ! time spent coupling

      end subroutine to_coupler

!=======================================================================
!BOP
!
! !IROUTINE: exit_coupler - exit from coupled/mpi environment
!
! !INTERFACE:
!
      subroutine exit_coupler
!
! !DESCRIPTION:
!
! Exit from coupled/MPI environment
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      include "mpif.h"         ! MPI library definitions

      integer (kind=int_kind) ::
     &  ierr  ! error flag

      if (my_task == master_task) then
         if (irbuf(cpl_fields_ibuf_stopnow) == 1) then
            write (nu_diag,*) '(ice) received final coupler msg'
         else
            write (nu_diag,*) '(ice) terminating before coupler'
            call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
         endif
      endif

      call cpl_interface_finalize (cpl_fields_icename)

      if (my_task == master_task) then
         write(nu_diag,*) '(ice) exit_coupler finished',my_task
      endif

      end subroutine exit_coupler

!=======================================================================

#endif                          ! coupled

      end module ice_coupling

!=======================================================================
