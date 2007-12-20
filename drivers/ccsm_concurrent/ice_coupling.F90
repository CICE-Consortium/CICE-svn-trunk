!=======================================================================
!
!BOP
!
! !MODULE: ice_coupling - message passing to and from the coupler
!
! !DESCRIPTION:
! (1) Message passing to and from the coupler
! (2) ESMF import/export state routines
! These are mutually exclusive.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! author: Elizabeth C. Hunke, LANL
!         Tony Craig, NCAR, Dec-30-2002, modified for cpl6
!         Philip W Jones, LANL, modified for ESMF import/export
!
! 2004: Block structure added by William Lipscomb
! 2005: Added ESMF import/export state routines
! 2006: Updated the code for newest CCSM cpl6 - D. Bailey
!
! !INTERFACE:
!
      module ice_coupling
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks
      use ice_boundary
      use ice_domain_size
      use ice_domain
      use ice_constants
      use ice_calendar
      use ice_grid
      use ice_state
      use ice_flux
      use ice_timers
      use ice_fileunits
      use ice_communicate, only: my_task, master_task
      use shr_sys_mod, only : shr_sys_flush
      use ice_restart, only : runtype
      use cpl_contract_mod
      use cpl_interface_mod
      use cpl_fields_mod

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: ice_coupling_setup, &
                init_cpl,           &
                to_coupler,         &
                from_coupler,       &
                exit_coupler

! !PUBLIC DATA MEMBERS:

      logical(kind=log_kind), public :: &
       &  l_CCSMcoupled     ! true if coupled to CCSM

!
!EOP
!BOC
!
      integer (kind=int_kind), dimension (cpl_fields_ibuf_total) :: &
         isbuf &
      ,  irbuf

      real (kind=dbl_kind), allocatable :: &
         sbuf(:,:)

      integer (kind=int_kind), save :: &
         nrecv, nsend

      real (kind=dbl_kind), allocatable :: &
         buffs(:,:) &
      ,  buffr(:,:)

      type(cpl_contract) :: &
         contractS &
      ,  contractR

      integer(kind=int_kind), save :: &
         nadv_i &
      ,  info_dbug

      integer(kind=int_kind), save :: index_i2c_Si_t        ! temperature
      integer(kind=int_kind), save :: index_i2c_Si_tref     ! 2m reference temperature
      integer(kind=int_kind), save :: index_i2c_Si_qref     ! 2m reference specific humidity
      integer(kind=int_kind), save :: index_i2c_Si_ifrac    ! fractional ice coverage
      integer(kind=int_kind), save :: index_i2c_Si_avsdr    ! albedo: visible, direct
      integer(kind=int_kind), save :: index_i2c_Si_anidr    ! albedo: near ir, direct
      integer(kind=int_kind), save :: index_i2c_Si_avsdf    ! albedo: visible, diffuse
      integer(kind=int_kind), save :: index_i2c_Si_anidf    ! albedo: near ir, diffuse
      integer(kind=int_kind), save :: index_i2c_index       ! global data compr index
      integer(kind=int_kind), save :: index_i2c_Faii_taux   ! wind stress, zonal
      integer(kind=int_kind), save :: index_i2c_Faii_tauy   ! wind stress, meridional
      integer(kind=int_kind), save :: index_i2c_Faii_lat    ! latent          heat flux
      integer(kind=int_kind), save :: index_i2c_Faii_sen    ! sensible        heat flux
      integer(kind=int_kind), save :: index_i2c_Faii_lwup   ! upward longwave heat flux
      integer(kind=int_kind), save :: index_i2c_Faii_evap   ! evaporation    water flux
      integer(kind=int_kind), save :: index_i2c_Faii_swnet  ! shortwave: net absorbed
      integer(kind=int_kind), save :: index_i2c_Fioi_swpen  ! net SW penetrating ice
      integer(kind=int_kind), save :: index_i2c_Fioi_melth  ! heat  flux from melting ice
      integer(kind=int_kind), save :: index_i2c_Fioi_meltw  ! water flux from melting ice
      integer(kind=int_kind), save :: index_i2c_Fioi_salt   ! salt  flux from melting ice
      integer(kind=int_kind), save :: index_i2c_Fioi_taux   ! ice/ocn stress, zonal
      integer(kind=int_kind), save :: index_i2c_Fioi_tauy   ! ice/ocn stress, meridional
      
      integer(kind=int_kind), save :: index_c2i_So_t        ! ocn temp
      integer(kind=int_kind), save :: index_c2i_So_s        ! ocn salinity
      integer(kind=int_kind), save :: index_c2i_So_u        ! ocn u velocity
      integer(kind=int_kind), save :: index_c2i_So_v        ! ocn v velocity
      integer(kind=int_kind), save :: index_c2i_So_dhdx     ! ocn surface slope, zonal
      integer(kind=int_kind), save :: index_c2i_So_dhdy     ! ocn surface slope, merid
      integer(kind=int_kind), save :: index_c2i_Sa_z        ! atm bottom layer height
      integer(kind=int_kind), save :: index_c2i_Sa_u        ! atm u velocity
      integer(kind=int_kind), save :: index_c2i_Sa_v        ! atm v velocity
      integer(kind=int_kind), save :: index_c2i_Sa_ptem     ! atm potential temp
      integer(kind=int_kind), save :: index_c2i_Sa_tbot     ! atm bottom temp
      integer(kind=int_kind), save :: index_c2i_Sa_shum     ! atm specfic humidity
      integer(kind=int_kind), save :: index_c2i_Sa_dens     ! atm air density
      integer(kind=int_kind), save :: index_c2i_Fioo_q      ! ocn freeze or melt heat
      integer(kind=int_kind), save :: index_c2i_Faxa_swndr  ! atm sw near-ir, direct
      integer(kind=int_kind), save :: index_c2i_Faxa_swvdr  ! atm sw visable, direct
      integer(kind=int_kind), save :: index_c2i_Faxa_swndf  ! atm sw near-ir, diffuse
      integer(kind=int_kind), save :: index_c2i_Faxa_swvdf  ! atm sw visable, diffuse
      integer(kind=int_kind), save :: index_c2i_Faxa_lwdn   ! long-wave down
      integer(kind=int_kind), save :: index_c2i_Faxc_rain   ! rain
      integer(kind=int_kind), save :: index_c2i_Faxc_snow   ! snow

      type (block) :: this_block
!
! EOC
!
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
      write(nu_diag,*) 'calling cpl_interface_init for model: ', &
           in_model_name,' ', trim(cpl_fields_icename)

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
! !USES:
!
      use ice_broadcast

! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer(kind=int_kind) :: i,j,n,iblk      ! local loop indices
      integer(kind=int_kind) :: ilo,ihi,jlo,jhi ! local loop indices
      integer(kind=int_kind) :: iyear     ! yyyy

      write(nu_diag,*) &
           '(ice_coupling,init_cpl) send initial msg. set contract'
      call shr_sys_flush(nu_diag)

      nadv_i = nint(secday/dt)

      isbuf                          = 0         ! default info-buffer value
      isbuf(cpl_fields_ibuf_cdate  ) = idate     ! initial date (coded: yyyymmdd)
      isbuf(cpl_fields_ibuf_sec    ) = sec       ! elapsed seconds into date
      isbuf(cpl_fields_ibuf_stopnow) = stop_now  ! stop now flag
      isbuf(cpl_fields_ibuf_userest) = 0         ! use model restart data initally
      isbuf(cpl_fields_ibuf_ncpl   ) = nadv_i    ! number of comms per day
      isbuf(cpl_fields_ibuf_lsize  ) &
         = (nx_block-2*nghost)*(ny_block-2*nghost)*nblocks
      isbuf(cpl_fields_ibuf_lisize ) = nx_block-2*nghost ! local size wrt i-index
      isbuf(cpl_fields_ibuf_ljsize ) = ny_block-2*nghost ! local size wrt j-index
      isbuf(cpl_fields_ibuf_gsize  ) = nx_global*ny_global ! size of global grid
      isbuf(cpl_fields_ibuf_gisize ) = nx_global  ! global size wrt i-index
      isbuf(cpl_fields_ibuf_gjsize ) = ny_global  ! global size wrt j-index
      isbuf(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
      isbuf(cpl_fields_ibuf_dead   ) = 0           ! not a dead model

      allocate(sbuf((nx_block-2*nghost)*(ny_block-2*nghost)*nblocks, &
         cpl_fields_grid_total))
      sbuf = -888.0_dbl_kind
      n=0
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo,jhi
         do i = ilo,ihi
            n=n+1
            sbuf(n,cpl_fields_grid_lon  ) = TLON(i,j,iblk)*rad_to_deg
            sbuf(n,cpl_fields_grid_lat  ) = TLAT(i,j,iblk)*rad_to_deg
            sbuf(n,cpl_fields_grid_area ) = tarea(i,j,iblk)/(radius*radius)
            sbuf(n,cpl_fields_grid_mask ) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
            sbuf(n,cpl_fields_grid_index) = rndex_global(i,j,iblk)
         enddo
         enddo
      enddo

      call cpl_interface_contractInit &
           (contractS, cpl_fields_icename, cpl_fields_cplname, &
            cpl_fields_i2c_fields, isbuf, sbuf)

      call cpl_interface_contractInit &
           (contractR, cpl_fields_icename, cpl_fields_cplname, &
            cpl_fields_c2i_fields, isbuf, sbuf)

      write(nu_diag,*) '(init_cpl) Initialized contracts with coupler'
      call shr_sys_flush(nu_diag)

      !-----------------------------------------------------------------
      ! Allocate memory for send and receive buffers
      !-----------------------------------------------------------------

      nsend = cpl_interface_contractNumatt(contractS)
      allocate(buffs((nx_block-2*nghost)*(ny_block-2*nghost)*nblocks, &
         nsend))

      nrecv = cpl_interface_contractNumatt(contractR)
      allocate(buffr((nx_block-2*nghost)*(ny_block-2*nghost)*nblocks, &
         nrecv))

      !-----------------------------------------------------------------
      ! Determine send indices
      !-----------------------------------------------------------------

      index_i2c_index      = cpl_interface_contractIndex(contractS,'index')
      index_i2c_Si_ifrac   = cpl_interface_contractIndex(contractS,'Si_ifrac')
      index_i2c_Si_t       = cpl_interface_contractIndex(contractS,'Si_t')
      index_i2c_Si_avsdr   = cpl_interface_contractIndex(contractS,'Si_avsdr')
      index_i2c_Si_anidr   = cpl_interface_contractIndex(contractS,'Si_anidr')
      index_i2c_Si_avsdf   = cpl_interface_contractIndex(contractS,'Si_avsdf')
      index_i2c_Si_anidf   = cpl_interface_contractIndex(contractS,'Si_anidf')
      index_i2c_Si_tref    = cpl_interface_contractIndex(contractS,'Si_tref')
      index_i2c_Si_qref    = cpl_interface_contractIndex(contractS,'Si_qref')
      index_i2c_Faii_taux  = cpl_interface_contractIndex(contractS,'Faii_taux')
      index_i2c_Faii_tauy  = cpl_interface_contractIndex(contractS,'Faii_tauy')
      index_i2c_Faii_lat   = cpl_interface_contractIndex(contractS,'Faii_lat')
      index_i2c_Faii_sen   = cpl_interface_contractIndex(contractS,'Faii_sen')
      index_i2c_Faii_lwup  = cpl_interface_contractIndex(contractS,'Faii_lwup')
      index_i2c_Faii_evap  = cpl_interface_contractIndex(contractS,'Faii_evap')
      index_i2c_Faii_swnet = cpl_interface_contractIndex(contractS,'Faii_swnet')
      index_i2c_Fioi_swpen = cpl_interface_contractIndex(contractS,'Fioi_swpen')
      index_i2c_Fioi_melth = cpl_interface_contractIndex(contractS,'Fioi_melth')
      index_i2c_Fioi_meltw = cpl_interface_contractIndex(contractS,'Fioi_meltw')
      index_i2c_Fioi_salt  = cpl_interface_contractIndex(contractS,'Fioi_salt')
      index_i2c_Fioi_taux  = cpl_interface_contractIndex(contractS,'Fioi_taux')
      index_i2c_Fioi_tauy  = cpl_interface_contractIndex(contractS,'Fioi_tauy')

      !-----------------------------------------------------------------
      ! Determine receive indices
      !-----------------------------------------------------------------

      index_c2i_So_t       = cpl_interface_contractIndex(contractR,'So_t')
      index_c2i_So_s       = cpl_interface_contractIndex(contractR,'So_s')
      index_c2i_So_u       = cpl_interface_contractIndex(contractR,'So_u')
      index_c2i_So_v       = cpl_interface_contractIndex(contractR,'So_v')
      index_c2i_Sa_z       = cpl_interface_contractIndex(contractR,'Sa_z')
      index_c2i_Sa_u       = cpl_interface_contractIndex(contractR,'Sa_u')
      index_c2i_Sa_v       = cpl_interface_contractIndex(contractR,'Sa_v')
      index_c2i_Sa_ptem    = cpl_interface_contractIndex(contractR,'Sa_ptem')
      index_c2i_Sa_tbot    = cpl_interface_contractIndex(contractR,'Sa_tbot')
      index_c2i_Sa_shum    = cpl_interface_contractIndex(contractR,'Sa_shum')
      index_c2i_Sa_dens    = cpl_interface_contractIndex(contractR,'Sa_dens')
      index_c2i_So_dhdx    = cpl_interface_contractIndex(contractR,'So_dhdx')
      index_c2i_So_dhdy    = cpl_interface_contractIndex(contractR,'So_dhdy')
      index_c2i_Fioo_q     = cpl_interface_contractIndex(contractR,'Fioo_q')
      index_c2i_Faxa_swvdr = cpl_interface_contractIndex(contractR,'Faxa_swvdr')
      index_c2i_Faxa_swndr = cpl_interface_contractIndex(contractR,'Faxa_swndr')
      index_c2i_Faxa_swvdf = cpl_interface_contractIndex(contractR,'Faxa_swvdf')
      index_c2i_Faxa_swndf = cpl_interface_contractIndex(contractR,'Faxa_swndf')
      index_c2i_Faxa_lwdn  = cpl_interface_contractIndex(contractR,'Faxa_lwdn')
      index_c2i_Faxc_rain  = cpl_interface_contractIndex(contractR,'Faxc_rain')
      index_c2i_Faxc_snow  = cpl_interface_contractIndex(contractR,'Faxc_snow')

      !-----------------------------------------------------------------
      ! Receive initial message from coupler.
      !-----------------------------------------------------------------

      call cpl_interface_ibufRecv(cpl_fields_cplname,irbuf)

      if (my_task==master_task) then
         write(nu_diag,*) &
              '(init_cpl) Received control buffer from coupler'
         call shr_sys_flush(nu_diag)

         if (trim(runtype)=='startup' .or. &
             trim(runtype)== 'hybrid') then
            idate = irbuf(cpl_fields_ibuf_cdate)
            write(nu_diag,*) '(init_cpl) idate from coupler = ',idate
            iyear   = (idate/10000)               ! integer year of basedate
            month = (idate-iyear*10000)/100       ! integer month of basedate
            mday  = idate-iyear*10000-month*100-1 ! day of year of basedate
            time  = ((iyear-1)*daycal(13)+daycal(month)+mday)*secday
            call calendar(time)                 ! recompute calendar info
            time_forc = time
            call shr_sys_flush(nu_diag)
         endif

      endif                     ! my_task==master_task

      call broadcast_scalar(idate,    master_task)
      call broadcast_scalar(time,     master_task)
      call broadcast_scalar(time_forc,master_task)

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
! !USES:
!
      use ice_global_reductions
      use ice_work, only: work1

! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i,j,n,n2,iblk   ! local loop indices
      integer(kind=int_kind)  :: ilo,ihi,jlo,jhi ! local loop indices
      integer (kind=int_kind) :: index           ! attribute vector property

      real (kind=dbl_kind) :: &
         gsum, workx, worky

      call ice_timer_start(timer_couple)  ! time spent coupling

      !-----------------------------------------------------------------
      ! Zero stuff while waiting, only filling in active cells.
      !-----------------------------------------------------------------

      zlvl   (:,:,:) = c0
      uatm   (:,:,:) = c0
      vatm   (:,:,:) = c0
      potT   (:,:,:) = c0
      Tair   (:,:,:) = c0
      Qa     (:,:,:) = c0
      rhoa   (:,:,:) = c0
      swvdr  (:,:,:) = c0
      swvdf  (:,:,:) = c0
      swidr  (:,:,:) = c0
      swidf  (:,:,:) = c0
      flw    (:,:,:) = c0
      frain  (:,:,:) = c0
      fsnow  (:,:,:) = c0
      sst    (:,:,:) = c0
      sss    (:,:,:) = c0
      uocn   (:,:,:) = c0
      vocn   (:,:,:) = c0
      ss_tltx(:,:,:) = c0
      ss_tlty(:,:,:) = c0
      frzmlt (:,:,:) = c0

      !-----------------------------------------------------------------
      ! recv input field msg
      !-----------------------------------------------------------------
     
      call ice_timer_start(timer_cplrecv)  ! time spent receiving

      call cpl_interface_contractRecv &
           (cpl_fields_cplname, contractR, irbuf, buffr)

      call ice_timer_stop(timer_cplrecv)
!     call ice_timer_start(17)  ! time spent cr-unpacking

      !--- unpack message
      n=0
      do iblk = 1, nblocks

       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
 
       do j = jlo,jhi
       do i = ilo,ihi

          n=n+1

          !--- ocn states--

          sst(i,j,iblk)  = buffr(n,index_c2i_So_t)
          sss(i,j,iblk)  = buffr(n,index_c2i_So_s)
          uocn(i,j,iblk) = buffr(n,index_c2i_So_u)
          vocn(i,j,iblk) = buffr(n,index_c2i_So_v)

          !--- atm states-

          zlvl (i,j,iblk) = buffr(n,index_c2i_Sa_z)
          uatm (i,j,iblk) = buffr(n,index_c2i_Sa_u)
          vatm (i,j,iblk) = buffr(n,index_c2i_Sa_v)
          potT (i,j,iblk) = buffr(n,index_c2i_Sa_ptem)
          Tair (i,j,iblk) = buffr(n,index_c2i_Sa_tbot)
          Qa   (i,j,iblk) = buffr(n,index_c2i_Sa_shum)
          rhoa (i,j,iblk) = buffr(n,index_c2i_Sa_dens)

          !--- ocn states--

          ss_tltx(i,j,iblk) = buffr(n,index_c2i_So_dhdx)
          ss_tlty(i,j,iblk) = buffr(n,index_c2i_So_dhdy)
          frzmlt (i,j,iblk) = buffr(n,index_c2i_Fioo_q)

          !--- atm fluxes--

          swvdr(i,j,iblk) = buffr(n,index_c2i_Faxa_swvdr)
          swidr(i,j,iblk) = buffr(n,index_c2i_Faxa_swndr)
          swvdf(i,j,iblk) = buffr(n,index_c2i_Faxa_swvdf)
          swidf(i,j,iblk) = buffr(n,index_c2i_Faxa_swndf)
          flw  (i,j,iblk) = buffr(n,index_c2i_Faxa_lwdn)
          frain(i,j,iblk) = buffr(n,index_c2i_Faxc_rain)
          fsnow(i,j,iblk) = buffr(n,index_c2i_Faxc_snow)

       end do
       end do
      end do

      call update_ghost_cells(sst    , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(sss    , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(uocn   , bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(vocn   , bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(zlvl   , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(uatm   , bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(vatm   , bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(potT   , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(Tair   , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(Qa     , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(rhoa   , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(ss_tltx, bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(ss_tlty, bndy_info, field_loc_center, &
                                                  field_type_vector)
      call update_ghost_cells(frzmlt , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(swvdr  , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(swidr  , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(swvdf  , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(swidf  , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(flw    , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(frain  , bndy_info, field_loc_center, &
                                                  field_type_scalar)
      call update_ghost_cells(fsnow  , bndy_info, field_loc_center, &
                                                  field_type_scalar)

!     call ice_timer_stop(17)  ! time spent cr-unpacking

      !-----------------------------------------------------------------
      ! broadcast dbug diagnostic level
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_infobug) >= 2 ) then
         if (my_task == master_task) write (nu_diag,*) &
              '(from_coupler) dbug level >= 2'
         info_dbug = 1
      endif

      !-----------------------------------------------------------------
      ! broadcast write_restart flag
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_resteod) == 1 .AND. new_day) then
         if (my_task == master_task) write (nu_diag,*) &
              '(from_coupler) received write restart signal'
         write_restart = 1
      endif

      !-----------------------------------------------------------------
      ! broadcast cpl_write_history flag
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_histeod) == 1 .AND. new_day) then
         if (my_task == master_task) write (nu_diag,*) &
              '(from_coupler) received write history signal'
         cpl_write_history = 1
      endif

      !-----------------------------------------------------------------
      ! broadcast stop_now flag
      !-----------------------------------------------------------------
      if (irbuf(cpl_fields_ibuf_stopnow) == 1) then
         if (my_task==master_task) write (nu_diag,*) &
              '(from_coupler) received terminate signal'
         stop_now = 1
      endif

      if (info_dbug == 1 .AND. stop_now /= 1) then

        do n=1,nrecv
           work1 = c0
           n2   = 0
           do iblk = 1, nblocks

              this_block = get_block(blocks_ice(iblk),iblk)         
              ilo = this_block%ilo
              ihi = this_block%ihi
              jlo = this_block%jlo
              jhi = this_block%jhi

              do j = jlo,jhi
              do i = ilo,ihi
                 n2 = n2 + 1 
                 if (hm(i,j,iblk) > p5) work1(i,j,iblk) = buffr(n2,n)
              enddo
              enddo
           enddo

           call update_ghost_cells (work1,            bndy_info, &
                                    field_loc_center, field_type_scalar)

           gsum = global_sum(work1, distrb_info, field_loc_center, &
                             tarea)

           if (my_task == master_task) then
              write (nu_diag,100) 'ice', 'recv', n, gsum
           endif
        enddo                   ! nrecv
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

      do iblk = 1, nblocks

       do j = 1,ny_block
       do i = 1,nx_block

          ! ocean
          workx      = uocn  (i,j,iblk) ! currents, m/s 
          worky      = vocn  (i,j,iblk)
          uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                         + worky*sin(ANGLET(i,j,iblk))
          vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                         - workx*sin(ANGLET(i,j,iblk))

          workx      = ss_tltx  (i,j,iblk)           ! sea sfc tilt, m/m
          worky      = ss_tlty  (i,j,iblk)
          ss_tltx(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid 
                            + worky*sin(ANGLET(i,j,iblk))
          ss_tlty(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                            - workx*sin(ANGLET(i,j,iblk))

          sst(i,j,iblk) = sst(i,j,iblk) - Tffresh       ! sea sfc temp (C)
          Tf (i,j,iblk) = -1.8_dbl_kind                 ! hardwired for NCOM
!         Tf (i,j,iblk) = -depressT*sss(i,j,iblk)       ! freezing temp (C)
!         Tf (i,j,iblk) = -depressT*max(sss(i,j,iblk),ice_ref_salinity)

       enddo
       enddo
      enddo

      ! Interpolate ocean dynamics variables from T-cell centers to 
      ! U-cell centers.

      call t2ugrid_vector(uocn)
      call t2ugrid_vector(vocn)
      call t2ugrid_vector(ss_tltx)
      call t2ugrid_vector(ss_tlty)

      ! Atmosphere variables are needed in T cell centers in
      ! subroutine stability and are interpolated to the U grid
      ! later as necessary.

      do iblk = 1, nblocks
       do j = 1, ny_block
       do i = 1, nx_block
         ! atmosphere
         workx      = uatm(i,j,iblk) ! wind velocity, m/s
         worky      = vatm(i,j,iblk) 
         uatm (i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & ! convert to POP grid
                         + worky*sin(ANGLET(i,j,iblk)) ! note uatm, vatm, wind
         vatm (i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & !  are on the T-grid here
                         - workx*sin(ANGLET(i,j,iblk))

         wind (i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
         fsw  (i,j,iblk) = swvdr(i,j,iblk) + swvdf(i,j,iblk) &
                         + swidr(i,j,iblk) + swidf(i,j,iblk)
       enddo
       enddo
      enddo

      time_forc=time

      call ice_timer_stop(timer_couple)   ! time spent coupling

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
      use ice_global_reductions
      use ice_work, only: work1

! !INPUT/OUTPUT PARAMETERS:
!
!EOP

      integer(kind=int_kind) :: i,j,n,n2,iblk    ! local loop indices
      integer(kind=int_kind) :: ilo,ihi,jlo,jhi ! local loop indices
      integer(kind=int_kind) :: index        ! attribute vector property

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         Tsrf  &      ! surface temperature
      ,  tauxa &     ! atmo/ice stress
      ,  tauya &
      ,  tauxo &       ! ice/ocean stress
      ,  tauyo &
      ,  ailohi        ! fractional ice area

      real (kind=dbl_kind) :: &
         gsum, workx, worky           ! tmps for converting grid

      logical :: flag
      flag=.false.

      call ice_timer_start(timer_couple)  ! time spent coupling

      do iblk = 1, nblocks
       do j = 1, ny_block
       do i = 1, nx_block

        ! ice fraction
        ailohi(i,j,iblk) = min(aice(i,j,iblk), c1)

        ! surface temperature
        Tsrf(i,j,iblk)  = Tffresh + trcr(i,j,1,iblk)             !Kelvin

        ! wind stress  (on POP T-grid:  convert to lat-lon)
        workx = strairxT(i,j,iblk)                             ! N/m^2
        worky = strairyT(i,j,iblk)                             ! N/m^2
        tauxa(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                        - worky*sin(ANGLET(i,j,iblk))
        tauya(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                        + workx*sin(ANGLET(i,j,iblk))

        ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
        workx = -strocnxT(i,j,iblk)                            ! N/m^2
        worky = -strocnyT(i,j,iblk)                            ! N/m^2
        tauxo(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                        - worky*sin(ANGLET(i,j,iblk))
        tauyo(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                        + workx*sin(ANGLET(i,j,iblk))

       enddo
       enddo
      enddo

      !--- set info buffer flags ---
      isbuf                          = 0       ! unused
      isbuf(cpl_fields_ibuf_stopnow) = 0       ! stop flag: 0 <=> able to continue
      isbuf(cpl_fields_ibuf_cdate)   = idate   ! model date, coded: yyyymmdd
      isbuf(cpl_fields_ibuf_sec)     = sec     ! elapsed seconds on model date
      isbuf(cpl_fields_ibuf_lisize)  = nx_block-2*nghost  ! local size wrt i-index
      isbuf(cpl_fields_ibuf_ljsize)  = ny_block-2*nghost  ! local size wrt j-index
      isbuf(cpl_fields_ibuf_ncpl)    = nadv_i  ! number of msg-pairs per day
      isbuf(cpl_fields_ibuf_lsize) &
         = (nx_block-2*nghost)*(ny_block-2*nghost)*nblocks 
      isbuf(cpl_fields_ibuf_dead)    = 0       ! not a dead model

!     call ice_timer_start(18)      ! Time spent packing

      !--- pack & send msg buffer ---

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
               flag = .true.
            endif
         end do
         end do
      end do
      if (flag) then
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) < c0 ) then
                 write(nu_diag,*) &
                  ' (ice) send: ERROR ailohi < 0.0 ',i,j,ailohi(i,j,iblk)
                 call shr_sys_flush(nu_diag)
               endif
            end do
            end do
         end do
      endif

      buffs(:,:)=spval_dbl
      n=0
      do iblk = 1, nblocks

       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       do j = jlo,jhi
       do i = ilo,ihi

         n=n+1

         !--- ice states
         buffs(n,index_i2c_Si_ifrac) = ailohi(i,j,iblk)     ! frac 

         !---  compression index
         buffs(n,index_i2c_index) = rndex_global(i,j,iblk)  ! global index

         if (tmask(i,j,iblk) .and. ailohi(i,j,iblk) > c0 ) then
            buffs(n,index_i2c_Si_t      ) = Tsrf(i,j,iblk)  ! temperature

            buffs(n,index_i2c_Si_avsdr  ) = alvdr(i,j,iblk) ! alb, visible/dir
            buffs(n,index_i2c_Si_anidr  ) = alidr(i,j,iblk) ! alb, near-ir/dir
            buffs(n,index_i2c_Si_avsdf  ) = alvdf(i,j,iblk) ! alb, visible/dif
            buffs(n,index_i2c_Si_anidf  ) = alidf(i,j,iblk) ! alb, near-ir/dif
            buffs(n,index_i2c_Si_tref   ) = Tref(i,j,iblk)  ! diagnostic: 2m ref temp
            buffs(n,index_i2c_Si_qref   ) = Qref(i,j,iblk)  ! diagnostic: 2m ref sp hum
            !--- a/i fluxes computed by ice
            buffs(n,index_i2c_Faii_taux ) = tauxa(i,j,iblk) ! stress: a/i zonal
            buffs(n,index_i2c_Faii_tauy ) = tauya(i,j,iblk) ! stress: a/i meridional
            buffs(n,index_i2c_Faii_lat  ) = flat(i,j,iblk)  ! latent heat flux
            buffs(n,index_i2c_Faii_sen  ) = fsens(i,j,iblk) ! sensible heat flux
            buffs(n,index_i2c_Faii_lwup ) = flwout(i,j,iblk)! upward longwave flux
            buffs(n,index_i2c_Faii_evap ) = evap(i,j,iblk)  ! evaporation h2o flux
            buffs(n,index_i2c_Faii_swnet) = fswabs(i,j,iblk)! sw net absorbed hf
            !--- i/o fluxes computed by ice
            buffs(n,index_i2c_Fioi_swpen) = fswthru(i,j,iblk) ! solar thru ice to ocn hf
            buffs(n,index_i2c_Fioi_melth) = fhocn(i,j,iblk) ! hf from melting
            buffs(n,index_i2c_Fioi_meltw) = fresh(i,j,iblk) ! h2o flux from melting
            buffs(n,index_i2c_Fioi_salt ) = fsalt(i,j,iblk) ! salt flux from melting
            buffs(n,index_i2c_Fioi_taux ) = tauxo(i,j,iblk) ! stress : i/o zonal
            buffs(n,index_i2c_Fioi_tauy ) = tauyo(i,j,iblk) ! stress : i/o meridional
         endif  ! tmask and ailohi > c0
       end do
       end do
      end do

!     call ice_timer_stop(18)      ! Time spent packing

      call ice_timer_start(timer_cplsend)     ! Time spent sending
      call cpl_interface_contractSend &
           (cpl_fields_cplname, contractS, isbuf, buffs)
      call ice_timer_stop(timer_cplsend)      ! Time spent sending

      !-----------------------------------------------------------------
      ! diagnostics
      !-----------------------------------------------------------------

      if (info_dbug == 1 .AND. stop_now /= 1) then

         do n=1,nsend
            work1 = c0
            n2 = 0
            do iblk = 1, nblocks

               this_block = get_block(blocks_ice(iblk),iblk)         
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

               do j = jlo,jhi
               do i = ilo,ihi
                  n2 = n2 + 1
                  if (ailohi(i,j,iblk) > c0 .and. &
                      ailohi(i,j,iblk) <= c1) then
                     work1(i,j,iblk) = buffs(n2,n)
                  endif
               enddo
               enddo
            enddo

            call update_ghost_cells (work1,            bndy_info, &
                                     field_loc_center, field_type_scalar)

            gsum = global_sum(work1, distrb_info, field_loc_center, &
                              tarea)

            if (my_task == master_task) then
               write (nu_diag,100) 'ice','send',n,gsum
            endif
         enddo
      endif
 100  format('comm_diag',1x,a3,1x,a4,1x,i3,es26.19)
      call shr_sys_flush(nu_diag)

      call ice_timer_stop(timer_couple)    ! time spent coupling

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

      integer (kind=int_kind) :: &
        ierr  ! error flag

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

!=======================================================================

      end module ice_coupling

!=======================================================================
