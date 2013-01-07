!=======================================================================
!
!BOP
!
! !MODULE: ice_meltpond_lvl - Meltpond parameterization
!
! !DESCRIPTION:
!
! This meltpond parameterization was developed for use with the delta-
! Eddington radiation scheme, and only affects the radiation budget in
! the model.  That is, although the pond volume is tracked, that liquid
! water is not used elsewhere in the model for mass budgets or other
! physical processes.
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors Elizabeth Hunke (LANL)
!         David Hebert (NRL Stennis)
!         Olivier Lecomte (Univ. Louvain)
!
! !INTERFACE:
!
      module ice_meltpond_lvl
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_fileunits
      use ice_read_write
      use ice_restart, only: lenstr, restart_dir, restart_file, &
                             pointer_file, runtype
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
!
!EOP
!
      implicit none

      logical (kind=log_kind) :: & 
         restart_pond_lvl, & ! if .true., read meltponds restart file
         snowinfil           ! if .true., adjust snow depth/area in dEdd
                             !            for infiltration of melt water

      character (len=char_len) :: &
         frzpnd              ! pond refreezing parameterization

      real (kind=dbl_kind) :: &
         dpscale, &          ! alter e-folding time scale for flushing 
         rfracmin, &         ! minimum retained fraction of meltwater
         rfracmax, &         ! maximum retained fraction of meltwater
         pndaspect, &        ! ratio of pond depth to pond fraction
         hs1                 ! tapering parameter for snow on pond ice

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks) :: &
         dhsn, &      ! depth difference for snow on sea ice and pond ice
         ffracn       ! fraction of fsurfn used to melt ipond

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_meltponds
!
! !DESCRIPTION:
!
!  Initialize melt ponds.
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_meltponds_lvl
!
! !USES:
!
      use ice_domain_size
      use ice_blocks
      use ice_domain
      use ice_flux
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      if (trim(runtype) == 'continue' .or. trim(runtype) == 'bering') &
           restart_pond_lvl = .true.

      if (restart_pond_lvl) then
         call read_restart_pond_lvl
      else
         trcrn(:,:,nt_apnd,:,:) = c0
         trcrn(:,:,nt_hpnd,:,:) = c0
         trcrn(:,:,nt_ipnd,:,:) = c0
         dhsn (:,:,:,:) = c0
      endif

      end subroutine init_meltponds_lvl

!=======================================================================
!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
!
      subroutine compute_ponds_lvl(nx_block,ny_block,   &
                                   ilo, ihi, jlo, jhi,  &
                                   rfrac, meltt, melts, &
                                   frain, Tair,  fsurfn,&
                                   dhs,   ffrac,        &
                                   aicen, vicen, vsnon, &
                                   trcrn)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_state, only: nt_Tsfc, nt_apnd, nt_hpnd, nt_ipnd, nt_alvl, nt_qice, nt_sice
      use ice_calendar, only: dt, istep
      use ice_domain_size, only: max_ntrcr
      use ice_itd, only: hi_min
      use ice_therm_shared, only: ktherm

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         rfrac, &    ! water fraction retained for melt ponds
         meltt, &
         melts, &
         frain, &
         Tair,  &    ! air temperature (K)
         fsurfn,&    ! atm-ice surface heat flux  (W/m2)
         aicen, &
         vicen, &
         vsnon

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_ntrcr), &
         intent(inout) :: &
         trcrn

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(out) :: &
         dhs, &    ! depth difference for snow on sea ice and pond ice
         ffrac     ! fraction of fsurfn over pond used to melt ipond

!     local temporary variables

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         volpn, &
         Tsfcn

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr) :: &
         qicen, &  ! ice layer enthalpy (J m-3)
         sicen     ! salinity (ppt)   

      real (kind=dbl_kind), dimension (nilyr) :: &
         Tmlt      ! melting temperature 

      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
         indxi, indxj     ! compressed indices for cells with ice melting

      integer (kind=int_kind) :: i,j,ij,icells,k

      real (kind=dbl_kind) :: &
         hi                     , & ! ice thickness (m)
         hs                     , & ! snow depth (m)
         dTs                    , & ! surface temperature diff for freeze-up (C)
         Tp                     , & ! pond freezing temperature (C)
         Ts                     , & ! surface air temperature (C)
         asnow                  , & ! area fraction of snow on ice
         apondn, &
         hpondn, &
         dvn   , &
         hlid, alid             , & ! refrozen lid thickness, area
         dhlid                  , & ! change in refrozen lid thickness
         bdt                    , & ! 2 kice dT dt / (rhoi Lfresh)
         alvl                   , & ! level ice fraction of ice area
         draft, deltah, pressure_head, perm, drain ! for permeability

      real (kind=dbl_kind), parameter :: &
         Td       = c2          , & ! temperature difference for freeze-up (C)
         rexp     = p01         , & ! pond contraction scaling
         viscosity= 1.79e-3_dbl_kind ! dynamic viscosity, kg/m/s

      !-----------------------------------------------------------------
      ! Initialize 
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         Tsfcn(i,j) = trcrn(i,j,nt_Tsfc)
         volpn(i,j) = trcrn(i,j,nt_hpnd) &
                    * aicen(i,j) * trcrn(i,j,nt_alvl) * trcrn(i,j,nt_apnd)
         ffrac(i,j) = c0
      enddo
      enddo
      do k = 1, nilyr
      do j = 1, ny_block
      do i = 1, nx_block
         qicen(i,j,k) = trcrn(i,j,nt_qice+k-1) 
         sicen(i,j,k) = trcrn(i,j,nt_sice+k-1) 
      enddo
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Identify grid cells where ponds can be
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (aicen(i,j)*trcrn(i,j,nt_alvl) > puny**2) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo                     ! i
      enddo                     ! j

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hi = vicen(i,j)/aicen(i,j)
         hs = vsnon(i,j)/aicen(i,j)
         alvl = trcrn(i,j,nt_alvl)

         if (hi < hi_min) then

         !--------------------------------------------------------------
         ! Remove ponds on thin ice
         !--------------------------------------------------------------
            apondn = c0
            hpondn = c0
            volpn (i,j) = c0
            hlid = c0

         else

            !-----------------------------------------------------------
            ! initialize pond area as fraction of ice
            !-----------------------------------------------------------
            apondn = trcrn(i,j,nt_apnd)*alvl

            !-----------------------------------------------------------
            ! update pond volume
            !-----------------------------------------------------------
            ! add melt water
            dvn = rfrac(i,j)/rhofresh*(meltt(i,j)*rhoi &
                +                      melts(i,j)*rhos &
                +                      frain(i,j)*  dt)*aicen(i,j)

            ! shrink pond volume under freezing conditions
            if (trim(frzpnd) == 'cesm') then
               Tp = Timelt - Td
               dTs = max(Tp - Tsfcn(i,j),c0)
               dvn = dvn - volpn(i,j) * (c1 - exp(rexp*dTs/Tp))

            elseif (trim(frzpnd) == 'hlid') then ! Stefan approximation
               ! assumes pond is fresh (freezing temperature = 0 C)
               ! and ice grows from existing pond ice
               hlid = trcrn(i,j,nt_ipnd)
               if (dvn == c0) then ! freeze pond
                  Ts = Tair(i,j) - Tffresh
                  if (Ts < c0) then
!                  if (Ts < -c2) then ! as in meltpond_cesm
                     bdt = -c2*Ts*kice*dt/(rhoi*Lfresh)
                     dhlid = p5*sqrt(bdt)                  ! open water freezing
                     if (hlid > dhlid) dhlid = p5*bdt/hlid ! existing ice
                     dhlid = min(dhlid, trcrn(i,j,nt_hpnd)*rhofresh/rhoi)
                     hlid = hlid + dhlid
                  endif
               else ! convert refrozen pond ice back to water
                  dhlid = max(fsurfn(i,j)*dt / (rhoi*Lfresh), c0) ! > 0
                  dhlid = -min(dhlid, hlid) ! < 0
                  hlid = max(hlid + dhlid, c0)
                  if (hs - dhs(i,j) < puny) then ! pond ice is snow-free
                     ffrac(i,j) = c1 ! fraction of fsurfn over pond used to melt ipond
                     if (fsurfn(i,j) > puny) &
                        ffrac(i,j) = min(-dhlid*rhoi*Lfresh/(dt*fsurfn(i,j)), c1)
                  endif
               endif
               alid = apondn * aicen(i,j)
               dvn = dvn - dhlid*alid*rhoi/rhofresh
            else ! Stefan approximation
                 ! assumes pond is fresh (freezing temperature = 0 C)
                 ! and ice grows from open water
               if (dvn == c0) then
                  Ts = Tair(i,j) - Tffresh
                  if (Ts < c0) then
                     dhlid = p5*sqrt(-c2*Ts*kice*dt/(rhoi*Lfresh))
                     alid = apondn * aicen(i,j)
                     dvn = dvn - dhlid*alid*rhoi/rhofresh
                  endif
               endif
            endif

            volpn(i,j) = volpn(i,j) + dvn

            !-----------------------------------------------------------
            ! update pond area and depth
            !-----------------------------------------------------------
            if (volpn(i,j) <= c0) then
               volpn(i,j) = c0
               apondn = c0
            endif

            if (apondn*aicen(i,j) > puny) then ! existing ponds
               apondn = max(c0, min(alvl, &
                            apondn + 0.5*dvn/(pndaspect*apondn*aicen(i,j))))
               hpondn = c0
               if (apondn > puny) &
                  hpondn = volpn(i,j)/(apondn*aicen(i,j))

            elseif (alvl*aicen(i,j) > c10*puny) then ! new ponds
               apondn = min (sqrt(volpn(i,j)/(pndaspect*aicen(i,j))), alvl)
               hpondn = pndaspect * apondn

            else           ! melt water runs off deformed ice      
               apondn = c0
               hpondn = c0
            endif
            apondn = max(apondn, c0)

            ! limit pond depth to maintain nonnegative freeboard
            hpondn = min(hpondn, ((rhow-rhoi)*hi - rhos*hs)/rhofresh)

            ! fraction of grid cell covered by ponds
            apondn = apondn * aicen(i,j)

            volpn(i,j) = hpondn*apondn
            if (volpn(i,j) <= c0) then
               volpn(i,j) = c0
               apondn = c0
               hpondn = c0
               hlid = c0
            endif

            !-----------------------------------------------------------
            ! drainage due to permeability (flushing)
            ! setting dpscale = 0 turns this off
            ! NOTE this uses the initial salinity and melting T profiles
            !-----------------------------------------------------------

            if (ktherm /= 2 .and. hpondn > c0 .and. dpscale > puny) then
               draft = (rhos*hs + rhoi*hi)/rhow + hpondn
               deltah = hpondn + hi - draft
               pressure_head = gravit * rhow * max(deltah, c0)
               Tmlt(:) = sicen(i,j,:) * depressT
               call brine_permeability(qicen(i,j,:), vicen(i,j), sicen(i,j,:), Tmlt, perm)
               drain = perm*pressure_head*dt / (viscosity*hi) * dpscale
               deltah = min(drain, hpondn)
               dvn = -deltah*apondn
               volpn(i,j) = volpn(i,j) + dvn
               apondn = max(c0, min(apondn &
                          + 0.5*dvn/(pndaspect*apondn), alvl*aicen(i,j)))
               hpondn = c0
               if (apondn > puny) hpondn = volpn(i,j)/apondn
            endif

         endif

         !-----------------------------------------------------------
         ! Reload tracer array
         !-----------------------------------------------------------

         trcrn(i,j,nt_hpnd) = hpondn
         trcrn(i,j,nt_apnd) = apondn / (aicen(i,j)*alvl)
         if (trim(frzpnd) == 'hlid') trcrn(i,j,nt_ipnd) = hlid

      enddo

      end subroutine compute_ponds_lvl

!=======================================================================
!BOP
!
! !ROUTINE: brine_permeability
!
! !DESCRIPTION:
!
! determine the liquid fraction of brine in the ice and the permeability
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine brine_permeability(qicen, vicen, salin, Tmlt, perm)
!
! !USES:
!
      use ice_therm_shared, only: calculate_Tin_from_qin
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qicen, &  ! enthalpy for each ice layer (J m-3)
         salin, &  ! salinity (ppt)   
         Tmlt      ! melting temperature (C)
    
      real (kind=dbl_kind), intent(in) :: &
         vicen     ! ice volume (m)
    
      real (kind=dbl_kind), intent(out) :: &
         perm      ! permeability (m^2)
!
!EOP
!
      real (kind=dbl_kind) ::   &
         Sbr       ! brine salinity

      real (kind=dbl_kind), dimension(nilyr) ::   &
         Tin, &    ! ice temperature (C)
         phi       ! liquid fraction

      integer (kind=int_kind) :: k
    
      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------

      do k = 1,nilyr
         Tin(k) = calculate_Tin_from_qin(qicen(k),Tmlt(k))
      enddo

      !-----------------------------------------------------------------
      ! brine salinity and liquid fraction
      !-----------------------------------------------------------------

      do k = 1,nilyr
         Sbr = c1/(1.e-3_dbl_kind - depressT/Tin(k)) ! Notz thesis eq 3.6
         phi(k) = salin(k)/Sbr ! liquid fraction
         if (phi(k) < 0.05) phi(k) = c0 ! impermeable
      enddo

      !-----------------------------------------------------------------
      ! permeability
      !-----------------------------------------------------------------

      perm = 3.0e-8_dbl_kind * (minval(phi))**3
    
      end subroutine brine_permeability
  
!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_pond - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_pond_lvl(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for restarting
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         iyear = nyr + year_init - 1
         imonth = month
         iday = mday
         
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.lvlpond.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_pond,filename,0)

      if (my_task == master_task) then
        write(nu_dump_pond) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      diag = .true.

      do n = 1, ncat
         call ice_write(nu_dump_pond,0, trcrn(:,:,nt_apnd,n,:),'ruf8',diag)
         call ice_write(nu_dump_pond,0, trcrn(:,:,nt_hpnd,n,:),'ruf8',diag)
         call ice_write(nu_dump_pond,0, trcrn(:,:,nt_ipnd,n,:),'ruf8',diag)
         call ice_write(nu_dump_pond,0,  dhsn(:,:,n,:),        'ruf8',diag)
         call ice_write(nu_dump_pond,0,ffracn(:,:,n,:),        'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_dump_pond)

      end subroutine write_restart_pond_lvl

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_pond - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_pond_lvl(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for a meltpond volume restart
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: &
         filename, filename0, string1, string2

      logical (kind=log_kind) :: &
         diag

      if (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         read(nu_rst_pointer,'(a)') filename0
         filename = trim(filename0)
         close(nu_rst_pointer)

         ! reconstruct path/file
         n = index(filename0,trim(restart_file))
         if (n == 0) call abort_ice('lvlpond restart: filename discrepancy')
         string1 = trim(filename0(1:n-1))
         string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
         write(filename,'(a,a,a,a)') &
            string1(1:lenstr(string1)), &
            restart_file(1:lenstr(restart_file)),'.lvlpond', &
            string2(1:lenstr(string2))
      endif ! master_task

      call ice_open(nu_restart_pond,filename,0)

      if (my_task == master_task) then
        read(nu_restart_pond) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      do n = 1, ncat
         call ice_read(nu_restart_pond,0, trcrn(:,:,nt_apnd,n,:),'ruf8',diag)
         call ice_read(nu_restart_pond,0, trcrn(:,:,nt_hpnd,n,:),'ruf8',diag)
         call ice_read(nu_restart_pond,0, trcrn(:,:,nt_ipnd,n,:),'ruf8',diag)
         call ice_read(nu_restart_pond,0,  dhsn(:,:,n,:),        'ruf8',diag)
         call ice_read(nu_restart_pond,0,ffracn(:,:,n,:),        'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_restart_pond)

      end subroutine read_restart_pond_lvl

!=======================================================================

      end module ice_meltpond_lvl

!=======================================================================
