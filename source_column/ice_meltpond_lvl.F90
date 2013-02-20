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
!
!EOP
!
      implicit none

      private
      public :: init_meltponds_lvl, compute_ponds_lvl

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
      subroutine init_meltponds_lvl(nx_block, ny_block, ncat, &
                                    apnd, hpnd, ipnd, dhsn)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
        integer(kind=int_kind), intent(in) :: &
             nx_block , &
             ny_block , &
             ncat

        real(kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
             intent(out) :: &
             apnd , & ! melt pond area fraction
             hpnd , & ! melt pond depth
             ipnd , & ! melt pond refrozen lid thickness
             dhsn     ! depth difference for snow on sea ice and pond ice
!
!EOP
!
      apnd(:,:,:) = c0
      hpnd(:,:,:) = c0
      ipnd(:,:,:) = c0
      dhsn(:,:,:) = c0

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
                                   dt,                  &
                                   dpscale, frzpnd,     &
                                   pndaspect,           &
                                   rfrac, meltt, melts, &
                                   frain, Tair,  fsurfn,&
                                   dhs,   ffrac,        &
                                   aicen, vicen, vsnon, &
                                   qicen, sicen,        &
                                   Tsfcn, alvl, &
                                   apnd,  hpnd, ipnd)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
      use ice_domain_size, only: max_ntrcr, nilyr
      use ice_itd, only: hi_min
      use ice_therm_shared, only: ktherm

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), intent(in) :: &
         dt,       & ! time step (s)  
         dpscale,  & ! alter e-folding time scale for flushing 
         pndaspect   ! ratio of pond depth to pond fraction

      character (len=char_len), intent(in) :: &
         frzpnd       ! pond refreezing parameterization

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         Tsfcn, &
         alvl,  &
         rfrac, &    ! water fraction retained for melt ponds
         meltt, &
         melts, &
         frain, &
         Tair,  &    ! air temperature (K)
         fsurfn,&    ! atm-ice surface heat flux  (W/m2)
         aicen, &
         vicen, &
         vsnon

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         apnd, hpnd, ipnd

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), intent(in) :: &
         qicen, &  ! ice layer enthalpy (J m-3)
         sicen     ! salinity (ppt)   

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         dhs       ! depth difference for snow on sea ice and pond ice

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(out) :: &
         ffrac     ! fraction of fsurfn over pond used to melt ipond

!     local temporary variables

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         volpn 

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
         alvl_tmp               , & ! level ice fraction of ice area
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
         volpn(i,j) = hpnd(i,j) &
                    * aicen(i,j) * alvl(i,j) * apnd(i,j)
         ffrac(i,j) = c0
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Identify grid cells where ponds can be
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (aicen(i,j)*alvl(i,j) > puny**2) then
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
         alvl_tmp = alvl(i,j)

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
            apondn = apnd(i,j)*alvl_tmp

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
               hlid = ipnd(i,j)
               if (dvn == c0) then ! freeze pond
                  Ts = Tair(i,j) - Tffresh
                  if (Ts < c0) then
!                  if (Ts < -c2) then ! as in meltpond_cesm
                     bdt = -c2*Ts*kice*dt/(rhoi*Lfresh)
                     dhlid = p5*sqrt(bdt)                  ! open water freezing
                     if (hlid > dhlid) dhlid = p5*bdt/hlid ! existing ice
                     dhlid = min(dhlid, hpnd(i,j)*rhofresh/rhoi)
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
               apondn = max(c0, min(alvl_tmp, &
                            apondn + 0.5*dvn/(pndaspect*apondn*aicen(i,j))))
               hpondn = c0
               if (apondn > puny) &
                  hpondn = volpn(i,j)/(apondn*aicen(i,j))

            elseif (alvl_tmp*aicen(i,j) > c10*puny) then ! new ponds
               apondn = min (sqrt(volpn(i,j)/(pndaspect*aicen(i,j))), alvl_tmp)
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
               Tmlt(:) = -sicen(i,j,:) * depressT
               call brine_permeability(qicen(i,j,:), vicen(i,j), sicen(i,j,:), Tmlt, perm)
               drain = perm*pressure_head*dt / (viscosity*hi) * dpscale
               deltah = min(drain, hpondn)
               dvn = -deltah*apondn
               volpn(i,j) = volpn(i,j) + dvn
               apondn = max(c0, min(apondn &
                          + 0.5*dvn/(pndaspect*apondn), alvl_tmp*aicen(i,j)))
               hpondn = c0
               if (apondn > puny) hpondn = volpn(i,j)/apondn
            endif

         endif

         !-----------------------------------------------------------
         ! Reload tracer array
         !-----------------------------------------------------------

         hpnd(i,j) = hpondn
         apnd(i,j) = apondn / (aicen(i,j)*alvl_tmp)
         if (trim(frzpnd) == 'hlid') ipnd(i,j) = hlid

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
      use ice_domain_size, only: nilyr
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

      end module ice_meltpond_lvl

!=======================================================================
