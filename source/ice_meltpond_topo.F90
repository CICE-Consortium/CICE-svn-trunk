!=======================================================================
!
!BOP
!
! !MODULE: ice_meltpond_topo - Meltpond parameterization
!
! !DESCRIPTION:
!
! Melt pond evolution based on the ice topography as inferred from
! the ice thickness distribution.  This code is based on (but differs
! from) that described in
!
! Flocco, D. and D. L. Feltham, 2007.  A continuum model of melt pond 
! evolution on Arctic sea ice.  J. Geophys. Res. 112, C08016, doi: 
! 10.1029/2006JC003836.
!
! Flocco, D., D. L. Feltham and A. K. Turner, 2010.  Incorporation of a
! physically based melt pond scheme into the sea ice component of a
! climate model.  J. Geophys. Res. 115, C08012, doi: 10.1029/2009JC005568.
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors Daniela Flocco (UCL)
!         Adrian Turner (UCL)
! 2010 ECH added module based on original code from Daniela Flocco, UCL
! 2012 DSCHR modifications
!
! !INTERFACE:
!
      module ice_meltpond_topo
!
! !USES:
!
      use ice_calendar, only: dt  
      use ice_communicate, only: my_task, master_task
      use ice_constants
      use ice_domain_size, only: max_ntrcr
      use ice_exit, only : abort_ice     
      use ice_flux
      use ice_fileunits
      use ice_itd, only: ilyr1, slyr1
      use ice_kinds_mod
      use ice_read_write
      use ice_restart, only: lenstr, restart_dir, restart_file, &
                             pointer_file, runtype
      use ice_state
      use ice_therm_vertical
      use ice_blocks
!
!EOP
!
      implicit none
      save

      logical (kind=log_kind) :: & 
         restart_pond_topo ! if .true., read meltponds restart file

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_meltponds_topo
!
! !DESCRIPTION:
!
!  Initialize melt ponds.
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_meltponds_topo
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
      if (trim(runtype) == 'continue') restart_pond_topo = .true.
      if (restart_pond_topo) then
         call read_restart_pond_topo
      else
         trcrn(:,:,nt_apnd,:,:) = c0
         trcrn(:,:,nt_hpnd,:,:) = c0
         trcrn(:,:,nt_ipnd,:,:) = c0
      endif ! .not restart_pond
        
      end subroutine init_meltponds_topo

!=======================================================================
!BOP
!
! !ROUTINE:   compute_ponds_topo
!
! !INTERFACE:
!
      subroutine compute_ponds_topo(nx_block,ny_block, &
                                    ilo, ihi, jlo, jhi,&
                                    aice,    aicen, &
                                    vice,    vicen, &
                                    vsno,    vsnon, &
                                    eicen,   esnon, &
                                    trcrn,          &
                                    potT,    meltt, &
                                    fsurf,   fpond)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aice, &    ! total ice area fraction
         vsno       ! total snow volume (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         vice, &    ! total ice volume (m)
         fpond      ! fresh water flux to ponds (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(in) :: &
         aicen, &   ! ice area fraction, per category
         vsnon      ! snow volume, per category (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         vicen      ! ice volume, per category (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(inout) :: &
         eicen      ! ice enthalpy, per category

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(inout) :: &
         esnon      ! snow enthalpy, per category

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(inout) :: &
         trcrn      ! ice tracer array

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         potT,  &   ! air potential temperature
         meltt, &   ! total surface meltwater flux
         fsurf      ! thermodynamic heat flux at ice/snow surface (W/m^2)
!
!EOP
!
      ! local variables
      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) :: &
         Tsfcn, & ! ice/snow surface temperature (C)
         volpn, & ! pond volume per unit area, per category (m)
         vuin     ! water-equivalent volume of ice lid on melt pond ('upper ice', m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) :: &
         apondn,& ! pond area fraction, per category
         hpondn   ! pond depth, per category (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         volp       ! total volume of pond, per unit area of pond (m)

      real (kind=dbl_kind) :: &
         hi,    & ! ice thickness (m)
         dHui,  & ! change in thickness of ice lid (m)
         omega,	& ! conduction
         dTice, & ! temperature difference across ice lid (C)
         dvice, & ! change in ice volume (m)
         Tavg,  & ! mean surface temperature across categories (C)
         Tp,    & ! pond freezing temperature (C)
         dvn      ! change in melt pond volume for fresh water budget
      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with ice melting

      integer (kind=int_kind) :: n,k,i,j,ij,icells,indxij ! loop indices

      integer (kind=int_kind), dimension (ncat) :: &
         kcells          ! cells where ice lid combines with vice

      real (kind=dbl_kind), dimension (nx_block*ny_block,nilyr) :: &
         qin             ! ice layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (nx_block*ny_block,nilyr+1) :: &
         zi1         , & ! depth of ice layer boundaries (m)
         zi2             ! adjusted depths, with equal hilyr (m)

      real (kind=dbl_kind) :: &
         dzi             ! ice layer thickness after growth/melting

      integer (kind=int_kind), dimension (nx_block*ny_block,ncat) :: &
         indxii, indxjj  ! i,j indices for kcells loop

      real (kind=dbl_kind), dimension (nx_block*ny_block) :: &
         hin         , & ! total ice thickness (m)
         hilyr           ! ice layer thickness (m)

      real (kind=dbl_kind), parameter :: &
         hicemin = p1           , & ! minimum ice thickness with ponds (m) 
         Td      = p15          , & ! temperature difference for freeze-up (C)
         rhoi_L  = Lfresh * rhoi, & ! (J/m^3)
         min_volp = 1.e-4_dbl_kind  ! minimum pond volume (m)

      !---------------------------------------------------------------
      ! initialize
      !---------------------------------------------------------------

      do j = 1, ny_block
         do i = 1, nx_block
            volp(i,j) = c0
         enddo
      enddo
      do n = 1, ncat
         do j = jlo, jhi
            do i = ilo, ihi
               ! load tracers
               volp(i,j) = volp(i,j) + trcrn(i,j,nt_hpnd,n) &
                                     * trcrn(i,j,nt_apnd,n) * aicen(i,j,n)
               Tsfcn(i,j,n) = trcrn (i,j,nt_Tsfc,n) !echmod - could get rid of this cp
               vuin (i,j,n) = trcrn(i,j,nt_ipnd,n) &
                            * trcrn(i,j,nt_apnd,n) * aicen(i,j,n)

               hpondn(i,j,n) = c0     ! pond depth, per category
               apondn(i,j,n) = c0     ! pond area,  per category
            enddo
         enddo
         indxii(:,n) = 0
         indxjj(:,n) = 0
         kcells  (n) = 0
      enddo

      ! The freezing temperature for meltponds is assumed slightly below 0C,
      ! as if meltponds had a little salt in them.  The salt budget is not
      ! altered for meltponds, but if it were then an actual pond freezing 
      ! temperature could be computed.

      Tp = Timelt - Td

      !-----------------------------------------------------------------
      ! Identify grid cells with ponds
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         hi = c0
         if (aice(i,j) > puny) hi = vice(i,j)/aice(i,j)
         if ( aice(i,j) > p01 .and. hi > hicemin .and. &
            volp(i,j) > min_volp*aice(i,j)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         else  ! remove ponds on thin ice
            fpond(i,j) = fpond(i,j) - volp(i,j)
            volpn(i,j,:) = c0
            vuin (i,j,:) = c0
            volp (i,j) = c0
         endif
      enddo                     ! i
      enddo                     ! j

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         !--------------------------------------------------------------
         ! calculate pond area and depth
         !--------------------------------------------------------------
         call pond_area(aice(i,j),vice(i,j),vsno(i,j), &
                   aicen(i,j,:), vicen(i,j,:), vsnon(i,j,:), eicen(i,j,:), &
                   volpn(i,j,:), volp(i,j), &
                   apondn(i,j,:),hpondn(i,j,:), dvn)

         fpond(i,j) = fpond(i,j) - dvn

         ! mean surface temperature
         Tavg = c0
         do n = 1, ncat
            Tavg = Tavg + Tsfcn(i,j,n)*aicen(i,j,n)
         enddo
         Tavg = Tavg / aice(i,j)

         do n = 1, ncat-1
                       
            if (vuin(i,j,n) > puny) then

         !----------------------------------------------------------------
         ! melting: floating upper ice layer melts in whole or part
         !----------------------------------------------------------------
!               if (Tsfcn(i,j,n) > Tp) then
               if (Tavg > Tp) then

                  dvice = min(meltt(i,j)*apondn(i,j,n), vuin(i,j,n))
                  if (dvice > puny) then
                     vuin (i,j,n) = vuin (i,j,n) - dvice
                     volpn(i,j,n) = volpn(i,j,n) + dvice
                     volp (i,j)   = volp (i,j)   + dvice
                     fpond(i,j)   = fpond(i,j)   + dvice
                       
                     if (vuin(i,j,n) < puny .and. volpn(i,j,n) > puny) then
                        ! ice lid melted and category is pond covered
                        volpn(i,j,n) = volpn(i,j,n) + vuin(i,j,n)
                        fpond(i,j)   = fpond(i,j)   + vuin(i,j,n)
                        vuin(i,j,n)  = c0
                     endif
                     hpondn(i,j,n) = volpn(i,j,n) / apondn(i,j,n)
                  endif

         !----------------------------------------------------------------
         ! freezing: existing upper ice layer grows
         !----------------------------------------------------------------
               else if (volpn(i,j,n) > puny) then ! Tavg <= Tp

                ! differential growth of base of surface floating ice layer
                  dTice = max(-Tavg, c0) ! > 0
                  omega = kice*DTice/rhoi_L
                  dHui = sqrt(omega*dt + p25*(vuin(i,j,n)/aicen(i,j,n))**2) &
                                       - p5 * vuin(i,j,n)/aicen(i,j,n)

                  dvice = min(dHui*apondn(i,j,n), volpn(i,j,n))   
                  if (dvice > puny) then
                     vuin (i,j,n) = vuin (i,j,n) + dvice
                     volpn(i,j,n) = volpn(i,j,n) - dvice
                     volp (i,j)   = volp (i,j)   - dvice
                     fpond(i,j)   = fpond(i,j)   - dvice
                     hpondn(i,j,n) = volpn(i,j,n) / apondn(i,j,n)
                  endif

               endif ! Tavg

         !----------------------------------------------------------------
         ! freezing: upper ice layer begins to form
         ! note: albedo does not change
         !----------------------------------------------------------------
            else ! vuin < puny
                    
               ! thickness of newly formed ice
               ! the surface temperature of a meltpond is the same as that
               ! of the ice underneath (0C), and the thermodynamic surface 
               ! flux is the same
               dHui = max(-fsurf(i,j)*dt/rhoi_L, c0)
               dvice = min(dHui*apondn(i,j,n), volpn(i,j,n))  
               if (dvice > puny) then
                  vuin (i,j,n) = dvice
                  volpn(i,j,n) = volpn(i,j,n) - dvice
                  volp (i,j)   = volp (i,j)   - dvice
                  fpond(i,j)   = fpond(i,j)   - dvice
                  hpondn(i,j,n)= volpn(i,j,n) / apondn(i,j,n)
               endif
                    
            endif  ! vuin

         enddo ! ncat

      enddo ! ij

      !---------------------------------------------------------------
      ! remove ice lid if there is no liquid pond
      ! vuin may be nonzero on category ncat due to dynamics
      !---------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         do n = 1, ncat
            if (aicen(i,j,n) > puny .and. volpn(i,j,n) < puny &
                                    .and. vuin (i,j,n) > puny) then
               kcells(n) = kcells(n) + 1
               indxij    = kcells(n)
               indxii(indxij,n) = i
               indxjj(indxij,n) = j
            endif
         enddo
      enddo                     ! i
      enddo                     ! j

      do n = 1, ncat

         if (kcells(n) > 0) then
         do ij = 1, kcells(n)
            i = indxii(ij,n)
            j = indxjj(ij,n)
            vuin(i,j,n) = c0
         enddo    ! ij
         endif

         ! reload tracers
         do j = jlo, jhi
            do i = ilo, ihi
               if (apondn(i,j,n) > puny) then
                  trcrn (i,j,nt_ipnd,n) = vuin(i,j,n) / apondn(i,j,n)
               else
                  vuin(i,j,n) = c0
                  trcrn (i,j,nt_ipnd,n) = c0
               endif
               if (aicen(i,j,n) > puny) then
                  trcrn (i,j,nt_apnd,n) = apondn(i,j,n) / aicen(i,j,n)
                  trcrn (i,j,nt_hpnd,n) = hpondn(i,j,n)
               else
                  trcrn (i,j,nt_apnd,n) = c0
                  trcrn (i,j,nt_hpnd,n) = c0
                  trcrn (i,j,nt_ipnd,n) = c0
               endif
            enddo ! i
         enddo    ! j

      enddo       ! n


 end subroutine compute_ponds_topo

!=======================================================================
!BOP
!
! !ROUTINE: pond_area
!
! !DESCRIPTION:
!
! Computes melt pond area, pond depth and melting rates
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine pond_area(aice,vice,vsno, &
                           aicen, vicen, vsnon, eicen, &
                           volpn, volp,  &
                           apondn,hpondn,dvolp)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
    
      real (kind=dbl_kind), intent(in) :: &
         aice,vice,vsno

      real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         aicen, vicen, vsnon

      real (kind=dbl_kind), dimension(ntilyr), intent(in) :: &
         eicen

      real (kind=dbl_kind), dimension(ncat), intent(inout) :: &
         volpn

      real (kind=dbl_kind), intent(inout) :: &
         volp, dvolp

      real (kind=dbl_kind), dimension(ncat), intent(out) :: &
         apondn, hpondn
!
!EOP
!
      integer (kind=int_kind) :: &
         n, ns,   &
	 m_index, &
         permflag

      real (kind=dbl_kind), dimension(ncat) :: &
         hicen, &
         hsnon, &
         asnon, &
         alfan, &
         betan, &
         cum_max_vol, &
         reduced_aicen        

      real (kind=dbl_kind), dimension(0:ncat) :: &
         cum_max_vol_tmp

      real (kind=dbl_kind) :: &
         hpond, &
         drain, &
         floe_weight, &
         pressure_head, &
         hsl_rel, &
         deltah, &
         perm, &
         apond

      real (kind=dbl_kind), parameter :: & 
         viscosity = 1.79e-3_dbl_kind  ! kinematic water viscosity in kg/m/s

 !-----------|
 !           |
 !           |-----------|
 !___________|___________|______________________________________sea-level
 !           |           |
 !           |           |---^--------|
 !           |           |   |        |
 !           |           |   |        |-----------|              |-------
 !           |           |   |alfan(n)|           |              |
 !           |           |   |        |           |--------------|
 !           |           |   |        |           |              |
 !---------------------------v-------------------------------------------
 !           |           |   ^        |           |              |
 !           |           |   |        |           |--------------|
 !           |           |   |betan(n)|           |              |
 !           |           |   |        |-----------|              |-------
 !           |           |   |        |
 !           |           |---v------- |
 !           |           |
 !           |-----------|
 !           |
 !-----------|
    
      !-------------------------------------------------------------------
      ! initialize
      !-------------------------------------------------------------------

      do n = 1, ncat

         apondn(n) = c0
         hpondn(n) = c0

         if (aicen(n) < puny)  then
            hicen(n) =  c0 
            hsnon(n) = c0
            reduced_aicen(n) = c0
         else
            hicen(n) = vicen(n) / aicen(n)
            hsnon(n) = vsnon(n) / aicen(n)
            reduced_aicen(n) = c1 ! n=ncat
            if (n < ncat) reduced_aicen(n) = aicen(n) &
                                 * (-0.024_dbl_kind*hicen(n) + 0.832_dbl_kind) 
            asnon(n) = reduced_aicen(n) 
         endif

! This choice for alfa and beta ignores hydrostatic equilibium of categories.
! Hydrostatic equilibium of the entire ITD is accounted for below, assuming
! a surface topography implied by alfa=0.6 and beta=0.4, and rigidity across all
! categories.  alfa and beta partition the ITD - they are areas not thicknesses!
! Multiplying by hicen, alfan and betan (below) are thus volumes per unit area.
! Here, alfa = 60% of the ice area (and since hice is constant in a category, 
! alfan = 60% of the ice volume) in each category lies above the reference line, 
! and 40% below. Note: p6 is an arbitrary choice, but alfa+beta=1 is required.

         alfan(n) = p6 * hicen(n)
         betan(n) = p4 * hicen(n)
       
         cum_max_vol(n)     = c0
         cum_max_vol_tmp(n) = c0
    
      enddo ! ncat

      cum_max_vol_tmp(0) = c0
      drain = c0
      dvolp = c0
    
      !--------------------------------------------------------------------------
      ! the maximum amount of water that can be contained up to each ice category
      !--------------------------------------------------------------------------
    
      do n = 1, ncat-1 ! last category can not hold any volume

         if (alfan(n+1) >= alfan(n) .and. alfan(n+1) > c0) then

            ! total volume in level including snow
            cum_max_vol_tmp(n) = cum_max_vol_tmp(n-1) + &
               (alfan(n+1) - alfan(n)) * sum(reduced_aicen(1:n)) 


            ! subtract snow solid volumes from lower categories in current level
            do ns = 1, n 
               cum_max_vol_tmp(n) = cum_max_vol_tmp(n) &
                  - rhos/rhow  * &    ! fraction of snow that is occupied by solid
                    asnon(ns)  * &    ! area of snow from that category
                    max(min(hsnon(ns)+alfan(ns)-alfan(n), alfan(n+1)-alfan(n)), c0)  
                                      ! thickness of snow from ns layer in n layer
            enddo

         else ! assume higher categories unoccupied
            cum_max_vol_tmp(n) = cum_max_vol_tmp(n-1)
         endif
         if (cum_max_vol_tmp(n) < c0) call abort_ice('negative melt pond volume')

      enddo
      cum_max_vol_tmp(ncat) = cum_max_vol_tmp(ncat-1)  ! last category holds no volume
      cum_max_vol  (1:ncat) = cum_max_vol_tmp(1:ncat)
    
      !----------------------------------------------------------------
      ! is there more meltwater than can be held in the floe?
      !----------------------------------------------------------------
      if (volp >= cum_max_vol(ncat)) then
         drain = volp - cum_max_vol(ncat) + puny
         volp = volp - drain
         dvolp = drain
         if (volp < puny) then
            dvolp = dvolp + volp
            volp = c0
         endif
      endif
    
      ! height and area corresponding to the remaining volume

      call calc_hpond(reduced_aicen, asnon, hsnon, alfan, &
           volp, cum_max_vol, hpond, m_index)
    
      do n=1, m_index
         hpondn(n) = hpond - alfan(n) + alfan(1)
         apondn(n) = reduced_aicen(n) 
      enddo
      apond = sum(apondn(1:m_index))
    
      !------------------------------------------------------------------------
      ! drainage due to ice permeability - Darcy's law
      !------------------------------------------------------------------------
    
      ! sea water level 
      floe_weight = (vsno*rhos + rhoi*vice + rhow*volp) / aice
      hsl_rel = floe_weight / rhow &
              - ((sum(betan(:)*aicen(:))/aice) + alfan(1))
    
      deltah = hpond - hsl_rel
      pressure_head = gravit * rhow * max(deltah, c0)

      ! drain if ice is permeable    
      permflag = 0
      if (pressure_head > c0) then
      do n = 1, ncat-1
         if (hicen(n) /= c0) then
            call permeability_phi(eicen(ilyr1(n):ilyr1(n)+nilyr-1),vicen(n),perm)
            if (perm > c0) permflag = 1
            drain = perm*apondn(n)*pressure_head*dt / (viscosity*hicen(n))
            dvolp = dvolp + min(drain, volp)
            volp = max(volp - drain, c0)
            if (volp < puny) then
               dvolp = dvolp + volp
               volp = c0
            endif
         endif
      enddo
 
      ! adjust melt pond dimensions
      if (permflag > 0) then
         ! recompute pond depth    
         call calc_hpond(reduced_aicen, asnon, hsnon, alfan, &
                         volp, cum_max_vol, hpond, m_index)
         do n=1, m_index
            hpondn(n) = hpond - alfan(n) + alfan(1)
            apondn(n) = reduced_aicen(n) 
         enddo
         apond = sum(apondn(1:m_index))
      endif
      endif ! pressure_head

      !------------------------------------------------------------------------
      ! total melt pond volume in category does not include snow volume
      ! snow in melt ponds is not melted
      !------------------------------------------------------------------------

      ! Calculate pond volume for lower categories
      do n=1,m_index-1
         volpn(n) = apondn(n) * hpondn(n) &
                  - (rhos/rhow) * asnon(n) * min(hsnon(n), hpondn(n))
      enddo

      ! Calculate pond volume for highest category = remaining pond volume
      if (m_index == 1) volpn(m_index) = volp
      if (m_index > 1) then
        if (volp > sum(volpn(1:m_index-1))) then
          volpn(m_index) = volp - sum(volpn(1:m_index-1))
        else
          volpn(m_index) = c0
          hpondn(m_index) = c0
          apondn(m_index) = c0
          ! If remaining pond volume is negative reduce pond volume of 
          ! lower category
          if (volp+puny < sum(volpn(1:m_index-1))) & 
            volpn(m_index-1) = volpn(m_index-1) - sum(volpn(1:m_index-1)) + &
                               volp
        endif
      endif

      do n=1,m_index
         if (apondn(n) > puny) then
             hpondn(n) = volpn(n) / apondn(n)
         else
            dvolp = dvolp + volpn(n)
            hpondn(n) = c0
            volpn(n) = c0
            apondn(n) = c0
         end if
      enddo
      do n = m_index+1, ncat
         hpondn(n) = c0
         apondn(n) = c0
         volpn (n) = c0
      enddo

      end subroutine pond_area
  
!=======================================================================
  
  subroutine calc_hpond(aicen, asnon, hsnon, alfan, &
       volp, cum_max_vol, &
       hpond, m_index)
    
    real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         aicen, &
         asnon, &
         hsnon, &
         alfan, &
         cum_max_vol
    
    real (kind=dbl_kind), intent(in) :: &
         volp
    
    real (kind=dbl_kind), intent(out) :: &
         hpond
    
    integer (kind=int_kind), intent(out) :: &
         m_index
    
    integer :: n, ns
    
    real (kind=dbl_kind), dimension(0:ncat+1) :: &
         hitl, &
         aicetl
    
    real (kind=dbl_kind) :: &
         rem_vol, &
         area, &
         vol, &
         tmp
    
    !----------------------------------------------------------------
    ! hpond is zero if volp is zero - have we fully drained? 
    !----------------------------------------------------------------
    
    if (volp < puny) then
       hpond = c0
       m_index = 0
    else
       
       !----------------------------------------------------------------
       ! Calculate the category where water fills up to 
       !----------------------------------------------------------------
       
       !----------|
       !          |
       !          |
       !          |----------|                                     -- --
       !__________|__________|_________________________________________ ^
       !          |          |             rem_vol     ^                | Semi-filled
       !          |          |----------|-- -- -- - ---|-- ---- -- -- --v layer
       !          |          |          |              |
       !          |          |          |              |hpond
       !          |          |          |----------|   |     |-------
       !          |          |          |          |   |     |
       !          |          |          |          |---v-----|
       !          |          | m_index  |          |         |
       !-------------------------------------------------------------
       
       m_index = 0  ! 1:m_index categories have water in them
       do n = 1, ncat
          if (volp <= cum_max_vol(n)) then
             m_index = n
             if (n == 1) then
                rem_vol = volp
             else
                rem_vol = volp - cum_max_vol(n-1)
             endif
             exit ! to break out of the loop
          endif
       enddo
       m_index = min(ncat-1, m_index)
       
       !----------------------------------------------------------------
       ! semi-filled layer may have m_index different snows in it
       !----------------------------------------------------------------
       
       !-----------------------------------------------------------  ^
       !                                                             |  alfan(m_index+1)
       !                                                             |
       !hitl(3)-->                             |----------|          |
       !hitl(2)-->                |------------| * * * * *|          |
       !hitl(1)-->     |----------|* * * * * * |* * * * * |          |
       !hitl(0)-->-------------------------------------------------  |  ^
       !                various snows from lower categories          |  |alfa(m_index)
       
       ! hitl - heights of the snow layers from thinner and current categories
       ! aicetl - area of each snow depth in this layer
       
       hitl(:) = c0
       aicetl(:) = c0
       do n = 1, m_index
          hitl(n)   = max(min(hsnon(n) + alfan(n) - alfan(m_index), &
                                 alfan(m_index+1) - alfan(m_index)), c0)
          aicetl(n) = asnon(n)
          
          aicetl(0) = aicetl(0) + (aicen(n) - asnon(n))
       enddo
       hitl(m_index+1) = alfan(m_index+1) - alfan(m_index)
       aicetl(m_index+1) = c0
       
       !----------------------------------------------------------------
       ! reorder array according to hitl 
       ! snow heights not necessarily in height order
       !----------------------------------------------------------------
       
       do ns = 1, m_index+1
          do n = 0, m_index - ns + 1
             if (hitl(n) > hitl(n+1)) then ! swap order
                tmp = hitl(n)
                hitl(n) = hitl(n+1)
                hitl(n+1) = tmp
                tmp = aicetl(n)
                aicetl(n) = aicetl(n+1)
                aicetl(n+1) = tmp
             endif
          enddo
       enddo
       
       !----------------------------------------------------------------
       ! divide semi-filled layer into set of sublayers each vertically homogenous
       !----------------------------------------------------------------
       
       !hitl(3)----------------------------------------------------------------
       !                                                       | * * * * * * * *  
       !                                                       |* * * * * * * * * 
       !hitl(2)----------------------------------------------------------------
       !                                    | * * * * * * * *  | * * * * * * * *  
       !                                    |* * * * * * * * * |* * * * * * * * * 
       !hitl(1)----------------------------------------------------------------
       !                 | * * * * * * * *  | * * * * * * * *  | * * * * * * * *  
       !                 |* * * * * * * * * |* * * * * * * * * |* * * * * * * * * 
       !hitl(0)----------------------------------------------------------------
       !    aicetl(0)         aicetl(1)           aicetl(2)          aicetl(3)            
       
       ! move up over layers incrementing volume
       do n = 1, m_index+1
          
          area = sum(aicetl(:)) - &                 ! total area of sub-layer
               (rhos/rhow) * sum(aicetl(n:ncat+1)) ! area of sub-layer occupied by snow
          
          vol = (hitl(n) - hitl(n-1)) * area      ! thickness of sub-layer times area
          
          if (vol >= rem_vol) then  ! have reached the sub-layer with the depth within
             hpond = rem_vol / area + hitl(n-1) + alfan(m_index) - alfan(1)
             exit
          else  ! still in sub-layer below the sub-layer with the depth
             rem_vol = rem_vol - vol
          endif
          
       enddo
       
    endif
    
  end subroutine calc_hpond
  
!=======================================================================
!BOP
!
! !ROUTINE: permeability_phi
!
! !DESCRIPTION:
!
! determine the liquid fraction of brine in the ice and the permeability
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine permeability_phi(eicen, vicen, perm)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(nilyr), intent(in) :: &
         eicen     ! energy of melting for each ice layer (J/m2)
    
      real (kind=dbl_kind), intent(in) :: &
         vicen     ! ice volume
    
      real (kind=dbl_kind), intent(out) :: &
         perm      ! permeability
!
!EOP
!
      real (kind=dbl_kind) ::   &
         Sbr, &    ! brine salinity
         qin       ! enthalpy

      real (kind=dbl_kind), dimension(nilyr) ::   &
         Tin, &    ! ice temperature
         phi       ! liquid fraction

      integer (kind=int_kind) :: k
    
      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------

      do k = 1,nilyr
         qin    = eicen(k)*real(nilyr,kind=dbl_kind) / vicen
         Tin(k) = calculate_Tin_from_qin(qin,Tmlt(k))
      enddo

      !-----------------------------------------------------------------
      ! brine salinity and liquid fraction
      !-----------------------------------------------------------------

      if (maxval(Tin) <= -c2) then

         do k = 1,nilyr
            Sbr = - 1.2_dbl_kind                 &
                  -21.8_dbl_kind     * Tin(k)    &
                  - 0.919_dbl_kind   * Tin(k)**2 &
                  - 0.01878_dbl_kind * Tin(k)**3
            phi(k) = salin(k)/Sbr ! liquid fraction
         enddo ! k
       
      else

         do k = 1,nilyr
            Sbr = -17.6_dbl_kind    * Tin(k)    &
                  - 0.389_dbl_kind  * Tin(k)**2 &
                  - 0.00362_dbl_kind* Tin(k)**3
            phi(k) = salin(k)/Sbr ! liquid fraction
         enddo

      endif
    
      !-----------------------------------------------------------------
      ! permeability
      !-----------------------------------------------------------------

      perm = 3.0e-08_dbl_kind * (minval(phi))**3
    
      end subroutine permeability_phi
  
!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_pond - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_pond_topo(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for restarting
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR
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
              restart_file(1:lenstr(restart_file)),'.volpn.', &
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
         call ice_write(nu_dump_pond,0,trcrn(:,:,nt_apnd,n,:),'ruf8',diag)
         call ice_write(nu_dump_pond,0,trcrn(:,:,nt_hpnd,n,:),'ruf8',diag)
         call ice_write(nu_dump_pond,0,trcrn(:,:,nt_ipnd,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_dump_pond)

      end subroutine write_restart_pond_topo

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_pond - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_pond_topo(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for a meltpond volume restart
!
! !REVISION HISTORY:
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR
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
         if (n == 0) call abort_ice('volpn restart: filename discrepancy')
         string1 = trim(filename0(1:n-1))
         string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
         write(filename,'(a,a,a,a)') &
            string1(1:lenstr(string1)), &
            restart_file(1:lenstr(restart_file)),'.volpn', &
            string2(1:lenstr(string2))
      endif ! master_task

      call ice_open(nu_restart_pond,filename,0)

      if (my_task == master_task) then
        read(nu_restart_pond) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      do n = 1, ncat
         call ice_read(nu_restart_pond,0,trcrn(:,:,nt_apnd,n,:),'ruf8',diag)
         call ice_read(nu_restart_pond,0,trcrn(:,:,nt_hpnd,n,:),'ruf8',diag)
         call ice_read(nu_restart_pond,0,trcrn(:,:,nt_ipnd,n,:),'ruf8',diag)
      enddo

      if (my_task == master_task) close(nu_restart_pond)

      end subroutine read_restart_pond_topo

!=======================================================================

      end module ice_meltpond_topo

!=======================================================================
