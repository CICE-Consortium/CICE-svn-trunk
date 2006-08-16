!=======================================================================
!BOP
!
! !MODULE: ice_therm_itd - thermo calculations after call to coupler
!
! !DESCRIPTION:
!
! Thermo calculations after call to coupler, related to ITD:
! ice thickness redistribution, lateral growth and melting.
!
! NOTE: The thermodynamic calculation is split in two for load balancing.
!       First ice_therm_vertical computes vertical growth rates and coupler
!       fluxes.  Then ice_therm_itd does thermodynamic calculations not
!       needed for coupling.
!       
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb.  
!
! !INTERFACE:
!
      module ice_therm_itd
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_domain_size
      use ice_constants
      use ice_fileunits
!
!EOP
!
      implicit none
      save

      real (kind=dbl_kind), parameter, private ::
     &   hfrazilmin = 0.05_dbl_kind ! min thickness of new frazil ice (m)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: linear_itd - ITD scheme that shifts ice among categories
!
! !INTERFACE:
!
      subroutine linear_itd (nx_block,    ny_block,
     &                       nghost,      trcr_depend,
     &                       aicen_init,  vicen_init,
     &                       aicen,       trcrn,
     &                       vicen,       vsnon,
     &                       eicen,       esnon,
     &                       aice,        aice0,
     &                       l_stop,
     &                       istop,       jstop)

! See Lipscomb, W. H.  Remapping the thickness distribution in sea
!     ice models. 2001, J. Geophys. Res., Vol 106, 13989--14000.
!
! Using the thermodynamic "velocities", interpolate to find the
! velocities in thickness space at the category boundaries, and
! compute the new locations of the boundaries.  Then for each
! category, compute the thickness distribution function,  g(h),
! between hL and hR, the left and right boundaries of the category.
! Assume g(h) is a linear polynomial that satisfies two conditions:
!
! (1) The ice area implied by g(h) equals aicen(n).
! (2) The ice volume implied by g(h) equals aicen(n)*hicen(n).
!
! Given g(h), at each boundary compute the ice area and volume lying
! between the original and new boundary locations.  Transfer area
! and volume across each boundary in the appropriate direction, thus
! restoring the original boundaries.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_itd, only: hin_max, hi_min, aggregate_area, shift_ice,
     &                   column_sum, column_conservation_check 
      use ice_work, only: worka
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::
     &   trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat),
     &   intent(in) ::
     &   aicen_init ! initial ice concentration (before vertical thermo)
     &,  vicen_init ! initial ice volume               (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(inout) ::
     &   aicen      ! ice concentration
     &,  vicen      ! volume per unit area of ice      (m)
     &,  vsnon      ! volume per unit area of snow     (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),
     &   intent(inout) ::
     &   trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),
     &   intent(inout) ::
     &   eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),
     &   intent(inout) ::
     &   esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   aice      ! concentration of ice
     &,  aice0     ! concentration of open water

      logical (kind=log_kind), intent(out) ::
     &   l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(out) ::
     &   istop, jstop    ! indices of grid cell where model aborts 
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j             ! horizontal indices
     &,  n, nd            ! category indices
     &,  k                ! ice layer index
     &,  ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      real (kind=dbl_kind) ::
     &   slope            ! rate of change of dhice with hice
     &,  dh0              ! change in ice thickness at h = 0
     &,  da0              ! area melting from category 1
     &,  damax            ! max allowed reduction in category 1 area
     &,  etamin, etamax   ! left and right limits of integration
     &,  x1               ! etamax - etamin
     &,  x2               ! (etamax^2 - etamin^2) / 2
     &,  x3               ! (etamax^3 - etamin^3) / 3
     &,  wk1, wk2         ! temporary variables

      real (kind=dbl_kind), dimension(nx_block,ny_block,0:ncat) ::
     &   hbnew            ! new boundary locations

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat) ::
     &   hicen            ! ice thickness for each cat     (m)
     &,  hicen_init       ! initial ice thickness for each cat (pre-thermo)
     &,  dhicen           ! thickness change for remapping (m)
     &,  g0               ! constant coefficient in g(h)
     &,  g1               ! linear coefficient in g(h)
     &,  hL               ! left end of range over which g(h) > 0
     &,  hR               ! right end of range over which g(h) > 0

      real (kind=dbl_kind), dimension(nx_block,ny_block) ::
     &   vice_init, vice_final  ! ice volume summed over categories
     &,  vsno_init, vsno_final  ! snow volume summed over categories
     &,  eice_init, eice_final  ! ice energy summed over categories
     &,  esno_init, esno_final  ! snow energy summed over categories

      ! NOTE: Third index of donor, daice, dvice should be ncat-1,
      !       except that compilers would have trouble when ncat = 1 
      integer (kind=int_kind), dimension(nx_block,ny_block,ncat) ::
     &   donor            ! donor category index

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat) ::
     &   daice            ! ice area transferred across boundary
     &,  dvice            ! ice volume transferred across boundary

      logical (kind=log_kind), dimension(nx_block,ny_block) ::
     &   remap_flag       ! remap ITD if remap_flag(i,j) is true

      character (len=char_len) ::
     &   fieldid           ! field identifier

      logical (kind=log_kind), parameter ::
     &   l_conservation_check = .true.   ! if true, check conservation
                                         ! (useful for debugging)

       integer (kind=int_kind) ::
     &   icells            ! number of grid cells with ice
     &,  ij                ! combined horizontal index

       integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &   indxi             ! compressed i/j indices
     &,  indxj

      logical (kind=log_kind) ::
     &   flag_changed      ! true if remapping is disabled somewhere


      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

      hin_max(ncat) = 999.9_dbl_kind ! arbitrary big number
      do n = 0, ncat
         hbnew(:,:,n) = hin_max(n)
      enddo

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      !-----------------------------------------------------------------
      ! Compute volume and energy sums that linear remapping should
      !  conserve.
      !-----------------------------------------------------------------

      if (l_conservation_check) then
         call column_sum (nx_block, ny_block,
     &                    nghost,   ncat,
     &                    vicen,    vice_init)
         call column_sum (nx_block, ny_block,
     &                    nghost,   ncat,
     &                    vsnon,    vsno_init)
         call column_sum (nx_block, ny_block,
     &                    nghost,   ntilyr,
     &                    eicen,    eice_init)
         call column_sum (nx_block, ny_block,
     &                    nghost,   ntslyr,
     &                    esnon,    esno_init)
      endif

      !-----------------------------------------------------------------
      ! Compute thickness change in each category.
      !-----------------------------------------------------------------

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block

            if (aicen_init(i,j,n) > puny) then
               hicen_init (i,j,n) = vicen_init(i,j,n) / 
     &                              aicen_init(i,j,n) 
            else
               hicen_init(i,j,n) = c0
            endif               ! aicen_init > puny

            if (aicen(i,j,n) > puny) then
               hicen (i,j,n) = vicen(i,j,n) / aicen(i,j,n) 
               dhicen(i,j,n) = hicen(i,j,n) - hicen_init(i,j,n)
            else
               hicen (i,j,n) = c0
               dhicen(i,j,n) = c0
            endif               ! aicen > puny

         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Compute fractional ice area in each grid cell.
      !-----------------------------------------------------------------
      call aggregate_area (nx_block, ny_block,
     &                     aicen,
     &                     aice,     aice0,
     &                     l_stop,
     &                     istop,    jstop)

      if (l_stop) return

      !-----------------------------------------------------------------
      ! Identify grid cells with ice and initialize remapping flag.
      ! Remapping is done wherever remap_flag = .true.
      ! In rare cases the category boundaries may shift too far for the
      !  remapping algorithm to work, and remap_flag is set to .false.
      ! In these cases the simpler 'rebin' subroutine will shift ice
      !  between categories if needed.
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo,jhi
      do i = ilo,ihi
         if (aice(i,j) > puny) then
            remap_flag(i,j) = .true.
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         else
            remap_flag(i,j) = .false.
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute new category boundaries, hbnew, based on changes in
      ! ice thickness from vertical thermodynamics.
      !-----------------------------------------------------------------

      do n = 1, ncat-1

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells       ! aice(i,j) > puny
            i = indxi(ij)
            j = indxj(ij)

            if (hicen_init(i,j,n)   > puny .and.
     &          hicen_init(i,j,n+1) > puny) then
                 ! interpolate between adjacent category growth rates
               slope = (dhicen(i,j,n+1) - dhicen(i,j,n)) /
     &                 (hicen_init(i,j,n+1) - hicen_init(i,j,n))
               hbnew(i,j,n) = hin_max(n) + dhicen(i,j,n)
     &                      + slope * (hin_max(n) - hicen_init(i,j,n))
            elseif (hicen_init(i,j,n) > puny) then ! hicen_init(n+1)=0
               hbnew(i,j,n) = hin_max(n) + dhicen(i,j,n)
            elseif (hicen_init(i,j,n+1) > puny) then ! hicen_init(n)=0
               hbnew(i,j,n) = hin_max(n) + dhicen(i,j,n+1)
            else
               hbnew(i,j,n) = hin_max(n)
            endif
         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Check that each boundary lies between adjacent values of hicen.
      ! If not, set remap_flag = .false.
      !-----------------------------------------------------------------

         flag_changed = .false.
         do ij = 1, icells       ! aice(i,j) > puny
            i = indxi(ij)
            j = indxj(ij)

            if (aicen(i,j,n) > puny .and.
     &          hicen(i,j,n) >= hbnew(i,j,n)) then
               remap_flag(i,j) = .false.
               flag_changed = .true.
            elseif (aicen(i,j,n+1) > puny .and.
     &              hicen(i,j,n+1) <= hbnew(i,j,n)) then
               remap_flag(i,j) = .false.
               flag_changed = .true.
            endif

      !-----------------------------------------------------------------
      ! Check that hbnew(n) lies between hin_max(n-1) and hin_max(n+1).
      ! If not, set remap_flag = .false.
      ! (In principle we could allow this, but it would make the code
      ! more complicated.)
      !-----------------------------------------------------------------

            if (hbnew(i,j,n) > hin_max(n+1)) then
               remap_flag(i,j) = .false.
               flag_changed = .true.
            endif

            if (hbnew(i,j,n) < hin_max(n-1)) then
               remap_flag(i,j) = .false.
               flag_changed = .true.
            endif

         enddo                  ! ij

      !-----------------------------------------------------------------
      ! Write diagnosis outputs if remap_flag was changed to false
      !-----------------------------------------------------------------

         if (flag_changed) then
            do ij = 1, icells   ! aice(i,j) > puny

               i = indxi(ij)
               j = indxj(ij)

               if (aicen(i,j,n) > puny .and.
     &             hicen(i,j,n) >= hbnew(i,j,n)) then
                  write(nu_diag,*) my_task,':',i,j,
     &                 'ITD: hicen(n) > hbnew(n)'
                  write(nu_diag,*) 'cat ',n
                  write(nu_diag,*) my_task,':',i,j,
     &                 'hicen(n) =', hicen(i,j,n)
                  write(nu_diag,*) my_task,':',i,j,
     &                 'hbnew(n) =', hbnew(i,j,n)
               elseif (aicen(i,j,n+1) > puny .and.
     &                 hicen(i,j,n+1) <= hbnew(i,j,n)) then
                  write(nu_diag,*) my_task,':',i,j,
     &                 'ITD: hicen(n+1) < hbnew(n)'
                  write(nu_diag,*) 'cat ',n
                  write(nu_diag,*) my_task,':',i,j,
     &                 'hicen(n+1) =', hicen(i,j,n+1)
                  write(nu_diag,*) my_task,':',i,j,
     &                 'hbnew(n) =', hbnew(i,j,n)
               endif

               if (hbnew(i,j,n) > hin_max(n+1)) then
                  write(nu_diag,*) my_task,':',i,j,
     &                 'ITD hbnew(n) > hin_max(n+1)'
                  write(nu_diag,*) 'cat ',n
                  write(nu_diag,*) my_task,':',i,j,
     &                 'hbnew(n) =', hbnew(i,j,n)
                  write(nu_diag,*) my_task,':',i,j,
     &                 'hin_max(n+1) =', hin_max(n+1)
               endif

               if (hbnew(i,j,n) < hin_max(n-1)) then
                  write(nu_diag,*) my_task,':',i,j,
     &                 'ITD: hbnew(n) < hin_max(n-1)'
                  write(nu_diag,*) 'cat ',n
                  write(nu_diag,*) my_task,':',i,j,
     &                 'hbnew(n) =', hbnew(i,j,n)
                  write(nu_diag,*) my_task,':',i,j,
     &                 'hin_max(n-1) =', hin_max(n-1)
               endif

            enddo               ! ij
         endif                  ! flag_changed

      enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Compute hbnew(ncat)
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells         ! aice(i,j) > puny
         i = indxi(ij)
         j = indxj(ij)
         if (aicen(i,j,ncat) > puny) then
            hbnew(i,j,ncat) = c3*hicen(i,j,ncat) - c2*hbnew(i,j,ncat-1)
         else
            hbnew(i,j,ncat) = hin_max(ncat)
         endif
         
         if (hbnew(i,j,ncat) < hin_max(ncat-1))
     &        hbnew(i,j,ncat) = hin_max(ncat-1)
      enddo

      !-----------------------------------------------------------------
      ! Identify cells where the ITD is to be remapped
      !-----------------------------------------------------------------

      icells = 0
      do j = jlo,jhi
      do i = ilo,ihi
         if (remap_flag(i,j)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo
      enddo


      !-----------------------------------------------------------------
      ! Compute g(h) for category 1 at start of time step
      ! (hicen = hicen_init)
      !-----------------------------------------------------------------

      worka(:,:) = hin_max(1)

      call fit_line(nx_block,       ny_block,
     &              nghost,
     &              aicen(:,:,1),   hicen_init(:,:,1),
     &              hbnew(:,:,0),   worka     (:,:),
     &              g0   (:,:,1),   g1        (:,:,1),
     &              hL   (:,:,1),   hR        (:,:,1),
     &              remap_flag)

      !-----------------------------------------------------------------
      ! Find area lost due to melting of thin (category 1) ice
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells    ! remap_flag = .true.
         i = indxi(ij)
         j = indxj(ij)

         if (aicen(i,j,1) > puny) then

            dh0 = dhicen(i,j,1)

            if (dh0 < c0) then   ! remove area from category 1
               dh0 = min(-dh0,hin_max(1))   ! dh0 --> |dh0|

      !-----------------------------------------------------------------
      ! Integrate g(1) from 0 to dh0 to estimate area melted
      !-----------------------------------------------------------------

               ! right integration limit (left limit = 0)
               etamax = min(dh0,hR(i,j,1)) - hL(i,j,1)

               if (etamax > c0) then
                  x1 = etamax
                  x2 = p5 * etamax*etamax
                  da0 = g1(i,j,1)*x2 + g0(i,j,1)*x1 ! ice area removed

               ! constrain new thickness <= hicen_init
                  damax = aicen(i,j,1)
     &                  * (c1-hicen(i,j,1)/hicen_init(i,j,1)) ! damax > 0
                  da0 = min (da0, damax)

               ! remove area, conserving volume
                  hicen(i,j,1) = hicen(i,j,1)
     &                         * aicen(i,j,1) / (aicen(i,j,1)-da0)
                  aicen(i,j,1) = aicen(i,j,1) - da0
               endif            ! etamax > 0

            else                ! dh0 >= 0
               hbnew(i,j,0) = min(dh0,hin_max(1))  ! shift hbnew(0) to right
            endif

         endif                  ! aicen(i,j,1) > puny
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Compute g(h) for each ice thickness category.
      !-----------------------------------------------------------------

      do n = 1, ncat
         call fit_line(nx_block,       ny_block,
     &                 nghost,
     &                 aicen(:,:,n),   hicen(:,:,n),
     &                 hbnew(:,:,n-1), hbnew(:,:,n),
     &                 g0   (:,:,n),   g1   (:,:,n),
     &                 hL   (:,:,n),   hR   (:,:,n),
     &                 remap_flag)

      enddo

      !-----------------------------------------------------------------
      ! Compute area and volume to be shifted across each boundary.
      !-----------------------------------------------------------------

      donor(:,:,:) = 0
      daice(:,:,:) = c0
      dvice(:,:,:) = c0

      do n = 1, ncat-1

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells   ! remap_flag = .true.
            i = indxi(ij)
            j = indxj(ij)

            if (hbnew(i,j,n) > hin_max(n)) then ! transfer from n to n+1

               ! left and right integration limits in eta space
               etamin = max(hin_max(n), hL(i,j,n)) - hL(i,j,n)
               etamax = min(hbnew(i,j,n), hR(i,j,n)) - hL(i,j,n)
               donor(i,j,n) = n

            else             ! hbnew(n) <= hin_max(n); transfer from n+1 to n

               ! left and right integration limits in eta space
               etamin = c0
               etamax = min(hin_max(n), hR(i,j,n+1)) - hL(i,j,n+1)
               donor(i,j,n) = n+1

            endif            ! hbnew(n) > hin_max(n)

            if (etamax > etamin) then
               x1  = etamax - etamin
               wk1 = etamin*etamin
               wk2 = etamax*etamax
               x2  = p5 * (wk2 - wk1)
               wk1 = wk1*etamin
               wk2 = wk2*etamax
               x3  = p333 * (wk2 - wk1)
               nd  = donor(i,j,n)
               daice(i,j,n) = g1(i,j,nd)*x2 + g0(i,j,nd)*x1
               dvice(i,j,n) = g1(i,j,nd)*x3 + g0(i,j,nd)*x2
     &                      + daice(i,j,n)*hL(i,j,nd)
            endif               ! etamax > etamin

            ! If daice or dvice is very small, shift no ice.

            nd = donor(i,j,n)

            if (daice(i,j,n) < aicen(i,j,nd)*puny) then
               daice(i,j,n) = c0
               dvice(i,j,n) = c0
               donor(i,j,n) = 0
            endif 

            if (dvice(i,j,n) < vicen(i,j,nd)*puny) then
               daice(i,j,n) = c0
               dvice(i,j,n) = c0
               donor(i,j,n) = 0
            endif

            ! If daice is close to aicen or dvice is close to vicen,
            ! shift entire category

            if (daice(i,j,n) > aicen(i,j,nd)*(c1-puny)) then
               daice(i,j,n) = aicen(i,j,nd)
               dvice(i,j,n) = vicen(i,j,nd)
            endif

            if (dvice(i,j,n) > vicen(i,j,nd)*(c1-puny)) then
               daice(i,j,n) = aicen(i,j,nd)
               dvice(i,j,n) = vicen(i,j,nd)
            endif

         enddo                  ! ij
      enddo                     ! boundaries, 1 to ncat-1

      !-----------------------------------------------------------------
      ! Shift ice between categories as necessary
      !-----------------------------------------------------------------

      call shift_ice (nx_block, ny_block,
     &                nghost,   trcr_depend,
     &                aicen,    trcrn,
     &                vicen,    vsnon,
     &                eicen,    esnon,
     &                hicen,    donor,
     &                daice,    dvice,
     &                l_stop,
     &                istop,    jstop)

      if (l_stop) return

      !-----------------------------------------------------------------
      ! Make sure hice(i,j,1) >= minimum ice thickness hi_min.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells          ! remap_flag = .true.
         i = indxi(ij)
         j = indxj(ij)
         if (hi_min > c0 .and.
     &        aicen(i,j,1) > puny .and. hicen(i,j,1) < hi_min) then
            aicen(i,j,1) = aicen(i,j,1) * hicen(i,j,1)/hi_min
            hicen(i,j,1) = hi_min
         endif
      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Update fractional ice area in each grid cell.
      !-----------------------------------------------------------------
      call aggregate_area (nx_block, ny_block,
     &                     aicen,
     &                     aice,     aice0,
     &                     l_stop,
     &                     istop,    jstop)

      if (l_stop) return

      !-----------------------------------------------------------------
      ! Check volume and energy conservation.
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         call column_sum (nx_block, ny_block,
     &                    nghost,   ncat,
     &                    vicen,    vice_final)
         fieldid = 'vice, ITD remap'
         call column_conservation_check (nx_block,  ny_block,
     &                                   nghost,    fieldid,
     &                                   vice_init, vice_final,
     &                                   puny,      l_stop,
     &                                   istop,     jstop)
         if (l_stop) return
         

         call column_sum (nx_block, ny_block,
     &                    nghost,   ncat,
     &                    vsnon,    vsno_final)
         fieldid = 'vsno, ITD remap'
         call column_conservation_check (nx_block,  ny_block,
     &                                   nghost,    fieldid,
     &                                   vsno_init, vsno_final,
     &                                   puny,      l_stop,
     &                                   istop,     jstop)
         if (l_stop) return

         call column_sum (nx_block, ny_block,
     &                    nghost,   ntilyr,
     &                    eicen,    eice_final)
         fieldid = 'eice, ITD remap'
         call column_conservation_check (nx_block,   ny_block,
     &                                   nghost,     fieldid,
     &                                   eice_init,  eice_final,
     &                                   puny*Lfresh*rhoi,
     &                                   l_stop,
     &                                   istop,     jstop)
         if (l_stop) return

         call column_sum (nx_block, ny_block,
     &                    nghost,   ntslyr,
     &                    esnon,    esno_final)
         fieldid = 'esno, ITD remap'
         call column_conservation_check (nx_block,   ny_block,
     &                                   nghost,     fieldid,
     &                                   esno_init,  esno_final,
     &                                   puny*Lfresh*rhos,
     &                                   l_stop,
     &                                   istop,     jstop)
         if (l_stop) return

      endif                     ! conservation check

      end subroutine linear_itd

!=======================================================================
!BOP
!
! !IROUTINE: fit_line - fit g(h) with a line using area, volume constraints
!
! !INTERFACE:
!

      subroutine fit_line (nx_block, ny_block,
     &                     nghost,
     &                     aicen,    hice,
     &                     hbL,      hbR,   
     &                     g0,       g1,
     &                     hL,       hR,
     &                     remap_flag)
!
! !DESCRIPTION:
!
! Fit g(h) with a line, satisfying area and volume constraints.
! To reduce roundoff errors caused by large values of g0 and g1,
! we actually compute g(eta), where eta = h - hL, and hL is the
! left boundary.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aicen           ! concentration of ice
     &,  hice            ! ice thickness

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   hbL, hbR        ! left and right category boundaries

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out)::
     &   g0, g1          ! coefficients in linear equation for g(eta)
     &,  hL              ! min value of range over which g(h) > 0
     &,  hR              ! max value of range over which g(h) > 0

      logical (kind=log_kind), dimension (nx_block,ny_block),
     &   intent(in) ::
     &   remap_flag
!
!EOP
!
      integer (kind=int_kind) ::
     &   i,j             ! horizontal indices
     &,  ilo,ihi,jlo,jhi ! beginning and end of physical domain

      real  (kind=dbl_kind) ::
     &   h13             ! hbL + 1/3 * (hbR - hbL)
     &,  h23             ! hbL + 2/3 * (hbR - hbL)
     &,  dhr             ! 1 / (hR - hL)
     &,  wk1, wk2        ! temporary variables

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      g0(:,:) = c0
      g1(:,:) = c0
      hL(:,:) = c0
      hR(:,:) = c0

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      !-----------------------------------------------------------------
      ! Compute g0, g1, hL, and hR for each category to be remapped.
      !-----------------------------------------------------------------

      do j = jlo,jhi
      do i = ilo,ihi

         if (remap_flag(i,j) .and. aicen(i,j) > puny
     &              .and. hbR(i,j) - hbL(i,j) > puny) then

         ! Initialize hL and hR

            hL(i,j) = hbL(i,j)
            hR(i,j) = hbR(i,j)

         ! Change hL or hR if hicen(n) falls outside central third of range

            h13 = p333 * (c2*hL(i,j) + hR(i,j))
            h23 = p333 * (hL(i,j) + c2*hR(i,j))
            if (hice(i,j) < h13) then
               hR(i,j) = c3*hice(i,j) - c2*hL(i,j)
            elseif (hice(i,j) > h23) then
               hL(i,j) = c3*hice(i,j) - c2*hR(i,j)
            endif

         ! Compute coefficients of g(eta) = g0 + g1*eta

            dhr = c1 / (hR(i,j) - hL(i,j))
            wk1 = c6 * aicen(i,j) * dhr
            wk2 = (hice(i,j) - hL(i,j)) * dhr
            g0(i,j) = wk1 * (p666 - wk2)
            g1(i,j) = c2*dhr * wk1 * (wk2 - p5)

         endif                  ! aicen > puny

      enddo                     ! i
      enddo                     ! j

      end subroutine fit_line

!=======================================================================
!BOP
!
! !ROUTINE: add_new_ice - add frazil ice to ice thickness distribution
!
! !DESCRIPTION:
!
! Given the volume of new ice grown in open water, compute its area
! and thickness and add it to the appropriate category or categories.
!
! NOTE: Usually all the new ice is added to category 1.  An exception is
!       made if there is no open water or if the new ice is too thick
!       for category 1, in which case ice is distributed evenly over the
!       entire cell.  Subroutine rebin should be called in case the ice
!       thickness lies outside category bounds after new ice formation.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine add_new_ice (nx_block,  ny_block,
     &                        nghost,
     &                        tmask,     dt,
     &                        aicen,     trcrn,
     &                        vicen,     eicen,
     &                        aice0,     aice,
     &                        frzmlt,    frazil,
     &                        frz_onset, yday,
     &                        l_stop,
     &                        istop,     jstop)
!
! !USES:
!
      use ice_itd, only: hin_max, ilyr1, column_sum,
     &                   column_conservation_check

! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      logical (kind=log_kind), dimension (nx_block,ny_block),
     &    intent(in) ::
     &   tmask     ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), intent(in) ::
     &   dt        ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aice      ! total concentration of ice
     &,  frzmlt    ! freezing/melting potential (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(inout) ::
     &   aicen     ! concentration of ice
     &,  vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),
     &   intent(inout) ::
     &   trcrn     ! ice tracers
                   ! 1: surface temperature

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),
     &   intent(inout) ::
     &   eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   aice0     ! concentration of open water
     &,  frazil    ! frazil ice growth        (m/step-->cm/day)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout), optional ::
     &   frz_onset ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(in), optional ::
     &   yday      ! day of year

      logical (kind=log_kind), intent(out) ::
     &   l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(out) ::
     &   istop, jstop    ! indices of grid cell where model aborts
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j             ! horizontal indices
     &,  n                ! ice category index
     &,  k                ! ice layer index
     &,  ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   ai0new           ! area of new ice added to cat 1
     &,  vi0new           ! volume of new ice added to cat 1
     &,  hsurp            ! thickness of new ice added to each cat
     &,  vlyr             ! ice layer volume

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   vice_init, vice_final  ! ice volume summed over categories

      real (kind=dbl_kind) ::
     &   fnew             ! heat flx to open water for new ice (W/m^2)
     &,  hi0new           ! thickness of new ice
     &,  hi0max           ! max allowed thickness of new ice
     &,  qi0(nilyr)       ! frazil ice enthalpy
     &,  qi0av            ! mean value of qi0 for new ice (J kg-1)
     &,  vsurp            ! volume of new ice added to each cat
     &,  area1            ! starting fractional area of existing ice
     &,  rnilyr           ! real(nilyr)
     &,  dfresh           ! change in fresh
     &,  dfsalt           ! change in fsalt

      integer (kind=int_kind) ::
     &   icells, jcells, kcells ! grid cell counters
     &,  ij                     ! combined i/j horizontal index

      integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &   indxi,  indxj          ! compressed i/j indices
     &,  indxi2, indxj2
     &,  indxi3, indxj3

      character (len=char_len) ::
     &   fieldid           ! field identifier

      l_stop = .false.
      istop = 0
      jstop = 0

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      if (ncat > 1) then
         hi0max = hin_max(1)*0.9_dbl_kind  ! not too close to boundary
      else
         hi0max = bignum                   ! big number
      endif

      ! initial ice volume in each grid cell
      call column_sum (nx_block, ny_block,
     &                 nghost,   ncat,
     &                 vicen,    vice_init)

      !-----------------------------------------------------------------
      ! Compute average enthalpy of new ice.
      !
      ! POP assumes new ice is fresh.  Otherwise, it would be better
      ! to do something like this:
      !  qi0(i,j,k) = -rhoi * (cp_ice*(Tmlt(k)-Tf(i,j))
      !             + Lfresh*(1.-Tmlt(k)/Tf(i,j)) - cp_ocn*Tmlt(k))
      !-----------------------------------------------------------------

      rnilyr = real(nilyr,kind=dbl_kind)
      qi0av = c0
      do k = 1, nilyr
         qi0(k) = -rhoi*Lfresh  ! note sign convention, qi < 0
         qi0av  = qi0av + qi0(k)
      enddo
      qi0av = qi0av/rnilyr

      !-----------------------------------------------------------------
      ! Identify ice/ocean grid points.
      !-----------------------------------------------------------------
      icells = 0
      do j = jlo, jhi
      do i = ilo, ihi
         if (tmask(i,j)) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo       ! i
      enddo       ! j

      !-----------------------------------------------------------------
      ! Compute the volume, area, and thickness of new ice.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         fnew = max (frzmlt(i,j), c0)   ! fnew > 0 iff frzmlt > 0
         vi0new(i,j) = -fnew*dt / qi0av ! note sign convention, qi < 0

         ! increment ice volume
         vice_init(i,j) = vice_init(i,j) + vi0new(i,j)

         ! history diagnostics
         frazil(i,j) = vi0new(i,j)

         if (present(frz_onset) .and. present(yday)) then
            if (frazil(i,j) > puny .and. frz_onset(i,j) < puny)
     &           frz_onset(i,j) = yday
         endif

      !-----------------------------------------------------------------
      ! Update fresh water and salt fluxes.
      !
      ! NOTE: POP assumes fresh water and salt flux due to frzmlt > 0
      !       is NOT included in fluxes fresh and fsalt.
      !-----------------------------------------------------------------

!!!         dfresh = -rhoi*vi0new(i,j)/dt  ! if POP had not already adjusted
                                           ! itself based on frzmlt
!!!         dfsalt = ice_ref_salinity*p001*dfresh

!!!         fresh(i,j)      = fresh(i,j)      + dfresh
!!!         fresh_hist(i,j) = fresh_hist(i,j) + dfresh
!!!         fsalt(i,j)      = fsalt(i,j)      + dfsalt
!!!         fsalt_hist(i,j) = fsalt_hist(i,j) + dfsalt

      !-----------------------------------------------------------------
      ! Decide how to distribute the new ice.
      !-----------------------------------------------------------------

         hsurp(i,j)  = c0
         ai0new(i,j) = c0

         if (vi0new(i,j) > c0) then

            ! new ice area and thickness
            ! hin_max(0) < new ice thickness < hin_max(1)
            if (aice0(i,j) > puny) then
               hi0new = max(vi0new(i,j)/aice0(i,j), hfrazilmin)
               if (hi0new > hi0max .and. aice0(i,j)+puny < c1) then
                  ! distribute excess volume over all categories (below)
                  hi0new = hi0max
                  ai0new(i,j) = aice0(i,j)
                  vsurp       = vi0new(i,j) - ai0new(i,j)*hi0new
                  hsurp(i,j)  = vsurp / aice(i,j)
                  vi0new(i,j) = ai0new(i,j)*hi0new
               else
                  ! put ice in a single category, with hsurp = 0
                  ai0new(i,j) = vi0new(i,j)/hi0new
               endif
            else                ! aice0 < puny
               hsurp(i,j) = vi0new(i,j)/aice(i,j) ! new thickness in each cat
               vi0new(i,j) = c0
            endif               ! aice0 > puny
         endif                  ! vi0new > puny

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! Identify grid cells receiving new ice.
      !-----------------------------------------------------------------
      jcells = 0
      kcells = 0

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         if (vi0new(i,j) > c0) then  ! add ice to category 1
            jcells = jcells + 1
            indxi2(jcells) = i
            indxj2(jcells) = j
         endif

         if (hsurp(i,j) > c0) then   ! add ice to all categories
            kcells = kcells + 1
            indxi3(kcells) = i
            indxj3(kcells) = j
         endif

      enddo

      !-----------------------------------------------------------------
      ! Distribute excess ice volume among ice categories by increasing
      ! ice thickness, leaving ice area unchanged.
      !
      ! NOTE: If new ice contains globally conserved tracers
      !       (e.g., isotopes from seawater), code must be added here.
      !-----------------------------------------------------------------

      do n = 1, ncat

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, kcells
            i = indxi3(ij)
            j = indxj3(ij)

            vicen(i,j,n) = vicen(i,j,n) + aicen(i,j,n)*hsurp(i,j)
            vlyr(i,j) = hsurp(i,j)/rnilyr * aicen(i,j,n)
         enddo                  ! ij

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, kcells
               i = indxi3(ij)
               j = indxj3(ij)

               eicen(i,j,ilyr1(n)+k-1) =
     &              eicen(i,j,ilyr1(n)+k-1) + qi0(k)*vlyr(i,j)
            enddo               ! ij
         enddo                  ! k

      enddo                     ! n

      !-----------------------------------------------------------------
      ! Combine new ice grown in open water with category 1 ice.
      ! Assume that vsnon and esnon are unchanged.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, jcells
         i = indxi2(ij)
         j = indxj2(ij)

         area1 = aicen(i,j,1)   ! save
         aicen(i,j,1) = aicen(i,j,1) + ai0new(i,j)
         aice0(i,j)   = aice0(i,j)   - ai0new(i,j)
         vicen(i,j,1) = vicen(i,j,1) + vi0new(i,j)
         trcrn(i,j,1,1) = (Tocnfrz*ai0new(i,j) + trcrn(i,j,1,1)*area1)
     &                / aicen(i,j,1)
         trcrn(i,j,1,1) = min (trcrn(i,j,1,1), c0)

      ! For other tracers, do something like this:
!         trcrn(i,j,99,1) = (tnew(i,j)*ai0new(i,j) + trcrn(i,j,99,1)*area1)
!     &                   / aicen(i,j,1)

         vlyr(i,j)    = vi0new(i,j) / rnilyr
      enddo                     ! ij

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, jcells
            i = indxi2(ij)
            j = indxj2(ij)
            eicen(i,j,k) = eicen(i,j,k) + qi0(k)*vlyr(i,j)
         enddo
      enddo

      call column_sum (nx_block, ny_block,
     &                 nghost,   ncat,
     &                 vicen,    vice_final)

      fieldid = 'vice, add_new_ice'
      call column_conservation_check (nx_block,  ny_block,
     &                                nghost,    fieldid,
     &                                vice_init, vice_final,
     &                                puny,      l_stop,
     &                                istop,     jstop)
      if (l_stop) return

      end subroutine add_new_ice

!=======================================================================
!BOP
!
! !ROUTINE: lateral_melt - melt ice laterally
!
! !DESCRIPTION:
!
! Given the fraction of ice melting laterally in each grid cell
!  (computed in subroutine frzmlt_bottom_lateral), melt ice.
!
! !REVISION HISTORY:
!
! author: C. M. Bitz, UW
! 2003:   Modified by William H. Lipscomb and Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine lateral_melt (nx_block,   ny_block,
     &                         nghost,     dt,
     &                         fresh,      fsalt,
     &                         fhocn,      fresh_hist,
     &                         fsalt_hist, fhocn_hist,
     &                         rside,      meltl,
     &                         aicen,      vicen,
     &                         vsnon,      eicen,
     &                         esnon)
!
! !USES:
!
      use ice_itd, only: ilyr1, slyr1
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), intent(in) ::
     &   dt        ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(inout) ::
     &   aicen       ! concentration of ice
     &,  vicen       ! volume per unit area of ice          (m)
     &,  vsnon       ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),
     &   intent(inout) ::
     &   eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),
     &   intent(inout) ::
     &   esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::
     &   rside       ! fraction of ice that melts laterally

      real (kind=dbl_kind), dimension(nx_block,ny_block),
     &   intent(inout) ::
     &   fresh       ! fresh water flux to ocean (kg/m^2/s)
     &,  fsalt       ! salt flux to ocean (kg/m^2/s)
     &,  fhocn       ! net heat flux to ocean (W/m^2)
     &,  fresh_hist  ! fresh water flux to ocean (kg/m^2/s)
     &,  fsalt_hist  ! salt flux to ocean (kg/m^2/s)
     &,  fhocn_hist  ! net heat flux to ocean (W/m^2)
     &,  meltl       ! lateral ice melt         (m/step-->cm/day)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j            ! horizontal indices
     &,  n               ! thickness category index
     &,  k               ! layer index
     &,  ilo,ihi,jlo,jhi ! beginning and end of physical domain
     &,  ij              ! horizontal index, combines i and j loops
     &,  icells          ! number of cells with aice > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) ::
     &   indxi, indxj    ! compressed indices for cells with aice > puny

      real (kind=dbl_kind) ::
     &   dfhocn      ! change in fhocn
     &,  dfresh      ! change in fresh
     &,  dfsalt      ! change in fsalt

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify grid cells with lateral melting.
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (rside(i,j) > c0) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo                  ! i
         enddo                  ! j

      !-----------------------------------------------------------------
      ! Melt the ice and increment fluxes.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            ! fluxes to coupler
            ! dfresh > 0, dfsalt > 0

            dfresh = (rhos*vsnon(i,j,n) + rhoi*vicen(i,j,n))
     &             * rside(i,j) / dt
            dfsalt = rhoi*vicen(i,j,n)*ice_ref_salinity*p001
     &             * rside(i,j) / dt

            fresh(i,j)      = fresh(i,j)      + dfresh
            fresh_hist(i,j) = fresh_hist(i,j) + dfresh

            fsalt(i,j)      = fsalt(i,j)      + dfsalt
            fsalt_hist(i,j) = fsalt_hist(i,j) + dfsalt

            ! history diagnostics
            meltl(i,j) = meltl(i,j) + vicen(i,j,n)*rside(i,j)

            ! state variables
            aicen(i,j,n) = aicen(i,j,n) * (c1 - rside(i,j))
            vicen(i,j,n) = vicen(i,j,n) * (c1 - rside(i,j))
            vsnon(i,j,n) = vsnon(i,j,n) * (c1 - rside(i,j))

         enddo                  ! ij

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               ! heat flux to coupler for ice melt (dfhocn < 0)

               dfhocn = eicen(i,j,ilyr1(n)+k-1)*rside(i,j) / dt
               fhocn(i,j)      = fhocn(i,j)      + dfhocn
               fhocn_hist(i,j) = fhocn_hist(i,j) + dfhocn

               ! ice energy
               eicen(i,j,ilyr1(n)+k-1) = eicen(i,j,ilyr1(n)+k-1)
     &                                 * (c1 - rside(i,j))
            enddo               ! ij
         enddo                  ! nilyr

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               ! heat flux to coupler for snow melt (dfhocn < 0)

               dfhocn = esnon(i,j,slyr1(n)+k-1)*rside(i,j) / dt
               fhocn(i,j)      = fhocn(i,j)      + dfhocn
               fhocn_hist(i,j) = fhocn_hist(i,j) + dfhocn

               ! snow energy
               esnon(i,j,slyr1(n)+k-1) = esnon(i,j,slyr1(n)+k-1)
     &                                 * (c1 - rside(i,j))
            enddo               ! ij
         enddo                  ! nslyr

      enddo  ! n

      end subroutine lateral_melt

!=======================================================================

      end module ice_therm_itd

!=======================================================================
