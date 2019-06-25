c $Id: ice_mechred_cice.F,v 1.16 2004/02/25 17:35:14 eclare Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_mechred_cice - mechanical redestribution and ice strength
!
! !DESCRIPTION:
!
! Ice mechanical redistribution (ridging) and strength computations 
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke (LANL)
!
! Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
!
! !INTERFACE:
!
      module ice_mechred_cice
!
! !USES:
! 
      use ice_model_size
      use ice_constants
      use ice_state
      use ice_itd
      use ice_grid
      use ice_fileunits
      use ice_domain
      use ice_calendar
      use ice_work, only:  worka
!
!EOP
!
      implicit none

!-----------------------------------------------------------------------
! Ridging parameters
!-----------------------------------------------------------------------

      integer (kind=int_kind) ::   ! defined in namelist
     &   kstrength        ! 1 for Rothrock 1975 pressure formulation
                          ! 0 for simple Hibler 1979 formulation
     &,  ksmooth = 0      ! 0 if ice strength is not smoothed
                          ! 1 if ice strength is smoothed
     
      real (kind=dbl_kind), parameter ::
     &   krdgmin = 1.1         ! min ridge thickness multiplier
     &,  Gstar = 0.15_dbl_kind ! max value of G(h) that participates
     &,  Hstar = 25._dbl_kind  ! determines mean thickness of ridged ice (m)
                               ! Flato & Hibler have Hstar = 100
     &,  Cf = 17._dbl_kind     ! ratio of ridging work to PE change in ridging
     &,  Cs = p25              ! fraction of shear energy contrbtng to ridging
     &,  Cp = p5 * gravit * (rhow-rhoi)*rhoi/rhow  ! proport const for PE
     &,  fsnowrdg = p5         ! snow fraction that survives in ridging

!-----------------------------------------------------------------------
!     Ridging diagnostic arrays for history files
!-----------------------------------------------------------------------
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::
     &   dardg1dt         ! rate of fractional area loss by ridging ice (1/s)
     &,  dardg2dt         ! rate of fractional area gain by new ridges (1/s)
     &,  dvirdgdt         ! rate of ice volume ridged (m/s)
     &,  opening          ! rate of opening due to divergence/shear (1/s)

!-----------------------------------------------------------------------
! Variables shared among ridging subroutines
!-----------------------------------------------------------------------
      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::     
     &   asum             ! sum of total ice and open water area
     &,  aksum            ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,0:ncat) ::     
     &   athorn           ! participation function; fraction of ridging/
                          !  closing associated w/ category n

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ncat) ::     
     &   hrmin            ! minimum ridge thickness
     &,  hrmax            ! maximum ridge thickness
     &,  krdg             ! mean ridge thickness/thickness of ridging ice 
                                       
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: ridge_ice - driver for mechanical redistribution
!
! !DESCRIPTION:
!
! Compute changes in the ice thickness distribution due to divergence
! and shear.
!
! See these references:
!
! Flato, G. M., and W. D. Hibler III, 1995: Ridging and strength
!  in modeling the thickness distribution of Arctic sea ice,
!  J. Geophys. Res., 100, 18,611-18,626.
!
! Hibler, W. D. III, 1980: Modeling a variable thickness sea ice
!  cover, Mon. Wea. Rev., 108, 1943-1973, 1980.
!
! Rothrock, D. A., 1975: The energetics of the plastic deformation of
!  pack ice by ridging, J. Geophys. Res., 80, 4514-4519.
!
! Thorndike, A. S., D. A. Rothrock, G. A. Maykut, and R. Colony, 
!  1975: The thickness distribution of sea ice, J. Geophys. Res., 
!  80, 4501-4513. 
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_ice (Delta, divu)
!
! !USES:
! 
      use ice_timers
      use ice_flux
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi), intent(in) ::
     &   Delta    ! (in the denominator of zeta, eta)             (1/s)
     &,  divu     ! strain rate I component, velocity divergence  (1/s)
!
!EOP
!
      integer (kind=int_kind), parameter ::
     &   nitermax = 20     ! max number of ridging iterations

      integer (kind=int_kind) ::
     &   i,j               ! horizontal indices
     &,  n                 ! thickness category index
     &,  niter             ! iteration counter

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi) ::
     &   closing_net       ! net rate at which area is removed    (1/s)
                           ! (ridging ice area - area of new ridges) / dt
     &,  divu_adv          ! divu as implied by transport scheme  (1/s)
     &,  opning            ! rate of opening due to divergence/shear
     &,  closing_gross     ! rate at which area removed, not counting
                           ! area of new ridges
     &,  msnow_mlt         ! mass of snow added to ocean (kg m-2)
     &,  esnow_mlt         ! energy needed to melt snow in ocean (J m-2)

      real (kind=dbl_kind) ::     
     &   w1                ! temporary variable
     &,  tmpfac            ! factor by which opening/closing rates are cut
     &,  dti               ! 1 / dt

      logical (kind=log_kind) ::
     &   iterate_ridging   ! if true, repeat the ridging

      real (kind=dbl_kind), parameter ::     
     &   big = 1.0e+8_dbl_kind

      logical (kind=log_kind) ::
     &   asum_error        ! flag for asum .ne. 1

      call ice_timer_start(6)  ! ridging 

      !-----------------------------------------------------------------
      ! Set hin_max(ncat) to a big value to ensure that all ridged ice 
      ! is thinner than hin_max(ncat).
      !-----------------------------------------------------------------
      hin_max(ncat) = big

      !-----------------------------------------------------------------
      ! Compute the ice strength, the thickness distribution of ridging 
      ! ice, and various quantities associated with the new ridged ice.
      !-----------------------------------------------------------------
      call ridge_prep

      do j = jlo, jhi
      do i = ilo, ihi

      !-----------------------------------------------------------------
      ! Initialize arrays.
      !-----------------------------------------------------------------

         msnow_mlt(i,j) = c0
         esnow_mlt(i,j) = c0
         dardg1dt(i,j) = c0
         dardg2dt(i,j) = c0
         dvirdgdt(i,j) = c0
         opening (i,j) = c0

      !-----------------------------------------------------------------
      ! Compute the net rate of closing due to convergence 
      ! and shear, based on Flato and Hibler (1995).
      ! 
      ! The energy dissipation rate is equal to the net closing rate
      ! times the ice strength.
      !
      ! NOTE: The NET closing rate is equal to the rate that open water 
      !  area is removed, plus the rate at which ice area is removed by 
      !  ridging, minus the rate at which area is added in new ridges.
      !  The GROSS closing rate is equal to the first two terms (open
      !  water closing and thin ice ridging) without the third term
      !  (thick, newly ridged ice).
      !-----------------------------------------------------------------

         closing_net(i,j) =
     &        Cs*p5*(Delta(i,j)-abs(divu(i,j))) - min(divu(i,j),c0)

      !-----------------------------------------------------------------
      ! Compute divu_adv, the divergence rate given by the transport/
      ! advection scheme, which may not be equal to divu as computed 
      ! from the velocity field.
      !
      ! If divu_adv < 0, make sure the closing rate is large enough
      ! to give asum = 1.0 after ridging.
      !-----------------------------------------------------------------

         divu_adv(i,j) = (c1-asum(i,j)) / dt  ! asum found in ridge_prep

         if (divu_adv(i,j) < c0) 
     &        closing_net(i,j) = max(closing_net(i,j), -divu_adv(i,j))

      !-----------------------------------------------------------------
      ! Compute the (non-negative) opening rate that will give 
      ! asum = 1.0 after ridging.
      !-----------------------------------------------------------------
         opning(i,j) = closing_net(i,j) + divu_adv(i,j)

      enddo
      enddo

      niter = 1                 ! iteration counter

  100 continue

      do j = jlo, jhi
      do i = ilo, ihi

      !-----------------------------------------------------------------
      ! Based on the ITD of ridging and ridged ice, convert the net
      !  closing rate to a gross closing rate.  
      ! NOTE: 0 < aksum <= 1
      !-----------------------------------------------------------------

         closing_gross(i,j) = closing_net(i,j) / aksum(i,j)

      !-----------------------------------------------------------------
      ! Reduce the closing rate if more than 100% of the open water 
      ! would be removed.  Reduce the opening rate proportionately.
      !-----------------------------------------------------------------

         if (aice0(i,j) > puny .and. athorn(i,j,0) > c0) then
            w1 = athorn(i,j,0) * closing_gross(i,j) * dt
            if (w1 > aice0(i,j)) then
               tmpfac = aice0(i,j) / w1
               closing_gross(i,j) = closing_gross(i,j) * tmpfac
               opning(i,j) = opning(i,j) * tmpfac
            endif
         endif

      enddo                     ! i
      enddo                     ! j

      !-----------------------------------------------------------------
      ! Reduce the closing rate if more than 100% of any ice category 
      ! would be removed.  Reduce the opening rate proportionately.
      !-----------------------------------------------------------------
      do n = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi

            if (aicen(i,j,n) > puny .and. athorn(i,j,n) > c0) then
               w1 = athorn(i,j,n) * closing_gross(i,j) * dt
               if (w1 > aicen(i,j,n)) then
                  tmpfac = aicen(i,j,n) / w1
                  closing_gross(i,j) = closing_gross(i,j) * tmpfac
                  opning(i,j) = opning(i,j) * tmpfac
               endif
            endif

         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Redistribute area, volume, and energy.
      !-----------------------------------------------------------------
      call ridge_shift (opning,    closing_gross, 
     &                  msnow_mlt, esnow_mlt)

      !-----------------------------------------------------------------
      ! Compute total area of ice plus open water after ridging.
      !-----------------------------------------------------------------
      call asum_ridging

      !-----------------------------------------------------------------
      ! Check whether asum = 1.  If not (because the closing and opening
      ! rates were reduced above), ridge again with new rates.
      !-----------------------------------------------------------------

      iterate_ridging = .false.

      do j = jlo, jhi
      do i = ilo, ihi
         if (abs(asum(i,j) - c1) < puny) then
            closing_net(i,j) = c0
            opning(i,j) = c0
         else
            iterate_ridging = .true.
            divu_adv(i,j) = (c1 - asum(i,j)) / dt
            closing_net(i,j) = max(c0, -divu_adv(i,j))
            opning(i,j) = max(c0, divu_adv(i,j))
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Repeat if necessary.
      ! NOTE: If strength smoothing is turned on, the ridging must be
      !       iterated globally because of the boundary update in the 
      !       smoothing.
      !-----------------------------------------------------------------

      niter = niter + 1

      if (iterate_ridging) then
         if (niter > nitermax) then
            write(nu_diag,*) 'istep1, my_task, nitermax =', 
     &                        istep1, my_task, nitermax
            call abort_ice('Exceeded max number of ridging iterations')
         endif
!         write(nu_diag,*) 'REPEAT RIDGING, istep1, my_task, niter =', 
!     &                                     istep1, my_task, niter
         call ridge_prep
         go to 100
      endif

      !-----------------------------------------------------------------
      ! Convert ridging rate diagnostics to correct units.
      !
      ! Update fresh water and heat fluxes due to snow melt.
      !-----------------------------------------------------------------

      dti = c1/dt

      asum_error = .false. 

      do j = jlo, jhi
      do i = ilo, ihi

         if (abs(asum(i,j) - c1) > puny) asum_error = .true.

         dardg1dt(i,j) = dardg1dt(i,j) * dti
         dardg2dt(i,j) = dardg2dt(i,j) * dti
         dvirdgdt(i,j) = dvirdgdt(i,j) * dti
         opening (i,j) = opening (i,j) * dti

         ! fresh water source for ocean
         fresh(i,j)      = fresh(i,j)      + msnow_mlt(i,j)*dti
         fresh_hist(i,j) = fresh_hist(i,j) + msnow_mlt(i,j)*dti
      
         ! heat sink for ocean
         fhnet(i,j)      = fhnet(i,j)      + esnow_mlt(i,j)*dti
         fhnet_hist(i,j) = fhnet_hist(i,j) + esnow_mlt(i,j)*dti

      enddo
      enddo

      if (asum_error) then
        do j = jlo, jhi
        do i = ilo, ihi

          if (abs(asum(i,j) - c1) > puny) then ! there is a bug
            write(nu_diag,*) ''
            write(nu_diag,*) 'Ridging error: total area =', asum(i,j)
            write(nu_diag,*) 'istep1, my_task, i, j:',
     &                  istep1, my_task, i, j
            write(nu_diag,*) 'n, aice, athorn:'
            write(nu_diag,*)  0, aice0(i,j), athorn(i,j,0)
            do n = 1, ncat
               write(nu_diag,*) n, aicen(i,j,n), athorn(i,j,n)
            enddo
            call abort_ice('ridging: total area must be <= 1')
          endif
        enddo
        enddo
      endif

      call ice_timer_stop(6)  ! ridging 

      end subroutine ridge_ice

!=======================================================================
!BOP
!
! !ROUTINE: ice_strength - compute ice strength
!
! !DESCRIPTION:
!
! Compute the strength of the ice pack, defined as the energy (J m-2) 
! dissipated per unit area removed from the ice pack under compression,
! and assumed proportional to the change in potential energy caused
! by ridging.  If ksmooth = 1, smooth the strength field to reduce noise.
!
! See Rothrock (1975) and Hibler (1980).
!
! For simpler strength parameterization, see this reference:
! Hibler, W. D. III, 1979: A dynamic-thermodynamic sea ice model,
!  J. Phys. Oceanog., 9, 817-846.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine ice_strength (kstrngth)
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   kstrngth    ! = 1 for Rothrock formulation, 0 for Hibler (1979)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i,j                 ! horizontal indices
     &,  n                   ! thickness category index

      real (kind=dbl_kind) ::     
     &   hi                  ! ice thickness (m)
     &,  w1                  ! temporary variable

      real (kind=dbl_kind), parameter ::     
     &   Pstar = 2.75e4_dbl_kind ! for Hibler strength formula
     &,  Chib  = 20._dbl_kind    ! for Hibler strength formula

      ! initialize
      strength = c0   

      !-----------------------------------------------------------------
      ! Compute thickness distribution of ridging and ridged ice.
      !-----------------------------------------------------------------
      call ridge_prep 

      if (kstrngth == 1) then

      !-----------------------------------------------------------------
      ! Compute ice strength based on change in potential energy,
      ! as in Rothrock (1975)
      !-----------------------------------------------------------------

         do n = 1, ncat
            do j = jlo, jhi
            do i = ilo, ihi

               if (aicen(i,j,n) > puny .and. athorn(i,j,n) > c0) then
                  hi = vicen(i,j,n) / aicen(i,j,n)
                  strength(i,j) = strength(i,j) - athorn(i,j,n) * hi**2 
                                ! PE loss from ridging ice
                  strength(i,j) = strength(i,j) 
     &                 + athorn(i,j,n)/krdg(i,j,n)
     &                 * p333 * (hrmax(i,j,n)**3 - hrmin(i,j,n)**3) 
     &                 / (hrmax(i,j,n)-hrmin(i,j,n))
                                ! PE gain from new ridge
               endif            ! aicen > puny

            enddo               ! i
            enddo               ! j
         enddo                  ! n

         do j = jlo, jhi
         do i = ilo, ihi

            strength(i,j) = Cf * Cp * strength(i,j) / aksum(i,j) 
                          ! Cp = (g/2)*(rhow-rhoi)*(rhoi/rhow)
                          ! Cf accounts for frictional dissipation

         enddo                  ! j
         enddo                  ! i

      else                      ! kstrngth ne 1:  Hibler (1979) form

      !-----------------------------------------------------------------
      ! Compute ice strength as in Hibler (1979)
      !-----------------------------------------------------------------
         do j = jlo, jhi
         do i = ilo, ihi
            strength(i,j) = Pstar*vice(i,j)*exp(-Chib*(c1-aice(i,j)))
         enddo                  ! j
         enddo                  ! i

      endif                     ! kstrngth

      !-----------------------------------------------------------------
      ! Smooth ice strength to reduce B-grid noise.
      ! Smoothing may be necessary to prevent large strength gradients 
      !  that lead to excessive ice velocities and thicknesses and may 
      !  crash the model (usually with a CFL violation in the horizontal
      !  transport).
      !-----------------------------------------------------------------

      if (ksmooth==1) then

         call bound(strength)

         do j = jlo, jhi
         do i = ilo, ihi
            if (asum(i,j) - aice0(i,j) > puny) then   ! ice exists
               worka(i,j) = (c4*strength(i,j) 
     &                        + strength(i-1,j) * hm(i-1,j)
     &                        + strength(i+1,j) * hm(i+1,j)
     &                        + strength(i,j-1) * hm(i,j-1)
     &                        + strength(i,j+1) * hm(i,j+1))
               w1 = c4 + hm(i-1,j) + hm(i+1,j)
     &                 + hm(i,j-1) + hm(i,j+1)
               worka(i,j) = worka(i,j) / w1
            else
               worka(i,j) = c0
            endif
         enddo
         enddo

         do j = jlo, jhi
         do i = ilo, ihi
            strength(i,j) = worka(i,j)
         enddo
         enddo

      endif

      end subroutine ice_strength

!=======================================================================
!BOP
!
! !ROUTINE: ridge_prep - preparation for ridging and strength calculations
!
! !DESCRIPTION:
!
! Compute the thickness distribution of the ice and open water 
! participating in ridging and of the resulting ridges.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_prep 
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) ::
     &   i,j              ! horizontal indices
     &,  n                ! thickness category index

      real (kind=dbl_kind), parameter ::     
     &   Gstari   = c1/Gstar    

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi,-1:ncat) ::     
     &   Gsum             ! Gsum(n) = sum of areas in categories 0 to n

      real (kind=dbl_kind) ::     
     &   hi               ! ice thickness for each cat (m)
     &,  hrmean           ! mean ridge thickness (m)

      !-----------------------------------------------------------------
      ! Zero out categories with very small areas
      !-----------------------------------------------------------------
      call zap_small_areas

      !-----------------------------------------------------------------
      ! Compute total area of ice plus open water.
      ! This is in general not equal to one because of divergence during
      !  transport.
      !-----------------------------------------------------------------
      call asum_ridging

      !-----------------------------------------------------------------
      ! Compute the cumulative thickness distribution function Gsum,
      ! where Gsum(n) is the fractional area in categories 0 to n.
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         worka(i,j)  = c1 / asum(i,j) ! used to normalize Gsum
         Gsum(i,j,-1) = c0
         Gsum(i,j,0) = aice0(i,j)*worka(i,j)
      enddo  
      enddo

      do n = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi
            Gsum(i,j,n) = Gsum(i,j,n-1) + aicen(i,j,n)*worka(i,j)
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute the participation function athorn; this is analogous to
      ! a(h) = b(h)g(h) as defined in Thorndike et al. (1975).
      !
      !              area lost from category n due to ridging/closing
      !  athorn(n) = ---------------------------------------------------
      !                  total area lost due to ridging/closing
      !
      ! Assume b(h) = (2/Gstar) * (1 - G(h)/Gstar). 
      ! The expressions for athorn are found by integrating b(h)g(h) between
      ! the category boundaries.
      !-----------------------------------------------------------------

      do n = 0, ncat
         do j = jlo, jhi
         do i = ilo, ihi

            if (Gsum(i,j,n) < Gstar) then
               athorn(i,j,n) = Gstari * (Gsum(i,j,n)-Gsum(i,j,n-1)) * 
     &              (c2 - (Gsum(i,j,n-1)+Gsum(i,j,n))*Gstari)
            elseif (Gsum(i,j,n-1) < Gstar) then
               athorn(i,j,n) = Gstari * (Gstar-Gsum(i,j,n-1)) * 
     &              (c2 - (Gsum(i,j,n-1)+Gstar)*Gstari)
            else
               athorn(i,j,n) = c0
            endif

         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Compute max and min ridged ice thickness for each ridging category.
      ! Assume ridged ice is uniformly distributed between hrmin and hrmax.
      ! 
      ! This parameterization is a modified version of Hibler (1980).
      ! The mean ridging thickness, hrmean, is proportional to hi^(0.5)
      !  and for very thick ridging ice must be >= krdgmin*hi
      !
      ! The minimum ridging thickness, hrmin, is equal to 2*hi 
      !  (i.e., rafting) and for very thick ridging ice is
      !  constrained by hrmin <= (hrmean + hi)/2.
      ! 
      ! The maximum ridging thickness, hrmax, is determined by
      !  hrmean and hrmin.
      !
      ! These modifications have the effect of reducing the ice strength
      ! (relative to the Hibler formulation) when very thick ice is
      ! ridging.
      !
      ! aksum = net area removed/ total area removed
      ! where total area removed = area of ice that ridges
      !         net area removed = total area removed - area of new ridges
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         aksum(i,j) = athorn(i,j,0)
      enddo
      enddo

      do n = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi

            if (aicen(i,j,n) > puny .and. athorn(i,j,n) > c0) then
               hi = vicen(i,j,n) / aicen(i,j,n)
               hrmean = max(sqrt(Hstar*hi), hi*krdgmin)
               hrmin(i,j,n) = min(2*hi, p5*(hrmean + hi))
               hrmax(i,j,n) = c2*hrmean - hrmin(i,j,n)
               krdg(i,j,n) = hrmean / hi
               aksum(i,j) = aksum(i,j) 
     &              + athorn(i,j,n) * (c1 - c1/krdg(i,j,n))
            else
               hrmin(i,j,n) = c0
               hrmax(i,j,n) = c0
               krdg (i,j,n) = c1
            endif

         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      end subroutine ridge_prep

!=======================================================================
!BOP
!
! !ROUTINE: ridge_shift - shift ridging ice among thickness categories
!
! !DESCRIPTION:
!
! Remove area, volume, and energy from each ridging category
! and add to thicker ice categories.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_shift (opning,    closing_gross, 
     &                        msnow_mlt, esnow_mlt)
!
! !USES:
!
      use ice_exit
! 
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(in) ::
     &   opning            ! rate of opening due to divergence/shear
     &,  closing_gross     ! rate at which area removed, not counting
                           ! area of new ridges

      real (kind=dbl_kind), dimension(ilo:ihi,jlo:jhi), intent(inout) ::
     &   msnow_mlt        ! mass of snow added to ocean (kg m-2)
     &,  esnow_mlt        ! energy needed to melt snow in ocean (J m-2)
!
!EOP
!
      integer (kind=int_kind) ::
     &   i,j               ! horizontal indices
     &,  n, n1, n2         ! thickness category indices
     &,  k                 ! ice layer index
     &,  ij                ! horizontal index, combines i and j loops
     &,  icells            ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(1:(ihi-ilo+1)*(jhi-jlo+1)) :: 
     &   indxi, indxj      ! compressed indices

      real (kind=dbl_kind), dimension (imt_local,jmt_local) ::     
     &   vice_init, vice_final  ! ice volume summed over categories
     &,  eice_init, eice_final  ! ice energy summed over layers

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ncat) ::     
     &   aicen_init        ! ice area before ridging
     &,  vicen_init        ! ice volume before ridging
     &,  vsnon_init        ! snow volume before ridging
     &,  esnon_init        ! snow energy before ridging

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,ntilay) ::     
     &   eicen_init        ! ice energy before ridging

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi) ::     
     &   afrac             ! fraction of category area ridged
     &,  ardg1             ! area of ice ridged
     &,  ardg2             ! area of new ridges
     &,  virdg             ! ice volume of ridging ice
     &,  vsrdg             ! snow volume of ridging ice
     &,  esrdg             ! snow energy of ridging ice
     &,  dhr               ! hrmax - hrmin
     &,  dhr2              ! hrmax^2 - hrmin^2
     &,  fvol              ! fraction of new ridge volume going to n2

      real (kind=dbl_kind), dimension (ilo:ihi,jlo:jhi,nilyr) ::     
     &   eirdg             ! ice energy of ridging ice

      real (kind=dbl_kind) ::     
     &   hL, hR            ! left and right limits of integration
     &,  farea             ! fraction of new ridge area going to n2

      character (len=char_len) ::
     &   fieldid           ! field identifier

      logical (kind=log_kind), parameter ::
     &   l_conservation_check = .true.  ! if true, check conservation 
                                        ! (useful for debugging)

      logical (kind=log_kind) ::
     &   neg_aice0         ! flag for aice0(i,j) < -puny
     &,  large_afrac       ! flag for afrac > 1

      !-----------------------------------------------------------------
      ! Compute quantities that ridging should conserve
      ! (not done for snow because snow may be dumped in ocean)
      !-----------------------------------------------------------------

      if (l_conservation_check) then
         call column_sum (ncat,   vicen, vice_init)
         call column_sum (ntilay, eicen, eice_init)
      endif

      !-----------------------------------------------------------------
      ! Compute change in open water area due to closing and opening.
      !-----------------------------------------------------------------

      neg_aice0 = .false.

      do j = jlo, jhi
      do i = ilo, ihi
         aice0(i,j) = aice0(i,j) 
     &              - athorn(i,j,0)*closing_gross(i,j)*dt
     &              + opning(i,j)*dt
         if (aice0(i,j) < -puny) then
            neg_aice0 = .true.
         elseif (aice0(i,j) < c0) then    ! roundoff error
            aice0(i,j) = c0
         endif
      enddo
      enddo

      if (neg_aice0) then       ! there is a bug
         do j = jlo, jhi
         do i = ilo, ihi
            if (aice0(i,j) < -puny) then 
               write (nu_diag,*) ''
               write (nu_diag,*) 'Ridging error: aice0 < 0'
               write (nu_diag,*) 'istep1, my_task, i, j:',
     &                            istep1, my_task, i, j
               write (nu_diag,*) 'aice0:', aice0(i,j)
               call abort_ice('ridging: aice0 must be >= 0')
            endif               ! aice0 < -puny
         enddo                  ! i
         enddo                  ! j
      endif                     ! neg_aice0

      !-----------------------------------------------------------------
      ! Save initial state variables
      !-----------------------------------------------------------------

      do n = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi
            aicen_init(i,j,n) = aicen(i,j,n)
            vicen_init(i,j,n) = vicen(i,j,n)
            vsnon_init(i,j,n) = vsnon(i,j,n)
            aicen_init(i,j,n) = aicen(i,j,n)
            esnon_init(i,j,n) = esnon(i,j,n)
         enddo
         enddo
      enddo
            
      do k = 1, ntilay
         do j = jlo, jhi
         do i = ilo, ihi
            eicen_init(i,j,k) = eicen(i,j,k)
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute the area, volume, and energy of ice ridging in each
      !  category, along with the area of the resulting ridge.
      !-----------------------------------------------------------------

      do n1 = 1, ncat

      !-----------------------------------------------------------------
      ! Identify grid cells with nonzero ridging
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen_init(i,j,n1) > puny .and. athorn(i,j,n1) > c0
     &           .and. closing_gross(i,j) > c0) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo                  ! i
         enddo                  ! j

         large_afrac = .false.

cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Compute area of ridging ice (ardg1) and of new ridge (ardg2)
      !-----------------------------------------------------------------

            ardg1(i,j) = athorn(i,j,n1)*closing_gross(i,j)*dt
            ardg2(i,j) = ardg1(i,j) / krdg(i,j,n1)

      !-----------------------------------------------------------------
      ! Compute ridging fraction (afrac), make sure afrac <=1 
      !-----------------------------------------------------------------

            afrac(i,j) = ardg1(i,j) / aicen_init(i,j,n1)
            if (afrac(i,j) > c1 + puny) then 
               large_afrac = .true.
            elseif (afrac(i,j) > c1) then  ! roundoff error
               afrac(i,j) = c1
            endif

      !-----------------------------------------------------------------
      ! Subtract area, volume, and energy from ridging category n1.
      ! (Ice energy in separate loop for vector friendliness)
      !-----------------------------------------------------------------

            virdg(i,j) = vicen_init(i,j,n1) * afrac(i,j)
            vsrdg(i,j) = vsnon_init(i,j,n1) * afrac(i,j)
            esrdg(i,j) = esnon_init(i,j,n1) * afrac(i,j)

            aicen(i,j,n1) = aicen(i,j,n1) - ardg1(i,j)
            vicen(i,j,n1) = vicen(i,j,n1) - virdg(i,j)
            vsnon(i,j,n1) = vsnon(i,j,n1) - vsrdg(i,j)
            esnon(i,j,n1) = esnon(i,j,n1) - esrdg(i,j)

      !-----------------------------------------------------------------
      ! Increment ridging diagnostics
      !-----------------------------------------------------------------

            dardg1dt(i,j) = dardg1dt(i,j) + ardg1(i,j)
            dardg2dt(i,j) = dardg2dt(i,j) + ardg2(i,j)
            dvirdgdt(i,j) = dvirdgdt(i,j) + virdg(i,j)
            opening(i,j)  = opening (i,j) + opning(i,j)*dt

      !-----------------------------------------------------------------
      !  Place part of the snow lost by ridging into the ocean. 
      !  Note that esnow_mlt < 0; the ocean must cool to melt snow.
      !  If the ocean temp = Tf already, new ice must grow.
      !  During the next time step, thermo_rates will determine whether
      !   the ocean cools or new ice grows.
      !-----------------------------------------------------------------
               
            msnow_mlt(i,j) = msnow_mlt(i,j) 
     &                     + rhos*vsrdg(i,j)*(c1-fsnowrdg)
            esnow_mlt(i,j) = esnow_mlt(i,j) 
     &                     + esrdg(i,j)*(c1-fsnowrdg)

      !-----------------------------------------------------------------
      ! Compute quantities used to apportion ice among categories
      ! in the n2 loop below
      !-----------------------------------------------------------------

            dhr(i,j)  = hrmax(i,j,n1) - hrmin(i,j,n1)
            dhr2(i,j) = hrmax(i,j,n1) * hrmax(i,j,n1) 
     &                - hrmin(i,j,n1) * hrmin(i,j,n1)

         enddo                  ! ij

         if (large_afrac) then  ! there is a bug
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (afrac(i,j) > c1 + puny) then 
                  write (nu_diag,*) ''
                  write (nu_diag,*) 'ardg > aicen'
                  write (nu_diag,*) 'istep1, my_task, i, j, n:',
     &                               istep1, my_task, i, j, n1
                  write (nu_diag,*) 'ardg, aicen_init:',
     &                 ardg1(i,j), aicen_init(i,j,n1)
                  call abort_ice ('ridging: ardg must be <= aicen')
               endif            ! afrac > 1 + puny
            enddo               ! if
         endif                  ! large_afrac


      !-----------------------------------------------------------------
      ! Subtract ice energy from ridging category n1. 
      !-----------------------------------------------------------------

         do k = 1, nilyr
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               
               eirdg(i,j,k) = eicen_init(i,j,ilyr1(n1)+k-1) * afrac(i,j)
               eicen(i,j,ilyr1(n1)+k-1) = eicen(i,j,ilyr1(n1)+k-1) 
     &                                  - eirdg(i,j,k)
            enddo
         enddo 

      !-----------------------------------------------------------------
      ! Add area, volume, and energy of new ridge to each category n2.
      !-----------------------------------------------------------------

         do n2 = 1, ncat
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

      !-----------------------------------------------------------------
      ! Compute the fraction of ridged ice area and volume going to 
      !  thickness category n2.
      ! Transfer area, volume, and energy accordingly.
      ! (Ice energy in separate loop for vector friendliness)
      !-----------------------------------------------------------------
               
               if (hrmin(i,j,n1) >= hin_max(n2) .or. 
     &             hrmax(i,j,n1) <= hin_max(n2-1)) then
                  hL = c0
                  hR = c0
               else
                  hL = max (hrmin(i,j,n1), hin_max(n2-1))
                  hR = min (hrmax(i,j,n1), hin_max(n2))
               endif

               ! fraction of ridged ice area and volume going to n2
               farea = (hR-hL) / dhr(i,j) 
               fvol(i,j) = (hR*hR - hL*hL) / dhr2(i,j)


               aicen(i,j,n2) = aicen(i,j,n2) + farea*ardg2(i,j)
               vicen(i,j,n2) = vicen(i,j,n2) + fvol(i,j)*virdg(i,j)
               vsnon(i,j,n2) = vsnon(i,j,n2) 
     &                       + fvol(i,j)*vsrdg(i,j)*fsnowrdg
               esnon(i,j,n2) = esnon(i,j,n2) 
     &                       + fvol(i,j)*esrdg(i,j)*fsnowrdg

            enddo               ! ij

      !-----------------------------------------------------------------
      ! Transfer ice energy to category n2 
      !-----------------------------------------------------------------
            do k = 1, nilyr
cdir$ ivdep      !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  eicen(i,j,ilyr1(n2)+k-1) = eicen(i,j,ilyr1(n2)+k-1) 
     &                                     + fvol(i,j)*eirdg(i,j,k)
               enddo            ! ij
            enddo               ! k

         enddo                  ! n2 (new ridges)            
      enddo                     ! n1 (ridging categories)

      !-----------------------------------------------------------------
      ! Check volume and energy conservation
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         call column_sum (ncat,   vicen, vice_final)
         fieldid = 'vice, ridging'
         call column_conservation_check (vice_init, vice_final, 
     &                                   puny,      fieldid)

         call column_sum (ntilay, eicen, eice_final)
         fieldid = 'eice, ridging'
         call column_conservation_check (eice_init, eice_final, 
     &                                   puny*Lfresh*rhoi, fieldid)

      endif

      end subroutine ridge_shift

!=======================================================================
!BOP
!
! !ROUTINE: asum_ridging - find total fractional area
!
! !DESCRIPTION:
!
! Find the total area of ice plus open water in each grid cell.
!
! This is similar to the aggregate_area subroutine except that the
! total area can be greater than 1, so the open water area is 
! included in the sum instead of being computed as a residual. 
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine asum_ridging
!
! !USES:
! 
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, n

      !-----------------------------------------------------------------
      ! open water
      !-----------------------------------------------------------------

      do j = jlo, jhi
      do i = ilo, ihi
         asum(i,j) = aice0(i,j)
      enddo
      enddo

      !-----------------------------------------------------------------
      ! ice categories
      !-----------------------------------------------------------------

      do n = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi
            asum(i,j) = asum(i,j) + aicen(i,j,n)
         enddo                  ! i
         enddo                  ! j
      enddo                     ! n 

      end subroutine asum_ridging
      
!=======================================================================

      end module ice_mechred_cice

!=======================================================================
