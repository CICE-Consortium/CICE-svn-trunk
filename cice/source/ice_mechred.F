c $Id: ice_mechred_cice.F,v 1.14 2004/02/09 17:58:46 lipscomb Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_mechred - driver for mechanical redestribution
!
! !DESCRIPTION:
!
! Driver for ice mechanical redistribution (ridging)
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
! Lipscomb, W. H., E. C. Hunke, W. Maslowski, and J. Jakacki, 2006: 
!  Ridging, strength, and stability in sea ice models, submitted 
!  to J. Geophys. Res. 
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
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: New options for participation and redistribution (WHL)
!
! !INTERFACE:
!
      module ice_mechred
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size
      use ice_constants
      use ice_fileunits
      use ice_itd, only: hin_max, ilyr1, slyr1, column_sum,
     &                   column_conservation_check, compute_tracers
!
!EOP
!
      implicit none
      save

!-----------------------------------------------------------------------
! Ridging parameters
!-----------------------------------------------------------------------

      integer (kind=int_kind) ::   ! defined in namelist 
     &   kstrength        ! 0 for simple Hibler (1979) formulation 
                          ! 1 for Rothrock (1975) pressure formulation 
     &,  krdg_partic      ! 0 for Thorndike et al. (1975) formulation 
                          ! 1 for exponential participation function 
     &,  krdg_redist      ! 0 for Hibler (1980) formulation 
                          ! 1 for exponential redistribution function 
 
      real (kind=dbl_kind), parameter :: 
     &   Cf = 17._dbl_kind   ! ratio of ridging work to PE change in ridging 
     &,  Cs = p25            ! fraction of shear energy contrbtng to ridging 
     &,  Cp = p5*gravit*(rhow-rhoi)*rhoi/rhow  ! proport const for PE 
     &,  fsnowrdg = p5           ! snow fraction that survives in ridging 
     &,  Gstar = 0.15_dbl_kind   ! max value of G(h) that participates 
                                 ! (krdg_partic = 0) 
     &,  astar = 0.05_dbl_kind   ! e-folding scale for G(h) participation 
                                 ! (krdg_partic = 1) 
     &,  maxraft = 1.0_dbl_kind  ! max value of hrmin - hi = max thickness 
                                 ! of ice that rafts (m) 
     &,  Hstar = 25._dbl_kind    ! determines mean thickness of ridged ice (m) 
                                 ! (krdg_redist = 0) 
                                 ! Flato & Hibler (1995) have Hstar = 100 
     &,  mu_rdg = 4.0_dbl_kind   ! gives e-folding scale of ridged ice (m^.5) 
                                 ! (krdg_redist = 1) 
     &,  Pstar = 2.75e4_dbl_kind ! constant in Hibler strength formula 
                                 ! (kstrength = 0) 
     &,  Cstar = 20._dbl_kind    ! constant in Hibler strength formula 
                                 ! (kstrength = 0) 

      logical (kind=log_kind), parameter ::
     &   l_conservation_check = .true.  ! if true, check conservation
                                        ! (useful for debugging)

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
! NOTE: This subroutine operates over a single block.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_ice (nx_block,    ny_block,
     &                      nghost,      dt,
!!     &                      Delt,        divu,
     &                      rdg_conv,    rdg_shear,
     &                      aicen,       trcrn,
     &                      vicen,       vsnon,
     &                      eicen,       esnon,
     &                      aice0,
     &                      trcr_depend, abort_flag,
     &                      dardg1dt,    dardg2dt,
     &                      dvirdgdt,    opening,
     &                      fresh,       fresh_hist,
     &                      fhocn,       fhocn_hist)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), intent(in) ::
     &   dt      ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::
!!!     &   Delt      ! (in the denominator of zeta, eta)            (1/s)
!!!     &,  divu      ! strain rate I component, velocity divergence (1/s)
     &   rdg_conv  ! normalized energy dissipation due to convergence (1/s)
     &,  rdg_shear ! normalized energy dissipation due to shear (1/s)
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(inout) ::
     &   aicen     ! concentration of ice
     &,  vicen     ! volume per unit area of ice          (m)
     &,  vsnon     ! volume per unit area of snow         (m)
 
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
     &   aice0     ! concentration of open water

      integer (kind=int_kind), dimension(ntrcr), intent(in) ::
     &   trcr_depend

      !  error flag
      logical (kind=log_kind), intent(out) ::
     &   abort_flag  ! if true, abort on return

      ! optional history fields
      real (kind=dbl_kind), dimension(nx_block,ny_block),
     &   intent(inout), optional ::
     &   dardg1dt    ! rate of fractional area loss by ridging ice (1/s)
     &,  dardg2dt    ! rate of fractional area gain by new ridges (1/s)
     &,  dvirdgdt    ! rate of ice volume ridged (m/s)
     &,  opening     ! rate of opening due to divergence/shear (1/s)
     &,  fresh       ! fresh water flux to ocean (kg/m^2/s)
     &,  fresh_hist  ! fresh water flux to ocean (kg/m^2/s)
     &,  fhocn       ! net heat flux to ocean (W/m^2)
     &,  fhocn_hist  ! net heat flux to ocean (W/m^2)
!
!EOP
!
      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   asum        ! sum of ice and open water area
     &,  aksum       ! ratio of area removed to area ridged
     &,  msnow_mlt   ! mass of snow added to ocean (kg m-2)
     &,  esnow_mlt   ! energy needed to melt snow in ocean (J m-2)
     &,  closing_net ! net rate at which area is removed    (1/s)
                     ! (ridging ice area - area of new ridges) / dt
     &,  divu_adv    ! divu as implied by transport scheme  (1/s)
     &,  opning      ! rate of opening due to divergence/shear
                     ! opning is a local variable;
                     ! opening is the history diagnostic variable
     &,  ardg1       ! fractional area loss by ridging ice
     &,  ardg2       ! fractional area gain by new ridges
     &,  virdg       ! ice volume ridged
     &,  aopen       ! area opening due to divergence/shear

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat) ::
     &   apartic          ! participation function; fraction of ridging
                          ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) ::
     &   hrmin            ! minimum ridge thickness
     &,  hrmax            ! maximum ridge thickness (krdg_redist = 0)
     &,  hrexp            ! ridge e-folding thickness (krdg_redist = 1) 
     &,  krdg             ! mean ridge thickness/thickness of ridging ice

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   vice_init, vice_final  ! ice volume summed over categories
     &,  vsno_init, vsno_final  ! snow volume summed over categories
     &,  eice_init, eice_final  ! ice energy summed over layers
     &,  esno_init, esno_final  ! snow energy summed over layers

      integer (kind=int_kind), parameter ::
     &   nitermax = 20    ! max number of ridging iterations

      integer (kind=int_kind) ::
     &   i,j              ! horizontal indices
     &,  n                ! thickness category index
     &,  niter            ! iteration counter
     &,  ilo,ihi,jlo,jhi  ! beg and end of physical domain

      real (kind=dbl_kind) ::
     &   dti              ! 1 / dt

      logical (kind=log_kind) ::
     &   iterate_ridging  ! if true, repeat the ridging
     &,  asum_error       ! flag for asum .ne. 1
     &,  conserv_err      ! flag for conservation error

      character (len=char_len) ::
     &   fieldid           ! field identifier

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      msnow_mlt(:,:) = c0
      esnow_mlt(:,:) = c0
      ardg1    (:,:) = c0
      ardg2    (:,:) = c0
      virdg    (:,:) = c0
      aopen    (:,:) = c0
      abort_flag = .false.

      !-----------------------------------------------------------------
      ! Compute area of ice plus open water before ridging.
      !-----------------------------------------------------------------
      call asum_ridging (nx_block,  ny_block,
     &                   nghost,    aicen,
     &                   aice0,     asum)

      !-----------------------------------------------------------------
      ! Compute the area opening and closing.
      !-----------------------------------------------------------------
      call ridge_prep (nx_block,  ny_block,
     &                 nghost,    dt,
!!     &                 Delt,      divu,
     &                 rdg_conv,  rdg_shear,
     &                 asum,      closing_net,
     &                 divu_adv,  opning)

      !-----------------------------------------------------------------
      ! Compute initial values of conserved quantities. 
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         call column_sum (nx_block,   ny_block,
     &                    nghost,     ncat,
     &                    vicen,      vice_init)

         call column_sum (nx_block,   ny_block,
     &                    nghost,     ncat,
     &                    vsnon,      vsno_init)

         call column_sum (nx_block,   ny_block,
     &                    nghost,     ntilyr,
     &                    eicen,      eice_init)

         call column_sum (nx_block,   ny_block,
     &                    nghost,     ntslyr,
     &                    esnon,      esno_init)

      endif            

      do niter = 1, nitermax

      !-----------------------------------------------------------------
      ! Compute the thickness distribution of ridging ice
      ! and various quantities associated with the new ridged ice.
      !-----------------------------------------------------------------

         call ridge_itd (nx_block,  ny_block,
     &                   nghost,     
     &                   aicen,     vicen,
     &                   aice0,     asum,
     &                   aksum,     apartic,
     &                   hrmin,     hrmax,
     &                   hrexp,     krdg)

      !-----------------------------------------------------------------
      ! Redistribute area, volume, and energy.
      !-----------------------------------------------------------------

         call ridge_shift (nx_block,  ny_block,
     &                     nghost,    dt,
     &                     aicen,     trcrn,
     &                     vicen,     vsnon,
     &                     eicen,     esnon,
     &                     aice0,     trcr_depend,
     &                     aksum,     apartic,
     &                     hrmin,     hrmax,
     &                     hrexp,     krdg,
     &                     closing_net, 
     &                     divu_adv,  opning,
     &                     ardg1,     ardg2,
     &                     virdg,     aopen,
     &                     msnow_mlt, esnow_mlt,
     &                     abort_flag)

      !-----------------------------------------------------------------
      ! Make sure the new area = 1.  If not (because the closing
      ! and opening rates were reduced above), prepare to ridge again
      ! with new rates.
      !-----------------------------------------------------------------

         call asum_ridging (nx_block,  ny_block,
     &                      nghost,    aicen,
     &                      aice0,     asum)

         call ridge_check (nx_block,   ny_block,
     &                     nghost,     dt,
     &                     asum,       closing_net,
     &                     divu_adv,   opning,
     &                     iterate_ridging)

      !-----------------------------------------------------------------
      ! If done, exit.  If not, prepare to ridge again.
      !-----------------------------------------------------------------

         if (iterate_ridging) then
            write(nu_diag,*) 'Repeat ridging, niter =', niter
         else
            exit
         endif

         if (niter == nitermax) then
            write(nu_diag,*) ' '
            write(nu_diag,*) 'Exceeded max number of ridging iterations'
            write(nu_diag,*) 'max =',nitermax
            abort_flag = .true.
         endif
            
      enddo                     ! niter

      !-----------------------------------------------------------------
      ! Compute final values of conserved quantities. 
      ! Check for conservation (allowing for snow thrown into ocean).
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         call column_sum (nx_block,   ny_block,
     &                    nghost,     ncat,
     &                    vicen,      vice_final)

         call column_sum (nx_block,   ny_block,
     &                    nghost,     ncat,
     &                    vsnon,      vsno_final)

         call column_sum (nx_block,   ny_block,
     &                    nghost,     ntilyr,
     &                    eicen,      eice_final)

         call column_sum (nx_block,   ny_block,
     &                    nghost,     ntslyr,
     &                    esnon,      esno_final)

         vsno_final(:,:) = vsno_final(:,:) + msnow_mlt(:,:)/rhos
         esno_final(:,:) = esno_final(:,:) + esnow_mlt(:,:)

         fieldid = 'vice, ridging'
         call column_conservation_check (nx_block,  ny_block,
     &                                   nghost,    fieldid,
     &                                   vice_init, vice_final,
     &                                   puny,      conserv_err)

         fieldid = 'vsno, ridging'
         call column_conservation_check (nx_block,  ny_block,
     &                                   nghost,    fieldid,
     &                                   vsno_init, vsno_final,
     &                                   puny,      conserv_err)

         fieldid = 'eice, ridging'
         call column_conservation_check (nx_block,  ny_block,
     &                                   nghost,    fieldid,
     &                                   eice_init, eice_final,
     &                                   puny*Lfresh*rhoi, 
     &                                   conserv_err)

         fieldid = 'esno, ridging'
         call column_conservation_check (nx_block,  ny_block,
     &                                   nghost,    fieldid,
     &                                   esno_init, esno_final,
     &                                   puny*Lfresh*rhos, 
     &                                   conserv_err)
         
         if (conserv_err) abort_flag = .true.

      endif                     ! l_conservation_check            

      !-----------------------------------------------------------------
      ! Compute ridging diagnostics.
      !-----------------------------------------------------------------

      dti = c1/dt

      if (present(dardg1dt))
     &     dardg1dt(:,:) = ardg1(:,:)*dti
      if (present(dardg2dt))
     &     dardg2dt(:,:) = ardg2(:,:)*dti
      if (present(dvirdgdt))
     &     dvirdgdt(:,:) = virdg(:,:)*dti
      if (present(opening))
     &     opening(:,:) = aopen(:,:)*dti

      !-----------------------------------------------------------------
      ! Update fresh water and heat fluxes due to snow melt.
      !-----------------------------------------------------------------

      if (present(fresh))
     &     fresh(:,:) = fresh(:,:) + msnow_mlt(:,:)*dti
      if (present(fresh_hist))
     &     fresh_hist(:,:) = fresh_hist(:,:) + msnow_mlt(:,:)*dti
      if (present(fhocn))
     &     fhocn(:,:) = fhocn(:,:) + esnow_mlt(:,:)*dti
      if (present(fhocn_hist))
     &     fhocn_hist(:,:) = fhocn_hist(:,:) + esnow_mlt(:,:)*dti

      !-----------------------------------------------------------------
      ! Check for fractional ice area > 1.
      !-----------------------------------------------------------------

      asum_error = .false.

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do j = jlo, jhi
      do i = ilo, ihi
         if (abs(asum(i,j) - c1) > puny) asum_error = .true.
      enddo
      enddo

      if (asum_error) then
        do j = jlo, jhi
        do i = ilo, ihi
          if (abs(asum(i,j) - c1) > puny) then ! insufficient ridging
            write(nu_diag,*) ' '
            write(nu_diag,*) 'Ridging error: total area > 1'
            write(nu_diag,*) 'i, j, area:', i, j, asum(i,j)
            write(nu_diag,*) 'n, aicen:'
            write(nu_diag,*)  0, aice0(i,j)
            do n = 1, ncat
               write(nu_diag,*) n, aicen(i,j,n)
            enddo
            abort_flag = .true.
          endif
        enddo
        enddo
      endif

      end subroutine ridge_ice

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
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine asum_ridging (nx_block, ny_block,
     &                         nghost,
     &                         aicen,    aice0,
     &                         asum)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(in) ::
     &   aicen          ! concentration of ice in each category

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aice0          ! concentration of open water

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out)::
     &   asum           ! sum of ice and open water area
!
!EOP
!
      integer (kind=int_kind) ::
     &   i, j, n

      !-----------------------------------------------------------------
      ! open water
      !-----------------------------------------------------------------
      asum(:,:) = aice0(:,:)

      !-----------------------------------------------------------------
      ! ice categories
      !-----------------------------------------------------------------

      do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do j = 1, ny_block
         do i = 1, nx_block
            asum(i,j) = asum(i,j) + aicen(i,j,n)
         enddo
         enddo
      enddo

      end subroutine asum_ridging

!=======================================================================
!BOP
!
! !ROUTINE: ridge_prep - preparation for ridging
!
! !DESCRIPTION: Initialize arrays, compute area of closing and opening
!
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_prep (nx_block,   ny_block,
     &                       nghost,     dt,
!!     &                       Delt,       divu,
     &                       rdg_conv,   rdg_shear,
     &                       asum,       closing_net,
     &                       divu_adv,   opning)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), intent(in) ::
     &   dt                  ! time step (s)

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::
!!     &   Delt      ! (in the denominator of zeta, eta)            (1/s)
!!     &,  divu      ! strain rate I component, velocity divergence (1/s)
     &   rdg_conv  ! normalized energy dissipation due to convergence (1/s)
     &,  rdg_shear ! normalized energy dissipation due to shear (1/s)
     &,  asum      ! sum of ice and open water area

      real (kind=dbl_kind), dimension(nx_block,ny_block),
     &   intent(out)::
     &   closing_net ! net rate at which area is removed    (1/s)
     &,  divu_adv  ! divu as implied by transport scheme  (1/s)
     &,  opning    ! rate of opening due to divergence/shear
!
!EOP
!
      real (kind=dbl_kind), parameter ::
     &   big = 1.0e+8_dbl_kind

      integer (kind=int_kind) ::
     &   i,j               ! horizontal indices
     &,  ilo,ihi,jlo,jhi   ! beginning and end of physical domain

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      closing_net(:,:) = c0
      divu_adv   (:,:) = c0
      opning     (:,:) = c0

      ! Set hin_max(ncat) to a big value to ensure that all ridged ice
      ! is thinner than hin_max(ncat).
      hin_max(ncat) = big

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do j = jlo, jhi
      do i = ilo, ihi

      !-----------------------------------------------------------------
      ! Compute the net rate of closing due to convergence
      ! and shear, based on Flato and Hibler (1995).
      !
      ! For the elliptical yield curve:
      !    rdg_conv  = -min (divu, 0)
      !    rdg_shear = (1/2) * (Delta - abs(divu))
      ! Note that the shear term also accounts for divergence.
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
!!     &        Cs*p5*(Delt(i,j)-abs(divu(i,j))) - min(divu(i,j),c0)
     &        Cs*rdg_shear(i,j) + rdg_conv(i,j)

      !-----------------------------------------------------------------
      ! Compute divu_adv, the divergence rate given by the transport/
      ! advection scheme, which may not be equal to divu as computed
      ! from the velocity field.
      !
      ! If divu_adv < 0, make sure the closing rate is large enough
      ! to give asum = 1.0 after ridging.
      !-----------------------------------------------------------------

         divu_adv(i,j) = (c1-asum(i,j)) / dt

         if (divu_adv(i,j) < c0)
     &        closing_net(i,j) = max(closing_net(i,j), -divu_adv(i,j))

      !-----------------------------------------------------------------
      ! Compute the (non-negative) opening rate that will give
      ! asum = 1.0 after ridging.
      !-----------------------------------------------------------------
         opning(i,j) = closing_net(i,j) + divu_adv(i,j)

      enddo
      enddo

      end subroutine ridge_prep

!=======================================================================
!BOP
!
! !ROUTINE: ridge_itd - thickness distribution of ridging and ridged ice
!
! !DESCRIPTION:
!
! Compute the thickness distribution of the ice and open water
! participating in ridging and of the resulting ridges.
!
! This version includes new options for ridging participation and
!  redistribution.
! The new participation scheme (krdg_partic = 1) improves stability
!  by increasing the time scale for large changes in ice strength.
! The new exponential redistribution function (krdg_redist = 1) improves 
!  agreement between ITDs of modeled and observed ridges.   
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! 2006: Changed subroutine name to ridge_itd
!       Added new options for ridging participation and redistribution.  
!
! !INTERFACE:
!
      subroutine ridge_itd (nx_block,    ny_block,
     &                      nghost,   
     &                      aicen,       vicen,
     &                      aice0,       asum,
     &                      aksum,       apartic,
     &                      hrmin,       hrmax,
     &                      hrexp,       krdg)
!
! !USES:
!
      use ice_work, only:  worka

! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(in) ::
     &   aicen     ! concentration of ice
     &,  vicen     ! volume per unit area of ice (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aice0     ! concentration of open water
     &,  asum      ! sum of ice and open water area

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out)::
     &   aksum            ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),
     &   intent(out) ::
     &   apartic          ! participation function; fraction of ridging
                          ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(out) ::
     &   hrmin            ! minimum ridge thickness
     &,  hrmax            ! maximum ridge thickness (krdg_redist = 0)
     &,  hrexp            ! ridge e-folding thickness (krdg_redist = 1) 
     &,  krdg             ! mean ridge thickness/thickness of ridging ice
!
!EOP
!
      integer (kind=int_kind) ::
     &   i,j              ! horizontal indices
     &,  n                ! thickness category index
     &,  ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      real (kind=dbl_kind), parameter ::
     &   Gstari   = c1/Gstar
     &,  astari   = c1/astar

      real (kind=dbl_kind), dimension(nx_block,ny_block,-1:ncat) ::
     &   Gsum             ! Gsum(n) = sum of areas in categories 0 to n

      real (kind=dbl_kind), dimension(nx_block,ny_block) ::
     &   expa0            ! exp(-aice0/asum)

      real (kind=dbl_kind) ::
     &   hi               ! ice thickness for each cat (m)
     &,  hieff            ! effective ice thickness (m) (krdg_redist = 2)
     &,  hrmean           ! mean ridge thickness (m)
     &,  xtmp             ! temporary variable

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      aksum(:,:) = c0

      do n = 0, ncat
         apartic(:,:,n) = c0
      enddo

      do n = 1, ncat
         hrmin (:,:,n) = c0
         hrmax (:,:,n) = c0
         hrexp (:,:,n) = c0
         krdg  (:,:,n) = c1
      enddo

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      !-----------------------------------------------------------------
      ! Compute the thickness distribution of ice participating in ridging.
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! First compute the cumulative thickness distribution function Gsum,
      !  where Gsum(n) is the fractional area in categories 0 to n.
      ! Ignore categories with very small areas.
      !-----------------------------------------------------------------

      Gsum(:,:,-1) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do j = jlo, jhi
      do i = ilo, ihi
         if (aice0(i,j) > puny) then
            Gsum(i,j,0) = aice0(i,j)
         else
            Gsum(i,j,0) = Gsum(i,j,-1)
         endif
      enddo  
      enddo

      do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n) > puny) then
               Gsum(i,j,n) = Gsum(i,j,n-1) + aicen(i,j,n)
            else
               Gsum(i,j,n) = Gsum(i,j,n-1)
            endif
         enddo
         enddo
      enddo

      ! normalize

      worka(:,:) = c1 / Gsum(:,:,ncat)
      do n = 0, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do j = jlo, jhi
         do i = ilo, ihi
            Gsum(i,j,n) = Gsum(i,j,n) * worka(i,j)
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute the participation function apartic; this is analogous to
      ! a(h) = b(h)g(h) as defined in Thorndike et al. (1975).
      !
      !              area lost from category n due to ridging/closing
      !  apartic(n) = ---------------------------------------------------
      !                  total area lost due to ridging/closing
      !
      !-----------------------------------------------------------------

      if (krdg_partic == 0) then  ! Thornike et al. 1975 formulation

      !-----------------------------------------------------------------
      ! Assume b(h) = (2/Gstar) * (1 - G(h)/Gstar).
      ! The expressions for apartic are found by integrating b(h)g(h) between
      ! the category boundaries.
      !-----------------------------------------------------------------

         do n = 0, ncat
            do j = jlo, jhi
            do i = ilo, ihi

               if (Gsum(i,j,n) < Gstar) then
                  apartic(i,j,n) = Gstari*(Gsum(i,j,n)-Gsum(i,j,n-1)) *
     &                 (c2 - (Gsum(i,j,n-1)+Gsum(i,j,n))*Gstari)
               elseif (Gsum(i,j,n-1) < Gstar) then
                  apartic(i,j,n) = Gstari * (Gstar-Gsum(i,j,n-1)) *
     &                 (c2 - (Gsum(i,j,n-1)+Gstar)*Gstari)
               endif

            enddo               ! i
            enddo               ! j
         enddo                  ! n

      elseif (krdg_partic==1) then   ! exponential dependence on G(h)

      !-----------------------------------------------------------------
      ! b(h) = exp(-G(h)/astar)
      ! apartic(n) = [exp(-G(n-1)/astar - exp(-G(n)/astar] / [1-exp(-1/astar)]. 
      ! The expression for apartic is found by integrating b(h)g(h)
      ! between the category boundaries.
      !-----------------------------------------------------------------

         ! precompute exponential terms using Gsum as work array

         xtmp = c1 / (c1 - exp(-astari))

         do n = -1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do j = jlo, jhi
            do i = ilo, ihi
               Gsum(i,j,n) = exp(-Gsum(i,j,n)*astari) * xtmp
            enddo               ! i
            enddo               ! j
         enddo                  ! n

         do n = 0, ncat
            do j = jlo, jhi
            do i = ilo, ihi
               apartic(i,j,n) = Gsum(i,j,n-1) - Gsum(i,j,n)
            enddo               ! i
            enddo               ! j
         enddo                  ! n

      endif                     ! krdg_partic

      !-----------------------------------------------------------------
      ! Compute variables related to ITD of ridged ice:
      ! 
      ! krdg = mean ridge thickness/ thickness of ridging ice
      ! hrmin  = min ridge thickness
      ! hrmax  = max ridge thickness (krdg_redist = 0)
      ! hrexp  = ridge e-folding scale (krdg_redist = 1)
      !----------------------------------------------------------------

      if (krdg_redist == 0) then  ! Hibler 1980 formulation

      !-----------------------------------------------------------------
      ! Assume ridged ice is uniformly distributed between hrmin and hrmax. 
      ! 
      ! This parameterization is a modified version of Hibler (1980). 
      ! In the original paper the min ridging thickness is hrmin = 2*hi, 
      !  and the max thickness is hrmax = 2*sqrt(hi*Hstar). 
      ! 
      ! Here the min thickness is hrmin = min(2*hi, hi+maxraft), 
      !  so thick ridging ice is not required to raft. 
      !
      !-----------------------------------------------------------------

         do n = 1, ncat
            do j = jlo, jhi
            do i = ilo, ihi

               if (aicen(i,j,n) > puny) then 
                  hi = vicen(i,j,n) / aicen(i,j,n) 
                  hrmin(i,j,n) = min(c2*hi, hi + maxraft) 
                  hrmax(i,j,n) = c2*sqrt(Hstar*hi) 
                  hrmax(i,j,n) = max(hrmax(i,j,n), hrmin(i,j,n)+puny) 
                  hrmean = p5 * (hrmin(i,j,n) + hrmax(i,j,n)) 
                  krdg(i,j,n) = hrmean / hi 
               endif 

            enddo               ! i
            enddo               ! j
         enddo                  ! n

      else               ! krdg_redist = 1; exponential redistribution
 
      !----------------------------------------------------------------- 
      ! The ridge ITD is a negative exponential: 
      ! 
      !  g(h) ~ exp[-(h-hrmin)/hrexp], h >= hrmin 
      ! 
      ! where hrmin is the minimum thickness of ridging ice and 
      ! hrexp is the e-folding thickness.
      ! 
      ! Here, assume as above that hrmin = min(2*hi, hi+maxraft).
      ! That is, the minimum ridge thickness results from rafting,
      !  unless the ice is thicker than maxraft.
      !
      ! Also, assume that hrexp = mu_rdg*sqrt(hi).
      ! The parameter mu_rdg is tuned to give e-folding scales mostly
      !  in the range 2-4 m as observed by upward-looking sonar.
      !
      ! Values of mu_rdg in the right column give ice strengths
      !  roughly equal to values of Hstar in the left column
      !  (within ~10 kN/m for typical ITDs):
      !
      !   Hstar     mu_rdg
      !
      !     25        3.0
      !     50        4.0
      !     75        5.0
      !    100        6.0
      !----------------------------------------------------------------- 

         do n = 1, ncat
            do j = jlo, jhi
            do i = ilo, ihi
               if (aicen(i,j,n) > puny) then
                  hi = vicen(i,j,n) / aicen(i,j,n)
                  hrmin(i,j,n) = min(c2*hi, hi + maxraft)
                  hrexp(i,j,n) = mu_rdg * sqrt(hi)
                  krdg(i,j,n) = (hrmin(i,j,n) + hrexp(i,j,n)) / hi
               endif
            enddo
            enddo
         enddo

      endif                     ! krdg_redist

      !----------------------------------------------------------------
      ! Compute aksum = net ice area removed / total area participating.
      ! For instance, if a unit area of ice with h = 1 participates in
      !  ridging to form a ridge with a = 1/3 and h = 3, then
      !  aksum = 1 - 1/3 = 2/3.
      !---------------------------------------------------------------- 

      do j = jlo, jhi
      do i = ilo, ihi
         aksum(i,j) = apartic(i,j,0)   ! area participating = area removed
      enddo
      enddo

      do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do j = jlo, jhi
         do i = ilo, ihi
            ! area participating > area removed
            aksum(i,j) = aksum(i,j) 
     &                 + apartic(i,j,n) * (c1 - c1/krdg(i,j,n)) 
         enddo
         enddo
      enddo

      end subroutine ridge_itd

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
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_shift (nx_block,    ny_block,
     &                        nghost,      dt,
     &                        aicen,       trcrn,
     &                        vicen,       vsnon,
     &                        eicen,       esnon,
     &                        aice0,       trcr_depend,     
     &                        aksum,       apartic,
     &                        hrmin,       hrmax, 
     &                        hrexp,       krdg,
     &                        closing_net,
     &                        divu_adv,    opning,
     &                        ardg1,       ardg2, 
     &                        virdg,       aopen,
     &                        msnow_mlt,   esnow_mlt,
     &                        abort_flag)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), intent(in) ::
     &   dt                  ! time step (s)

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::
     &   trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(inout) ::
     &   aicen     ! concentration of ice
     &,  vicen     ! volume per unit area of ice          (m)
     &,  vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),
     &   intent(inout) ::
     &   trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr),
     &   intent(inout) ::
     &   eicen     ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr),
     &   intent(inout) ::
     &   esnon     ! energy of melting for each snow layer (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aksum             ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),
     &   intent(in) ::
     &   apartic          ! participation function; fraction of ridging
                          ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(in) ::
     &   hrmin            ! minimum ridge thickness
     &,  hrmax            ! maximum ridge thickness (krdg_redist = 0)
     &,  hrexp            ! ridge e-folding thickness (krdg_redist = 1) 
     &,  krdg             ! mean ridge thickness/thickness of ridging ice

      real (kind=dbl_kind), dimension(nx_block,ny_block),
     &   intent(inout) ::
     &   ardg1       ! fractional area loss by ridging ice
     &,  ardg2       ! fractional area gain by new ridges
     &,  virdg       ! ice volume ridged (m)
     &,  aopen       ! area opened due to divergence/shear
     &,  closing_net ! net rate at which area is removed    (1/s)
     &,  divu_adv    ! divu as implied by transport scheme  (1/s)
     &,  opning      ! rate of opening due to divergence/shear (1/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block),
     &   intent(inout) ::
     &   msnow_mlt   ! mass of snow added to ocean (kg m-2)
     &,  esnow_mlt   ! energy needed to melt snow in ocean (J m-2)

      logical (kind=log_kind), intent(inout) ::
     &   abort_flag  ! if true, abort on return
!
!EOP
!
      integer (kind=int_kind) ::
     &   i,j               ! horizontal indices
     &,  n, nr             ! thickness category indices
     &,  k                 ! ice layer index
     &,  it                ! tracer index
     &,  ij                ! horizontal index, combines i and j loops
     &,  icells            ! number of cells with aicen > puny
     &,  ilo,ihi,jlo,jhi   ! beginning and end of physical domain

      integer (kind=int_kind), dimension (nx_block*ny_block) ::
     &   indxi, indxj      ! compressed indices

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) ::
     &   aicen_init        ! ice area before ridging
     &,  vicen_init        ! ice volume before ridging
     &,  vsnon_init        ! snow volume before ridging

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr) ::
     &   eicen_init        ! ice energy before ridging

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr) ::
     &   esnon_init        ! snow energy before ridging

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr,ncat) ::
     &   atrcrn            ! aicen*trcrn

      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   closing_gross     ! rate at which area removed, not counting
                           ! area of new ridges
     &,  afrac             ! fraction of category area ridged
     &,  ardg1n            ! area of ice ridged
     &,  ardg2n            ! area of new ridges
     &,  virdgn            ! ridging ice volume
     &,  vsrdgn            ! ridging snow volume
     &,  dhr               ! hrmax - hrmin
     &,  dhr2              ! hrmax^2 - hrmin^2
     &,  farea             ! fraction of new ridge area going to nr
     &,  fvol              ! fraction of new ridge volume going to nr

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr) ::
     &   eirdgn            ! ridging ice energy

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) ::
     &   esrdgn            ! ridging snow energy

      real (kind=dbl_kind) ::
     &   hi1               ! thickness of ridging ice
     &,  hexp              ! ridge e-folding thickness
     &,  hL, hR            ! left and right limits of integration
     &,  expL, expR        ! exponentials involving hL, hR
     &,  tmpfac            ! factor by which opening/closing rates are cut
     &,  wk1               ! work variable

      logical (kind=log_kind) ::
     &   neg_aice0         ! flag for aice0(i,j) < -puny
     &,  large_ardg        ! flag for afrac > 1

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      !-----------------------------------------------------------------
      ! Define variables equal to aicen*trcrn, vicen*trcrn, vsnon*trcrn
      !-----------------------------------------------------------------

      do n = 1, ncat
         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then ! ice area tracer
               atrcrn(:,:,it,n) = aicen(:,:,n)*trcrn(:,:,it,n)
            elseif (trcr_depend(it) == 1) then  ! ice volume tracer
               atrcrn(:,:,it,n) = vicen(:,:,n)*trcrn(:,:,it,n)
            elseif (trcr_depend(it) == 2) then  ! snow volume tracer
               atrcrn(:,:,it,n) = vsnon(:,:,n)*trcrn(:,:,it,n)
            endif
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Based on the ITD of ridging and ridged ice, convert the net
      !  closing rate to a gross closing rate.
      ! NOTE: 0 < aksum <= 1
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do j = jlo, jhi
      do i = ilo, ihi

         closing_gross(i,j) = closing_net(i,j) / aksum(i,j)

      !-----------------------------------------------------------------
      ! Reduce the closing rate if more than 100% of the open water
      ! would be removed.  Reduce the opening rate proportionately.
      !-----------------------------------------------------------------

         if (apartic(i,j,0) > c0) then
            wk1 = apartic(i,j,0) * closing_gross(i,j) * dt
            if (wk1 > aice0(i,j)) then
               tmpfac = aice0(i,j) / wk1
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
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do j = jlo, jhi
         do i = ilo, ihi

            if (aicen(i,j,n) > puny .and. apartic(i,j,n) > c0) then
               wk1 = apartic(i,j,n) * closing_gross(i,j) * dt
               if (wk1 > aicen(i,j,n)) then
                  tmpfac = aicen(i,j,n) / wk1
                  closing_gross(i,j) = closing_gross(i,j) * tmpfac
                  opning(i,j) = opning(i,j) * tmpfac
               endif
            endif

         enddo                  ! i
         enddo                  ! j
      enddo                     ! n

      !-----------------------------------------------------------------
      ! Compute change in open water area due to closing and opening.
      !-----------------------------------------------------------------

      neg_aice0 = .false.

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do j = jlo, jhi
      do i = ilo, ihi
         aice0(i,j) = aice0(i,j)
     &              - apartic(i,j,0)*closing_gross(i,j)*dt
     &              + opning(i,j)*dt
         if (aice0(i,j) < -puny) then
            neg_aice0 = .true.
         elseif (aice0(i,j) < c0) then    ! roundoff error
            aice0(i,j) = c0
         endif

         aopen(i,j) = aopen(i,j) + opning(i,j)*dt  ! optional diagnostic

      enddo
      enddo

      if (neg_aice0) then       ! there is a bug
         do j = jlo, jhi
         do i = ilo, ihi
            if (aice0(i,j) < -puny) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'Ridging error: aice0 < 0'
               write (nu_diag,*) 'i, j, aice0:', i, j, aice0(i,j)
            endif               ! aice0 < -puny
         enddo                  ! i
         enddo                  ! j
         abort_flag = .true.
      endif                     ! neg_aice0

      !-----------------------------------------------------------------
      ! Save initial state variables
      !-----------------------------------------------------------------

      aicen_init(:,:,:) = aicen(:,:,:)
      vicen_init(:,:,:) = vicen(:,:,:)
      vsnon_init(:,:,:) = vsnon(:,:,:)
      eicen_init(:,:,:) = eicen(:,:,:)
      esnon_init(:,:,:) = esnon(:,:,:)

      !-----------------------------------------------------------------
      ! Compute the area, volume, and energy of ice ridging in each
      !  category, along with the area of the resulting ridge.
      !-----------------------------------------------------------------

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify grid cells with nonzero ridging
      !-----------------------------------------------------------------

         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen_init(i,j,n) > puny .and. apartic(i,j,n) > c0
     &           .and. closing_gross(i,j) > c0) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo                  ! i
         enddo                  ! j

         large_ardg = .false.

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Compute area of ridging ice (ardg1n) and of new ridge (ardg2n).
      ! Make sure ridging fraction <=1.  (Roundoff errors can give
      !  ardg1 slightly greater than aicen.) 
      !-----------------------------------------------------------------

            ardg1n(i,j) = apartic(i,j,n)*closing_gross(i,j)*dt

            if (ardg1n(i,j) > aicen_init(i,j,n) + puny) then
               large_ardg = .true.
            else
               ardg1n(i,j) = min(aicen_init(i,j,n), ardg1n(i,j))
            endif

            ardg2n(i,j) = ardg1n(i,j) / krdg(i,j,n)
            afrac(i,j) = ardg1n(i,j) / aicen_init(i,j,n)

      !-----------------------------------------------------------------
      ! Subtract area, volume, and energy from ridging category n.
      ! Note: Tracer values are unchanged.
      !-----------------------------------------------------------------

            virdgn(i,j) = vicen_init(i,j,n) * afrac(i,j)
            vsrdgn(i,j) = vsnon_init(i,j,n) * afrac(i,j)

            aicen(i,j,n) = aicen(i,j,n) - ardg1n(i,j)
            vicen(i,j,n) = vicen(i,j,n) - virdgn(i,j)
            vsnon(i,j,n) = vsnon(i,j,n) - vsrdgn(i,j)

      !-----------------------------------------------------------------
      ! Increment ridging diagnostics
      !-----------------------------------------------------------------

            ardg1(i,j) = ardg1(i,j) + ardg1n(i,j)
            ardg2(i,j) = ardg2(i,j) + ardg2n(i,j)
            virdg(i,j) = virdg(i,j) + virdgn(i,j)

      !-----------------------------------------------------------------
      !  Place part of the snow lost by ridging into the ocean.
      !-----------------------------------------------------------------

            msnow_mlt(i,j) = msnow_mlt(i,j)
     &                     + rhos*vsrdgn(i,j)*(c1-fsnowrdg)

      !-----------------------------------------------------------------
      ! Compute quantities used to apportion ice among categories
      ! in the nr loop below
      !-----------------------------------------------------------------

            dhr(i,j)  = hrmax(i,j,n) - hrmin(i,j,n)
            dhr2(i,j) = hrmax(i,j,n) * hrmax(i,j,n)
     &                - hrmin(i,j,n) * hrmin(i,j,n)

         enddo                  ! ij

         if (large_ardg) then  ! there is a bug
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (ardg1(i,j) > aicen_init(i,j,n) + puny) then 
                  write (nu_diag,*) ' '
                  write (nu_diag,*) 'Ridging error: ardg > aicen'
                  write (nu_diag,*) 'i, j, n:', i, j, n
                  write (nu_diag,*) 'ardg, aicen:',
     &                               ardg1n(i,j), aicen_init(i,j,n)
                  abort_flag = .true.
               endif            ! afrac > 1 + puny
            enddo               ! ij
         endif                  ! large_ardg

      !-----------------------------------------------------------------
      ! Subtract ice energy from ridging category n.
      !-----------------------------------------------------------------

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               eirdgn(i,j,k) = eicen_init(i,j,ilyr1(n)+k-1) * afrac(i,j)
               eicen(i,j,ilyr1(n)+k-1) = eicen (i,j,ilyr1(n)+k-1)
     &                                 - eirdgn(i,j,k)
            enddo
         enddo

      !-----------------------------------------------------------------
      ! Subtract snow energy from ridging category n.
      ! Increment energy needed to melt snow in ocean.
      ! Note that esnow_mlt < 0; the ocean must cool to melt snow.
      !-----------------------------------------------------------------

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               esrdgn(i,j,k) = esnon_init(i,j,slyr1(n)+k-1) * afrac(i,j)
               esnon(i,j,slyr1(n)+k-1) = esnon (i,j,slyr1(n)+k-1)
     &                                 - esrdgn(i,j,k)
               esnow_mlt(i,j) = esnow_mlt(i,j)
     &                        + esrdgn(i,j,k)*(c1-fsnowrdg)
           enddo
         enddo

      !-----------------------------------------------------------------
      ! Subtract area- and volume-weighted tracers from category n.
      !-----------------------------------------------------------------

         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then ! ice area tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcrn(i,j,it,n) = atrcrn(i,j,it,n)
     &                             - ardg1n(i,j)*trcrn(i,j,it,n)
               enddo

            elseif (trcr_depend(it) == 1) then ! ice volume tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcrn(i,j,it,n) = atrcrn(i,j,it,n)
     &                             - virdgn(i,j)*trcrn(i,j,it,n)
               enddo

               elseif (trcr_depend(it) == 2) then ! snow volume tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  atrcrn(i,j,it,n) = atrcrn(i,j,it,n)
     &                             - vsrdgn(i,j)*trcrn(i,j,it,n)
               enddo
            endif               ! trcr_depend
         enddo                  ! ntrcr


      !-----------------------------------------------------------------
      ! Add area, volume, and energy of new ridge to each category nr.
      !-----------------------------------------------------------------

         do nr = 1, ncat

            if (krdg_redist == 0) then ! Hibler 1980 formulation

               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

      !-----------------------------------------------------------------
      ! Compute the fraction of ridged ice area and volume going to
      !  thickness category nr.
      !-----------------------------------------------------------------

                  if (hrmin(i,j,n) >= hin_max(nr) .or.
     &                hrmax(i,j,n) <= hin_max(nr-1)) then
                     hL = c0
                     hR = c0
                  else
                     hL = max (hrmin(i,j,n), hin_max(nr-1))
                     hR = min (hrmax(i,j,n), hin_max(nr))
                  endif

                  farea(i,j) = (hR-hL) / dhr(i,j)
                  fvol (i,j) = (hR*hR - hL*hL) / dhr2(i,j)

               enddo            ! ij

            else         ! krdg_redist = 1; 2005 exponential formulation

      !-----------------------------------------------------------------
      ! Compute the fraction of ridged ice area and volume going to
      !  thickness category nr.
      !-----------------------------------------------------------------

               if (nr < ncat) then

                  do ij = 1, icells
                     i = indxi(ij)
                     j = indxj(ij)

                     hi1  = hrmin(i,j,n)
                     hexp = hrexp(i,j,n)

                     if (hi1 >= hin_max(nr)) then
                        farea(i,j) = c0
                        fvol (i,j) = c0
                     else
                        hL = max (hi1, hin_max(nr-1))
                        hR = hin_max(nr)
                        expL = exp(-(hL-hi1)/hexp)
                        expR = exp(-(hR-hi1)/hexp)
                        farea(i,j) = expL - expR
                        fvol (i,j) = ((hL + hexp)*expL 
     &                              - (hR + hexp)*expR) / (hi1 + hexp)
                     endif
                  enddo         ! ij

               else             ! nr = ncat

                  do ij = 1, icells
                     i = indxi(ij)
                     j = indxj(ij)

                     hi1  = hrmin(i,j,n)
                     hexp = hrexp(i,j,n)

                     hL = max (hi1, hin_max(nr-1))
                     expL = exp(-(hL-hi1)/hexp)
                     farea(i,j) = expL
                     fvol (i,j) = (hL + hexp)*expL / (hi1 + hexp)

                  enddo

               endif            ! nr < ncat

            endif               ! krdg_redist

      !-----------------------------------------------------------------
      ! Transfer ice area, ice volume, and snow volume to category nr.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               aicen(i,j,nr) = aicen(i,j,nr) + farea(i,j)*ardg2n(i,j)
               vicen(i,j,nr) = vicen(i,j,nr) + fvol(i,j) *virdgn(i,j)
               vsnon(i,j,nr) = vsnon(i,j,nr)
     &                       + fvol(i,j)*vsrdgn(i,j)*fsnowrdg
            enddo

      !-----------------------------------------------------------------
      ! Transfer ice energy to category nr
      !-----------------------------------------------------------------
            do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  eicen(i,j,ilyr1(nr)+k-1) = eicen(i,j,ilyr1(nr)+k-1)
     &                                     + fvol(i,j)*eirdgn(i,j,k)
               enddo            ! ij
            enddo               ! k

      !-----------------------------------------------------------------
      ! Transfer snow energy to category nr
      !-----------------------------------------------------------------
            do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  esnon(i,j,slyr1(nr)+k-1) = esnon(i,j,slyr1(nr)+k-1)
     &                              + fvol(i,j)*esrdgn(i,j,k)*fsnowrdg
               enddo            ! ij
            enddo               ! k

      !-----------------------------------------------------------------
      ! Transfer area-weighted and volume-weighted tracers to category nr.
      ! Note: The global sum aicen*trcrn of ice area tracers 
      !       (trcr_depend = 0) is not conserved by ridging.
      !       However, ridging conserves the global sum of volume
      !       tracers (trcr_depend = 1 or 2).
      !-----------------------------------------------------------------

            do it = 1, ntrcr
               if (trcr_depend(it) == 0) then  ! ice area tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, icells
                     i = indxi(ij)
                     j = indxj(ij)
                     atrcrn(i,j,it,nr) = atrcrn(i,j,it,nr)
     &                          + farea(i,j)*ardg2n(i,j)*trcrn(i,j,it,n)

                  enddo
               elseif (trcr_depend(it) == 1) then ! ice volume tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, icells
                     i = indxi(ij)
                     j = indxj(ij)
                     atrcrn(i,j,it,nr) = atrcrn(i,j,it,nr)
     &                           + fvol(i,j)*virdgn(i,j)*trcrn(i,j,it,n)

                  enddo
               elseif (trcr_depend(it) == 2) then ! snow volume tracer
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
                  do ij = 1, icells
                     i = indxi(ij)
                     j = indxj(ij)
                     atrcrn(i,j,it,nr) = atrcrn(i,j,it,nr)
     &                  + fvol(i,j)*vsrdgn(i,j)*fsnowrdg*trcrn(i,j,it,n)

                  enddo
               endif            ! trcr_depend
            enddo               ! ntrcr

         enddo                  ! nr (new ridges)
      enddo                     ! n (ridging categories)


      !-----------------------------------------------------------------
      ! Compute new tracers
      !-----------------------------------------------------------------

      do n = 1, ncat
         call compute_tracers (nx_block,       ny_block,
     &                         trcr_depend,
     &                         atrcrn(:,:,:,n), aicen(:,:,  n),
     &                         vicen (:,:,  n), vsnon(:,:,  n),
     &                         trcrn(:,:,:,n))
      enddo

      end subroutine ridge_shift

!=======================================================================
!BOP
!
! !ROUTINE: ridge_check - check for ice area > 1
!
! !DESCRIPTION: Make sure ice area <=1.  If not, prepare to repeat ridging.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine ridge_check (nx_block,  ny_block,
     &                        nghost,    dt,
     &                        asum,      closing_net,
     &                        divu_adv,  opning,
     &                        iterate_ridging)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), intent(in) ::
     &   dt                  ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   asum             ! sum of ice and open water area

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(inout) ::
     &   closing_net ! net rate at which area is removed    (1/s)
     &,  divu_adv    ! divu as implied by transport scheme  (1/s)
     &,  opning      ! rate of opening due to divergence/shear

      logical (kind=log_kind), intent(out) ::
     &   iterate_ridging      ! if true, repeat the ridging
!
!EOP
!
      integer (kind=int_kind) ::
     &   i,j              ! horizontal indices
     &,  ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      iterate_ridging = .false.

      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      do j = jlo, jhi
      do i = ilo, ihi
         if (abs(asum(i,j) - c1) < puny) then
            closing_net(i,j) = c0
            opning     (i,j) = c0
         else
            iterate_ridging = .true.
            divu_adv(i,j) = (c1 - asum(i,j)) / dt
            closing_net(i,j) = max(c0, -divu_adv(i,j))
            opning(i,j) = max(c0, divu_adv(i,j))
         endif
      enddo
      enddo

      end subroutine ridge_check

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
! by ridging.
!
! See Rothrock (1975) and Hibler (1980).
!
! For simpler strength parameterization, see this reference:
! Hibler, W. D. III, 1979: A dynamic-thermodynamic sea ice model,
!  J. Phys. Oceanog., 9, 817-846.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      subroutine ice_strength (nx_block, ny_block,
     &                         nghost,
     &                         aice,     vice, 
     &                         aice0,    aicen, 
     &                         vicen,    strength)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) ::
     &   nx_block, ny_block  ! block dimensions
     &,  nghost              ! number of ghost cells

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) ::
     &   aice       ! concentration of ice
     &,  vice       ! volume per unit area of ice  (m)
     &,  aice0      ! concentration of open water

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),
     &   intent(in) ::
     &   aicen      ! concentration of ice
     &,  vicen      ! volume per unit area of ice  (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block),
     &   intent(out) ::
     &   strength   ! ice strength (N/m)
!
!EOP
!
! LOCAL VARIABLES
!
      real (kind=dbl_kind), dimension (nx_block,ny_block) ::
     &   asum             ! sum of ice and open water area
     &,  aksum            ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat) ::
     &   apartic          ! participation function; fraction of ridging
                          ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) ::
     &   hrmin            ! minimum ridge thickness
     &,  hrmax            ! maximum ridge thickness (krdg_redist = 0)
     &,  hrexp            ! ridge e-folding thickness (krdg_redist = 1) 
     &,  krdg             ! mean ridge thickness/thickness of ridging ice

      integer (kind=int_kind) ::
     &   i,j                 ! horizontal indices
     &,  n                   ! thickness category index
     &,  ilo,ihi,jlo,jhi     ! beg and end of physical domain

      real (kind=dbl_kind) ::
     &   hi                  ! ice thickness (m)
     &,  h2rdg               ! mean value of h^2 for new ridge
     &,  dh2rdg              ! change in mean value of h^2 per unit area
                             ! consumed by ridging 

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      ilo = 1 + nghost
      ihi = nx_block - nghost
      jlo = 1 + nghost
      jhi = ny_block - nghost

      strength(:,:) = c0

      if (kstrength == 1) then  ! Rothrock '75 formulation

      !-----------------------------------------------------------------
      ! Compute thickness distribution of ridging and ridged ice.
      !-----------------------------------------------------------------

         call asum_ridging (nx_block, ny_block,
     &                      nghost,
     &                      aicen,    aice0,
     &                      asum)

         call ridge_itd (nx_block,    ny_block,
     &                   nghost,   
     &                   aicen,       vicen,
     &                   aice0,       asum,
     &                   aksum,       apartic, 
     &                   hrmin,       hrmax,
     &                   hrexp,       krdg)

      !-----------------------------------------------------------------
      ! Compute ice strength based on change in potential energy,
      ! as in Rothrock (1975)
      !-----------------------------------------------------------------

         if (krdg_redist==0) then ! Hibler 1980 formulation

            do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do j = jlo, jhi
               do i = ilo, ihi
                  if (aicen(i,j,n) > puny .and. apartic(i,j,n) > c0)then
                     hi = vicen(i,j,n) / aicen(i,j,n)
                     h2rdg = p333 * (hrmax(i,j,n)**3 - hrmin(i,j,n)**3) 
     &                            / (hrmax(i,j,n) - hrmin(i,j,n)) 
                     dh2rdg = -hi*hi + h2rdg/krdg(i,j,n)
                     strength(i,j) = strength(i,j)
     &                             + apartic(i,j,n) * dh2rdg
                  endif         ! aicen > puny
               enddo            ! i
               enddo            ! j
            enddo               ! n

         elseif (krdg_redist==1) then ! exponential formulation

            do n = 1, ncat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do j = jlo, jhi
               do i = ilo, ihi
                  if (aicen(i,j,n) > puny .and. apartic(i,j,n) > c0)then
                     hi = vicen(i,j,n) / aicen(i,j,n)
                     h2rdg =    hrmin(i,j,n)*hrmin(i,j,n) 
     &                     + c2*hrmin(i,j,n)*hrexp(i,j,n)
     &                     + c2*hrexp(i,j,n)*hrexp(i,j,n)
                     dh2rdg = -hi*hi + h2rdg/krdg(i,j,n)
                     strength(i,j) = strength(i,j)
     &                             + apartic(i,j,n) * dh2rdg
                  endif
               enddo            ! i
               enddo            ! j
            enddo               ! n

         endif                  ! krdg_redist

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do j = jlo, jhi
         do i = ilo, ihi
            strength(i,j) = Cf * Cp * strength(i,j) / aksum(i,j)
                          ! Cp = (g/2)*(rhow-rhoi)*(rhoi/rhow)
                          ! Cf accounts for frictional dissipation
         enddo                  ! j
         enddo                  ! i

      else                      ! kstrength /= 1:  Hibler (1979) form

      !-----------------------------------------------------------------
      ! Compute ice strength as in Hibler (1979)
      !-----------------------------------------------------------------
         do j = jlo, jhi
         do i = ilo, ihi
            strength(i,j) = Pstar*vice(i,j)*exp(-Cstar*(c1-aice(i,j)))
         enddo                  ! j
         enddo                  ! i

      endif                     ! kstrength

      end subroutine ice_strength

!=======================================================================

      end module ice_mechred

!=======================================================================
