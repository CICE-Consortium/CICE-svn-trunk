!=========================================================================
!BOP
!
! !MODULE: ice_therm_bl99 - Bitz and Lipscomb 1999 thermodynamics
!
! !DESCRIPTION:
!
! Update ice and snow internal temperatures
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW
!          Elizabeth C. Hunke, LANL
!
! 2012: Split from ice_therm_vertical.F90
!
! !INTERFACE:
!
      module ice_therm_bl99
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size, only: nilyr, nslyr, max_ntrcr, n_aero, ncat
      use ice_constants
      use ice_fileunits, only: nu_diag
      use ice_therm_shared, only: conduct, calc_Tsfc, ferrmax, l_brine, hfrazilmin
!
!EOP
!
      implicit none
      save

      private
      public :: surface_fluxes, add_new_ice_bl99, temperature_changes

      real (kind=dbl_kind), parameter, public :: &
         hs_min = 1.e-4_dbl_kind    ! min snow thickness for computing Tsno (m)

      real (kind=dbl_kind), parameter :: &
         betak   = 0.13_dbl_kind, & ! constant in formula for k (W m-1 ppt-1)
         kimin   = 0.10_dbl_kind    ! min conductivity of saline ice (W m-1 deg-1)


!=======================================================================

      contains

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
      subroutine add_new_ice_bl99 (nx_block,  ny_block,   &
                              ntrcr,     icells,     &
                              indxi,     indxj,      &
                              tmask,     dt,         &
                              aicen,     trcrn,      &
                              vicen,                 &
                              aice0,     aice,       &
                              frzmlt,    frazil,     &
                              frz_onset, yday,       &
                              update_ocn_f,          &
                              fresh,     fsalt,      &
                              Tf,        sss,        &
                              salinz,    l_stop,     &
                              istop,     jstop,      &
                              fsice, &
                              flux_bio,  nbltrcr, &
                              ocean_bio)
!
! !USES:
!
      use ice_itd, only: hin_max, column_sum, &
                         column_conservation_check
      use ice_state, only: nt_Tsfc, nt_iage, nt_FY, nt_alvl, nt_vlvl, nt_aero, &
                           nt_sice, nt_qice, &
                           nt_apnd, tr_pond_cesm, tr_pond_lvl, tr_pond_topo, &
                           tr_iage, tr_FY, tr_lvl, tr_aero
      use ice_therm_shared, only: solve_Sin
      use ice_zbgc, only: add_new_ice_bgc
      use ice_zbgc_public, only: initbio_frac, rhosi

! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         icells                ! number of ice/ocean grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj          ! compressed i/j indices

      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
         tmask     ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aice  , & ! total concentration of ice
         frzmlt, & ! freezing/melting potential (W/m^2)
         Tf    , & ! freezing temperature (C)
         sss       ! sea surface salinity (ppt)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         aicen , & ! concentration of ice
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), &
         intent(inout) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         aice0     , & ! concentration of open water
         frazil    , & ! frazil ice growth        (m/step-->cm/day)
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt         ! salt flux to ocean (kg/m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout), optional :: &
         frz_onset ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(in), optional :: &
         yday      ! day of year

      real (kind=dbl_kind), dimension(nx_block,ny_block,nilyr+1), intent(in) :: &
         salinz     ! initial salinity profile
      
      logical (kind=log_kind), intent(in) :: &
         update_ocn_f ! if true, update fresh water and salt fluxes

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop    ! indices of grid cell where model aborts

      ! BGC
      integer (kind=int_kind), intent(in) :: &
         nbltrcr         ! number of biology tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         fsice      ! salt flux to ocean from prognostic S (kg/m^2/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbltrcr), &
         intent(inout) :: &
         flux_bio   ! tracer flux to ocean from biology (mmol/m^2/s) 
        
       real (kind=dbl_kind), dimension (nx_block,ny_block,nbltrcr), &
         intent(in) :: &
         ocean_bio   ! ocean concentration of biological tracer
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j         , & ! horizontal indices
         n            , & ! ice category index
         k            , & ! ice layer index
         it               ! aerosol tracer index

      real (kind=dbl_kind), dimension (icells) :: &
         ai0new       , & ! area of new ice added to cat 1
         vi0new       , & ! volume of new ice added to cat 1
         hsurp            ! thickness of new ice added to each cat

      real (kind=dbl_kind), dimension (icells) :: &
         vice1        , & ! starting volume of existing ice
         vice_init, vice_final  ! ice volume summed over categories

      real (kind=dbl_kind), dimension (icells,nilyr) :: &
         qi0              ! frazil ice enthalpy

      real (kind=dbl_kind) :: &
         fnew         , & ! heat flx to open water for new ice (W/m^2)
         hi0new       , & ! thickness of new ice
         hi0max       , & ! max allowed thickness of new ice
         qi0av        , & ! mean value of qi0 for new ice (J kg-1)
         vsurp        , & ! volume of new ice added to each cat
         vtmp         , & ! total volume of new and old ice
         area1        , & ! starting fractional area of existing ice
         alvl         , & ! starting level ice area
         rnilyr       , & ! real(nilyr)
         dfresh       , & ! change in fresh
         dfsalt           ! change in fsalt

      integer (kind=int_kind) :: &
         jcells, kcells     , & ! grid cell counters
         ij, m                  ! combined i/j horizontal indices

      integer (kind=int_kind), dimension (icells) :: &
         indxij2,indxij3    , & ! compressed i/j indices
         indxi2, indxj2     , &
         indxi3, indxj3

      character (len=char_len) :: &
         fieldid           ! field identifier

      ! BGC
      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat) :: &
         aicen_init, &    ! fractional area of ice
         vicen_init       ! volume per unit area of ice (m)

      real (kind=dbl_kind), dimension (icells) :: &
         vi0_init         ! volume of new ice

      l_stop = .false.
      istop = 0
      jstop = 0

      jcells = 0
      kcells = 0

      if (ncat > 1) then
         hi0max = hin_max(1)*0.9_dbl_kind  ! not too close to boundary
      else
         hi0max = bignum                   ! big number
      endif

      ! initial ice volume in each grid cell
      call column_sum (nx_block, ny_block,       &
                       icells,   indxi,   indxj, &
                       ncat,                     &
                       vicen,    vice_init)

      ! for bgc
      aicen_init(:,:,:) = aicen(:,:,:)
      vicen_init(:,:,:) = vicen(:,:,:)

      !-----------------------------------------------------------------
      ! Compute average enthalpy of new ice.
      !
      ! POP assumes new ice is fresh.  Otherwise, it would be better
      ! to do something like this:
      !  qi0(i,j,k) = -rhoi * (cp_ice*(Tmlt(k)-Tf(i,j))
      !             + Lfresh*(1.-Tmlt(k)/Tf(i,j)) - cp_ocn*Tmlt(k))
      !
      ! Compute the volume, area, and thickness of new ice.
      !-----------------------------------------------------------------

      rnilyr = real(nilyr,kind=dbl_kind)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         qi0av = c0
         do k = 1, nilyr
          qi0(ij,k) = -rhoi*Lfresh  ! note sign convention, qi < 0
          qi0av  = qi0av + qi0(ij,k)
         enddo
         qi0av = qi0av/rnilyr

         fnew = max (frzmlt(i,j), c0)   ! fnew > 0 iff frzmlt > 0
         vi0new(ij) = -fnew*dt / qi0av ! note sign convention, qi < 0
         vi0_init(ij) = vi0new(ij) ! for bgc

         ! increment ice volume
         vice_init(ij) = vice_init(ij) + vi0new(ij)

         ! history diagnostics
         frazil(i,j) = vi0new(ij)

         if (present(frz_onset) .and. present(yday)) then
            if (frazil(i,j) > puny .and. frz_onset(i,j) < puny) &
                 frz_onset(i,j) = yday
         endif

      !-----------------------------------------------------------------
      ! Update fresh water and salt fluxes.
      !
      ! NOTE: POP assumes fresh water and salt flux due to frzmlt > 0
      !       is NOT included in fluxes fresh and fsalt.
      !-----------------------------------------------------------------

         if (update_ocn_f) then
            dfresh = -rhoi*vi0new(ij)/dt 
            dfsalt = ice_ref_salinity*p001*dfresh

            fresh(i,j)      = fresh(i,j)      + dfresh
            fsalt(i,j)      = fsalt(i,j)      + dfsalt

            ! bgc
            fsice(i,j)      = fsice(i,j)  &
                  - rhosi*vi0new(ij)/dt*p001*sss(i,j)*salt_loss
            do k = 1, nbltrcr  ! only correct for dissolved tracers
               flux_bio(i,j,k) = flux_bio(i,j,k) &
                  - vi0new(ij)/dt*ocean_bio(i,j,k)*initbio_frac
            enddo
         endif

      !-----------------------------------------------------------------
      ! Decide how to distribute the new ice.
      !-----------------------------------------------------------------

         hsurp(ij)  = c0
         ai0new(ij) = c0

         if (vi0new(ij) > c0) then

            ! new ice area and thickness
            ! hin_max(0) < new ice thickness < hin_max(1)
            if (aice0(i,j) > puny) then
               hi0new = max(vi0new(ij)/aice0(i,j), hfrazilmin)
               if (hi0new > hi0max .and. aice0(i,j)+puny < c1) then
                  ! distribute excess volume over all categories (below)
                  hi0new = hi0max
                  ai0new(ij) = aice0(i,j)
                  vsurp      = vi0new(ij) - ai0new(ij)*hi0new
                  hsurp(ij)  = vsurp / aice(i,j)
                  vi0new(ij) = ai0new(ij)*hi0new
               else
                  ! put ice in a single category, with hsurp = 0
                  ai0new(ij) = vi0new(ij)/hi0new
               endif
            else                ! aice0 < puny
               hsurp(ij) = vi0new(ij)/aice(i,j) ! new thickness in each cat
               vi0new(ij) = c0
            endif               ! aice0 > puny
         endif                  ! vi0new > puny

      !-----------------------------------------------------------------
      ! Identify grid cells receiving new ice.
      !-----------------------------------------------------------------

         i = indxi(ij)
         j = indxj(ij)

         if (vi0new(ij) > c0) then  ! add ice to category 1
            jcells = jcells + 1
            indxi2(jcells) = i
            indxj2(jcells) = j
            indxij2(jcells) = ij
         endif

         if (hsurp(ij) > c0) then   ! add ice to all categories
            kcells = kcells + 1
            indxi3(kcells) = i
            indxj3(kcells) = j
            indxij3(kcells) = ij
         endif

      enddo                     ! ij

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
            m = indxij3(ij)

            vsurp = hsurp(m) * aicen(i,j,n)

            ! update ice age due to freezing (new ice age = dt)
            vtmp = vicen(i,j,n) + vsurp
            if (tr_iage .and. vtmp > puny) &
                trcrn(i,j,nt_iage,n) = &
               (trcrn(i,j,nt_iage,n)*vicen(i,j,n) + dt*vsurp) / vtmp

            if (tr_lvl .and. vicen(i,j,n) > puny) then
                trcrn(i,j,nt_vlvl,n) = &
               (trcrn(i,j,nt_vlvl,n)*vicen(i,j,n) + &
                trcrn(i,j,nt_alvl,n)*vsurp) / vtmp
            endif

            if (tr_aero) then
               do it = 1, n_aero
                  trcrn(i,j,nt_aero+2+4*(it-1),n) = &
                  trcrn(i,j,nt_aero+2+4*(it-1),n)*vicen(i,j,n) / vtmp
                  trcrn(i,j,nt_aero+3+4*(it-1),n) = &
                  trcrn(i,j,nt_aero+3+4*(it-1),n)*vicen(i,j,n) / vtmp
               enddo
            endif

            ! update category volumes
            vicen(i,j,n) = vtmp
         enddo                  ! ij

         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, kcells
               i = indxi3(ij)
               j = indxj3(ij)
               m = indxij3(ij)

               ! factor of nilyr cancels out
               vsurp = hsurp(m) * aicen(i,j,n)  ! note - save this above?
               vtmp = vicen(i,j,n) - vsurp      ! vicen is the new volume
               if (vicen(i,j,n) > c0) then
                  ! enthalpy
                  trcrn(i,j,nt_qice+k-1,n) = &
                 (trcrn(i,j,nt_qice+k-1,n)*vtmp + qi0(ij,k)*vsurp) / vicen(i,j,n)
                  if (.not. solve_Sin) then !bgc
                     ! Otherwise this is done for nt_bgc_S and mapped to nt_sice
                     trcrn(i,j,nt_sice+k-1,n) = &
                    (trcrn(i,j,nt_sice+k-1,n)*vtmp + salinz(i,j,k)*vsurp) / vicen(i,j,n) 
                  endif
               endif
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
         m = indxij2(ij)

         area1        = aicen(i,j,1)   ! save
         vice1(ij)    = vicen(i,j,1)   ! save
         aicen(i,j,1) = aicen(i,j,1) + ai0new(m)
         aice0(i,j)   = aice0(i,j)   - ai0new(m)
         vicen(i,j,1) = vicen(i,j,1) + vi0new(m)

         trcrn(i,j,nt_Tsfc,1) = &
            (trcrn(i,j,nt_Tsfc,1)*area1 + Tf(i,j)*ai0new(m))/aicen(i,j,1)
!bgc    ! But Tf(i,j) is not consistent with salinz(1) and qi0(1)
!bgc    ! trcrn(i,j,nt_Tsfc,1) = &
!bgc    !     (trcrn(i,j,nt_Tsfc,1)*area1 + Ti0(ij)*ai0new(m))/aicen(i,j,1)
         trcrn(i,j,nt_Tsfc,1) = min (trcrn(i,j,nt_Tsfc,1), c0)

         if (tr_FY) then
            trcrn(i,j,nt_FY,1) = &
           (trcrn(i,j,nt_FY,1)*area1 + ai0new(m))/aicen(i,j,1)
            trcrn(i,j,nt_FY,1) = min(trcrn(i,j,nt_FY,1), c1)
         endif

         if (vicen(i,j,1) > puny) then
            if (tr_iage) &
               trcrn(i,j,nt_iage,1) = &
              (trcrn(i,j,nt_iage,1)*vice1(ij) + dt*vi0new(m))/vicen(i,j,1)

            if (tr_aero) then
               do it = 1, n_aero
                  trcrn(i,j,nt_aero+2+4*(it-1),1) = &
                  trcrn(i,j,nt_aero+2+4*(it-1),1)*vice1(ij)/vicen(i,j,1)
                  trcrn(i,j,nt_aero+3+4*(it-1),1) = &
                  trcrn(i,j,nt_aero+3+4*(it-1),1)*vice1(ij)/vicen(i,j,1)
               enddo
            endif

            if (tr_lvl) then
                alvl = trcrn(i,j,nt_alvl,1) !echmod
                trcrn(i,j,nt_alvl,1) = &
               (trcrn(i,j,nt_alvl,1)*area1 + ai0new(m))/aicen(i,j,1)
                trcrn(i,j,nt_vlvl,1) = &
               (trcrn(i,j,nt_vlvl,1)*vice1(ij) + vi0new(m))/vicen(i,j,1)
            endif

            if (tr_pond_cesm .or. tr_pond_topo) then
               trcrn(i,j,nt_apnd,1) = &
               trcrn(i,j,nt_apnd,1)*area1/aicen(i,j,1)
            elseif (tr_pond_lvl) then
               if (trcrn(i,j,nt_alvl,1) > puny) then
                  trcrn(i,j,nt_apnd,1) = &
                  trcrn(i,j,nt_apnd,1) * alvl*area1 &
                                       / (trcrn(i,j,nt_alvl,1)*aicen(i,j,1))
               endif
            endif
         endif
      enddo                     ! ij

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, jcells
            i = indxi2(ij)
            j = indxj2(ij)
            m = indxij2(ij)

            if (vicen(i,j,1) > c0) then
            ! factor of nilyr cancels out
               ! enthalpy
               trcrn(i,j,nt_qice+k-1,1) = &
              (trcrn(i,j,nt_qice+k-1,1)*vice1(ij) &
                             + qi0(ij,k)*vi0new(m))/vicen(i,j,1) 
              
               ! salinity
              if (.NOT. solve_Sin) & !bgc
               trcrn(i,j,nt_sice+k-1,1) = &
              (trcrn(i,j,nt_sice+k-1,1)*vice1(ij) &
                            + salinz(i,j,k)*vi0new(m))/vicen(i,j,1)
            endif    !vicen
         enddo       !ij
      enddo          !k

      call column_sum (nx_block, ny_block,       &
                       icells,   indxi,   indxj, &
                       ncat,                     &
                       vicen,    vice_final)

      fieldid = 'vice, add_new_ice'
      call column_conservation_check (nx_block,  ny_block,      &
                                      icells,   indxi,   indxj, &
                                      fieldid,                  &
                                      vice_init, vice_final,    &
                                      puny,      l_stop,        &
                                      istop,     jstop)
      if (l_stop) return

      !-----------------------------------------------------------------
      ! Biogeochemistry
      !-----------------------------------------------------------------     
      call add_new_ice_bgc (nx_block,  ny_block,             &
                           icells,     jcells,     kcells,   &
                           indxi,      indxj,                &
                           indxi2,     indxj2,     indxij2,  &
                           indxi3,     indxj3,     indxij3,  &
                           aicen_init, vicen_init, vi0_init, &
                           aicen,      vicen,      vi0new,   &
                           ntrcr,      trcrn,      nbltrcr,  &
                           sss,        ocean_bio,            &
                           hsurp,      &
                           l_stop,     istop,      jstop)

      end subroutine add_new_ice_bl99

!=======================================================================
!BOP
!
! !ROUTINE: temperature_changes  - new vertical temperature profile
!
! !DESCRIPTION:
!
! Compute new surface temperature and internal ice and snow
! temperatures.  Include effects of salinity on sea ice heat
! capacity in a way that conserves energy (Bitz and Lipscomb, 1999).
!
! New temperatures are computed iteratively by solving a tridiagonal
! system of equations; heat capacity is updated with each iteration.
! Finite differencing is backward implicit.
!
! See Bitz, C.M., and W.H. Lipscomb, 1999:
! An energy-conserving thermodynamic model of sea ice,
! J. Geophys. Res., 104, 15,669-15,677.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine temperature_changes (nx_block, ny_block, &
                                      my_task,  istep1,   &
                                      dt,       icells,   & 
                                      indxi,    indxj,    &
                                      rhoa,     flw,      &
                                      potT,     Qa,       &
                                      shcoef,   lhcoef,   &
                                      fswsfc,   fswint,   &
                                      fswthrun, Sswabs,   &
                                      Iswabs,             &
                                      hilyr,    hslyr,    &
                                      qin,      Tin,      &
                                      qsn,      Tsn,      &
                                      Sin,                &
                                      Tsf,      Tbot,     &
                                      fsensn,   flatn,    &
                                      fswabsn,  flwoutn,  &
                                      fsurfn,             &
                                      fcondtopn,fcondbot, &
                                      einit,    l_stop,   &
                                      istop,    jstop,    &
                                      hin,      einex)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         my_task     , & ! task number (diagnostic only)
         istep1      , & ! time step index (diagnostic only)
         icells          ! number of cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         Tbot            ! ice bottom surface temperature (deg C)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         fswint      , & ! SW absorbed in ice interior below surface (W m-2)
         fswthrun        ! SW through ice to ocean         (W m-2)

      real (kind=dbl_kind), dimension (icells), intent(inout) :: &
         einex           ! excess energy from dqmat to ocean

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         einit       , & ! initial energy of melting (J m-2)
         hin             ! initial ice thickness

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(inout) :: &
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(inout) :: &
         Iswabs          ! SW radiation absorbed in ice layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout):: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fswabsn     , & ! shortwave absorbed by ice (W m-2)
         flwoutn         ! upward LW at surface (W m-2)

      real (kind=dbl_kind), dimension (icells), intent(out):: &
         fcondbot        ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), dimension (icells), &
         intent(inout):: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3)
         Tin         , & ! internal ice layer temperatures
         Sin             ! internal ice layer salinities

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(inout) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         Tsn             ! internal snow layer temperatures

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      integer (kind=int_kind), intent(inout) :: &
         istop, jstop    ! i and j indices of cell where model fails
!
!EOP
!
      integer (kind=int_kind), parameter :: &
         nitermax = 100, & ! max number of iterations in temperature solver
         nmat = nslyr + nilyr + 1  ! matrix dimension

      real (kind=dbl_kind), parameter :: &
         Tsf_errmax = 5.e-4_dbl_kind ! max allowed error in Tsf
                                     ! recommend Tsf_errmax < 0.01 K

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         k           , & ! ice layer index
         niter           ! iteration counter in temperature solver

      integer (kind=int_kind) :: &
         isolve          ! number of cells with temps not converged

      integer (kind=int_kind), dimension (icells) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), dimension (icells) :: &
         indxij          ! compressed 1D index for cells not converged

      logical (kind=log_kind), dimension (icells) :: &
         l_snow      , & ! true if snow temperatures are computed
         l_cold          ! true if surface temperature is computed

      real (kind=dbl_kind), dimension (:), allocatable :: &
         Tsf_start   , & ! Tsf at start of iteration
         dTsf        , & ! Tsf - Tsf_start
         dTi1        , & ! Ti1(1) - Tin_start(1)
         dfsurf_dT   , & ! derivative of fsurf wrt Tsf
         avg_Tsi     , & ! = 1. if new snow/ice temps avg'd w/starting temps
         enew            ! new energy of melting after temp change (J m-2)

      real (kind=dbl_kind), dimension (icells) :: &
         dTsf_prev   , & ! dTsf from previous iteration
         dTi1_prev   , & ! dTi1 from previous iteration
         dfsens_dT   , & ! deriv of fsens wrt Tsf (W m-2 deg-1)
         dflat_dT    , & ! deriv of flat wrt Tsf (W m-2 deg-1)
         dflwout_dT  , & ! deriv of flwout wrt Tsf (W m-2 deg-1)
         dt_rhoi_hlyr, & ! dt/(rhoi*hilyr)
         ferr            ! energy conservation error (W m-2)

      real (kind=dbl_kind), dimension (icells,nilyr) :: &
         Tin_init    , & ! Tin at beginning of time step
         Tin_start   , & ! Tin at start of iteration
         dTmat       , & ! Tin - matrix solution before limiting
         dqmat       , & ! associated enthalpy difference
         Tmlts           ! melting temp, -depressT * salinity

      real (kind=dbl_kind), dimension (icells,nslyr) :: &
         Tsn_init    , & ! Tsn at beginning of time step
         Tsn_start   , & ! Tsn at start of iteration
         etas            ! dt / (rho * cp * h) for snow layers

      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         etai        , & ! dt / (rho * cp * h) for ice layers
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs         , & ! rhs of tri-diagonal matrix equation
         Tmat            ! matrix output temperatures

      real (kind=dbl_kind), dimension(icells,nilyr+nslyr+1):: &
         kh              ! effective conductivity at interfaces (W m-2 deg-1)

      real (kind=dbl_kind) :: &
         ci          , & ! specific heat of sea ice (J kg-1 deg-1)
         avg_Tsf     , & ! = 1. if Tsf averaged w/Tsf_start, else = 0.
         Iswabs_tmp  , & ! energy to melt through fraction frac of layer
         Sswabs_tmp  , & ! same for snow
         dswabs      , & ! difference in swabs and swabs_tmp
         frac        , & ! fraction of layer that can be melted through
         dTemp           ! minimum temperature difference for absorption

      logical (kind=log_kind), dimension (icells) :: &
         converged      ! = true when local solution has converged

      logical (kind=log_kind) :: &
         all_converged  ! = true when all cells have converged

      logical (kind=log_kind) , dimension (icells,nilyr) :: &
         reduce_kh        ! reduce conductivity when T exceeds Tmlt

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      all_converged   = .false.

      do ij = 1, icells

         converged (ij) = .false.
         l_snow    (ij) = .false.
         l_cold    (ij) = .true.
         fcondbot  (ij) = c0
         dTsf_prev (ij) = c0
         dTi1_prev (ij) = c0
         dfsens_dT (ij) = c0
         dflat_dT  (ij) = c0
         dflwout_dT(ij) = c0  
         einex    (ij) = c0
         dt_rhoi_hlyr(ij) = dt / (rhoi*hilyr(ij))  ! hilyr > 0
         if (hslyr(ij) > hs_min/real(nslyr,kind=dbl_kind)) &
            l_snow(ij) = .true.
      enddo                     ! ij

      do k = 1, nslyr
         do ij = 1, icells
            Tsn_init (ij,k) = Tsn(ij,k) ! beginning of time step
            Tsn_start(ij,k) = Tsn(ij,k) ! beginning of iteration
            if (l_snow(ij)) then
               etas(ij,k) = dt/(rhos*cp_ice*hslyr(ij))
            else
               etas(ij,k) = c0
            endif
         enddo                  ! ij
      enddo                     ! k

      do k = 1, nilyr
         do ij = 1, icells
            Tin_init (ij,k) =  Tin(ij,k)   ! beginning of time step
            Tin_start(ij,k) =  Tin(ij,k)   ! beginning of iteration
            Tmlts    (ij,k) = -Sin(ij,k) * depressT
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute thermal conductivity at interfaces (held fixed during
      !  subsequent iterations).
      ! Ice and snow interfaces are combined into one array (kh) to
      !  simplify the logic.
      !-----------------------------------------------------------------

      call conductivity (nx_block, ny_block,         &
                         l_snow,   icells,           &
                         indxi,    indxj,    indxij, &
                         hilyr,    hslyr,            &
                         Tin,      kh,       Sin)

      !-----------------------------------------------------------------
      ! Check for excessive absorbed solar radiation that may result in
      ! temperature overshoots. Convergence is particularly difficult
      ! if the starting temperature is already very close to the melting 
      ! temperature and extra energy is added.   In that case, or if the
      ! amount of energy absorbed is greater than the amount needed to
      ! melt through a given fraction of a layer, we put the extra 
      ! energy into the surface.
      ! NOTE: This option is not available if the atmosphere model
      !       has already computed fsurf.  (Unless we adjust fsurf here)
      !-----------------------------------------------------------------
!mclaren: Should there be an if calc_Tsfc statement here then?? 

      frac = 0.9
      dTemp = 0.02_dbl_kind
      do k = 1, nilyr
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            Iswabs_tmp = c0 ! all Iswabs is moved into fswsfc
            if (Tin_init(ij,k) <= Tmlts(ij,k) - dTemp) then
               if (l_brine) then
                  ci = cp_ice - Lfresh * Tmlts(ij,k) / (Tin_init(ij,k)**2)
                  Iswabs_tmp = min(Iswabs(i,j,k), &
                     frac*(Tmlts(ij,k)-Tin_init(ij,k))*ci/dt_rhoi_hlyr(ij))
               else
                  ci = cp_ice
                  Iswabs_tmp = min(Iswabs(i,j,k), &
                     frac*(       -Tin_init(ij,k))*ci/dt_rhoi_hlyr(ij))
               endif
            endif
            if (Iswabs_tmp < puny) Iswabs_tmp = c0

            dswabs = min(Iswabs(i,j,k) - Iswabs_tmp, fswint(i,j))

            fswsfc(i,j)   = fswsfc(i,j) + dswabs
            fswint(i,j)   = fswint(i,j) - dswabs
            Iswabs(i,j,k) = Iswabs_tmp

         enddo
      enddo

      do k = 1, nslyr
         do ij = 1, icells
            if (l_snow(ij)) then
               i = indxi(ij)
               j = indxj(ij)

               Sswabs_tmp = c0
               if (Tsn_init(ij,k) <= -dTemp) then
                  Sswabs_tmp = min(Sswabs(i,j,k), &
                          -frac*Tsn_init(ij,k)/etas(ij,k))
               endif
               if (Sswabs_tmp < puny) Sswabs_tmp = c0

               dswabs = min(Sswabs(i,j,k) - Sswabs_tmp, fswint(i,j))

               fswsfc(i,j)   = fswsfc(i,j) + dswabs
               fswint(i,j)   = fswint(i,j) - dswabs
               Sswabs(i,j,k) = Sswabs_tmp

            endif
         enddo
      enddo

!lipscomb - This could be done in the shortwave module instead.
!           (Change zerolayer routine also)
      ! absorbed shortwave flux for coupler

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         fswabsn(i,j) = fswsfc(i,j) + fswint(i,j) + fswthrun(i,j)
      enddo

      !-----------------------------------------------------------------
      ! Solve for new temperatures.
      ! Iterate until temperatures converge with minimal energy error.
      !-----------------------------------------------------------------

      do niter = 1, nitermax

      !-----------------------------------------------------------------
      ! Identify cells, if any, where calculation has not converged.
      !-----------------------------------------------------------------

         if (all_converged) then  ! thermo calculation is done
            exit
         else                     ! identify cells not yet converged
            isolve = 0
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               if (.not.converged(ij)) then
                  isolve = isolve + 1
                  indxii(isolve) = i
                  indxjj(isolve) = j
                  indxij(isolve) = ij
               endif
            enddo               ! ij
         endif

      !-----------------------------------------------------------------
      ! Allocate and initialize
      !-----------------------------------------------------------------

         allocate(   sbdiag(isolve,nilyr+nslyr+1))
         allocate(     diag(isolve,nilyr+nslyr+1))
         allocate(   spdiag(isolve,nilyr+nslyr+1))
         allocate(      rhs(isolve,nilyr+nslyr+1))
         allocate(     Tmat(isolve,nilyr+nslyr+1))
         allocate(     etai(isolve,nilyr))
         allocate(Tsf_start(isolve))
         allocate(     dTsf(isolve))
         allocate(dfsurf_dT(isolve))
         allocate(  avg_Tsi(isolve))
         allocate(     enew(isolve))
         allocate(     dTi1(isolve))

         all_converged = .true.

         do ij = 1, isolve
            m = indxij(ij)
            converged(m)  = .true.
            dfsurf_dT(ij) = c0
            avg_Tsi  (ij) = c0
            enew     (ij) = c0
            einex    (m) = c0
         enddo

      !-----------------------------------------------------------------
      ! Update specific heat of ice layers.
      ! To ensure energy conservation, the specific heat is a function of
      ! both the starting temperature and the (latest guess for) the
      ! final temperature.
      !-----------------------------------------------------------------

         do k = 1, nilyr
            do ij = 1, isolve
               m = indxij(ij)
               i = indxii(ij)
               j = indxjj(ij)

               if (l_brine) then
                  ci = cp_ice - Lfresh*Tmlts(m,k) /  &
                                (Tin(m,k)*Tin_init(m,k))
               else
                  ci = cp_ice
               endif
               etai(ij,k) = dt_rhoi_hlyr(m) / ci

            enddo
         enddo

         if (calc_Tsfc) then

      !-----------------------------------------------------------------
      ! Update radiative and turbulent fluxes and their derivatives
      ! with respect to Tsf.
      !-----------------------------------------------------------------

          call surface_fluxes(nx_block,    ny_block,          &
                              isolve,      icells,            &
                              indxii,      indxjj,    indxij, &
                              Tsf,         fswsfc,            &
                              rhoa,        flw,               &
                              potT,        Qa,                &
                              shcoef,      lhcoef,            &
                              flwoutn,     fsensn,            &
                              flatn,       fsurfn,            &
                              dflwout_dT,  dfsens_dT,         &
                              dflat_dT,    dfsurf_dT)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Compute conductive flux at top surface, fcondtopn.
      ! If fsurfn < fcondtopn and Tsf = 0, then reset Tsf to slightly less
      !  than zero (but not less than -puny).
      !-----------------------------------------------------------------

            if (l_snow(m)) then
               fcondtopn(i,j) = kh(m,1) * (Tsf(m) - Tsn(m,1))
            else
               fcondtopn(i,j) = kh(m,1+nslyr) * (Tsf(m) - Tin(m,1))
            endif

            if (fsurfn(i,j) < fcondtopn(i,j)) &
                 Tsf(m) = min (Tsf(m), -puny)

      !-----------------------------------------------------------------
      ! Save surface temperature at start of iteration
      !-----------------------------------------------------------------

            Tsf_start(ij) = Tsf(m)

            if (Tsf(m) <= -puny) then
               l_cold(m) = .true.
            else
               l_cold(m) = .false.
            endif
          enddo                  ! ij

      !-----------------------------------------------------------------
      ! Compute elements of tridiagonal matrix.
      !-----------------------------------------------------------------

            call get_matrix_elements_calc_Tsfc &
                                  (nx_block, ny_block,         &
                                   isolve,   icells,           &
                                   indxii,   indxjj,   indxij, &
                                   l_snow,   l_cold,           &
                                   Tsf,      Tbot,             &
                                   fsurfn,   dfsurf_dT,        &
                                   Tin_init, Tsn_init,         &
                                   kh,       Sswabs,           &
                                   Iswabs,                     &
                                   etai,     etas,             &
                                   sbdiag,   diag,             &
                                   spdiag,   rhs)

         else
            call get_matrix_elements_know_Tsfc &
                                  (nx_block, ny_block,         &
                                   isolve,   icells,           &
                                   indxii,   indxjj,   indxij, &
                                   l_snow,   Tbot,             &
                                   Tin_init, Tsn_init,         &
                                   kh,       Sswabs,           &
                                   Iswabs,                     &
                                   etai,     etas,             &
                                   sbdiag,   diag,             &
                                   spdiag,   rhs,              &
                                   fcondtopn)
         endif  ! calc_Tsfc

      !-----------------------------------------------------------------
      ! Solve tridiagonal matrix to obtain the new temperatures.
      !-----------------------------------------------------------------

         call tridiag_solver (nx_block, ny_block, &
                              isolve,   icells,   &
                              indxii,   indxjj,   &
                              nmat,     sbdiag,   &
                              diag,     spdiag,   &
                              rhs,      Tmat)

      !-----------------------------------------------------------------
      ! Determine whether the computation has converged to an acceptable
      ! solution.  Five conditions must be satisfied:
      !
      !    (1) Tsf <= 0 C.
      !    (2) Tsf is not oscillating; i.e., if both dTsf(niter) and
      !        dTsf(niter-1) have magnitudes greater than puny, then
      !        dTsf(niter)/dTsf(niter-1) cannot be a negative number
      !        with magnitude greater than 0.5.  
      !    (3) abs(dTsf) < Tsf_errmax
      !    (4) If Tsf = 0 C, then the downward turbulent/radiative 
      !        flux, fsurfn, must be greater than or equal to the downward
      !        conductive flux, fcondtopn.
      !    (5) The net energy added to the ice per unit time must equal 
      !        the net change in internal ice energy per unit time,
      !        within the prescribed error ferrmax.
      !
      ! For briny ice (the standard case), Tsn and Tin are limited
      !  to prevent them from exceeding their melting temperatures.
      !  (Note that the specific heat formula for briny ice assumes
      !  that T < Tmlt.)  
      ! For fresh ice there is no limiting, since there are cases
      !  when the only convergent solution has Tsn > 0 and/or Tin > 0.
      !  Above-zero temperatures are then reset to zero (with melting 
      !  to conserve energy) in the thickness_changes subroutine.
      !-----------------------------------------------------------------

         if (calc_Tsfc) then

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          do ij = 1, isolve
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Reload Tsf from matrix solution
      !-----------------------------------------------------------------

            if (l_cold(m)) then
               if (l_snow(m)) then
                  Tsf(m) = Tmat(ij,1)
               else
                  Tsf(m) = Tmat(ij,1+nslyr)
               endif
            else                ! melting surface
               Tsf(m) = c0
            endif

      !-----------------------------------------------------------------
      ! Initialize convergence flag (true until proven false), dTsf,
      !  and temperature-averaging coefficients.
      ! Average only if test 1 or 2 fails.
      ! Initialize energy.
      !-----------------------------------------------------------------

            dTsf(ij) = Tsf(m) - Tsf_start(ij)
            avg_Tsf  = c0

      !-----------------------------------------------------------------
      ! Condition 1: check for Tsf > 0
      ! If Tsf > 0, set Tsf = 0, then average Tsn and Tin to force
      ! internal temps below their melting temps.
      !-----------------------------------------------------------------

            if (Tsf(m) > puny) then
               Tsf(m) = c0
               dTsf(ij) = -Tsf_start(ij)
               if (l_brine) avg_Tsi(ij) = c1   ! avg with starting temp
               converged(m) = .false.
               all_converged = .false.

      !-----------------------------------------------------------------
      ! Condition 2: check for oscillating Tsf
      ! If oscillating, average all temps to increase rate of convergence.
      !-----------------------------------------------------------------

            elseif (niter > 1 &                ! condition (2)
              .and. Tsf_start(ij) <= -puny &
              .and. abs(dTsf(ij)) > puny &
              .and. abs(dTsf_prev(m)) > puny &
              .and. -dTsf(ij)/(dTsf_prev(m)+puny*puny) > p5) then

               if (l_brine) then ! average with starting temp
                  avg_Tsf  = c1    
                  avg_Tsi(ij) = c1
               endif
               dTsf(ij) = p5 * dTsf(ij)
               converged(m) = .false.
               all_converged = .false.
            endif

!!!            dTsf_prev(m) = dTsf(ij)

      !-----------------------------------------------------------------
      ! If condition 2 failed, average new surface temperature with
      !  starting value.
      !-----------------------------------------------------------------
            Tsf(m)  = Tsf(m) &
                      + avg_Tsf * p5 * (Tsf_start(ij) - Tsf(m))

          enddo  ! ij

         endif   ! calc_Tsfc

         do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, isolve
               m = indxij(ij)

      !-----------------------------------------------------------------
      ! Reload Tsn from matrix solution
      !-----------------------------------------------------------------

               if (l_snow(m)) then
                  Tsn(m,k) = Tmat(ij,k+1)
               else
                  Tsn(m,k) = c0
               endif
               if (l_brine) Tsn(m,k) = min(Tsn(m,k), c0)

      !-----------------------------------------------------------------
      ! If condition 1 or 2 failed, average new snow layer
      !  temperatures with their starting values.
      !-----------------------------------------------------------------
               Tsn(m,k) = Tsn(m,k) &
                         + avg_Tsi(ij)*p5*(Tsn_start(m,k)-Tsn(m,k))

      !-----------------------------------------------------------------
      ! Compute qsn and increment new energy.
      !-----------------------------------------------------------------
               qsn(m,k) = -rhos * (Lfresh - cp_ice*Tsn(m,k))
               enew(ij) = enew(ij) + hslyr(m) * qsn(m,k)

               Tsn_start(m,k) = Tsn(m,k) ! for next iteration

            enddo               ! ij
         enddo                  ! nslyr

         dTmat(:,:) = c0
         dqmat(:,:) = c0
         reduce_kh(:,:) = .false.
         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, isolve
               m = indxij(ij)

      !-----------------------------------------------------------------
      ! Reload Tin from matrix solution
      !-----------------------------------------------------------------

               Tin(m,k) = Tmat(ij,k+1+nslyr)

               if (l_brine .and. Tin(m,k) > Tmlts(m,k) - puny) then
                  dTmat(m,k) = Tin(m,k) - Tmlts(m,k)
!echmod: return this energy to the ocean
                  dqmat(m,k) = rhoi * dTmat(m,k) &
                             * (cp_ice - Lfresh * Tmlts(m,k)/Tin(m,k)**2)
! use this for the case that Tmlt changes by an amount dTmlt=Tmltnew-Tmlt(k)
!                             + rhoi * dTmlt &
!                             * (cp_ocn - cp_ice + Lfresh/Tin(m,k))
                  Tin(m,k) = Tmlts(m,k)
                  reduce_kh(m,k) = .true.
               endif

      !-----------------------------------------------------------------
      ! Condition 2b: check for oscillating Tin(1)
      ! If oscillating, average all ice temps to increase rate of convergence.
      !-----------------------------------------------------------------

               if (k==1 .and. .not.calc_Tsfc) then
                  dTi1(ij) = Tin(m,k) - Tin_start(m,k)

                  if (niter > 1 &                    ! condition 2b    
                      .and. abs(dTi1(ij)) > puny &
                      .and. abs(dTi1_prev(m)) > puny &
                      .and. -dTi1(ij)/(dTi1_prev(m)+puny*puny) > p5) then

                     if (l_brine) avg_Tsi(ij) = c1
                     dTi1(ij) = p5 * dTi1(ij)
                     converged(m) = .false.
                     all_converged = .false.
                  endif
                  dTi1_prev(m) = dTi1(ij)
               endif   ! k = 1 .and. calc_Tsfc = F

      !-----------------------------------------------------------------
      ! If condition 1 or 2 failed, average new ice layer
      !  temperatures with their starting values.
      !-----------------------------------------------------------------
               Tin(m,k) = Tin(m,k) &
                         + avg_Tsi(ij)*p5*(Tin_start(m,k)-Tin(m,k))

      !-----------------------------------------------------------------
      ! Compute qin and increment new energy.
      !-----------------------------------------------------------------
               if (l_brine) then
                  qin(m,k) = -rhoi * (cp_ice*(Tmlts(m,k)-Tin(m,k)) &
                                      + Lfresh*(c1-Tmlts(m,k)/Tin(m,k)) &
                                      - cp_ocn*Tmlts(m,k))
               else
                  qin(m,k) = -rhoi * (-cp_ice*Tin(m,k) + Lfresh)
               endif
               enew(ij) = enew(ij) + hilyr(m) * (qin(m,k) + dqmat(m,k))
               einex(m) = einex(m) + hilyr(m) * dqmat(m,k)

               Tin_start(m,k) = Tin(m,k) ! for next iteration

            enddo               ! ij
         enddo                  ! nilyr

         if (calc_Tsfc) then

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
          do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Condition 3: check for large change in Tsf
      !-----------------------------------------------------------------

            if (abs(dTsf(ij)) > Tsf_errmax) then
               converged(m) = .false.
               all_converged = .false.
            endif

      !-----------------------------------------------------------------
      ! Condition 4: check for fsurfn < fcondtopn with Tsf > 0
      !-----------------------------------------------------------------

            fsurfn(i,j) = fsurfn(i,j) + dTsf(ij)*dfsurf_dT(ij)
            if (l_snow(m)) then
               fcondtopn(i,j) = kh(m,1) * (Tsf(m)-Tsn(m,1))
            else
               fcondtopn(i,j) = kh(m,1+nslyr) * (Tsf(m)-Tin(m,1))
            endif

            if (Tsf(m) > -puny .and. fsurfn(i,j) < fcondtopn(i,j)) then
               converged(m) = .false.
               all_converged = .false.
            endif

            dTsf_prev(m) = dTsf(ij)

          enddo                  ! ij
         endif                   ! calc_Tsfc

      !-----------------------------------------------------------------
      ! Condition 5: check for energy conservation error
      ! Change in internal ice energy should equal net energy input.
      !-----------------------------------------------------------------

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

            fcondbot(m) = kh(m,1+nslyr+nilyr) * &
                          (Tin(m,nilyr) - Tbot(i,j))

            ! Flux extra energy out of the ice
            fcondbot(m) = fcondbot(m) + einex(m)/dt    

            ferr(m) = abs( (enew(ij)-einit(m))/dt &
                    - (fcondtopn(i,j) - fcondbot(m) + fswint(i,j)) )

            ! factor of 0.9 allows for roundoff errors later
            if (ferr(m) > 0.9_dbl_kind*ferrmax) then         ! condition (5)

               converged(m) = .false.
               all_converged = .false.

                  ! reduce conductivity for next iteration
               do k = 1, nilyr
                  if (reduce_kh(m,k) .and. dqmat(m,k) > c0) then
                     frac = max(0.5*(c1-ferr(m)/abs(fcondtopn(i,j)-fcondbot(m))),p1)
!                     frac = p1
                     kh(m,k+nslyr+1) = kh(m,k+nslyr+1) * frac
                     kh(m,k+nslyr)   = kh(m,k+nslyr+1)
                  endif
               enddo

            endif               ! ferr 
         enddo                  ! ij

         deallocate(sbdiag)
         deallocate(diag)
         deallocate(spdiag)
         deallocate(rhs)
         deallocate(Tmat)
         deallocate(etai)
         deallocate(Tsf_start)
         deallocate(dTsf)
         deallocate(dfsurf_dT)
         deallocate(avg_Tsi)
         deallocate(enew)
         deallocate(dTi1)

      enddo                     ! temperature iteration niter

      if (.not.all_converged) then

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

      !-----------------------------------------------------------------
      ! Check for convergence failures.
      !-----------------------------------------------------------------
            if (.not.converged(ij)) then
               write(nu_diag,*) 'Thermo iteration does not converge,', &
                                'istep1, my_task, i, j:', &
                                 istep1, my_task, i, j
               write(nu_diag,*) 'Ice thickness:',  hilyr(ij)*nilyr
               write(nu_diag,*) 'Snow thickness:', hslyr(ij)*nslyr
               write(nu_diag,*) 'dTsf, Tsf_errmax:',dTsf_prev(ij), &
                                 Tsf_errmax
               write(nu_diag,*) 'Tsf:', Tsf(ij)
               write(nu_diag,*) 'fsurf:', fsurfn(i,j)
               write(nu_diag,*) 'fcondtop, fcondbot, fswint', &
                                 fcondtopn(i,j), fcondbot(ij), fswint(i,j)
               write(nu_diag,*) 'fswsfc, fswthrun', &
                                 fswsfc(i,j), fswthrun(i,j)
               write(nu_diag,*) 'Iswabs',(Iswabs(i,j,k),k=1,nilyr)
               write(nu_diag,*) 'Flux conservation error =', ferr(ij)
               write(nu_diag,*) 'Initial snow temperatures:'
               write(nu_diag,*) (Tsn_init(ij,k),k=1,nslyr)
               write(nu_diag,*) 'Initial ice temperatures:'
               write(nu_diag,*) (Tin_init(ij,k),k=1,nilyr)
               write(nu_diag,*) 'Matrix ice temperature diff:'
               write(nu_diag,*) (dTmat(ij,k),k=1,nilyr)
               write(nu_diag,*) 'dqmat*hilyr/dt:'
               write(nu_diag,*) (hilyr(ij)*dqmat(ij,k)/dt,k=1,nilyr)
               write(nu_diag,*) 'Final snow temperatures:'
               write(nu_diag,*) (Tsn(ij,k),k=1,nslyr)
               write(nu_diag,*) 'Matrix ice temperature diff:'
               write(nu_diag,*) (dTmat(ij,k),k=1,nilyr)
               write(nu_diag,*) 'dqmat*hilyr/dt:'
               write(nu_diag,*) (hilyr(ij)*dqmat(ij,k)/dt,k=1,nilyr)
               write(nu_diag,*) 'Final ice temperatures:'
               write(nu_diag,*) (Tin(ij,k),k=1,nilyr)
               write(nu_diag,*) 'Ice melting temperatures:'
               write(nu_diag,*) (Tmlts(ij,k),k=1,nilyr)
               write(nu_diag,*) 'Ice bottom temperature:', Tbot(i,j)
               write(nu_diag,*) 'dT initial:'
               write(nu_diag,*) (Tmlts(ij,k)-Tin_init(ij,k),k=1,nilyr)
               write(nu_diag,*) 'dT final:'
               write(nu_diag,*) (Tmlts(ij,k)-Tin(ij,k),k=1,nilyr)
               write(nu_diag,*) 'Sin'
               write(nu_diag,*) (Sin(ij,k),k=1,nilyr)
               l_stop = .true.
               istop = i
               jstop = j
               return
            endif
         enddo                  ! ij
      endif                     ! all_converged

      if (calc_Tsfc) then
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            ! update fluxes that depend on Tsf
            flwoutn(i,j) = flwoutn(i,j) + dTsf_prev(ij) * dflwout_dT(ij)
            fsensn(i,j)  = fsensn(i,j)  + dTsf_prev(ij) * dfsens_dT(ij)
            flatn(i,j)   = flatn(i,j)   + dTsf_prev(ij) * dflat_dT(ij)

         enddo                     ! ij
      endif                        ! calc_Tsfc

      end subroutine temperature_changes

!=======================================================================
!BOP
!
! !ROUTINE: conductivity - compute ice thermal conductivity
!
! !DESCRIPTION:
!
! Compute thermal conductivity at interfaces (held fixed during
!  the subsequent iteration).
!
! NOTE: Ice conductivity must be >= kimin
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine conductivity (nx_block, ny_block,         &
                               l_snow,   icells,           &
                               indxi,    indxj,    indxij, &
                               hilyr,    hslyr,            &
                               Tin,      kh,       Sin)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells          ! number of cells with aicen > puny

      logical (kind=log_kind), dimension(icells), &
         intent(in) :: &
         l_snow          ! true if snow temperatures are computed

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      integer (kind=int_kind), dimension (icells), &
         intent(in) :: &
         indxij          ! compressed 1D index for cells not converged

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         hilyr       , & ! ice layer thickness (same for all ice layers)
         hslyr           ! snow layer thickness (same for all snow layers)

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(in) :: &
         Tin         , & ! internal ice layer temperatures
         Sin             ! internal ice layer salinities

      real (kind=dbl_kind), dimension (icells,nilyr+nslyr+1), &
         intent(out) :: &
         kh              ! effective conductivity at interfaces (W m-2 deg-1)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         k               ! vertical index

      real (kind=dbl_kind), dimension (icells,nilyr) :: &
         kilyr           ! thermal cond at ice layer midpoints (W m-1 deg-1)

      real (kind=dbl_kind), dimension (icells,nslyr) :: &
         kslyr           ! thermal cond at snow layer midpoints (W m-1 deg-1)

      ! interior snow layers (simple for now, but may be fancier later)
      do k = 1, nslyr
         do ij = 1, icells
            kslyr(ij,k) = ksno
         enddo
      enddo                     ! nslyr

      ! interior ice layers
      if (conduct == 'MU71') then
         ! Maykut and Untersteiner 1971 form (with Wettlaufer 1991 constants)
         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               kilyr(ij,k) = kice + betak*Sin(ij,k)/min(-puny,Tin(ij,k))
               kilyr(ij,k) = max (kilyr(ij,k), kimin)
            enddo
         enddo                     ! nilyr
      else
         ! Pringle et al JGR 2007 'bubbly brine'
         do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               kilyr(ij,k) = (2.11_dbl_kind - 0.011_dbl_kind*Tin(ij,k) &
                            + 0.09_dbl_kind*Sin(ij,k)/min(-puny,Tin(ij,k))) &
                            * rhoi / 917._dbl_kind
               kilyr(ij,k) = max (kilyr(ij,k), kimin)
            enddo
         enddo                     ! nilyr
      endif ! conductivity

      ! top snow interface, top and bottom ice interfaces
      do ij = 1, icells
         ! top of snow layer; top surface of top ice layer
         if (l_snow(ij)) then
            kh(ij,1)       = c2 * kslyr(ij,1) / hslyr(ij)
            kh(ij,1+nslyr) = c2 * kslyr(ij,nslyr) * kilyr(ij,1) / &
                             ( kslyr(ij,nslyr)*hilyr(ij) +  &
                               kilyr(ij,1    )*hslyr(ij) )
         else
            kh(ij,1)       = c0
            kh(ij,1+nslyr) = c2 * kilyr(ij,1) / hilyr(ij)
         endif

         ! bottom surface of bottom ice layer
         kh(ij,1+nslyr+nilyr) = c2 * kilyr(ij,nilyr) / hilyr(ij)

      enddo                     ! ij

      ! interior snow interfaces

      if (nslyr > 1) then
         do k = 2, nslyr
            do ij = 1, icells
               if (l_snow(ij)) then
                  kh(ij,k) = c2 * kslyr(ij,k-1) * kslyr(ij,k) / &
                            ((kslyr(ij,k-1) + kslyr(ij,k))*hslyr(ij))
               else
                  kh(ij,k) = c0
               endif
            enddo                  ! ij
         enddo                     ! nilyr
      endif ! nslyr > 1

      ! interior ice interfaces
      do k = 2, nilyr
         do ij = 1, icells
            kh(ij,k+nslyr) = c2 * kilyr(ij,k-1) * kilyr(ij,k) / &
                            ((kilyr(ij,k-1) + kilyr(ij,k))*hilyr(ij))
         enddo                  ! ij
      enddo                     ! nilyr

      end subroutine conductivity

!=======================================================================
!BOP
!
! !ROUTINE: surface_fluxes - surface radiative and turbulent fluxes
!
! !DESCRIPTION:
!
! Compute radiative and turbulent fluxes and their derivatives
! with respect to Tsf.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine surface_fluxes (nx_block,   ny_block,          &
                                 isolve,     icells,            &
                                 indxii,     indxjj,    indxij, &
                                 Tsf,        fswsfc,            &
                                 rhoa,       flw,               &
                                 potT,       Qa,                &
                                 shcoef,     lhcoef,            &
                                 flwoutn,    fsensn,            &
                                 flatn,      fsurfn,            &
                                 dflwout_dT, dfsens_dT,         &
                                 dflat_dT,   dfsurf_dT)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         isolve            , & ! number of cells with temps not converged
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension(icells), &
         intent(in) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), dimension (icells) :: &
         indxij          ! compressed 1D index for cells not converged

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsurfn          ! net flux to top surface, excluding fcondtopn

      real (kind=dbl_kind), dimension (icells), &
         intent(inout) :: &
         dfsens_dT   , & ! deriv of fsens wrt Tsf (W m-2 deg-1)
         dflat_dT    , & ! deriv of flat wrt Tsf (W m-2 deg-1)
         dflwout_dT      ! deriv of flwout wrt Tsf (W m-2 deg-1)

      real (kind=dbl_kind), dimension (isolve), &
         intent(inout) :: &
         dfsurf_dT       ! derivative of fsurfn wrt Tsf
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         TsfK        , & ! ice/snow surface temperature (K)
         Qsfc        , & ! saturated surface specific humidity (kg/kg)
         dQsfcdT     , & ! derivative of Qsfc wrt surface temperature
         qsat        , & ! the saturation humidity of air (kg/m^3)
         flwdabs     , & ! downward longwave absorbed heat flx (W/m^2)
         tmpvar          ! 1/TsfK

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, isolve
         i = indxii(ij)         ! NOTE: not indxi and indxj
         j = indxjj(ij)
         m = indxij(ij)

         ! ice surface temperature in Kelvin
         TsfK = Tsf(m) + Tffresh
         tmpvar = c1/TsfK

         ! saturation humidity
         qsat    = qqqice * exp(-TTTice*tmpvar)
         Qsfc    = qsat / rhoa(i,j)
         dQsfcdT = TTTice * tmpvar*tmpvar * Qsfc

         ! longwave radiative flux
         flwdabs =  emissivity * flw(i,j)
         flwoutn(i,j) = -emissivity * stefan_boltzmann * TsfK**4

         ! downward latent and sensible heat fluxes
         fsensn(i,j) = shcoef(i,j) * (potT(i,j) - TsfK)
         flatn(i,j)  = lhcoef(i,j) * (Qa(i,j) - Qsfc)

         ! derivatives wrt surface temp
         dflwout_dT(m) = - emissivity*stefan_boltzmann * c4*TsfK**3
         dfsens_dT(m)  = - shcoef(i,j)
         dflat_dT(m)   = - lhcoef(i,j) * dQsfcdT

         fsurfn(i,j) = fswsfc(i,j) + flwdabs + flwoutn(i,j) &
                     + fsensn(i,j) + flatn(i,j)
         dfsurf_dT(ij) = dflwout_dT(m) &
                         + dfsens_dT(m) + dflat_dT(m)

      enddo                     ! ij

      end subroutine surface_fluxes

!=======================================================================
!BOP
!
! !ROUTINE: get_matrix_elements - compute tridiagonal matrix elements
!
! !DESCRIPTION:
!
! Compute terms in tridiagonal matrix that will be solved to find
!  the new vertical temperature profile
! This routine is for the case in which Tsfc is being computed.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
!
! March 2004 by William H. Lipscomb for multiple snow layers
! April 2008 by E. C. Hunke, divided into two routines based on calc_Tsfc 
!
! !INTERFACE:
!
      subroutine get_matrix_elements_calc_Tsfc &
                                     (nx_block, ny_block,         &
                                      isolve,   icells,           &
                                      indxii,   indxjj,   indxij, &
                                      l_snow,   l_cold,           &
                                      Tsf,      Tbot,             &
                                      fsurfn,   dfsurf_dT,        &
                                      Tin_init, Tsn_init,         &
                                      kh,       Sswabs,           &
                                      Iswabs,                     &
                                      etai,     etas,             &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         isolve            , & ! number of cells with temps not converged
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(icells), &
         intent(in) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), dimension (icells), &
         intent(in) :: &
         indxij          ! compressed 1D index for cells not converged

      logical (kind=log_kind), dimension (icells), &
         intent(in) :: &
         l_snow      , & ! true if snow temperatures are computed
         l_cold          ! true if surface temperature is computed

      real (kind=dbl_kind), dimension (icells), intent(in) :: &
         Tsf             ! ice/snow top surface temp (deg C)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn (W/m^2)
         Tbot            ! ice bottom surface temperature (deg C)

      real (kind=dbl_kind), dimension (isolve), intent(in) :: &
         dfsurf_dT       ! derivative of fsurf wrt Tsf

      real (kind=dbl_kind), dimension (isolve,nilyr), &
         intent(in) :: &
         etai            ! dt / (rho*cp*h) for ice layers

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(in) :: &
         Tin_init        ! ice temp at beginning of time step

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(in) :: &
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(in) :: &
         Iswabs          ! absorbed SW flux in ice layers

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(in) :: &
         etas        , & ! dt / (rho*cp*h) for snow layers
         Tsn_init        ! snow temp at beginning of time step
                         ! Note: no absorbed SW in snow layers

      real (kind=dbl_kind), dimension (icells,nslyr+nilyr+1), &
         intent(in) :: &
         kh              ! effective conductivity at layer interfaces

      real (kind=dbl_kind), dimension (isolve,nslyr+nilyr+1), &
         intent(inout) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         k, ks, ki, kr   ! vertical indices and row counters

      !-----------------------------------------------------------------
      ! Initialize matrix elements.
      ! Note: When we do not need to solve for the surface or snow
      !       temperature, we solve dummy equations with solution T = 0.
      !       Ice layers are fully initialized below.
      !-----------------------------------------------------------------

      do k = 1, nslyr+1
         do ij = 1, isolve
            sbdiag(ij,k) = c0
            diag  (ij,k) = c1
            spdiag(ij,k) = c0
            rhs   (ij,k) = c0
         enddo
      enddo
            
      !-----------------------------------------------------------------
      ! Compute matrix elements
      !
      ! Four possible cases to solve:
      !   (1) Cold surface (Tsf < 0), snow present
      !   (2) Melting surface (Tsf = 0), snow present
      !   (3) Cold surface (Tsf < 0), no snow
      !   (4) Melting surface (Tsf = 0), no snow
      !-----------------------------------------------------------------

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

      !-----------------------------------------------------------------
      ! Tsf equation for case of cold surface (with or without snow)
      !-----------------------------------------------------------------
            if (l_cold(m)) then
               if (l_snow(m)) then
                  k = 1
               else                ! no snow
                  k = 1 + nslyr
               endif
               kr = k

               sbdiag(ij,kr) = c0
               diag  (ij,kr) = dfsurf_dT(ij) - kh(m,k)
               spdiag(ij,kr) = kh(m,k)
               rhs   (ij,kr) = dfsurf_dT(ij)*Tsf(m) - fsurfn(i,j)
            endif                  ! l_cold

      !-----------------------------------------------------------------
      ! top snow layer
      !-----------------------------------------------------------------
!           k = 1
!           kr = 2

            if (l_snow(m)) then
               if (l_cold(m)) then
                  sbdiag(ij,2) = -etas(m,1) * kh(m,1)
                  spdiag(ij,2) = -etas(m,1) * kh(m,2)
                  diag  (ij,2) = c1 &
                                + etas(m,1) * (kh(m,1) + kh(m,2))
                  rhs   (ij,2) = Tsn_init(m,1) &
                                + etas(m,1) * Sswabs(i,j,1)
               else                ! melting surface
                  sbdiag(ij,2) = c0
                  spdiag(ij,2) = -etas(m,1) * kh(m,2)
                  diag  (ij,2) = c1 &
                                + etas(m,1) * (kh(m,1) + kh(m,2))
                  rhs   (ij,2) = Tsn_init(m,1) &
                                + etas(m,1)*kh(m,1)*Tsf(m) &
                                + etas(m,1) * Sswabs(i,j,1)
               endif               ! l_cold
            endif                  ! l_snow

         enddo                    ! ij

      !-----------------------------------------------------------------
      ! remaining snow layers
      !-----------------------------------------------------------------

      if (nslyr > 1) then

         do k = 2, nslyr
            kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m)) then
                  sbdiag(ij,kr) = -etas(m,k) * kh(m,k)
                  spdiag(ij,kr) = -etas(m,k) * kh(m,k+1)
                  diag  (ij,kr) = c1 &
                               + etas(m,k) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tsn_init(m,k) &
                               + etas(m,k) * Sswabs(i,j,k)
               endif
            enddo               ! ij
         enddo                  ! nslyr

      endif                     ! nslyr > 1


      if (nilyr > 1) then

      !-----------------------------------------------------------------
      ! top ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m) .or. l_cold(m)) then
                  sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
                  spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
                  diag  (ij,kr) = c1 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki) &
                                 + etai(ij,ki)*Iswabs(i,j,ki)
               else    ! no snow, warm surface
                  sbdiag(ij,kr) = c0
                  spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
                  diag  (ij,kr) = c1 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki) &
                                 + etai(ij,ki)*Iswabs(i,j,ki) &
                                 + etai(ij,ki)*kh(m,k)*Tsf(m)
               endif
            enddo    ! ij

      !-----------------------------------------------------------------
      ! bottom ice layer
      !-----------------------------------------------------------------

         ki = nilyr
         k  = ki + nslyr
         kr = k + 1
      
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)
            sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
            spdiag(ij,kr) = c0
            diag  (ij,kr) = c1  &
                           + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
            rhs   (ij,kr) = Tin_init(m,ki) &
                           + etai(ij,ki)*Iswabs(i,j,ki) &
                           + etai(ij,ki)*kh(m,k+1)*Tbot(i,j)
         enddo                   ! ij
      
      else         ! nilyr = 1

      !-----------------------------------------------------------------
      ! single ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m) .or. l_cold(m)) then
                  sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
                  spdiag(ij,kr) = c0
                  diag  (ij,kr) = c1                                 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki)                     &
                                 + etai(ij,ki) * Iswabs(i,j,ki)      &
                                 + etai(ij,ki) * kh(m,k+1)*Tbot(i,j)
               else   ! no snow, warm surface
                  sbdiag(ij,kr) = c0
                  spdiag(ij,kr) = c0
                  diag  (ij,kr) = c1                                 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki)                     &
                                 + etai(ij,ki) * Iswabs(i,j,ki)      &
                                 + etai(ij,ki) * kh(m,k)*Tsf(m)      &
                                 + etai(ij,ki) * kh(m,k+1)*Tbot(i,j)
               endif
            enddo                  ! ij
 
      endif        ! nilyr > 1

      !-----------------------------------------------------------------
      ! interior ice layers
      !-----------------------------------------------------------------

      do ki = 2, nilyr-1
           
         k  = ki + nslyr
         kr = k + 1
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

            sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
            spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
            diag  (ij,kr) = c1 &
                           + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
            rhs   (ij,kr) = Tin_init(m,ki) &
                           + etai(ij,ki)*Iswabs(i,j,ki)
         enddo                  ! ij
      enddo                     ! nilyr

      end subroutine get_matrix_elements_calc_Tsfc

!=======================================================================
!BOP
!
! !ROUTINE: get_matrix_elements - compute tridiagonal matrix elements
!
! !DESCRIPTION:
!
! Compute terms in tridiagonal matrix that will be solved to find
!  the new vertical temperature profile
! This routine is for the case in which Tsfc is already known.
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
!
! March 2004 by William H. Lipscomb for multiple snow layers
! April 2008 by E. C. Hunke, divided into two routines based on calc_Tsfc 
!
! !INTERFACE:
!
      subroutine get_matrix_elements_know_Tsfc &
                                     (nx_block, ny_block,         &
                                      isolve,   icells,           &
                                      indxii,   indxjj,   indxij, &
                                      l_snow,   Tbot,             &
                                      Tin_init, Tsn_init,         &
                                      kh,       Sswabs,           &
                                      Iswabs,                     &
                                      etai,     etas,             &
                                      sbdiag,   diag,             &
                                      spdiag,   rhs,              &
                                      fcondtopn)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         isolve            , & ! number of cells with temps not converged
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(icells), &
         intent(in) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), dimension (icells), &
         intent(in) :: &
         indxij          ! compressed 1D index for cells not converged

      logical (kind=log_kind), dimension (icells), &
         intent(in) :: &
         l_snow          ! true if snow temperatures are computed

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tbot            ! ice bottom surface temperature (deg C)

      real (kind=dbl_kind), dimension (isolve,nilyr), &
         intent(in) :: &
         etai            ! dt / (rho*cp*h) for ice layers

      real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(in) :: &
         Tin_init        ! ice temp at beginning of time step

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(in) :: &
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(in) :: &
         Iswabs          ! absorbed SW flux in ice layers

      real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(in) :: &
         etas        , & ! dt / (rho*cp*h) for snow layers
         Tsn_init        ! snow temp at beginning of time step
                         ! Note: no absorbed SW in snow layers

      real (kind=dbl_kind), dimension (icells,nslyr+nilyr+1), &
         intent(in) :: &
         kh              ! effective conductivity at layer interfaces

      real (kind=dbl_kind), dimension (isolve,nslyr+nilyr+1), &
         intent(inout) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in),  &
         optional :: &
         fcondtopn       ! conductive flux at top sfc, positive down (W/m^2)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij, m       , & ! horizontal indices, combine i and j loops
         k, ks, ki, kr   ! vertical indices and row counters

      !-----------------------------------------------------------------
      ! Initialize matrix elements.
      ! Note: When we do not need to solve for the surface or snow
      !       temperature, we solve dummy equations with solution T = 0.
      !       Ice layers are fully initialized below.
      !-----------------------------------------------------------------

      do k = 1, nslyr+1
         do ij = 1, isolve
            sbdiag(ij,k) = c0
            diag  (ij,k) = c1
            spdiag(ij,k) = c0
            rhs   (ij,k) = c0
         enddo
      enddo
            
      !-----------------------------------------------------------------
      ! Compute matrix elements
      !
      ! Four possible cases to solve:
      !   (1) Cold surface (Tsf < 0), snow present
      !   (2) Melting surface (Tsf = 0), snow present
      !   (3) Cold surface (Tsf < 0), no snow
      !   (4) Melting surface (Tsf = 0), no snow
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! top snow layer
      !-----------------------------------------------------------------
!        k = 1
!        kr = 2

         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

            if (l_snow(m)) then
               sbdiag(ij,2) = c0
               spdiag(ij,2) = -etas(m,1) * kh(m,2)
               diag  (ij,2) = c1                                 &
                             + etas(m,1) * kh(m,2)
               rhs   (ij,2) = Tsn_init(m,1)                      &
                             + etas(m,1) * Sswabs(i,j,1)         &
                             + etas(m,1) * fcondtopn(i,j)
            endif   ! l_snow
         enddo   ! ij

      !-----------------------------------------------------------------
      ! remaining snow layers
      !-----------------------------------------------------------------

      if (nslyr > 1) then

         do k = 2, nslyr
            kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m)) then
                  sbdiag(ij,kr) = -etas(m,k) * kh(m,k)
                  spdiag(ij,kr) = -etas(m,k) * kh(m,k+1)
                  diag  (ij,kr) = c1 &
                               + etas(m,k) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tsn_init(m,k) &
                               + etas(m,k) * Sswabs(i,j,k)
               endif
            enddo               ! ij
         enddo                  ! nslyr

      endif                     ! nslyr > 1


      if (nilyr > 1) then

      !-----------------------------------------------------------------
      ! top ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m)) then

                  sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
                  spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
                  diag  (ij,kr) = c1                                &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki)                    &
                                 + etai(ij,ki) * Iswabs(i,j,ki)
               else                  
                  sbdiag(ij,kr) = c0
                  spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
                  diag  (ij,kr) = c1                                &
                                 + etai(ij,ki) * kh(m,k+1)
                  rhs   (ij,kr) = Tin_init(m,ki)                    &
                                 + etai(ij,ki) * Iswabs(i,j,ki)       &
                                 + etai(ij,ki) * fcondtopn(i,j)
               endif  ! l_snow
            enddo   ! ij

      !-----------------------------------------------------------------
      ! bottom ice layer
      !-----------------------------------------------------------------

         ki = nilyr
         k  = ki + nslyr
         kr = k + 1
      
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)
            sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
            spdiag(ij,kr) = c0
            diag  (ij,kr) = c1  &
                           + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
            rhs   (ij,kr) = Tin_init(m,ki) &
                           + etai(ij,ki)*Iswabs(i,j,ki) &
                           + etai(ij,ki)*kh(m,k+1)*Tbot(i,j)

         enddo                   ! ij
      
      else         ! nilyr = 1

      !-----------------------------------------------------------------
      ! single ice layer
      !-----------------------------------------------------------------

         ki = 1
         k  = ki + nslyr
         kr = k + 1

            do ij = 1, isolve
               i = indxii(ij)
               j = indxjj(ij)
               m = indxij(ij)

               if (l_snow(m)) then
                  sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
                  spdiag(ij,kr) = c0
                  diag  (ij,kr) = c1                                 &
                                 + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
                  rhs   (ij,kr) = Tin_init(m,ki)                     &
                                 + etai(ij,ki) * Iswabs(i,j,ki)      &
                                 + etai(ij,ki) * kh(m,k+1)*Tbot(i,j)
               else
                  sbdiag(ij,kr) = c0
                  spdiag(ij,kr) = c0
                  diag  (ij,kr) = c1                                 &
                                 + etai(ij,ki) * kh(m,k+1)
                  rhs   (ij,kr) = Tin_init(m,ki)                     &
                                 + etai(ij,ki) * Iswabs(i,j,ki)      &
                                 + etai(ij,ki) * fcondtopn(i,j)      &
                                 + etai(ij,ki) * kh(m,k+1)*Tbot(i,j)
               endif
            enddo                     ! ij

      endif        ! nilyr > 1

      !-----------------------------------------------------------------
      ! interior ice layers
      !-----------------------------------------------------------------

      do ki = 2, nilyr-1
           
         k  = ki + nslyr
         kr = k + 1
         do ij = 1, isolve
            i = indxii(ij)
            j = indxjj(ij)
            m = indxij(ij)

            sbdiag(ij,kr) = -etai(ij,ki) * kh(m,k)
            spdiag(ij,kr) = -etai(ij,ki) * kh(m,k+1)
            diag  (ij,kr) = c1 &
                           + etai(ij,ki) * (kh(m,k) + kh(m,k+1))
            rhs   (ij,kr) = Tin_init(m,ki) &
                           + etai(ij,ki)*Iswabs(i,j,ki)

         enddo                  ! ij
      enddo                     ! nilyr

      end subroutine get_matrix_elements_know_Tsfc

!=======================================================================
!BOP
!
! !ROUTINE: tridiag_solver - tridiagonal matrix solver
!
! !DESCRIPTION:
!
! Tridiagonal matrix solver--used to solve the implicit vertical heat
! equation in ice and snow
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine tridiag_solver (nx_block, ny_block, &
                                 isolve,   icells,   &
                                 indxii,   indxjj,   &
                                 nmat,     sbdiag,   &
                                 diag,     spdiag,   &
                                 rhs,      xout)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         isolve            , & ! number of cells with temps not converged
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(icells), &
         intent(in) :: &
         indxii, indxjj  ! compressed indices for cells not converged

      integer (kind=int_kind), intent(in) :: &
         nmat            ! matrix dimension

      real (kind=dbl_kind), dimension (isolve,nmat), &
           intent(in) :: &
         sbdiag      , & ! sub-diagonal matrix elements
         diag        , & ! diagonal matrix elements
         spdiag      , & ! super-diagonal matrix elements
         rhs             ! rhs of tri-diagonal matrix eqn.

      real (kind=dbl_kind), dimension (isolve,nmat), &
           intent(inout) :: &
         xout            ! solution vector
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k               ! row counter

      real (kind=dbl_kind), dimension (isolve) :: &
         wbeta           ! temporary matrix variable

      real (kind=dbl_kind), dimension(isolve,nilyr+nslyr+1):: &
         wgamma          ! temporary matrix variable

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, isolve
         wbeta(ij) = diag(ij,1)
         xout(ij,1) = rhs(ij,1) / wbeta(ij)
      enddo                     ! ij

      do k = 2, nmat
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            wgamma(ij,k) = spdiag(ij,k-1) / wbeta(ij)
            wbeta(ij) = diag(ij,k) - sbdiag(ij,k)*wgamma(ij,k)
            xout(ij,k) = (rhs(ij,k) - sbdiag(ij,k)*xout(ij,k-1)) &
                         / wbeta(ij)
         enddo                  ! ij
      enddo                     ! k

      do k = nmat-1, 1, -1
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, isolve
            xout(ij,k) = xout(ij,k) - wgamma(ij,k+1)*xout(ij,k+1)
         enddo                  ! ij
      enddo                     ! k

      end subroutine tridiag_solver

!=======================================================================

      end module ice_therm_bl99

!=======================================================================
