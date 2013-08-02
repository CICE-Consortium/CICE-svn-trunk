#define debug 0
#define diag 0

#define flushing 1
#define flushing_down 1
#define flushing_up 0
 
#define load_save_state 0

module ice_therm_mushy

  use ice_kinds_mod
  use ice_constants
  use ice_domain_size, only: ncat, nilyr, nslyr, max_ntrcr, n_aero
  use ice_calendar, only: istep, istep1, dt
  use ice_state, only: nt_apnd, nt_hpnd, tr_pond
  use ice_exit, only: abort_ice
  use ice_communicate, only: my_task, master_task
  use ice_therm_shared, only: ferrmax, hfrazilmin
  use ice_fileunits, only: nu_diag

  implicit none

  private
  public :: temperature_mush, &
            temperature_snow, &
            liquid_fraction, &
            temperature_mush_liquid_fraction, &
            liquidus_brine_salinity_mush, &
            liquidus_temperature_mush, &
            temperature_changes_salinity, &
            enthalpy_mush, &
            enthalpy_of_melting, &
            add_new_ice_mushy, &
            permeability
  
  !-----------------------------------------------------------------
  ! namelist parameters
  !-----------------------------------------------------------------
  real(kind=dbl_kind), public :: &
       a_rapid_mode      , & ! channel radius for rapid drainage mode (m)
       Rac_rapid_mode    , & ! critical Rayleigh number for rapid drainage mode
       aspect_rapid_mode , & ! aspect ratio for rapid drainage mode (larger is wider)
       dSdt_slow_mode    , & ! slow mode drainage strength (m s-1 K-1)
       phi_c_slow_mode   , & ! critical liquid fraction porosity cutoff for slow mode
       phi_i_mushy           ! liquid fraction of congelation ice

  !-----------------------------------------------------------------
  ! physical parameters for passing to residual function through solver
  !-----------------------------------------------------------------

  integer(kind=int_kind) :: &
       g_i , & ! current grid i-index (for diagnostics)
       g_j , & ! current grid j-index  (for diagnostics)
       g_n     ! current category (for diagnostics)

  !-----------------------------------------------------------------
  ! Constants for Liquidus relation from Assur (1958)
  !-----------------------------------------------------------------

  ! liquidus relation - higher temperature region
  real(kind=dbl_kind), parameter :: &
       az1_liq = -18.48_dbl_kind, &
       bz1_liq =    0.0_dbl_kind

  ! liquidus relation - lower temperature region
  real(kind=dbl_kind), parameter :: &
       az2_liq = -10.3085_dbl_kind, &
       bz2_liq =     62.4_dbl_kind

  ! liquidus break
  real(kind=dbl_kind), parameter :: &
       Tb_liq = -7.6362968855167352_dbl_kind, & ! temperature of liquidus break
       Sb_liq =  123.66702800276086_dbl_kind    ! salinity of liquidus break

  ! basic liquidus relation constants
  real(kind=dbl_kind), parameter :: &
       az1p_liq = az1_liq / 1000.0_dbl_kind, &
       bz1p_liq = bz1_liq / 1000.0_dbl_kind, &
       az2p_liq = az2_liq / 1000.0_dbl_kind, &
       bz2p_liq = bz2_liq / 1000.0_dbl_kind
  
  ! quadratic constants - higher temperature region
  real(kind=dbl_kind), parameter :: &
       AS1_liq = az1p_liq * (rhow * cp_ocn - rhoi * cp_ice)                                   , &
       AC1_liq = rhoi * cp_ice * az1_liq                                                      , &
       BS1_liq = (c1 + bz1p_liq) * (rhow * cp_ocn - rhoi * cp_ice) + rhoi * Lfresh * az1p_liq , &
       BQ1_liq = -az1_liq                                                                     , &
       BC1_liq = rhoi * cp_ice * bz1_liq - rhoi * Lfresh * az1_liq                            , &
       CS1_liq = rhoi * Lfresh * (c1 + bz1p_liq)                                              , &
       CQ1_liq = -bz1_liq                                                                     , &
       CC1_liq = -rhoi * Lfresh * bz1_liq
  
  ! quadratic constants - lower temperature region
  real(kind=dbl_kind), parameter :: &
       AS2_liq = az2p_liq * (rhow * cp_ocn - rhoi * cp_ice)                                   , &
       AC2_liq = rhoi * cp_ice * az2_liq                                                      , &
       BS2_liq = (c1 + bz2p_liq) * (rhow * cp_ocn - rhoi * cp_ice) + rhoi * Lfresh * az2p_liq , &
       BQ2_liq = -az2_liq                                                                     , &
       BC2_liq = rhoi * cp_ice * bz2_liq - rhoi * Lfresh * az2_liq                            , &
       CS2_liq = rhoi * Lfresh * (c1 + bz2p_liq)                                              , &
       CQ2_liq = -bz2_liq                                                                     , &
       CC2_liq = -rhoi * Lfresh * bz2_liq
  
  ! break enthalpy constants
  real(kind=dbl_kind), parameter :: &
       D_liq = ((c1 + az1p_liq * Tb_liq + bz1p_liq) / (az1_liq * Tb_liq + bz1_liq)) * &
               ((cp_ocn * rhow - cp_ice * rhoi) * Tb_liq + Lfresh * rhoi), &
       E_liq = cp_ice * rhoi * Tb_liq - Lfresh * rhoi
  
  ! just fully melted enthapy constants
  real(kind=dbl_kind), parameter :: &
       F1_liq = (-1000.0_dbl_kind * cp_ocn * rhow) / az1_liq , &
       G1_liq = -1000.0_dbl_kind                             , &
       H1_liq = (-bz1_liq * cp_ocn * rhow) / az1_liq         , &
       F2_liq = (-1000.0_dbl_kind * cp_ocn * rhow) / az2_liq , &
       G2_liq = -1000.0_dbl_kind                             , &
       H2_liq = (-bz2_liq * cp_ocn * rhow) / az2_liq
  
  ! warmer than fully melted constants
  real(kind=dbl_kind), parameter :: &
       I_liq = c1 / (cp_ocn * rhow)

  ! temperature to brine salinity
  real(kind=dbl_kind), parameter :: &
       J1_liq = bz1_liq / az1_liq         , &
       K1_liq = c1 / 1000.0_dbl_kind      , &
       L1_liq = (c1 + bz1p_liq) / az1_liq , &
       J2_liq = bz2_liq / az2_liq         , &
       K2_liq = c1 / 1000.0_dbl_kind      , &
       L2_liq = (c1 + bz2p_liq) / az2_liq

  ! brine salinity to temperature
  real(kind=dbl_kind), parameter :: &
       M1_liq = az1_liq            , &
       N1_liq = -az1p_liq          , &
       O1_liq = -bz1_liq / az1_liq , &
       M2_liq = az2_liq            , &
       N2_liq = -az2p_liq          , &
       O2_liq = -bz2_liq / az2_liq

 !-----------------------------------------------------------------
 ! Other parameters
 !-----------------------------------------------------------------

  real(kind=dbl_kind), parameter :: &
       serrmax = 1e-11_dbl_kind  ! max allowed salt flux error (ppt m s-1)

  real(kind=dbl_kind), parameter :: &
       dTsf_errmax = 5.e-4_dbl_kind  ! max allowed error in dTsf

  real(kind=dbl_kind), parameter :: &
       dTsf_errcon = dTsf_errmax  ! max allowed error in dTsf for consistency

  real(kind=dbl_kind), parameter :: &
       ferrcon = ferrmax ! max allowed surface flux error for consistency

  real(kind=dbl_kind), parameter :: &
       hs_min = 1.e-4_dbl_kind  ! min snow thickness for computing Tsno (m)

!=======================================================================

contains

!=======================================================================

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
      subroutine add_new_ice_mushy (nx_block,  ny_block,   &
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
                              phi_init,              &
                              dSin0_frazil,          &
                              l_stop,                &
                              istop,     jstop)
!
! !USES:
!
      use ice_itd, only: hin_max, column_sum, &
                         column_conservation_check 
      use ice_therm_itd, only: update_vertical_tracers
      use ice_state, only: nt_Tsfc, nt_iage, nt_FY, nt_alvl, nt_vlvl, nt_aero, &
                           nt_sice, nt_qice, &
                           nt_apnd, tr_pond_cesm, tr_pond_lvl, tr_pond_topo, &
                           tr_iage, tr_FY, tr_lvl, tr_aero

! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         icells                ! number of ice/ocean grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi,  indxj         ! compressed i/j indices

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

      real (kind=dbl_kind), intent(in) :: &
         phi_init     , & ! initial frazil liquid fraction
         dSin0_frazil     ! initial frazil bulk salinity reduction from sss

      logical (kind=log_kind), intent(in) :: &
         update_ocn_f ! if true, update fresh water and salt fluxes

      logical (kind=log_kind), intent(out) :: &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(out) :: &
         istop, jstop    ! indices of grid cell where model aborts
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
         hsurp        , & ! thickness of new ice added to each cat
         vlyr             ! ice layer volume

      real (kind=dbl_kind), dimension (icells) :: &
         vice1        , & ! starting volume of existing ice
         vice_init, vice_final  ! ice volume summed over categories

      real (kind=dbl_kind) :: &
         fnew         , & ! heat flx to open water for new ice (W/m^2)
         hi0new       , & ! thickness of new ice
         hi0max       , & ! max allowed thickness of new ice
         qi0(nilyr)   , & ! frazil ice enthalpy
         qi0av        , & ! mean value of qi0 for new ice (J kg-1)
         Si0(nilyr)   , & ! frazil ice salinity
         vsurp        , & ! volume of new ice added to each cat
         vtmp         , & ! total volume of new and old ice
         area1        , & ! starting fractional area of existing ice
         alvl         , & ! starting level ice area
         rnilyr       , & ! real(nilyr)
         dfresh       , & ! change in fresh
         dfsalt       , & ! change in fsalt
         Ti               ! frazil temperature
      
      real (kind=dbl_kind), dimension (icells) :: &
         qi0new       , & ! frazil ice enthalpy
         Si0new           ! frazil ice bulk salinity

      integer (kind=int_kind) :: &
         jcells, kcells     , & ! grid cell counters
         ij, m                  ! combined i/j horizontal indices

      integer (kind=int_kind), dimension (icells) :: &
         indxij2,  indxij3  , & ! compressed i/j indices
         indxi2, indxj2     , &
         indxi3, indxj3

      character (len=char_len) :: &
         fieldid           ! field identifier

      l_stop = .false.
      istop = 0
      jstop = 0

      jcells = 0
      kcells = 0

      rnilyr = real(nilyr,kind=dbl_kind)


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

      !-----------------------------------------------------------------
      ! Compute average enthalpy of new ice.
      !
      ! POP assumes new ice is fresh.  Otherwise, it would be better
      ! to do something like this:
      !  qi0(i,j,k) = -rhoi * (cp_ice*(Tmlt(k)-Tf(i,j))
      !             + Lfresh*(1.-Tmlt(k)/Tf(i,j)) - cp_ocn*Tmlt(k))

!echmod: This is inconsistent.  For now, mimic constant salinity profile
!echmod: in tracer nt_sice.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            
            Si0new(ij) = sss(i,j) - dSin0_frazil
            
            Ti = liquidus_temperature_mush(Si0new(ij) / phi_init)
            
            qi0new(ij) = enthalpy_mush(Ti, Si0new(ij))
            
         enddo ! ij

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
         vi0new(ij) = -fnew*dt / qi0new(ij) ! note sign convention, qi < 0

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
                  vsurp       = vi0new(ij) - ai0new(ij)*hi0new
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
            vlyr(m) = vsurp/rnilyr

         enddo                  ! ij

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, kcells
               i = indxi3(ij)
               j = indxj3(ij)
               m = indxij3(ij)
               
               vsurp = hsurp(m) * aicen(i,j,n)  ! note - save this above?
               vtmp = vicen(i,j,n) - vsurp      ! vicen is the new volume
               if (vicen(i,j,n) > c0) then
                  
                  call update_vertical_tracers(trcrn(i,j,nt_qice:nt_qice+nilyr-1,n), vtmp, vicen(i,j,n), qi0new(m))
                  call update_vertical_tracers(trcrn(i,j,nt_sice:nt_sice+nilyr-1,n), vtmp, vicen(i,j,n), Si0new(m))
                  
               endif

            enddo               ! ij
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
                alvl = trcrn(i,j,nt_alvl,1)
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

         vlyr(m)    = vi0new(m) / rnilyr
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
                                 + qi0new(m)*vi0new(m))/vicen(i,j,1)
                  ! salinity
                  trcrn(i,j,nt_sice+k-1,1) = &
                 (trcrn(i,j,nt_sice+k-1,1)*vice1(ij) &
                                 + Si0new(m)*vi0new(m))/vicen(i,j,1)
               endif
            enddo
         enddo

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

      end subroutine add_new_ice_mushy

!=======================================================================

  subroutine temperature_changes_salinity(nx_block, ny_block, &
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
                                          trcrn,              &
                                          Tsf,      Tbot,     &
                                          sss,                &
                                          fsensn,   flatn,    &
                                          fswabsn,  flwoutn,  &
                                          fsurfn,             &
                                          fcondtopn,fcondbot, &
                                          fadvocn,  snoice,   &
                                          einit,    l_stop,   &
                                          istop,    jstop)

    ! solve the changes in enthalpy and bulk salinity for the mushy vertical thermodynamics

    integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         my_task     , & ! task number (diagnostic only)
         istep1      , & ! time step index (diagnostic only)
         icells          ! number of cells with aicen > puny
    
    real (kind=dbl_kind), intent(in) :: &
         dt              ! time step (s)
    
    integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny
    
    real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         Tbot        , & ! ice bottom surface temperature (deg C)
         sss             ! sea surface salinity (PSU)

    real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         fswint      , & ! SW absorbed in ice interior below surface (W m-2)
         fswthrun        ! SW through ice to ocean         (W m-2)
    
    real (kind=dbl_kind), dimension (icells), intent(inout) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr           ! snow layer thickness (m)

    real (kind=dbl_kind), dimension (icells), intent(in) :: &
         einit           ! initial energy of melting (J m-2)

    real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(inout) :: &
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)
    
    real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(inout) :: &
         Iswabs          ! SW radiation absorbed in ice layers (W m-2)
    
    real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fswabsn     , & ! shortwave absorbed by ice (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fadvocn         ! advection heat flux to ocean
    
    real (kind=dbl_kind), dimension (icells), intent(out):: &
         fcondbot        ! downward cond flux at bottom surface (W m-2)
    
    real (kind=dbl_kind), dimension (icells), &
         intent(inout):: &
         Tsf             ! ice/snow surface temperature (C)
    
    real (kind=dbl_kind), dimension (icells,nilyr), &
         intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3)
         Tin         , & ! internal ice layer temperatures
         Sin             ! internal ice layer salinities

    real (kind=dbl_kind), dimension (icells,nslyr), &
         intent(inout) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         Tsn             ! internal snow layer temperatures

    real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr), &
         intent(inout) :: &
         trcrn
    
    real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         snoice          ! snow-ice formation       (m/step-->cm/day) 

    logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model
    
    integer (kind=int_kind), intent(inout) :: &
         istop, jstop    ! i and j indices of cell where model fails
 
    real(kind=dbl_kind) :: &
         trc_hpnd, &     ! pond depth tracer
         trc_apnd        ! pond area tracer

    integer(kind=int_kind) :: &
         ij, &           ! icells index
         i,  &           ! i index
         j               ! j index
    
    ! loop over cells
    do ij = 1, icells
       i = indxi(ij)
       j = indxj(ij)
       
       g_i = i
       g_j = j
       g_n = 1
       
       if (tr_pond) then
          trc_hpnd = trcrn(i,j,nt_hpnd)
          trc_apnd = trcrn(i,j,nt_apnd)
       else
          trc_hpnd = c0
          trc_apnd = c0
       endif
          
       !call temperature_changes_column(istep1,        dt,            & 
       !                                rhoa(i,j),     flw(i,j),      &
       !                                potT(i,j),     Qa(i,j),       &
       !                                shcoef(i,j),   lhcoef(i,j),   &
       !                                fswsfc(i,j),   fswint(i,j),   &
       !                                fswthrun(i,j),                &
       !                                Sswabs(i,j,:), Iswabs(i,j,:), &
       !                                hilyr(ij),     hslyr(ij),     &
       !                                qin(ij,:),     Tin(ij,:),     &
       !                                qsn(ij,:),     Tsn(ij,:),     &
       !                                Sin(ij,:),                    &
       !                                trc_hpnd,      trc_apnd,      &
       !                                Tsf(ij),       Tbot(i,j),     &
       !                                sss(i,j),                     &
       !                                fsensn(i,j),   flatn(i,j),    &
       !                                fswabsn(i,j),  flwoutn(i,j),  &
       !                                fsurfn(i,j),                  &
       !                                fcondtopn(i,j),fcondbot(ij),  &
       !                                fadvocn(i,j),  snoice(i,j),   &
       !                                einit(ij),     l_stop)


       ! comment out above call and comment in this call for debugging
       call temperature_changes_column_debug(istep1,        dt,            & 
                                             rhoa(i,j),     flw(i,j),      &
                                             potT(i,j),     Qa(i,j),       &
                                             shcoef(i,j),   lhcoef(i,j),   &
                                             fswsfc(i,j),   fswint(i,j),   &
                                             fswthrun(i,j),                &
                                             Sswabs(i,j,:), Iswabs(i,j,:), &
                                             hilyr(ij),     hslyr(ij),     &
                                             qin(ij,:),     Tin(ij,:),     &
                                             qsn(ij,:),     Tsn(ij,:),     &
                                             Sin(ij,:),                    &
                                             trc_hpnd,      trc_apnd,      &
                                             Tsf(ij),       Tbot(i,j),     &
                                             sss(i,j),                     &
                                             fsensn(i,j),   flatn(i,j),    &
                                             fswabsn(i,j),  flwoutn(i,j),  &
                                             fsurfn(i,j),                  &
                                             fcondtopn(i,j),fcondbot(ij),  &
                                             fadvocn(i,j),  snoice(i,j),   &
                                             einit(ij),     l_stop)

       if (tr_pond) then
          trcrn(i,j,nt_hpnd) = trc_hpnd
          trcrn(i,j,nt_apnd) = trc_apnd
       endif
       
       if (l_stop) then
          istop = i
          jstop = j
          if (my_task == master_task) then
             write(nu_diag,*) "ice_therm_mushy solver failure: istep, n, i, j:", istep, ncat, i, j
          endif
          call abort_ice("ice_therm_mushy solver failure")
       endif
       
    enddo ! ij

  end subroutine temperature_changes_salinity

!=======================================================================
! Solver with separation of cold and melting phases
!=======================================================================

  subroutine temperature_changes_column(istep1,   dt,       & 
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
                                        hpond,    apond,    &
                                        Tsf,      Tbot,     &
                                        sss,                &
                                        fsensn,   flatn,    &
                                        fswabsn,  flwoutn,  &
                                        fsurfn,             &
                                        fcondtop, fcondbot, &
                                        fadvheat, snoice,   &
                                        einit_old, lstop)

    ! solve the enthalpy and bulk salinity of the ice for a single column

    integer (kind=int_kind), intent(in) :: &
         istep1          ! time step index (diagnostic only)
    
    real (kind=dbl_kind), intent(in) :: &
         dt              ! time step (s)
    
    real (kind=dbl_kind), intent(in) :: &
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         Tbot        , & ! ice bottom surfce temperature (deg C)
         sss             ! sea surface salinity (PSU)
         
    real (kind=dbl_kind), intent(inout) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         fswint      , & ! SW absorbed in ice interior below surface (W m-2)
         fswthrun        ! SW through ice to ocean  (W m-2)
    
    real (kind=dbl_kind), intent(inout) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr           ! snow layer thickness (m)

    real (kind=dbl_kind), intent(in) :: &
         einit_old       ! initial energy of melting (J m-2)
    
    real (kind=dbl_kind), dimension (nslyr), &
         intent(inout) :: &
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)
    
    real (kind=dbl_kind), dimension (nilyr), &
         intent(inout) :: &
         Iswabs          ! SW radiation absorbed in ice layers (W m-2)
    
    real (kind=dbl_kind), intent(inout):: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtop    , & ! downward cond flux at top surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fswabsn     , & ! shortwave absorbed by ice (W m-2)
         flwoutn         ! upward LW at surface (W m-2)
    
    real (kind=dbl_kind), intent(out):: &
         fcondbot    , & ! downward cond flux at bottom surface (W m-2)
         fadvheat    , & ! flow of heat to ocean due to advection (W m-2)
         snoice          ! snow ice formation

    real (kind=dbl_kind), intent(inout):: &
         Tsf         , & ! ice/snow surface temperature (C)
         hpond       , & ! melt pond depth (m)
         apond           ! melt pond area
    
    real (kind=dbl_kind), dimension (nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3)
         Tin         , & ! internal ice layer temperatures
         Sin             ! internal ice layer salinities

    real (kind=dbl_kind), dimension (nslyr), intent(inout) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         Tsn             ! internal snow layer temperatures
    
    logical (kind=log_kind), intent(inout) :: &
         lstop           ! solver failure flag 
    
    ! local variables
    real(kind=dbl_kind), dimension(1:nilyr) :: &
         qin0        , & ! ice layer enthalpy (J m-3) at start of timestep
         Tin0        , & ! internal ice layer temperatures (C) at start of timestep
         Sin0        , & ! internal ice layer salinities (ppt) at start of timestep
         phi         , & ! liquid fraction
         km          , & ! ice conductivity (W m-1 K-1)
         dSdt            ! gravity drainage desalination rate for slow mode (ppt s-1)

    real(kind=dbl_kind), dimension(1:nilyr+1) :: &
         Sbr         , & ! brine salinity (ppt)
         qbr             ! brine enthalpy (J m-3)

    real(kind=dbl_kind), dimension(0:nilyr) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         qsn0        , & ! snow layer enthalpy (J m-3) at start of timestep
         Tsn0        , & ! internal snow layer temperatures (C) at start of timestep
         ks              ! snow conductivity (W m-1 K-1)

    real(kind=dbl_kind) :: &
         Tsf0        , & ! ice/snow surface temperature (C) at start of timestep
         hin         , & ! ice thickness (m)
         hsn         , & ! snow thickness (m)
         hslyr_min   , & ! minimum snow layer thickness (m)
         w           , & ! vertical flushing Darcy velocity (m/s)
         wdown       , & ! downwards vertical flushing Darcy velocity (m/s)
         wup         , & ! upwards vertical flushing Darcy velocity (m/s)
         w_scaling   , & !
         einit       , & ! initial total energy (J)
         efinal      , & ! final total energy (J)
         sinit       , & ! initial total salt content (ppt m)
         sfinal      , & ! final total salt content
         fadvsalt    , & ! heat flux to ocean from brine advection (W m-2)
         snowice_en  , & ! energy added as snowice
         snowice_st  , & ! salt added as snowice
         dhhead      , & ! hydraulic head (m)
         w_up_max    , & 
         qocn        , & ! ocean brine enthalpy (J m-3)
         qpond       , & ! melt pond brine enthalpy (J m-3)
         Spond           ! melt pond salinity (ppt)

    integer(kind=int_kind) :: &
         k               ! ice/snow layer index

    logical(kind=log_kind) :: &
         lsnow           ! snow presence: T: has snow, F: no snow

    lstop = .false.
    fadvsalt = c0
    snoice = c0

#if defined oned
    if (hilyr <= c0) return
#endif

    Tsf0 = Tsf
    qsn0 = qsn
    qin0 = qin
    Sin0 = Sin
    Tsn0 = Tsn
    Tin0 = Tin

    Spond = 0.0_dbl_kind
    qpond = enthalpy_brine(c0)

    hslyr_min = hs_min / real(nslyr, dbl_kind)

    lsnow = (hslyr > hslyr_min)

    hin = hilyr * real(nilyr,dbl_kind)

    qocn = enthalpy_brine(Tbot)

    if (lsnow) then
       hsn = hslyr * real(nslyr,dbl_kind)
    else
       hsn = c0
    endif

    do k = 1, nilyr
       phi(k) = liquid_fraction(temperature_mush(qin(k),Sin(k)),Sin(k))
    enddo ! k

#if flushing == 1
    ! calculate vertical bulk darcy flow
    call flushing_velocity(Tin,    Sin,   &
                           phi,           &
                           hin,    hsn,   &
                           hilyr,         &
                           hpond,  apond, & 
                           dt,     w,     &
                           wdown,  wup,   &
                           dhhead, w_up_max)

!    open(55,file='history/w.txt',position='append')
!    write(55,*) istep, abs(max(w,c0)), abs(min(w,c0))
!    close(55)
#else
    w = c0
    wdown = c0
    wup = c0
#endif

    ! calculate quantities related to drainage
    call explicit_flow_velocities(Sin,    qin,    &
                                  Tin,    Tsf,    &
                                  Tbot,   q,      &
                                  dSdt,   Sbr,    &
                                  qbr,    dt,     &
                                  sss,    qocn,   &
                                  hilyr,  hin)
    !q = c0 ; dSdt = c0
    
#if diag == 1
    call diagnose_salt_flux(q, dSdt)
#endif

    ! calculate the conductivities
    call conductivity_mush_array(qin0, Sin0, km)

    if (lsnow) then
       ! case with snow

       ! calculate the snow conductivities
       call conductivity_snow_array(ks)

       ! run the two stage solver
       call two_stage_solver_snow(Tsf,         Tsf0,       &
                                  qsn,         qsn0,       &
                                  qin,         qin0,       &
                                  Sin,         Sin0,       &
                                  Tsn,         Tsn0,       &
                                  Tin,         Tin0,       &
                                  phi,         Tbot,       &
                                  km,          ks,         &
                                  q,           dSdt,       &
                                  w,                       &
                                  fswint,      fswsfc,     &
                                  rhoa,        flw,        &
                                  potT,        Qa,         &
                                  shcoef,      lhcoef,     &
                                  Iswabs,      Sswabs,     &
                                  qpond,       qocn,       &
                                  Spond,       sss,        &
                                  hilyr,       hslyr,      &
                                  fcondtop,    fcondbot,   &
                                  fadvheat,                &
                                  flwoutn,     fsensn,     &
                                  flatn,       fsurfn,     &
                                  lstop)

       ! given the updated enthalpy and bulk salinity calculate other quantities
       do k = 1, nslyr
          Tsn(k) = temperature_snow(qsn(k))
       enddo ! k

       do k = 1, nilyr
          Tin(k) = temperature_mush_liquid_fraction(qin(k), phi(k))
          Sbr(k) = liquidus_brine_salinity_mush(Tin(k)) 
          qbr(k) = enthalpy_brine(Tin(k))
       enddo ! k
          
    else
       ! case without snow

       ! run the two stage solver
       call two_stage_solver_nosnow(Tsf,         Tsf0,       &
                                    qsn,         qsn0,       &
                                    qin,         qin0,       &
                                    Sin,         Sin0,       &
                                    Tsn,         Tsn0,       &
                                    Tin,         Tin0,       &
                                    phi,         Tbot,       &
                                    km,          ks,         &
                                    q,           dSdt,       &
                                    w,                       &
                                    fswint,      fswsfc,     &
                                    rhoa,        flw,        &
                                    potT,        Qa,         &
                                    shcoef,      lhcoef,     &
                                    Iswabs,      Sswabs,     &
                                    qpond,       qocn,       &
                                    Spond,       sss,        &
                                    hilyr,       hslyr,      &
                                    fcondtop,    fcondbot,   &
                                    fadvheat,                &
                                    flwoutn,     fsensn,     &
                                    flatn,       fsurfn,     &
                                    lstop)

       ! given the updated enthalpy and bulk salinity calculate other quantities
       do k = 1, nilyr
          Tin(k) = temperature_mush_liquid_fraction(qin(k), phi(k))
          Sbr(k) = liquidus_brine_salinity_mush(Tin(k)) 
          qbr(k) = enthalpy_brine(Tin(k))
       enddo ! k
       
    endif

    if (lstop) then
       return
    end if

#if flushing == 1

    call flushing_advection_flow_restriction(w_scaling,   &
                                             wdown,  wup, &
                                             Sin0,   Sbr, &
                                             sss,    Spond, &
                                             dt)

    w = w * w_scaling

    ! now perform final flushing related changes to the ice
    ! drain ponds from flushing 
    call flush_pond(w, hin, hpond, apond, dt)
       
    ! flood snow ice
   ! call flood_ice(w,          dt,       &
   !                dhhead,     w_up_max, &
   !                hsn,        hin,      &
   !                hslyr,      hilyr,    & 
   !                qsn,        qin,      &
   !                phi,                  &
   !                Sin,        Sbr,      &
   !                qbr,        snoice,   &
   !                snowice_en, snowice_st)

    fadvheat = fadvheat - snowice_en

#endif

  end subroutine temperature_changes_column

!=======================================================================

  subroutine two_stage_solver_snow(Tsf,         Tsf0,       &
                                   qsn,         qsn0,       &
                                   qin,         qin0,       &
                                   Sin,         Sin0,       &
                                   Tsn,         Tsn0,       &
                                   Tin,         Tin0,       &
                                   phi,         Tbot,       &
                                   km,          ks,         &
                                   q,           dSdt,       &
                                   w,                       &
                                   fswint,      fswsfc,     &
                                   rhoa,        flw,        &
                                   potT,        Qa,         &
                                   shcoef,      lhcoef,     &
                                   Iswabs,      Sswabs,     &
                                   qpond,       qocn,       &
                                   Spond,       sss,        &
                                   hilyr,       hslyr,      &
                                   fcondtop,    fcondbot,   &
                                   fadvheat,                &
                                   flwoutn,     fsensn,     &
                                   flatn,       fsurfn,     &
                                   lstop)

    ! solve the vertical temperature and salt change for case with snow
    ! 1) determine what type of surface condition existed previously - cold or melting
    ! 2) solve the system assuming this condition persists
    ! 3) check the consistency of the surface condition of the solution
    ! 4) If the surface condition is inconsistent resolve for the other surface condition
    ! 5) If neither solution is consistent the resolve the inconsistency

    real(kind=dbl_kind), intent(inout) :: &
         Tsf             ! snow surface temperature (C)

    real(kind=dbl_kind), intent(out) :: &
         fcondtop    , & ! downward cond flux at top surface (W m-2)
         fcondbot    , & ! downward cond flux at bottom surface (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtop
         fadvheat        ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), intent(in) :: &
         Tsf0            ! snow surface temperature (C) at beginning of timestep

    real(kind=dbl_kind), dimension(1:nslyr), intent(inout) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         Tsn             ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         qsn0        , & ! snow layer enthalpy (J m-3) at beginning of timestep
         Tsn0        , & ! snow layer temperature (C) at beginning of timestep
         ks          , & ! snow conductivity (W m-1 K-1)
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), dimension(1:nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         Sin         , & ! ice layer bulk salinity (ppt)
         Tin         , & ! ice layer temperature (C)
         phi             ! ice layer liquid fraction

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         qin0        , & ! ice layer enthalpy (J m-3) at beginning of timestep
         Sin0        , & ! ice layer bulk salinity (ppt) at beginning of timestep
         Tin0        , & ! ice layer temperature (C) at beginning of timestep
         km          , & ! ice conductivity (W m-1 K-1)
         Iswabs      , & ! SW radiation absorbed in ice layers (W m-2)
         dSdt            ! gravity drainage desalination rate for slow mode (ppt s-1)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         Tbot        , & ! ice bottom surfce temperature (deg C)
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         fswint      , & ! SW absorbed in ice interior below surface (W m-2)
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         w           , & ! vertical flushing Darcy velocity (m/s)
         qpond       , & ! melt pond brine enthalpy (J m-3)
         qocn        , & ! ocean brine enthalpy (J m-3)
         Spond       , & ! melt pond salinity (ppt)
         sss             ! sea surface salinity (PSU)

    logical(kind=log_kind), intent(inout) :: &
         lstop           ! solver failure flag

    ! determine if surface is initially cold or melting
    if (Tsf < c0) then

       ! initially cold

#if debug == 1
       write(*,*) istep1, "Two stage snow cold: 1"
#endif

       ! solve the system for cold and snow
       call picard_solver(.true.,  .true.,    &
                          Tsf,      qsn,      &
                          qin,      Sin,      &
                          Tin,      Tsn,      &
                          phi,                &
                          hilyr,    hslyr,    &
                          km,       ks,       &
                          Iswabs,   Sswabs,   &
                          Tbot,               &
                          fswint,   fswsfc,   &
                          rhoa,     flw,      &
                          potT,     Qa,       &
                          shcoef,   lhcoef,   &
                          fcondtop, fcondbot, &
                          fadvheat,           &
                          flwoutn,  fsensn,   &
                          flatn,    fsurfn,   &
                          qpond,    qocn,     &
                          Spond,    sss,      &
                          q,        dSdt,     &
                          w,        lstop)

       ! halt if solver failed
       if (lstop) return

#if debug == 1
       write(*,*) istep1, "Two stage snow: 1", Tsf, Tsf < c0
#endif

       ! check if solution is consistent - surface should still be cold
       if (Tsf < dTsf_errcon) then

          ! solution is consistent - have solution so finish
          return

       else
          
          ! solution is inconsistent - surface is warmer than melting 
          ! resolve assuming surface is melting

#if debug == 1
          write(*,*) istep1, "Two stage snow melt: 2"
#endif

          ! reset the solution to initial values
          Tsf = c0
          Tsf = Tsf0
          qsn = qsn0
          qin = qin0
          Sin = Sin0

          ! solve the system for melting and snow             
          call picard_solver(.true.,   .false.,  &
                             Tsf,      qsn,      &
                             qin,      Sin,      &
                             Tin,      Tsn,      &
                             phi,                &
                             hilyr,    hslyr,    &
                             km,       ks,       &
                             Iswabs,   Sswabs,   &
                             Tbot,               &
                             fswint,   fswsfc,   &
                             rhoa,     flw,      &
                             potT,     Qa,       &
                             shcoef,   lhcoef,   &
                             fcondtop, fcondbot, &
                             fadvheat,           &
                             flwoutn,  fsensn,   &
                             flatn,    fsurfn,   &
                             qpond,    qocn,     &
                             Spond,    sss,      &
                             q,        dSdt,     &
                             w,        lstop)

          ! halt if solver failed
          if (lstop) return

#if debug == 1
          write(*,*) istep1, "Two stage snow: 2", fsurfn, fcondtop, fsurfn > fcondtop
#endif

          ! check if solution is consistent - surface conductive heat flux should be less than incoming surface heat flux
          if (fcondtop - fsurfn < ferrcon) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent

             write(*,*) istep1, "inconsistency! 1" ; stop

          endif ! surface flux consistency

       endif ! surface temperature consistency

    else

       ! initially melting

#if debug == 1
       write(*,*) istep1, "Two stage snow melt: 3"
#endif

       Tsf = c0

       ! solve the system for melting and snow
       call picard_solver(.true.,   .false.,  &
                          Tsf,      qsn,      &
                          qin,      Sin,      &
                          Tin,      Tsn,      &
                          phi,                &
                          hilyr,    hslyr,    &
                          km,       ks,       &
                          Iswabs,   Sswabs,   &
                          Tbot,               &
                          fswint,   fswsfc,   &
                          rhoa,     flw,      &
                          potT,     Qa,       &
                          shcoef,   lhcoef,   &
                          fcondtop, fcondbot, &
                          fadvheat,           &
                          flwoutn,  fsensn,   &
                          flatn,    fsurfn,   &
                          qpond,    qocn,     &
                          Spond,    sss,      &
                          q,        dSdt,     &
                          w,        lstop)

       ! halt if solver failed
       if (lstop) return
       
#if debug == 1
       write(*,*) istep1, "Two stage snow: 3", fsurfn, fcondtop, fsurfn > fcondtop
#endif
       
       ! check if solution is consistent - surface conductive heat flux should be less than incoming surface heat flux
       if (fcondtop - fsurfn < ferrcon) then

          ! solution is consistent - have solution so finish
          return

       else

          ! solution is inconsistent - resolve assuming other surface condition
          ! assume surface is cold

#if debug == 1
          write(*,*) istep1, "Two stage snow cold: 4"
#endif

          ! reset the solution to initial values
          Tsf = Tsf0
          qsn = qsn0
          qin = qin0
          Sin = Sin0

          ! solve the system for cold and snow      
          call picard_solver(.true.,   .true.,   &       
                             Tsf,      qsn,      &
                             qin,      Sin,      &
                             Tin,      Tsn,      &
                             phi,                &
                             hilyr,    hslyr,    & 
                             km,       ks,       &
                             Iswabs,   Sswabs,   &
                             Tbot,               &
                             fswint,   fswsfc,   &
                             rhoa,     flw,      &
                             potT,     Qa,       &
                             shcoef,   lhcoef,   &
                             fcondtop, fcondbot, &
                             fadvheat,           &
                             flwoutn,  fsensn,   &
                             flatn,    fsurfn,   &
                             qpond,    qocn,     &
                             Spond,    sss,      &
                             q,        dSdt,     &
                             w,        lstop)

          ! halt if solver failed
          if (lstop) return

#if debug == 1
          write(*,*) istep1, "Two stage snow: 4", Tsf, Tsf < c0
#endif

          ! check if solution is consistent - surface should be cold
          if (Tsf < dTsf_errcon) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent
             ! failed to find a solution so need to refine solutions until consistency found
             
             write(*,*) istep1, "inconsistency! 2" ; stop

          endif ! surface temperature consistency
          
       endif ! surface flux consistency

    endif

  end subroutine two_stage_solver_snow

!=======================================================================

  subroutine two_stage_solver_nosnow(Tsf,         Tsf0,       &
                                     qsn,         qsn0,       &
                                     qin,         qin0,       &
                                     Sin,         Sin0,       &
                                     Tsn,         Tsn0,       &
                                     Tin,         Tin0,       &
                                     phi,         Tbot,       &
                                     km,          ks,         &
                                     q,           dSdt,       &
                                     w,                       &
                                     fswint,      fswsfc,     &
                                     rhoa,        flw,        &
                                     potT,        Qa,         &
                                     shcoef,      lhcoef,     &
                                     Iswabs,      Sswabs,     &
                                     qpond,       qocn,       &
                                     Spond,       sss,        &
                                     hilyr,       hslyr,      &
                                     fcondtop,    fcondbot,   &
                                     fadvheat,                &
                                     flwoutn,     fsensn,     &
                                     flatn,       fsurfn,     &
                                     lstop)
    
    ! solve the vertical temperature and salt change for case with no snow
    ! 1) determine what type of surface condition existed previously - cold or melting
    ! 2) solve the system assuming this condition persists
    ! 3) check the consistency of the surface condition of the solution
    ! 4) If the surface condition is inconsistent resolve for the other surface condition
    ! 5) If neither solution is consistent the resolve the inconsistency

    real(kind=dbl_kind), intent(inout) :: &
         Tsf             ! ice surface temperature (C)

    real(kind=dbl_kind), intent(out) :: &
         fcondtop    , & ! downward cond flux at top surface (W m-2)
         fcondbot    , & ! downward cond flux at bottom surface (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtop
         fadvheat        ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), intent(in) :: &
         Tsf0            ! ice surface temperature (C) at beginning of timestep

    real(kind=dbl_kind), dimension(1:nslyr), intent(inout) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         Tsn             ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         qsn0        , & ! snow layer enthalpy (J m-3) at beginning of timestep
         Tsn0        , & ! snow layer temperature (C) at beginning of timestep
         ks          , & ! snow conductivity (W m-1 K-1)
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), dimension(1:nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         Sin         , & ! ice layer bulk salinity (ppt)
         Tin         , & ! ice layer temperature (C)
         phi             ! ice layer liquid fraction

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         qin0        , & ! ice layer enthalpy (J m-3) at beginning of timestep
         Sin0        , & ! ice layer bulk salinity (ppt) at beginning of timestep
         Tin0        , & ! ice layer temperature (C) at beginning of timestep
         km          , & ! ice conductivity (W m-1 K-1)
         Iswabs      , & ! SW radiation absorbed in ice layers (W m-2)
         dSdt            ! gravity drainage desalination rate for slow mode (ppt s-1)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         Tbot        , & ! ice bottom surfce temperature (deg C)
         fswint      , & ! SW absorbed in ice interior below surface (W m-2)
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         w           , & ! vertical flushing Darcy velocity (m/s)
         qpond       , & ! melt pond brine enthalpy (J m-3)
         qocn        , & ! ocean brine enthalpy (J m-3)
         Spond       , & ! melt pond salinity (ppt)
         sss             ! sea surface salinity (PSU)

    logical, intent(inout) :: &
         lstop           ! solver failure flag

    real(kind=dbl_kind) :: &
         Tmlt            ! upper ice layer melting temperature (C)

    ! initial surface melting temperature
    Tmlt = liquidus_temperature_mush(Sin(1))

    ! determine if surface is initially cold or melting
    if (Tsf < Tmlt) then

       ! initially cold

#if debug == 1
       write(*,*) istep1, "Two stage nosnow cold: 1"
#endif

       ! solve the system for cold and no snow          
       call picard_solver(.false.,  .true.,   &
                          Tsf,      qsn,      &
                          qin,      Sin,      &
                          Tin,      Tsn,      &
                          phi,                &
                          hilyr,    hslyr,    &
                          km,       ks,       &
                          Iswabs,   Sswabs,   &
                          Tbot,               &
                          fswint,   fswsfc,   &
                          rhoa,     flw,      &
                          potT,     Qa,       &
                          shcoef,   lhcoef,   &
                          fcondtop, fcondbot, &
                          fadvheat,           &
                          flwoutn,  fsensn,   &
                          flatn,    fsurfn,   &
                          qpond,    qocn,     &
                          Spond,    sss,      &
                          q,        dSdt,     &
                          w,        lstop)

       ! halt if solver failed
       if (lstop) return

       ! surface melting temperature
       Tmlt = liquidus_temperature_mush(Sin(1))

#if debug == 1
       write(*,*) istep1, "Two stage nosnow: 1", Tsf, Tmlt, Tsf < Tmlt, Tsf - Tmlt
#endif

       ! check if solution is consistent - surface should still be cold
       if (Tsf < Tmlt + dTsf_errcon) then

          ! solution is consistent - have solution so finish
          return

       else
          ! solution is inconsistent - surface is warmer than melting 
          ! resolve assuming surface is melting

#if debug == 1
          write(*,*) istep1, "Two stage nosnow melt: 2: "
#endif

          ! reset the solution to initial values
          Tsf = liquidus_temperature_mush(Sin(1))
          !Tsf = Tsf0
          qin = qin0
          Sin = Sin0

          ! solve the system for melt and no snow
          call picard_solver(.false.,  .false.,  &
                             Tsf,      qsn,      &
                             qin,      Sin,      &
                             Tin,      Tsn,      &
                             phi,                &
                             hilyr,    hslyr,    &
                             km,       ks,       &
                             Iswabs,   Sswabs,   &
                             Tbot,               &
                             fswint,   fswsfc,   &
                             rhoa,     flw,      & 
                             potT,     Qa,       &
                             shcoef,   lhcoef,   &
                             fcondtop, fcondbot, &
                             fadvheat,           &
                             flwoutn,  fsensn,   &
                             flatn,    fsurfn,   &
                             qpond,    qocn,     &
                             Spond,    sss,      &
                             q,        dSdt,     &
                             w,        lstop)

          ! halt if solver failed
          if (lstop) return

#if debug == 1
          write(*,*) istep1, "Two stage nosnow: 2", fsurfn, fcondtop, fsurfn > fcondtop, fsurfn - fcondtop
#endif

          ! check if solution is consistent - surface conductive heat flux should be less than incoming surface heat flux
          if (fcondtop - fsurfn < ferrcon) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent

             write(*,*) istep1, "inconsistency! 3" ; stop

          endif

       endif

    else
       ! initially melting

#if debug == 1
       write(*,*) istep1, "Two stage nosnow melt: 3"
#endif

       ! solve the system for melt and no snow
       
       Tsf = liquidus_temperature_mush(Sin(1))
       
       call picard_solver(.false.,  .false.,  &
                          Tsf,      qsn,      &
                          qin,      Sin,      &
                          Tin,      Tsn,      &
                          phi,                &
                          hilyr,    hslyr,    &
                          km,       ks,       &
                          Iswabs,   Sswabs,   &
                          Tbot,               &
                          fswint,   fswsfc,   &
                          rhoa,     flw,      &
                          potT,     Qa,       &
                          shcoef,   lhcoef,   &
                          fcondtop, fcondbot, &
                          fadvheat,           &
                          flwoutn,  fsensn,   &
                          flatn,    fsurfn,   &
                          qpond,    qocn,     &
                          Spond,    sss,      &
                          q,        dSdt,     &
                          w,        lstop)

       ! halt if solver failed
       if (lstop) return

#if debug == 1
       write(*,*) istep1, "Two stage nosnow: 3", fsurfn, fcondtop, fsurfn > fcondtop
#endif

       ! check if solution is consistent - surface conductive heat flux should be less than incoming surface heat flux
       if (fcondtop - fsurfn < ferrcon) then

          ! solution is consistent - have solution so finish
          return

       else

          ! solution is inconsistent - resolve assuming other surface condition
          ! assume surface is cold

#if debug == 1
          write(*,*) istep1, "Two stage nosnow cold: 4"
#endif

          ! reset the solution to initial values
          Tsf = Tsf0
          qin = qin0
          Sin = Sin0

          ! solve the system for cold and no snow
          call picard_solver(.false.,  .true.,   &
                             Tsf,      qsn,      &
                             qin,      Sin,      &
                             Tin,      Tsn,      &
                             phi,                &
                             hilyr,    hslyr,    &
                             km,       ks,       &
                             Iswabs,   Sswabs,   &
                             Tbot,               &
                             fswint,   fswsfc,   &
                             rhoa,     flw,      &
                             potT,     Qa,       &
                             shcoef,   lhcoef,   &
                             fcondtop, fcondbot, &
                             fadvheat,           &
                             flwoutn,  fsensn,   &
                             flatn,    fsurfn,   &
                             qpond,    qocn,     &
                             Spond,    sss,      &
                             q,        dSdt,     &
                             w,        lstop)

          ! halt if solver failed
          if (lstop) return

#if debug == 1
          write(*,*) istep1, "Two stage nosnow: 4", Tsf, Tsf < c0
#endif

          ! surface melting temperature
          Tmlt = liquidus_temperature_mush(Sin(1))

          ! check if solution is consistent - surface should be cold
          if (Tsf < Tmlt + dTsf_errcon) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent

             write(*,*) istep1, "inconsistency! 4" ; stop

          endif
          
       endif

    endif

  end subroutine two_stage_solver_nosnow

!=======================================================================
! Picard/TDMA based solver
!=======================================================================

  subroutine prep_picard(lsnow, qsn,    &
                         qin,   Sin,    &
                         hilyr, hslyr,  &
                         km,    ks,     &
                         Tin,   Tsn,    &
                         Sbr,   phi,    &
                         dxp,   kcstar, &
                         einit)
    
    logical, intent(in) :: &
         lsnow      ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qin    , & ! ice layer enthalpy (J m-3)
         Sin    , & ! ice layer bulk salinity (ppt)
         km         ! ice conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         qsn    , & ! snow layer enthalpy (J m-3)
         ks         ! snow conductivity (W m-1 K-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr  , & ! ice layer thickness (m)
         hslyr      ! snow layer thickness (m)

    real(kind=dbl_kind), dimension(nilyr), intent(out) :: &
         Tin    , & ! ice layer temperature (C)
         Sbr    , & ! ice layer brine salinity (ppt)
         phi        ! ice layer liquid fraction

    real(kind=dbl_kind), dimension(nslyr), intent(out) :: &
         Tsn        ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(out) :: &
         dxp    , & ! distances between grid points (m)
         kcstar     ! interface conductivities  (W m-1 K-1)

    real(kind=dbl_kind), intent(out) :: &
         einit      ! initial total energy (J)

    integer :: &
         k          ! vertical layer index

    ! calculate initial ice temperatures
    do k = 1, nilyr
       Tin(k) = temperature_mush(qin(k), Sin(k))
       Sbr(k) = liquidus_brine_salinity_mush(Tin(k))
       phi(k) = liquid_fraction(Tin(k), Sin(k))
    enddo ! k

    if (lsnow) then

       do k = 1, nslyr
          Tsn(k) = temperature_snow(qsn(k))
       enddo ! k

    endif ! lsnow

    ! interface distances
    call calc_intercell_thickness(lsnow, hilyr, hslyr, dxp)

    ! interface conductivities
    call calc_intercell_conductivity(lsnow, km, ks, hilyr, hslyr, kcstar)

    ! total energy content
    call total_energy_content(lsnow,        &
                              qin,   qsn,   &
                              hilyr, hslyr, &
                              einit)

  end subroutine prep_picard

!=======================================================================

  subroutine picard_solver(lsnow,    lcold,    &
                           Tsf,      qsn,      &
                           qin,      Sin,      &
                           Tin,      Tsn,      &
                           phi,                &
                           hilyr,    hslyr,    &
                           km,       ks,       &
                           Iswabs,   Sswabs,   &
                           Tbot,               &
                           fswint,   fswsfc,   &
                           rhoa,     flw,      &
                           potT,     Qa,       &
                           shcoef,   lhcoef,   &
                           fcondtop, fcondbot, &
                           fadvheat,           &
                           flwoutn,  fsensn,   &
                           flatn,    fsurfn,   &
                           qpond,    qocn,     &
                           Spond,    sss,      &
                           q,        dSdt,     &
                           w,        lstop)

    use ice_therm_shared, only: surface_heat_flux, dsurface_heat_flux_dTsf

    logical, intent(in) :: &
         lsnow        , & ! snow presence: T: has snow, F: no snow
         lcold            ! surface cold: T: surface is cold, F: surface is melting

    real(kind=dbl_kind), intent(inout) :: &
         Tsf              ! snow surface temperature (C)

    real(kind=dbl_kind), intent(out) :: &
         fcondtop     , & ! downward cond flux at top surface (W m-2)
         fcondbot     , & ! downward cond flux at bottom surface (W m-2)
         fadvheat         ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), dimension(nilyr), intent(inout) :: &
         qin          , & ! ice layer enthalpy (J m-3)
         Sin          , & ! ice layer bulk salinity (ppt)
         Tin          , & ! ice layer temperature (C)
         phi              ! ice layer liquid fraction

    real(kind=dbl_kind), dimension(nslyr), intent(inout) :: &
         qsn          , & ! snow layer enthalpy (J m-3)
         Tsn              ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         km           , & ! ice conductivity (W m-1 K-1)
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         dSdt             ! gravity drainage desalination rate for slow mode (ppt s-1)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         ks           , & ! snow conductivity (W m-1 K-1)
         Sswabs           ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), intent(out) :: &
         flwoutn      , & ! upward LW at surface (W m-2)
         fsensn       , & ! surface downward sensible heat (W m-2)
         flatn        , & ! surface downward latent heat (W m-2)
         fsurfn           ! net flux to top surface, excluding fcondtop

    real(kind=dbl_kind), intent(in) :: &
         hilyr        , & ! ice layer thickness (m)
         hslyr        , & ! snow layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         fswint       , & ! SW absorbed in ice interior below surface (W m-2)
         fswsfc       , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa         , & ! air density (kg/m^3)
         flw          , & ! incoming longwave radiation (W/m^2)
         potT         , & ! air potential temperature (K)
         Qa           , & ! specific humidity (kg/kg)
         shcoef       , & ! transfer coefficient for sensible heat
         lhcoef       , & ! transfer coefficient for latent heat
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         Spond        , & ! melt pond salinity (ppt)
         sss          , & ! sea surface salinity (ppt)
         w                ! vertical flushing Darcy velocity (m/s)
      
    logical(kind=log_kind), intent(inout) :: &
         lstop            ! solver failure flag 

    real(kind=dbl_kind), dimension(nilyr) :: &
         Sbr          , & ! ice layer brine salinity (ppt)
         qbr          , & ! ice layer brine enthalpy (J m-3)
         Tin0         , & ! ice layer temperature (C) at start of timestep
         qin0         , & ! ice layer enthalpy (J m-3) at start of timestep
         Sin0         , & ! ice layer bulk salinity (ppt) at start of timestep
         Tin_prev         ! ice layer temperature at previous iteration
    
    real(kind=dbl_kind), dimension(nslyr) :: &
         qsn0         , & ! snow layer enthalpy (J m-3) at start of timestep
         Tsn0         , & ! snow layer temperature (C) at start of timestep
         Tsn_prev         ! snow layer temperature at previous iteration

    real(kind=dbl_kind), dimension(nslyr+nilyr+1) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    real(kind=dbl_kind) :: &
         Tsf0         , & ! snow surface temperature (C) at start of timestep
         dfsurfn_dTsf , & ! derivative of net flux to top surface, excluding fcondtopn
         Tsf_prev     , & ! snow surface temperature at previous iteration
         einit        , & ! initial total energy (J)
         fadvheat_nit     ! heat to ocean due to advection (W m-2) during iteration

    logical :: &
         lconverged       ! has Picard solver converged?

    integer :: &
         nit          , & ! Picard iteration count
         k                ! vertical layer index

    integer, parameter :: &
         nit_max = 100    ! maximum number of Picard iterations

    real(kind=dbl_kind) :: tmpflux

    lconverged = .false.

    ! prepare quantities for picard iteration
    call prep_picard(lsnow, qsn,    &
                     qin,   Sin,    &
                     hilyr, hslyr,  &
                     km,    ks,     &
                     Tin,   Tsn,    &
                     Sbr,   phi,    &
                     dxp,   kcstar, &
                     einit)

    Tsf0 = Tsf
    qin0 = qin
    qsn0 = qsn
    Tin0 = Tin
    Tsn0 = Tsn
    Sin0 = Sin

    ! set prev variables
    Tsf_prev = Tsf
    Tsn_prev = Tsn
    Tin_prev = Tin

    ! picard iteration
    picard: do nit = 1, nit_max

       ! surface heat flux
       call surface_heat_flux(Tsf,     fswsfc, &
                              rhoa,    flw,    &
                              potT,    Qa,     &
                              shcoef,  lhcoef, &
                              flwoutn, fsensn, &
                              flatn,   fsurfn)

       ! derivative of heat flux with respect to surface temperature
       call dsurface_heat_flux_dTsf(Tsf,     fswsfc, &
                                    rhoa,    flw,    &
                                    potT,    Qa,     &
                                    shcoef,  lhcoef, &
                                    dfsurfn_dTsf, &
                                    tmpflux, tmpflux, tmpflux)

       ! tridiagonal solve of new temperatures
       call solve_heat_conduction(lsnow,     lcold,        &
                                  Tsf,       Tbot,         &
                                  qin0,      qsn0,         &
                                  phi,       dt,           &
                                  qpond,     qocn,         &
                                  q,         w,            &
                                  hilyr,     hslyr,        &
                                  dxp,       kcstar,       &
                                  Iswabs,    Sswabs,       &
                                  fsurfn,    dfsurfn_dTsf, &
                                  Tin,       Tsn,nit)

       ! update brine enthalpy
       call picard_updates_enthalpy(Tin, qbr)

       ! drainage fluxes
       call picard_drainage_fluxes(fadvheat_nit, q,    &
                                   qbr,          qocn, &
                                   hilyr)

       ! flushing fluxes
       call picard_flushing_fluxes(fadvheat_nit, w, &
                                   qbr,             &
                                   qocn, qpond)

       ! perform convergence check
       call check_picard_convergence(lsnow,                &
                                     lconverged, nit,      & 
                                     Tsf,        Tsf_prev, &
                                     Tin,        Tin_prev, &
                                     Tsn,        Tsn_prev, &
                                     phi,        Tbot,     &
                                     qin,        qsn,      &
                                     km,         ks,       &
                                     hilyr,      hslyr,    &
                                     fswint,               &
                                     einit,                &
                                     fcondtop,   fcondbot, &
                                     fadvheat_nit)

       if (lconverged) exit

       Tsf_prev = Tsf
       Tsn_prev = Tsn
       Tin_prev = Tin

    enddo picard

    fadvheat = fadvheat_nit

    ! update the picard iterants
    call picard_updates(Tin, &
                        Sbr, qbr)
    
    ! solve for the salinity
    call solve_salinity(Sin,   Sbr,   &
                        Spond, sss,   &
                        q,     dSdt,  &
                        w,     hilyr, &
                        dt)

    ! final surface heat flux
    call surface_heat_flux(Tsf,     fswsfc, &
                           rhoa,    flw,    &
                           potT,    Qa,     &
                           shcoef,  lhcoef, &
                           flwoutn, fsensn, &
                           flatn,   fsurfn)

    ! if not converged
    if (.not. lconverged) then

       call picard_nonconvergence(Tsf0, Tsf, &
                                  Tsn0, Tsn, &
                                  Tin0, Tin, &
                                  Sin0, Sin, &
                                  qsn0, qsn, &
                                  qin0, qin, &
                                  phi)
       lstop = .true.

    endif

  end subroutine picard_solver

!=======================================================================

  subroutine picard_nonconvergence(Tsf0, Tsf, &
                                   Tsn0, Tsn, &
                                   Tin0, Tin, &
                                   Sin0, Sin, &
                                   qsn0, qsn, &
                                   qin0, qin, &
                                   phi)

    real(kind=dbl_kind), intent(in) :: &
         Tsf0 , & ! snow surface temperature (C) at beginning of timestep
         Tsf      ! snow surface temperature (C)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         Tsn0 , & ! snow layer temperature (C) at beginning of timestep
         Tsn  , & ! snow layer temperature (C)
         qsn0 , &
         qsn

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Tin0 , & ! ice layer temperature (C)
         Tin  , & ! ice layer temperature (C)
         Sin0 , & ! ice layer bulk salinity (ppt)
         Sin  , & ! ice layer bulk salinity (ppt)
         phi  , & ! ice layer liquid fraction
         qin0 , &
         qin


    integer :: &
         k        ! vertical layer index

    write(nu_diag,*) istep1, "convergence failed!"

    write(nu_diag,*) 0, Tsf0, Tsf
    
    do k = 1, nslyr
       write(nu_diag,*) k, Tsn0(k), Tsn(k), qsn0(k)
    enddo ! k          
    
    do k = 1, nilyr
       write(nu_diag,*) k, Tin0(k), Tin(k), Sin0(k), Sin(k), phi(k), qin0(k)
    enddo ! k

  end subroutine picard_nonconvergence

!=======================================================================
  
  subroutine check_picard_convergence(lsnow,                &
                                      lconverged, nit,      & 
                                      Tsf,        Tsf_prev, &
                                      Tin,        Tin_prev, &
                                      Tsn,        Tsn_prev, &
                                      phi,        Tbot,     &
                                      qin,        qsn,      &
                                      km,         ks,       &
                                      hilyr,      hslyr,    &
                                      fswint,               &
                                      einit,                &
                                      fcondtop,   fcondbot, &
                                      fadvheat)

    logical, intent(inout) :: &
         lconverged   ! has Picard solver converged?

    logical, intent(in) :: &
         lsnow        ! snow presence: T: has snow, F: no snow

    integer, intent(in) :: &
         nit          ! Picard iteration count

    real(kind=dbl_kind), intent(in) :: &
         Tsf      , & ! snow surface temperature (C)
         Tsf_prev , & ! snow surface temperature at previous iteration
         hilyr    , & ! ice layer thickness (m)
         hslyr    , & ! snow layer thickness (m)
         fswint   , & ! SW absorbed in ice interior below surface (W m-2)
         einit    , & ! initial total energy (J)
         Tbot     , & ! ice bottom surfce temperature (deg C)
         fadvheat     ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Tin      , & ! ice layer temperature (C)
         Tin_prev , & ! ice layer temperature at previous iteration
         phi      , & ! ice layer liquid fraction
         km           ! ice conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(nilyr), intent(out) :: &
         qin          ! ice layer enthalpy (J m-3)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         Tsn      , & ! snow layer temperature (C)
         Tsn_prev , & ! snow layer temperature at previous iteration
         ks           ! snow conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(nslyr), intent(out) :: &
         qsn          ! snow layer enthalpy (J m-3)

    real(kind=dbl_kind), intent(out) :: &    
         fcondtop , & ! downward cond flux at top surface (W m-2)
         fcondbot     ! downward cond flux at bottom surface (W m-2)

    real(kind=dbl_kind) :: &
         ferr     , & ! energy flux error
         efinal   , & ! initial total energy (J) at iteration
         dTsn     , & ! change in snow temperature (C) between iterations
         dTin     , & ! change in ice temperature (C) between iterations
         dTsf         ! change in surface temperature (C) between iterations
         
    real(kind=dbl_kind), parameter :: &
         dTmp_errmax = 5.0e-4_dbl_kind ! max allowed chnage in temperature changes between iterations

    call picard_final(lsnow,    &
                      qin, qsn, &
                      Tin, Tsn, &
                      phi)

    call total_energy_content(lsnow,         &
                              qin,    qsn,   &
                              hilyr,  hslyr, &
                              efinal)

    call maximum_variables_changes(lsnow,               &
                                   Tsf, Tsf_prev, dTsf, &
                                   Tsn, Tsn_prev, dTsn, & 
                                   Tin, Tin_prev, dTin)


    fcondbot = c2 * km(nilyr) * ((Tin(nilyr) - Tbot) / hilyr)

    if (lsnow) then
       fcondtop = c2 * ks(1) * ((Tsf - Tsn(1)) / hslyr)
    else
       fcondtop = c2 * km(1) * ((Tsf - Tin(1)) / hilyr)
    endif

    ferr = (efinal - einit) / dt - (fcondtop - fcondbot + fswint - fadvheat)

    lconverged = (dTsf < dTmp_errmax .and. &
                  dTsn < dTmp_errmax .and. &
                  dTin < dTmp_errmax .and. &
                  abs(ferr) < 0.9_dbl_kind*ferrmax)

    !write(53,*) nit, dTsf < dTmp_errmax, dTsn < dTmp_errmax, dTin < dTmp_errmax, abs(ferr) < 0.9_dbl_kind*ferrmax

    !write(54,*)  ferr, (efinal - einit) / dt, (fcondtop - fcondbot + fswint - fadvheat), &
    !     fcondtop, fcondbot, fswint, fadvheat

  end subroutine check_picard_convergence

!=======================================================================

  subroutine picard_drainage_fluxes(fadvheat, q,    &
                                    qbr,      qocn, &
                                    hilyr)

    real(kind=dbl_kind), intent(out) :: &
         fadvheat ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q        ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qbr      ! ice layer brine enthalpy (J m-3)

    real(kind=dbl_kind), intent(in) :: &
         qocn , & ! ocean brine enthalpy (J m-3)
         hilyr    ! ice layer thickness (m)

    integer :: &
         k        ! vertical layer index

    fadvheat = c0

    ! calculate fluxes from base upwards
    do k = 1, nilyr-1

       fadvheat = fadvheat - q(k) * (qbr(k+1) - qbr(k))

    enddo ! k

    k = nilyr

    fadvheat = fadvheat - q(k) * (qocn - qbr(k))

  end subroutine picard_drainage_fluxes

!=======================================================================

  subroutine picard_flushing_fluxes(fadvheat, w,   &
                                    qbr,           &
                                    qocn,     qpond)

   real(kind=dbl_kind), intent(inout) :: &
         fadvheat  ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qbr       ! ice layer brine enthalpy (J m-3)

    real(kind=dbl_kind), intent(in) :: &
         w     , & ! vertical flushing Darcy velocity (m/s)
         qocn  , & ! ocean brine enthalpy (J m-3)
         qpond     ! melt pond brine enthalpy (J m-3)
    
    real(kind=dbl_kind) :: &
         wdown , & ! downwards vertical flushing Darcy velocity (m/s)
         wup       ! upwards vertical flushing Darcy velocity (m/s)

    wdown = -min(w,c0)
    wup   =  max(w,c0)

    fadvheat = fadvheat - wup   * (qocn       - qbr(1)) &
                        + wdown * (qbr(nilyr) - qpond)

  end subroutine picard_flushing_fluxes

!=======================================================================

  subroutine maximum_variables_changes(lsnow,               &
                                       Tsf, Tsf_prev, dTsf, &
                                       Tsn, Tsn_prev, dTsn, & 
                                       Tin, Tin_prev, dTin)

    logical, intent(in) :: &
         lsnow        ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), intent(in) :: &
         Tsf      , & ! snow surface temperature (C)
         Tsf_prev     ! snow surface temperature at previous iteration

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         Tsn      , & ! snow layer temperature (C)
         Tsn_prev     ! snow layer temperature at previous iteration

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Tin      , & ! ice layer temperature (C)
         Tin_prev     ! ice layer temperature at previous iteration

    real(kind=dbl_kind), intent(out) :: &
         dTsf     , & ! change in surface temperature (C) between iterations
         dTsn     , & ! change in snow temperature (C) between iterations
         dTin         ! change in surface temperature (C) between iterations

    integer :: &
         k            ! vertical layer index
    
    dTsf = abs(Tsf - Tsf_prev)

    if (lsnow) then
       dTsn = maxval(abs(Tsn - Tsn_prev))
    else ! lsnow
       dTsn = c0 
    endif ! lsnow

    dTin = maxval(abs(Tin - Tin_prev))

  end subroutine maximum_variables_changes

!=======================================================================

  subroutine total_energy_content(lsnow,         &
                                  qin,    qsn,   &
                                  hilyr,  hslyr, &
                                  energy)
    logical, intent(in) :: &
         lsnow     ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qin       ! ice layer enthalpy (J m-3)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         qsn       ! snow layer enthalpy (J m-3)

    real(kind=dbl_kind), intent(in) :: &
         hilyr , & ! ice layer thickness (m)
         hslyr     ! snow layer thickness (m)

    real(kind=dbl_kind), intent(out) :: &    
         energy    ! total energy of ice and snow

    integer :: &
         k         ! vertical layer index

    energy = c0
    
    if (lsnow) then
       
       do k = 1, nslyr

          energy = energy + hslyr * qsn(k)

       enddo ! k

    endif ! lsnow

    do k = 1, nilyr

       energy = energy + hilyr * qin(k)

    enddo ! k

  end subroutine total_energy_content

!=======================================================================

  subroutine picard_updates(Tin, &
                            Sbr, qbr)

    ! update brine salinity and liquid fraction based on new temperatures

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Tin     ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(nilyr), intent(inout) :: &
         Sbr , & ! ice layer brine salinity (ppt)
         qbr     ! ice layer brine enthalpy (J m-3)

    integer :: &
         k       ! vertical layer index

    do k = 1, nilyr

       Sbr(k) = liquidus_brine_salinity_mush(Tin(k))
       qbr(k) = enthalpy_brine(Tin(k))

    enddo ! k

  end subroutine picard_updates

!=======================================================================

  subroutine picard_updates_enthalpy(Tin, qbr)

    ! update brine salinity and liquid fraction based on new temperatures

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Tin ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(nilyr), intent(inout) :: &
         qbr ! ice layer brine enthalpy (J m-3)

    integer :: &
         k   ! vertical layer index

    do k = 1, nilyr

       qbr(k) = enthalpy_brine(Tin(k))

    enddo ! k

  end subroutine picard_updates_enthalpy

!=======================================================================

  subroutine picard_final(lsnow,    &
                          qin, qsn, &
                          Tin, Tsn, &
                          phi)

    logical, intent(in) :: &
         lsnow   ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), dimension(nilyr), intent(out) :: &
         qin     ! ice layer enthalpy (J m-3)

    real(kind=dbl_kind), dimension(nslyr), intent(out) :: &
         qsn     ! snow layer enthalpy (J m-3)

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Tin , & ! ice layer temperature (C)
         phi     ! ice layer liquid fraction

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         Tsn     ! snow layer temperature (C)

    integer :: &
         k       ! vertical layer index

    do k = 1, nilyr

       qin(k) = enthalpy_mush_liquid_fraction(Tin(k), phi(k))

    enddo ! k

    if (lsnow) then

       do k = 1, nslyr
          
          qsn(k) = enthalpy_snow(Tsn(k))
          
       enddo ! k

    endif ! lsnow

  end subroutine picard_final

!=======================================================================

  subroutine calc_intercell_thickness(lsnow, hilyr, hslyr, dxp)

    logical, intent(in) :: &
         lsnow     ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), intent(in) :: &
         hilyr , & ! ice layer thickness (m)
         hslyr     ! snow layer thickness (m)

    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(out) :: &
         dxp       ! distances between grid points (m)

    integer :: &
         l         ! vertical index

    if (lsnow) then

       dxp(1) = hslyr / c2
       
       do l = 2, nslyr

          dxp(l) = hslyr

       enddo ! l

       dxp(nslyr+1) = (hilyr + hslyr) / c2

       do l = nslyr+2, nilyr+nslyr

          dxp(l) = hilyr

       enddo ! l

       dxp(nilyr+nslyr+1) = hilyr / c2

    else ! lsnow

       dxp(1) = hilyr / c2

       do l = 2, nilyr

          dxp(l) = hilyr

       enddo ! l

       dxp(nilyr+1) = hilyr / c2

       do l = nilyr+2, nilyr+nslyr+1

          dxp(l) = c0

       enddo ! l

    endif ! lsnow

  end subroutine calc_intercell_thickness

!=======================================================================

  subroutine calc_intercell_conductivity(lsnow,        &
                                         km,    ks,    &
                                         hilyr, hslyr, &
                                         kcstar)

    logical, intent(in) :: &
         lsnow      ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         km         ! ice conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         ks         ! snow conductivity (W m-1 K-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr  , & ! ice layer thickness (m)
         hslyr      ! snow layer thickness (m)

    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(out) :: &
         kcstar     ! interface conductivities (W m-1 K-1)

    real(kind=dbl_kind) :: &
         fe         ! distance fraction at interface

    integer :: &
         k, &       ! vertical layer index
         l          ! vertical index

    if (lsnow) then

       kcstar(1) = ks(1)
       
       do l = 2, nslyr

          k = l
          kcstar(l) = (c2 * ks(k) * ks(k-1)) / (ks(k) + ks(k-1))

       enddo ! l

       fe = hilyr / (hilyr + hslyr)
       kcstar(nslyr+1) = c1 / ((c1 - fe) / ks(nslyr) + fe / km(1))

       do k = 2, nilyr

          l = k + nslyr
          kcstar(l) = (c2 * km(k) * km(k-1)) / (km(k) + km(k-1))

       enddo ! k

       kcstar(nilyr+nslyr+1) = km(nilyr)

    else ! lsnow

       kcstar(1) = km(1)

       do k = 2, nilyr
          
          l = k
          kcstar(l) = (c2 * km(k) * km(k-1)) / (km(k) + km(k-1))

       enddo ! k

       kcstar(nilyr+1) = km(nilyr)

       do l = nilyr+2, nilyr+nslyr+1

          kcstar(l) = c0

       enddo ! l

    endif ! lsnow

  end subroutine calc_intercell_conductivity

!=======================================================================
  
  subroutine solve_heat_conduction(lsnow,  lcold,        &
                                   Tsf,    Tbot,         &
                                   qin0,   qsn0,         &
                                   phi,    dt,           &
                                   qpond,  qocn,         &
                                   q,      w,            &
                                   hilyr,  hslyr,        &
                                   dxp,    kcstar,       &
                                   Iswabs, Sswabs,       &
                                   fsurfn, dfsurfn_dTsf, &
                                   Tin,    Tsn,nit)

    logical, intent(in) :: &
         lsnow        , & ! snow presence: T: has snow, F: no snow
         lcold            ! surface cold: T: surface is cold, F: surface is melting

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qin0         , & ! ice layer enthalpy (J m-3) at beggining of timestep
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         phi              ! ice layer liquid fraction

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         qsn0         , & ! snow layer enthalpy (J m-3) at start of timestep
         Sswabs           ! SW radiation absorbed in snow layers (W m-2)
    
    real(kind=dbl_kind), intent(inout) :: &
         Tsf              ! snow surface temperature (C)

    real(kind=dbl_kind), intent(in) :: &
         dt           , & ! timestep (s)
         hilyr        , & ! ice layer thickness (m)
         hslyr        , & ! snow layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         w            , & ! vertical flushing Darcy velocity (m/s)
         fsurfn       , & ! net flux to top surface, excluding fcondtop
         dfsurfn_dTsf     ! derivative of net flux to top surface, excluding fcondtopn

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(in) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    real(kind=dbl_kind), dimension(nilyr), intent(inout) :: &
         Tin              ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(nslyr), intent(out) :: &
         Tsn              ! snow layer temperature (C)

    integer, intent(in) :: &
         nit              ! Picard iteration count

    real(kind=dbl_kind), dimension(nilyr+nslyr+1) :: &
         Ap           , & ! diagonal of tridiagonal matrix
         As           , & ! lower off-diagonal of tridiagonal matrix
         An           , & ! upper off-diagonal of tridiagonal matrix
         b            , & ! right hand side of matrix solve
         T                ! ice and snow temperatures

    real(kind=dbl_kind) :: &
         wup          , & ! upwards vertical flushing Darcy velocity (m/s)
         wdown            ! downwards vertical flushing Darcy velocity (m/s)

    integer :: &
         k            , & ! vertical layer index
         l            , & ! vertical index
         nyn              ! matrix size

    wdown = -min(w,c0)
    wup   =  max(w,c0)

    ! set up matrix and right hand side - snow
    if (lsnow) then

       if (lcold) then

          call matrix_elements_snow_cold(Ap, As, An, b, nyn, &
                                         Tsf,    Tbot,         &
                                         qin0,   qsn0,         &
                                         Tin,                  &
                                         qpond,  qocn,         &
                                         phi,    q,            &
                                         wup,    wdown,        &
                                         hilyr,  hslyr,        &
                                         dxp,    kcstar,       &
                                         Iswabs, Sswabs,       &
                                         fsurfn, dfsurfn_dTsf, &
                                         dt)

       else ! lcold

          call matrix_elements_snow_melt(Ap, As, An, b, nyn, &
                                         Tsf,    Tbot,         &
                                         qin0,   qsn0,         &
                                         Tin,                  &
                                         qpond,  qocn,         &
                                         phi,    q,            &
                                         wup,    wdown,        &
                                         hilyr,  hslyr,        &
                                         dxp,    kcstar,       &
                                         Iswabs, Sswabs,       &
                                         fsurfn, dfsurfn_dTsf, &
                                         dt)

       endif ! lcold

    else ! lsnow

       if (lcold) then

          call matrix_elements_nosnow_cold(Ap, As, An, b, nyn, &
                                           Tsf,    Tbot,         &
                                           qin0,   qsn0,         &
                                           Tin,                  &
                                           qpond,  qocn,         &
                                           phi,    q,            &
                                           wup,    wdown,        &
                                           hilyr,  hslyr,        &
                                           dxp,    kcstar,       &
                                           Iswabs, Sswabs,       &
                                           fsurfn, dfsurfn_dTsf, &
                                           dt)

       else ! lcold

          call matrix_elements_nosnow_melt(Ap, As, An, b, nyn, &
                                           Tsf,    Tbot,         &
                                           qin0,   qsn0,         &
                                           Tin,                  &
                                           qpond,  qocn,         &
                                           phi,    q,            &
                                           wup,    wdown,        &
                                           hilyr,  hslyr,        &
                                           dxp,    kcstar,       &
                                           Iswabs, Sswabs,       &
                                           fsurfn, dfsurfn_dTsf, &
                                           dt)

       endif ! lcold

    endif ! lsnow

    ! tridiag to get new temperatures
    call tdma_solve_sparse(An(1:nyn), Ap(1:nyn), As(1:nyn), b(1:nyn), T(1:nyn), nyn)

    call update_temperatures(lsnow, lcold, &
                             T,     Tsf,   &
                             Tin,   Tsn)

  end subroutine solve_heat_conduction

!=======================================================================

  subroutine update_temperatures(lsnow, lcold, &
                                 T,     Tsf,   &
                                 Tin,   Tsn)

    logical, intent(in) :: &
         lsnow , & ! snow presence: T: has snow, F: no snow
         lcold     ! surface cold: T: surface is cold, F: surface is melting

    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(in) :: &
         T         ! matrix solution vector

    real(kind=dbl_kind), intent(inout) :: &
         Tsf       ! snow surface temperature (C)

    real(kind=dbl_kind), dimension(nilyr), intent(inout) :: &
         Tin       ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(nslyr), intent(inout) :: &
         Tsn       ! snow layer temperature (C)

    integer :: &
         l     , & ! vertical index
         k         ! vertical layer index

    if (lsnow) then
       
       if (lcold) then

          Tsf = T(1)

          do k = 1, nslyr
             l = k + 1
             Tsn(k) = T(l)
          enddo ! k
          
          do k = 1, nilyr
             l = k + nslyr + 1
             Tin(k) = T(l)
          enddo ! k

       else ! lcold
          
          do k = 1, nslyr
             l = k
             Tsn(k) = T(l)
          enddo ! k
          
          do k = 1, nilyr
             l = k + nslyr
             Tin(k) = T(l)
          enddo ! k

       endif ! lcold

    else ! lsnow
       
       if (lcold) then

          Tsf = T(1)
          
          do k = 1, nilyr
             l = k + 1
             Tin(k) = T(l)
          enddo ! k

       else ! lcold

          do k = 1, nilyr
             l = k
             Tin(k) = T(l)
          enddo ! k

       endif ! lcold

    endif ! lsnow

  end subroutine update_temperatures

!=======================================================================
  
  subroutine matrix_elements_nosnow_melt(Ap, As, An, b, nyn, &
                                         Tsf,    Tbot,         &
                                         qin0,   qsn0,         &
                                         Tin,                  &
                                         qpond,  qocn,         &
                                         phi,    q,            &
                                         wup,    wdown,        &
                                         hilyr,  hslyr,        &
                                         dxp,    kcstar,       &
                                         Iswabs, Sswabs,       &
                                         fsurfn, dfsurfn_dTsf, &
                                         dt)
     
    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(out) :: &
         Ap           , & ! diagonal of tridiagonal matrix
         As           , & ! lower off-diagonal of tridiagonal matrix
         An           , & ! upper off-diagonal of tridiagonal matrix
         b                ! right hand side of matrix solve

    integer, intent(out) :: &
         nyn              ! matrix size

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qin0         , & ! ice layer enthalpy (J m-3) at beggining of timestep
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         phi          , & ! ice layer liquid fraction
         Tin              ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         qsn0         , & ! snow layer enthalpy (J m-3) at start of timestep
         Sswabs           ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), intent(in) :: &
         Tsf          , & ! snow surface temperature (C)
         dt           , & ! timestep (s)
         hilyr        , & ! ice layer thickness (m)
         hslyr        , & ! snow layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         wup          , & ! upwards vertical flushing Darcy velocity (m/s)
         wdown        , & ! downwards vertical flushing Darcy velocity (m/s)
         fsurfn       , & ! net flux to top surface, excluding fcondtop
         dfsurfn_dTsf     ! derivative of net flux to top surface, excluding fcondtopn

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(in) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    integer :: &
         k            , & ! vertical layer index
         l                ! vertical index
    
    ! surface layer
    k = 1
    l = k
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(k+1) / dxp(k+1) + &
             kcstar(k)   / dxp(k) + &
             q(k)  * cp_ocn * rhow + &
             wdown * cp_ocn * rhow + &
             wup   * cp_ocn * rhow
    As(l) = -kcstar(k+1) / dxp(k+1) - &
             q(k)  * cp_ocn * rhow - &
             wup   * cp_ocn * rhow
    An(l) = c0
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k) + &
            (kcstar(k) / dxp(k)) * Tsf + &
             wdown * qpond

    ! interior ice layers
    do k = 2, nilyr-1
          
       l = k
          
       Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
                kcstar(k+1) / dxp(k+1) + &
                kcstar(k)   / dxp(k) + &
                q(k)  * cp_ocn * rhow + &
                wdown * cp_ocn * rhow + &
                wup   * cp_ocn * rhow
       As(l) = -kcstar(k+1) / dxp(k+1) - &
                q(k)  * cp_ocn * rhow - &
                wup   * cp_ocn * rhow
       An(l) = -kcstar(k)   / dxp(k) - &
                wdown * cp_ocn * rhow
       b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k)

    enddo ! k
    
    ! bottom layer
    k = nilyr
    l = k
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(k+1) / dxp(k+1) + &
             kcstar(k)   / dxp(k) + &
             q(k)  * cp_ocn * rhow + &
             wdown * cp_ocn * rhow + &
             wup   * cp_ocn * rhow
    As(l) = c0
    An(l) = -kcstar(k) / dxp(k) - &
            wdown * cp_ocn * rhow
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k) + &
            (kcstar(k+1) * Tbot) / dxp(k+1) + &
            q(k) * qocn + &
            wup  * qocn
    
    nyn = nilyr

  end subroutine matrix_elements_nosnow_melt

!=======================================================================

  subroutine matrix_elements_nosnow_cold(Ap, As, An, b, nyn, &
                                         Tsf,    Tbot,         &
                                         qin0,   qsn0,         &
                                         Tin,                  &
                                         qpond,  qocn,         &
                                         phi,    q,            &
                                         wup,    wdown,        &
                                         hilyr,  hslyr,        &
                                         dxp,    kcstar,       &
                                         Iswabs, Sswabs,       &
                                         fsurfn, dfsurfn_dTsf, &
                                         dt)
     
    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(out) :: &
         Ap           , & ! diagonal of tridiagonal matrix
         As           , & ! lower off-diagonal of tridiagonal matrix
         An           , & ! upper off-diagonal of tridiagonal matrix
         b                ! right hand side of matrix solve

    integer, intent(out) :: &
         nyn              ! matrix size

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qin0         , & ! ice layer enthalpy (J m-3) at beggining of timestep
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         phi          , & ! ice layer liquid fraction
         Tin              ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         qsn0         , & ! snow layer enthalpy (J m-3) at start of timestep
         Sswabs           ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), intent(in) :: &
         Tsf          , & ! snow surface temperature (C)
         dt           , & ! timestep (s)
         hilyr        , & ! ice layer thickness (m)
         hslyr        , & ! snow layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         wup          , & ! upwards vertical flushing Darcy velocity (m/s)
         wdown        , & ! downwards vertical flushing Darcy velocity (m/s)
         fsurfn       , & ! net flux to top surface, excluding fcondtop
         dfsurfn_dTsf     ! derivative of net flux to top surface, excluding fcondtopn

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(in) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    integer :: &
         k            , & ! vertical layer index
         l                ! vertical index

    ! surface temperature
    l = 1
    Ap(l) = dfsurfn_dTsf - kcstar(1) / dxp(1)
    As(l) = kcstar(1) / dxp(1)
    An(l) = c0
    b (l) = dfsurfn_dTsf * Tsf - fsurfn

    ! surface layer
    k = 1
    l = k + 1
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(k+1) / dxp(k+1) + &
             kcstar(k)   / dxp(k) + &
             q(k)  * cp_ocn * rhow + &
             wdown * cp_ocn * rhow + &
             wup   * cp_ocn * rhow
    As(l) = -kcstar(k+1) / dxp(k+1) - &
             q(k)  * cp_ocn * rhow - &
             wup   * cp_ocn * rhow
    An(l) = -kcstar(k)   / dxp(k)
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k) + &
             wdown * qpond

    ! interior ice layers
    do k = 2, nilyr-1
          
       l = k + 1
          
       Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
                kcstar(k+1) / dxp(k+1) + &
                kcstar(k)   / dxp(k) + &
                q(k)  * cp_ocn * rhow + &
                wdown * cp_ocn * rhow + &
                wup   * cp_ocn * rhow
       As(l) = -kcstar(k+1) / dxp(k+1) - &
                q(k)  * cp_ocn * rhow - &
                wup   * cp_ocn * rhow
       An(l) = -kcstar(k)   / dxp(k) - &
                wdown * cp_ocn * rhow
       b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k)

    enddo ! k
    
    ! bottom layer
    k = nilyr
    l = k + 1
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(k+1) / dxp(k+1) + &
             kcstar(k)   / dxp(k) + &
             q(k)  * cp_ocn * rhow + &
             wdown * cp_ocn * rhow + &
             wup   * cp_ocn * rhow
    As(l) = c0
    An(l) = -kcstar(k) / dxp(k) - &
            wdown * cp_ocn * rhow
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k) + &
            (kcstar(k+1) * Tbot) / dxp(k+1) + &
            q(k) * qocn + &
            wup  * qocn
    
    nyn = nilyr + 1

  end subroutine matrix_elements_nosnow_cold

!=======================================================================

  subroutine matrix_elements_snow_melt(Ap, As, An, b, nyn, &
                                       Tsf,    Tbot,         &
                                       qin0,   qsn0,         &
                                       Tin,                  &
                                       qpond,  qocn,         &
                                       phi,    q,            &
                                       wup,    wdown,        &
                                       hilyr,  hslyr,        &
                                       dxp,    kcstar,       &
                                       Iswabs, Sswabs,       &
                                       fsurfn, dfsurfn_dTsf, &
                                       dt)
     
    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(out) :: &
         Ap           , & ! diagonal of tridiagonal matrix
         As           , & ! lower off-diagonal of tridiagonal matrix
         An           , & ! upper off-diagonal of tridiagonal matrix
         b                ! right hand side of matrix solve

    integer, intent(out) :: &
         nyn              ! matrix size

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qin0         , & ! ice layer enthalpy (J m-3) at beggining of timestep
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         phi          , & ! ice layer liquid fraction
         Tin              ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         qsn0         , & ! snow layer enthalpy (J m-3) at start of timestep
         Sswabs           ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), intent(in) :: &
         Tsf          , & ! snow surface temperature (C)
         dt           , & ! timestep (s)
         hilyr        , & ! ice layer thickness (m)
         hslyr        , & ! snow layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         wup          , & ! upwards vertical flushing Darcy velocity (m/s)
         wdown        , & ! downwards vertical flushing Darcy velocity (m/s)
         fsurfn       , & ! net flux to top surface, excluding fcondtop
         dfsurfn_dTsf     ! derivative of net flux to top surface, excluding fcondtopn

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(in) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    integer :: &
         k            , & ! vertical layer index
         l                ! vertical index

    ! surface layer
    k = 1
    l = k
    
    Ap(l) = ((rhos * cp_ice) / dt) * hslyr + &
             kcstar(l+1) / dxp(l+1) + &
             kcstar(l)   / dxp(l)
    As(l) = -kcstar(l+1) / dxp(l+1)
    An(l) = c0
    b (l) = ((rhos * Lfresh + qsn0(k)) / dt) * hslyr + Sswabs(k) + &
            (kcstar(l) * Tsf) / dxp(l)

    ! interior snow layers
    do k = 2, nslyr
       
       l = k

       Ap(l) = ((rhos * cp_ice) / dt) * hslyr + &
                kcstar(l+1) / dxp(l+1) + &
                kcstar(l)   / dxp(l)
       As(l) = -kcstar(l+1) / dxp(l+1)
       An(l) = -kcstar(l)   / dxp(l)
       b (l) = ((rhos * Lfresh + qsn0(k)) / dt) * hslyr + Sswabs(k)
          
    enddo ! k
    
    ! top ice layer
    k = 1
    l = nslyr + k
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(l+1) / dxp(l+1) + &
             kcstar(l)   / dxp(l) + &
             q(k)  * cp_ocn * rhow + &
             wdown * cp_ocn * rhow + &
             wup   * cp_ocn * rhow
    As(l) = -kcstar(l+1) / dxp(l+1) - &
             q(k)  * cp_ocn * rhow - &
             wup   * cp_ocn * rhow
    An(l) = -kcstar(l)   / dxp(l)
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k) + &
             wdown * qpond

    ! interior ice layers
    do k = 2, nilyr-1
          
       l = nslyr + k
       
       Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
                kcstar(l+1) / dxp(l+1) + &
                kcstar(l)   / dxp(l) + &
                q(k)  * cp_ocn * rhow + &
                wdown * cp_ocn * rhow + &
                wup   * cp_ocn * rhow
       As(l) = -kcstar(l+1) / dxp(l+1) - &
                q(k)  * cp_ocn * rhow - &
                wup   * cp_ocn * rhow
       An(l) = -kcstar(l)   / dxp(l) - &
                wdown * cp_ocn * rhow
       b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k)
          
    enddo ! k

    ! bottom layer
    k = nilyr
    l = nilyr + nslyr
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(l+1) / dxp(l+1) + &
             kcstar(l)   / dxp(l) + &
             q(k)  * cp_ocn * rhow + &
             wdown * cp_ocn * rhow + &
             wup   * cp_ocn * rhow
    As(l) = c0
    An(l) = -kcstar(l)   / dxp(l) - &
             wdown * cp_ocn * rhow
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k) + &
            (kcstar(l+1) * Tbot) / dxp(l+1) + &
             q(k) * qocn + &
             wup  * qocn
       
    nyn = nilyr + nslyr

  end subroutine matrix_elements_snow_melt

!=======================================================================

  subroutine matrix_elements_snow_cold(Ap, As, An, b, nyn, &
                                       Tsf,    Tbot,         &
                                       qin0,   qsn0,         &
                                       Tin,                  &
                                       qpond,  qocn,         &
                                       phi,    q,            &
                                       wup,    wdown,        &
                                       hilyr,  hslyr,        &
                                       dxp,    kcstar,       &
                                       Iswabs, Sswabs,       &
                                       fsurfn, dfsurfn_dTsf, &
                                       dt)
     
    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(out) :: &
         Ap           , & ! diagonal of tridiagonal matrix
         As           , & ! lower off-diagonal of tridiagonal matrix
         An           , & ! upper off-diagonal of tridiagonal matrix
         b                ! right hand side of matrix solve

    integer, intent(out) :: &
         nyn              ! matrix size

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qin0         , & ! ice layer enthalpy (J m-3) at beggining of timestep
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         phi          , & ! ice layer liquid fraction
         Tin              ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         qsn0         , & ! snow layer enthalpy (J m-3) at start of timestep
         Sswabs           ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), intent(in) :: &
         Tsf          , & ! snow surface temperature (C)
         dt           , & ! timestep (s)
         hilyr        , & ! ice layer thickness (m)
         hslyr        , & ! snow layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         wup          , & ! upwards vertical flushing Darcy velocity (m/s)
         wdown        , & ! downwards vertical flushing Darcy velocity (m/s)
         fsurfn       , & ! net flux to top surface, excluding fcondtop
         dfsurfn_dTsf     ! derivative of net flux to top surface, excluding fcondtopn

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(nilyr+nslyr+1), intent(in) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    integer :: &
         k            , & ! vertical layer index
         l            , & ! matrix index
         m                ! vertical index

    ! surface temperature
    l = 1
    Ap(l) = dfsurfn_dTsf - kcstar(1) / dxp(1)
    As(l) = kcstar(1) / dxp(1)
    An(l) = c0
    b (l) = dfsurfn_dTsf * Tsf - fsurfn

    ! surface layer
    k = 1
    l = k + 1
    m = 1
    
    Ap(l) = ((rhos * cp_ice) / dt) * hslyr + &
             kcstar(m+1) / dxp(m+1) + &
             kcstar(m)   / dxp(m)
    As(l) = -kcstar(m+1) / dxp(m+1)
    An(l) = -kcstar(m)   / dxp(m)
    b (l) = ((rhos * Lfresh + qsn0(k)) / dt) * hslyr + Sswabs(k)

    ! interior snow layers
    do k = 2, nslyr
       
       l = k + 1
       m = k

       Ap(l) = ((rhos * cp_ice) / dt) * hslyr + &
                kcstar(m+1) / dxp(m+1) + &
                kcstar(m)   / dxp(m)
       As(l) = -kcstar(m+1) / dxp(m+1)
       An(l) = -kcstar(m)   / dxp(m)
       b (l) = ((rhos * Lfresh + qsn0(k)) / dt) * hslyr + Sswabs(k)

    enddo ! k
    
    ! top ice layer
    k = 1
    l = nslyr + k + 1
    m = k + nslyr
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(m+1) / dxp(m+1) + &
             kcstar(m)   / dxp(m) + &
             q(k)  * cp_ocn * rhow + &
             wdown * cp_ocn * rhow + &
             wup   * cp_ocn * rhow
    As(l) = -kcstar(m+1) / dxp(m+1) - &
             q(k)  * cp_ocn * rhow - &
             wup   * cp_ocn * rhow
    An(l) = -kcstar(m)   / dxp(m)
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k) + &
              wdown * qpond

    ! interior ice layers
    do k = 2, nilyr-1
          
       l = nslyr + k + 1
       m = k + nslyr
       
       Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
                kcstar(m+1) / dxp(m+1) + &
                kcstar(m)   / dxp(m) + &
                q(k)  * cp_ocn * rhow + &
                wdown * cp_ocn * rhow + &
                wup   * cp_ocn * rhow
       As(l) = -kcstar(m+1) / dxp(m+1) - &
                q(k)  * cp_ocn * rhow - &
                wup   * cp_ocn * rhow
       An(l) = -kcstar(m)   / dxp(m) - &
                wdown * cp_ocn * rhow
       b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k)

    enddo ! k

    ! bottom layer
    k = nilyr
    l = nilyr + nslyr + 1
    m = k + nslyr
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(m+1) / dxp(m+1) + &
             kcstar(m)   / dxp(m) + &
             q(k)  * cp_ocn * rhow + &
             wdown * cp_ocn * rhow + &
             wup   * cp_ocn * rhow
    As(l) = c0
    An(l) = -kcstar(m)   / dxp(m) - &
             wdown * cp_ocn * rhow
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + qin0(k)) / dt) * hilyr + Iswabs(k) + &
            (kcstar(m+1) * Tbot) / dxp(m+1) + &
             q(k) * qocn + &
             wup  * qocn

    nyn = nilyr + nslyr + 1

  end subroutine matrix_elements_snow_cold

!=======================================================================

  subroutine solve_salinity(Sin,   Sbr,   &
                            Spond, sss,   &
                            q,     dSdt,  &
                            w,     hilyr, &
                            dt)

    real(kind=dbl_kind), dimension(nilyr), intent(inout) :: &
         Sin       ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Sbr   , & ! ice layer brine salinity (ppt)
         dSdt      ! gravity drainage desalination rate for slow mode (ppt s-1)
         
    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q         ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         Spond , & ! melt pond salinity (ppt)
         sss   , & ! sea surface salinity (ppt)
         w     , & ! vertical flushing Darcy velocity (m/s)
         hilyr , & ! ice layer thickness (m)
         dt        ! timestep (s)

    real(kind=dbl_kind) :: &
         wdown , & ! downwards vertical flushing Darcy velocity (m/s)
         wup       ! upwards vertical flushing Darcy velocity (m/s)

    integer :: &
         k         ! vertical layer index

    real(kind=dbl_kind), parameter :: &
         S_min = 0.01_dbl_kind
    


    real(kind=dbl_kind), dimension(nilyr) :: &
         Sin0

    Sin0 = Sin


    wdown = -min(w,c0)
    wup   =  max(w,c0)

    k = 1
    Sin(k) = Sin(k) + max(S_min - Sin(k), &
         ((q(k)  * (Sbr(k+1) - Sbr(k))) / hilyr + &
          dSdt(k)                               + &
          (wdown * (Spond    - Sbr(k))) / hilyr + &
          (wup   * (Sbr(k+1) - Sbr(k))) / hilyr) * dt)

    do k = 2, nilyr-1

       Sin(k) = Sin(k) + max(S_min - Sin(k), &
            ((q(k)  * (Sbr(k+1) - Sbr(k))) / hilyr + &
             dSdt(k)                               + &
             (wdown * (Sbr(k-1) - Sbr(k))) / hilyr + &
             (wup   * (Sbr(k+1) - Sbr(k))) / hilyr) * dt)

    enddo ! k

    k = nilyr
    Sin(k) = Sin(k) + max(S_min - Sin(k), &
         ((q(k)  * (sss      - Sbr(k))) / hilyr + &
          dSdt(k)                               + &
          (wdown * (Sbr(k-1) - Sbr(k))) / hilyr + &
          (wup   * (sss      - Sbr(k))) / hilyr) * dt)


    if (minval(Sin) < c0) then


         write(*,*) (q(k)  * (Sbr(k+1) - Sbr(k))) / hilyr, &
          dSdt(k)                               , &
          (wdown * (Spond    - Sbr(k))) / hilyr , &
          (wup   * (Sbr(k+1) - Sbr(k))) / hilyr

       do k = 1, nilyr
       
          write(*,*) k, Sin(k), Sin0(k)

       enddo

       stop

    endif

  end subroutine solve_salinity

!=======================================================================

  subroutine solve_salinity_diffusion(Sin,   Sbr, &
                                      Spond, sss, &
                                      kappa, w,   &
                                      hilyr, dt)

    real(kind=dbl_kind), dimension(nilyr), intent(inout) :: &
         Sin       ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Sbr       ! ice layer brine salinity (ppt)

    real(kind=dbl_kind), dimension(nilyr+1), intent(in) :: &
         kappa     ! salt diffusivity

    real(kind=dbl_kind), intent(in) :: &
         Spond , & ! melt pond salinity (ppt)
         sss   , & ! sea surface salinity (ppt)
         w     , & ! vertical flushing Darcy velocity (m/s)
         hilyr , & ! ice layer thickness (m)
         dt        ! timestep (s)

    real(kind=dbl_kind) :: &
         wdown , & ! downwards vertical flushing Darcy velocity (m/s)
         wup       ! upwards vertical flushing Darcy velocity (m/s)

    integer :: &
         k         ! vertical layer index

    wdown = -min(w,c0)
    wup   =  max(w,c0)

    k = 1
    Sin(k) = Sin(k) + ( &
         (kappa(k+1) * (Sbr(k+1) - Sbr(k)) - kappa(k) * (Sbr(k) - Spond)) / hilyr**2 + &
         (wdown * (Spond    - Sbr(k))) / hilyr + &
         (wup   * (Sbr(k+1) - Sbr(k))) / hilyr) * dt

    do k = 2, nilyr-1
       
       Sin(k) = Sin(k) + ( &
            (kappa(k+1) * (Sbr(k+1) - Sbr(k)) - kappa(k) * (Sbr(k) - Sbr(k-1))) / hilyr**2 + &
            (wdown * (Sbr(k-1) - Sbr(k))) / hilyr + &
            (wup   * (Sbr(k+1) - Sbr(k))) / hilyr) * dt
       
    enddo ! k
    
    k = nilyr
    Sin(k) = Sin(k) + ( &
         (kappa(k+1) * (sss - Sbr(k)) - kappa(k) * (Sbr(k) - Sbr(k-1))) / hilyr**2 + &
         (wdown * (Sbr(k-1) - Sbr(k))) / hilyr + &
         (wup   * (sss      - Sbr(k))) / hilyr) * dt
    
  end subroutine solve_salinity_diffusion

!=======================================================================

  subroutine tdma_solve_sparse(a, b, c, d, x, n)
    
    ! perform a tri-diagonal solve with TDMA using a sparse tridiagoinal matrix

    real(kind=dbl_kind), dimension(1:n), intent(in) :: &
         a  , & ! matrix lower off-diagonal
         b  , & ! matrix diagonal
         c  , & ! matrix upper off-diagonal
         d      ! right hand side vector
 
    real(kind=dbl_kind), dimension(1:n), intent(out) :: &
         x      ! solution vector

    integer(kind=int_kind), intent(in) :: &
         n      ! matrix size
    
    real(kind=dbl_kind), dimension(nilyr+nslyr+1) :: &
         cp , & ! modified upper off-diagonal vector
         dp     ! modified right hand side vector

    integer(kind=int_kind) :: &
         i      ! vector index
    
    ! forward sweep
    cp(1) = c(1) / b(1)
    do i = 2, n-1
       cp(i) = c(i) / (b(i) - cp(i-1)*a(i))
    enddo
    
    dp(1) = d(1) / b(1)
    do i = 2, n
       dp(i) = (d(i) - dp(i-1)*a(i)) / (b(i) - cp(i-1)*a(i))
    enddo

    ! back substitution
    x(n) = dp(n)
    do i = n-1,1,-1
       x(i) = dp(i) - cp(i)*x(i+1)
    enddo

  end subroutine tdma_solve_sparse

!=======================================================================
! Effect of salinity
!=======================================================================

  subroutine diagnose_salt_flux(q, dSdt)

    ! write out diagnostics related to gravity drainage

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q    ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         dSdt ! gravity drainage desalination rate for slow mode (ppt s-1)

    integer(kind=int_kind) :: &
         k    ! ice layer index

    open(11,file='history/salt_flux.txt',position='append')

    write(11,*) istep, 0, q(0)
    do k = 1, nilyr
       write(11,*) istep, k, q(k), dSdt(k)
    enddo ! k
    write(11,*) istep, k, c0, c0, c0

    close(11)

  end subroutine diagnose_salt_flux

!=======================================================================

  function permeability(phi) result(perm)

    ! given the liquid fraction calculate the permeability
    ! See Golden et al. 2007

    real(kind=dbl_kind), intent(in) :: &
         phi                  ! liquid fraction

    real(kind=dbl_kind) :: &
         perm                 ! permeability (m2)
    
    real(kind=dbl_kind), parameter :: &
         phic = 0.05_dbl_kind ! critical liquid fraction for impermeability

    !perm = 3.0e-8_dbl_kind * phi**3

    perm = 3.0e-8_dbl_kind * max(phi - phic, c0)**3

    !perm = 3.0e-8_dbl_kind * max(phi - phic, c0)**2

  end function permeability

!=======================================================================

  subroutine explicit_flow_velocities(Sin,   qin,  &
                                      Tin,   Tsf,  &
                                      Tbot,  q,    &
                                      dSdt,  Sbr,  &
                                      qbr,   dt,   &
                                      sss,   qocn, &
                                      hilyr, hin)

    ! calculate the rapid gravity drainage mode Darcy velocity and the
    ! slow mode drainage rate

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         Sin , &   ! ice layer bulk salinity (ppt)
         qin , &   ! ice layer enthalpy (J m-3)
         Tin       ! ice layer temperature (C)

    real(kind=dbl_kind), intent(in) :: &
         Tsf   , & ! ice/snow surface temperature (C)
         Tbot  , & ! ice bottom temperature (C)
         dt    , & ! time step (s)
         sss   , & ! sea surface salinty (ppt)
         qocn  , & ! ocean enthalpy (J m-3)
         hilyr , & ! ice layer thickness (m)
         hin       ! ice thickness (m)

    real(kind=dbl_kind), dimension(0:nilyr), intent(out) :: &
         q         ! rapid mode upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(1:nilyr), intent(out) :: &
         dSdt      ! slow mode drainage rate (ppt s-1)

    real(kind=dbl_kind), dimension(1:nilyr+1), intent(out) :: &
         Sbr , &   ! ice layer brine salinity (ppt)
         qbr       ! ice layer brine enthalpy (J m-3)

    real(kind=dbl_kind), parameter :: &
         kappal        = 8.824e-8_dbl_kind                 , & ! heat diffusivity of liquid
         viscosity_dyn = 1.79e-3_dbl_kind                  , & ! dynamic viscosity of brine
         ra_constants  = gravit / (viscosity_dyn * kappal) , & ! Rayleigh number constants
         fracmax       = 0.2_dbl_kind                      , & ! limiting advective fraction of layer
         Sin_min       = 0.1_dbl_kind                      , & ! minimum bulk salinity (ppt)
         safety_factor = 10.0_dbl_kind                         ! safety factor for getting negative salinities

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         phi           ! ice layer liquid fraction

    real(kind=dbl_kind), dimension(0:nilyr) :: &
         rho           ! ice layer brine density (kg m-3)

    real(kind=dbl_kind) :: &
         rho_ocn   , & ! ocean density (kg m-3)
         perm_min  , & ! minimum permeability from layer to ocean (m2)
         perm_harm , & ! harmonic mean of permeability from layer to ocean (m2)
         rho_sum   , & ! sum of the brine densities from layer to ocean (kg m-3)
         rho_pipe  , & ! density of the brine in the channel (kg m-3)
         z         , & ! distance to layer from top surface (m)
         perm      , & ! ice layer permeability (m2)
         drho      , & ! brine density difference between layer and ocean (kg m-3)
         Ra        , & ! local mush Rayleigh number
         rn        , & ! real value of number of layers considered
         L         , & ! thickness of the layers considered (m)
         dx        , & ! horizontal size of convective flow (m)
         dx2       , & ! square of the horizontal size of convective flow (m2)
         Am        , & ! A parameter for mush
         Bm        , & ! B parameter for mush
         Ap        , & ! A parameter for channel
         Bp        , & ! B parameter for channel
         qlimit    , & ! limit to vertical Darcy flow for numerical stability
         dS_guess  , & ! expected bulk salinity without limits
         alpha         ! desalination limiting factor

    integer(kind=int_kind) :: &
         k             ! ice layer index

    ! initial downward sweep - determine derived physical quantities
    do k = 1, nilyr
       
       Sbr(k) = liquidus_brine_salinity_mush(Tin(k))
       phi(k) = liquid_fraction(Tin(k), Sin(k))
       qbr(k) = enthalpy_brine(Tin(k))
       rho(k) = density_brine(Sbr(k))

    enddo ! k

    rho(0) = rho(1)

    ! ocean conditions
    Sbr(nilyr+1) = sss
    qbr(nilyr+1) = qocn
    rho_ocn = density_brine(sss)

    ! initialize accumulated quantities
    perm_min = bignum
    perm_harm = c0
    rho_sum = c0

    ! limit to q for numerical stability
    qlimit = (fracmax * hilyr) / dt

    ! no flow through ice top surface
    q(0) = c0

    ! first iterate over layers going up
    do k = nilyr, 1, -1

       ! vertical position from ice top surface
       z = ((real(k, dbl_kind) - p5) / real(nilyr, dbl_kind)) * hin

       ! permeabilities
       perm = permeability(phi(k))
       perm_min = min(perm_min,perm)
       perm_harm = perm_harm + (c1 / max(perm,1.0e-30_dbl_kind))

       ! densities
       rho_sum = rho_sum + rho(k)
       !rho_pipe = rho(k)
       rho_pipe = p5 * (rho(k) + rho(k-1))
       drho = max(rho(k) - rho_ocn, c0)

       ! mush Rayleigh number 
       Ra = drho * (hin-z) * perm_min * ra_constants

       ! height of mush layer to layer k
       rn = real(nilyr-k+1,dbl_kind)
       L = rn * hilyr

       ! horizontal size of convection
       dx = L * c2 * aspect_rapid_mode
       dx2 = dx**2

       ! determine vertical Darcy flow
       Am = (dx2 * rn) / (viscosity_dyn * perm_harm)
       Bm = (-gravit * rho_sum) / rn

       Ap = (pi * a_rapid_mode**4) / (c8 * viscosity_dyn)
       Bp = -rho_pipe * gravit

       q(k) = max((Am / dx2) * ((-Ap*Bp - Am*Bm) / (Am + Ap) + Bm), 1.0e-30_dbl_kind)

       ! modify by Rayleigh number and advection limit
       q(k) = min(q(k) * (max(Ra - Rac_rapid_mode, c0) / (Ra+puny)), qlimit)

       ! late stage drainage
       dSdt(k) = dSdt_slow_mode * (max((Sin(k) - phi_c_slow_mode*Sbr(k)), c0) * max((Tbot - Tsf), c0)) / (hin + 0.001_dbl_kind)

       dSdt(k) = max(dSdt(k), (-Sin(k) * 0.5_dbl_kind) / dt)

       ! restrict flows to prevent too much salt loss
       dS_guess = (((q(k) * (Sbr(k+1) - Sbr(k))) / hilyr + dSdt(k)) * dt) * safety_factor

       if (abs(dS_guess) < puny) then
          alpha = c1
       else
          alpha = (Sin_min - Sin(k)) / dS_guess
       endif

       if (alpha < c0 .or. alpha > c1) alpha = c1

       q(k)    = q(k)    * alpha
       dSdt(k) = dSdt(k) * alpha

    enddo ! k

  end subroutine explicit_flow_velocities

!=======================================================================
! Flushing
!=======================================================================

  subroutine flushing_velocity(Tin,    Sin,   &
                               phi,           &
                               hin,    hsn,   &
                               hilyr,         &
                               hpond,  apond, &
                               dt,     w,     &
                               wdown,  wup,   &
                               dhhead, w_up_max)
   
    ! calculate the vertical flushing Darcy velocity
    ! negative - downward flushing
    ! positive - upward flushing

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Tin       , & ! ice layer temperature (C)
         Sin       , & ! ice layer bulk salinity (ppt)
         phi           ! ice layer liquid fraction

    real(kind=dbl_kind), intent(in) :: &
         hilyr     , & ! ice layer thickness (m)
         hpond     , & ! melt pond thickness (m)
         apond     , & ! melt pond area (-)
         hsn       , & ! snow thickness (m)
         hin       , & ! ice thickness (m)
         dt            ! time step (s)

    real(kind=dbl_kind), intent(out) :: &
         w         , & ! vertical flushing Darcy flow rate (m s-1)
         wdown     , & ! downwards vertical flushing Darcy velocity (m/s)
         wup       , & ! upwards vertical flushing Darcy velocity (m/s)
         dhhead    , & ! hydraulic head (m)
         w_up_max      ! maximum upward flushing Darcy flow rate (m s-1) 

    real(kind=dbl_kind), parameter :: &
         viscosity_dyn = 1.79e-3_dbl_kind, & ! dynamic viscosity
         advection_limit = 0.005_dbl_kind      ! limit to fraction of brine in any layer that can be advected 

    real(kind=dbl_kind) :: &
         perm       , & ! ice layer permeability (m2)
         Mice       , & ! mass of ice (kg m-2)
         perm_harm  , & ! harmonic mean of ice permeability (m2)
         hocn       , & ! ocean surface height above ice base (m)
         hbrine     , & ! brine surface height above ice base (m)
         w_down_max , & ! maximum downward flushing Darcy flow rate (m s-1) 
         phi_min    , & ! minimum porosity in the mush
         wlimit     , & ! limit to w to avoid advecting all brine in layer
         rtmp

    integer(kind=int_kind) :: &
         k             ! ice layer index

    ! initialize
    Mice = c0
    perm_harm = c0
    phi_min = c1

    do k = 1, nilyr

       ! liquid fraction
       !phi = liquid_fraction(Tin(k), Sin(k))
       phi_min = min(phi_min,phi(k))

       ! permeability
       perm = permeability(phi(k))
       !perm = 3.0e-8 * phi(k)**3
       !xperm = 3.0e-8 * (0.05_dbl_kind)**3! * 100.0_dbl_kind

       ! ice mass
       Mice = Mice + hilyr * (phi(k) * rhow + (c1 - phi(k)) * rhoi)

       ! permeability harmonic mean
       perm_harm = perm_harm + c1 / (perm + 1e-30_dbl_kind)

    enddo ! k

    perm_harm = real(nilyr,dbl_kind) / perm_harm 

    ! calculate ocean surface height above bottom of ice
    hocn = (Mice + hpond * apond * rhow + hsn * rhos) / rhow

    ! calculate brine height above bottom of ice
    hbrine = hin + hpond

    ! hydraulic head
    !dhhead = hbrine - hocn
    !if (hbrine == hin) then
    !   dhhead = min(dhhead,c0)
    !endif

    if (hocn > hbrine) then
       ! snowice formation
       dhhead = hbrine - hocn
    else
       ! pond formation
       dhhead = hpond
    endif

    ! darcy flow through ice
    w = -(perm_harm * rhow * gravit * (dhhead / hin)) / viscosity_dyn

    ! maximum down flow to drain pond
    w_down_max = -(hpond * apond) / dt
    
    ! maximum up flow to flood to sea level
    w_up_max = max(((Mice + hsn * rhos - hin * rhow) * (rhoi - rhos)) / (dt * rhos * rhow), c0)

    !open(44,file='./history/wup.txt',position='append')
    !write(44,*) istep, w, w_up_max, min(max(w,w_down_max),w_up_max), hocn, hbrine, hin, hpond, dhhead, Mice, perm_harm
    !close(44)

    !open(44,file='./history/wdn.txt',position='append')
    !write(44,*) istep, w, abs(min(w,c0)), w_down_max, hpond, apond, &
    !     min(max(w,w_down_max),w_up_max), perm_harm, hbrine, hocn, dhhead, hbrine - hocn, hin
    !close(44)

    ! limit flow
    w = min(max(w,w_down_max),w_up_max)
    !w = c0!

    ! limit amount of brine that can be advected out of any particular layer
    wlimit = (advection_limit * phi_min * hilyr) / dt

    !open(44,file='./history/wup2.txt',position='append')
    !write(44,*) istep, wlimit, w, w * max(min(abs(wlimit/w),c1),c0)
    !close(44)

    if (abs(w) > puny) then
       w = w * max(min(abs(wlimit/w),c1),c0)
    else
       w = c0
    endif

#if flushing_up == 0
    ! no flooding
    w = min(w,c0)
#endif

#if flushing_down == 0
    ! no flushing
    w = max(w,c0)
#endif

    wdown = -min(w,c0) / hilyr
    wup   =  max(w,c0) / hilyr

  end subroutine flushing_velocity

!=======================================================================

  subroutine flushing_advection_flow_restriction(w_scaling,     &
                                                 wdown,  wup,   &
                                                 Sin0,   Sbr,   &
                                                 sss,    Spond, &
                                                 dt)

    real(kind=dbl_kind), intent(out) :: &
         w_scaling

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Sin0, &
         Sbr

    real(kind=dbl_kind), intent(in) :: &
         wdown, & ! downwards vertical flushing Darcy velocity (m/s)
         wup, &! upwards vertical flushing Darcy velocity (m/s)
         sss, &
         Spond, & ! melt pond salinity (ppt)
         dt

    real(kind=dbl_kind) :: &
         deltaS

    integer :: &
         k

    real(kind=dbl_kind), parameter :: &
         safety_factor = 0.1_dbl_kind

    ! calculate salinity change and scaling factor
    w_scaling = c1

    k = 1
    deltaS = wdown * (Spond    - Sbr(k)) + &
             wup   * (Sbr(k+1) - Sbr(k)) 

    w_scaling = min(w_scaling,c1/(max(deltaS*dt / (-safety_factor * Sin0(k)),c1)))

    do k = 2, nilyr-1

       deltaS = wdown * (Sbr(k-1) - Sbr(k)) + &
                wup   * (Sbr(k+1) - Sbr(k))

       w_scaling = min(w_scaling,c1/(max(deltaS*dt / (-safety_factor * Sin0(k)),c1)))

    enddo ! k

    k = nilyr
    deltaS = wdown * (Sbr(k-1) - Sbr(k)) + &
             wup   * (sss      - Sbr(k))

    w_scaling = min(w_scaling,c1/(max(deltaS*dt / (-safety_factor * Sin0(k)),c1)))
    !w_scaling = 1.0_dbl_kind

  end subroutine flushing_advection_flow_restriction

!=======================================================================

  subroutine flush_pond(w, hin, hpond, apond, dt)

    ! given a flushing velocity drain the meltponds

    real(kind=dbl_kind), intent(in) :: &
         w     , & ! vertical flushing Darcy flow rate (m s-1)
         hin   , & ! ice thickness (m)
         apond , & ! melt pond area (-)
         dt        ! time step (s)

    real(kind=dbl_kind), intent(inout) :: &
         hpond     ! melt pond thickness (m)
    
    real(kind=dbl_kind), dimension(1:nilyr) ::&
         phi       ! ice layer liquid fraction

    real(kind=dbl_kind), parameter :: &
         lambda_pond = c1 / (10.0_dbl_kind * 24.0_dbl_kind * 3600.0_dbl_kind), &
         hpond0 = 0.01_dbl_kind

    if (apond > c0 .and. hpond > c0) then

       hpond = hpond + (min(w,c0) * dt) / apond
       
       hpond = max(hpond, c0)

       !open(55,file="history/flush_pond.txt",position='append')
       !write(55,*) istep, w, abs(min(w,c0)), (min(w,c0) * dt), (min(w,c0) * dt) / apond, hpond
       !close(55)

       ! exponential decay of pond
       !write(65,*) istep1, hpond, hpond - lambda_pond * dt * (hpond + hpond0), - lambda_pond * dt * (hpond + hpond0)

       hpond = hpond - lambda_pond * dt * (hpond + hpond0)

       hpond = max(hpond, c0)

    endif

  end subroutine flush_pond

 !=======================================================================

  subroutine flood_ice(w,      dt,       &
                       dhhead, w_up_max, &
                       hsn,    hin,      &
                       hslyr,  hilyr,    & 
                       qsn,    qin,      &
                       phi,              &
                       Sin,    Sbr,      &
                       qbr,    snoice,   &
                       eadded, sadded)

    ! given upwards flushing brine flow calculate amount of snow ice and
    ! convert snow to ice with appropriate properties

    real(kind=dbl_kind), intent(in) :: &
         w           , & ! vertical flushing Darcy flow rate (m s-1)
         dt          , & ! time step (s)
         hsn         , & ! snow thickness (m)
         hin         , & ! ice thickness (m)
         dhhead      , & ! hydraulic head (m)
         w_up_max

    real(kind=dbl_kind), dimension(nslyr), intent(inout) :: &
         qsn             ! snow layer enthalpy (J m-2)

    real(kind=dbl_kind), dimension(nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-2)
         Sin         , & ! ice layer bulk salinity (ppt)
         phi

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Sbr         , & ! ice layer brine salinity (ppt)
         qbr             ! ice layer brine enthalpy (J m-2)

    real(kind=dbl_kind), intent(inout) :: &
         hslyr       , & ! ice layer thickness (m)
         hilyr           ! snow layer thickness (m)

    real(kind=dbl_kind), intent(out) :: &
         snoice      , & ! snow ice formation
         eadded      , & ! energy added by forming sea ice 
         sadded          ! salt added by forming sea ice

    real(kind=dbl_kind) :: &
         hin2        , & ! new ice thickness (m)
         hsn2        , & ! new snow thickness (m)
         hilyr2      , & ! new ice layer thickness (m)
         hslyr2      , & ! new snow layer thickness (m)
         dh          , & ! thickness of snowice formation (m)
         phi_snowice , & ! liquid fraction of new snow ice
         Sin_snowice , & ! bulk salinity of new snowice (ppt)
         qin_snowice , & ! ice enthalpy of new snowice (J m-2)
         qsn_snowice , & ! snow enthalpy of snow thats becoming snowice (J m-2)
         w_lateral   , &
         w_combined, &
         w_lateral_limit

    real(kind=dbl_kind), parameter :: &
         dhhead_lateral = -0.1_dbl_kind

    integer :: k
    real(kind=dbl_kind) :: einit, sinit
    real(kind=dbl_kind) :: efinal, sfinal
    real(kind=dbl_kind) :: euncon, suncon

    eadded = c0
    sadded = c0
    snoice = c0

    !write(61,*) istep1, w, dhhead, w_up_max, hsn

    if (w > c0 .or. dhhead < dhhead_lateral) then 

       ! initial energy and salt
       !einit = c0
       !sinit = c0
       
       !do k = 1, nslyr
       !   einit = einit + qsn(k) * hslyr
       !enddo ! k
       
       !do k = 1, nilyr
       !   einit = einit + qin(k) * hilyr
       !   sinit = sinit + Sin(k) * hilyr
       !enddo ! k
      
       ! lateral flow
       w_lateral = w_up_max * (abs(min(max(dhhead - dhhead_lateral,-hsn),c0)) / hsn)
       !w_lateral = c0

       ! combined flow effective rate
       w_combined = w + w_lateral

       ! limit combined flow rate 
       w_combined = min(w_combined, w_up_max)

       ! sea ice fraction of newly formed snow ice
       phi_snowice = (c1 - rhos / rhoi)

       ! calculate thickness of new ice added
       dh = (w_combined * dt) / phi_snowice
       
       ! enthalpy of snow that becomes snowice
       call enthalpy_snow_snowice(dh, hsn, qsn, qsn_snowice)
       
       ! change thicknesses
       hin2 = hin + dh
       hsn2 = hsn - dh
       
       hilyr2 = hin2 / real(nilyr,dbl_kind)
       hslyr2 = hsn2 / real(nslyr,dbl_kind)
       
       ! properties of new snow ice
       Sin_snowice = phi_snowice * Sbr(1)
       qin_snowice = phi_snowice * qbr(1) + qsn_snowice
       
       ! change snow properties
       call update_vertical_tracers_snow(qsn, hslyr, hslyr2)
       
       ! change ice properties
       call update_vertical_tracers_ice(qin, hilyr, hilyr2, &
                                        hin, hin2, qin_snowice)
       call update_vertical_tracers_ice(Sin, hilyr, hilyr2, &
                                        hin, hin2, Sin_snowice)
       call update_vertical_tracers_ice(phi, hilyr, hilyr2, &
                                        hin, hin2, phi_snowice)

       ! change thicknesses
       hilyr = hilyr2
       hslyr = hslyr2
       snoice = dh
       
       ! final energy and salt
       !efinal = c0
       !sfinal = c0
       
       !do k = 1, nslyr
       !   efinal = efinal + qsn(k) * hslyr
       !enddo ! k
       
       !do k = 1, nilyr
       !   efinal = efinal + qin(k) * hilyr
       !   sfinal = sfinal + Sin(k) * hilyr
       !enddo ! k
       
       eadded = w_combined * qbr(1)
       sadded = w_combined * Sbr(1)
       
       !euncon = efinal - einit - eadded * dt
       !suncon = sfinal - sinit - sadded * dt
       
       !open(11,file='./history/flood.txt',position='append')
       !write(11,*) istep, einit, efinal, eadded, euncon, euncon / einit, &
       !                   sinit, sfinal, sadded, suncon, suncon / sinit
       !close(11)

    endif
    
  end subroutine flood_ice

!=======================================================================

  subroutine enthalpy_snow_snowice(dh, hsn, qsn, qsn_snowice)

    ! determine enthalpy of the snow being converted to snow ice
    
    real(kind=dbl_kind), intent(in) :: &
         dh      , & ! thickness of new snowice formation (m)
         hsn         ! initial snow thickness

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         qsn         ! snow layer enthalpy (J m-2)

    real(kind=dbl_kind), intent(out) :: &
         qsn_snowice ! enthalpy of snow becoming snowice (J m-2)

    real(kind=dbl_kind) :: &
         rnlyr       ! real value of number of snow layers turning to snowice

    integer(kind=int_kind) :: &
         nlyr    , & ! no of snow layers involved in snow ice formation
         k           ! snow layer index

    ! snow depth and snow layers affected by snowice formation
    rnlyr = (dh / hsn) * nslyr
    nlyr = min(floor(rnlyr),nslyr-1)

    qsn_snowice = c0

    ! loop over full snow layers affected
    do k = nslyr, nslyr-nlyr+1, -1

       qsn_snowice = qsn_snowice + qsn(k) / rnlyr

    enddo ! k

    ! partially converted snow layer
    qsn_snowice = qsn_snowice + &
                  ((rnlyr - real(nlyr,dbl_kind)) / rnlyr) * qsn(nslyr)

    qsn_snowice = qsn_snowice

    !open(11,file='./history/qsn_snowice.txt',position='append')
    !write(11,*) istep, dh, hsn, rnlyr, nlyr, (rnlyr - real(nlyr,dbl_kind)), qsn_snowice, qsn(nslyr), qsn_snowice - qsn(nslyr)
    !close(11)

  end subroutine enthalpy_snow_snowice

!=======================================================================
  
  subroutine update_vertical_tracers_snow(trc, hlyr1, hlyr2)

    ! given some snow ice formation regrid snow layers
    
    real(kind=dbl_kind), dimension(1:nslyr), intent(inout) :: &
         trc         ! vertical tracer

    real(kind=dbl_kind), intent(in) :: &
         hlyr1   , & ! old cell thickness
         hlyr2       ! new cell thickness

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         trc2        ! temporary array for updated tracer
    
    ! vertical indexes for old and new grid
    integer(kind=int_kind) :: &
         k1      , & ! vertical index for old grid
         k2          ! vertical index for new grid
    
    real(kind=dbl_kind) :: &
         z1a     , & ! lower boundary of old cell
         z1b     , & ! upper boundary of old cell
         z2a     , & ! lower boundary of new cell
         z2b     , & ! upper boundary of new cell
         overlap     ! overlap between old and new cell
    
    ! loop over new grid cells
    do k2 = 1, nslyr
       
       ! initialize new tracer
       trc2(k2) = c0
       
       ! calculate upper and lower boundary of new cell
       z2a = (k2 - 1) * hlyr2
       z2b = k2       * hlyr2
       
       ! loop over old grid cells
       do k1 = 1, nslyr
          
          ! calculate upper and lower boundary of old cell
          z1a = (k1 - 1) * hlyr1
          z1b = k1       * hlyr1
          
          ! calculate overlap between old and new cell
          overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)
          
          ! aggregate old grid cell contribution to new cell
          trc2(k2) = trc2(k2) + overlap * trc(k1)
          
       enddo ! k1

       ! renormalize new grid cell
       trc2(k2) = trc2(k2) / hlyr2
       
    enddo ! k2
    
    ! update vertical tracer array with the adjusted tracer
    trc = trc2

  end subroutine update_vertical_tracers_snow

!=======================================================================
  
  subroutine update_vertical_tracers_ice(trc, hlyr1, hlyr2, &
                                         h1, h2, trc0)

    ! given some snow ice formation regrid ice layers

    real(kind=dbl_kind), dimension(1:nilyr), intent(inout) :: &
         trc         ! vertical tracer
    
    real(kind=dbl_kind), intent(in) :: &
         hlyr1 , &   ! old cell thickness
         hlyr2 , &   ! new cell thickness
         h1    , &   ! old total thickness
         h2    , &   ! new total thickness
         trc0        ! tracer value of added snow ice on ice top
    
    real(kind=dbl_kind), dimension(1:nilyr) :: &
         trc2        ! temporary array for updated tracer
    
    integer(kind=int_kind) :: &
         k1 , &      ! vertical indexes for old grid
         k2          ! vertical indexes for new grid
    
    real(kind=dbl_kind) :: &
         z1a     , & ! lower boundary of old cell
         z1b     , & ! upper boundary of old cell
         z2a     , & ! lower boundary of new cell
         z2b     , & ! upper boundary of new cell
         overlap     ! overlap between old and new cell
    
    ! loop over new grid cells
    do k2 = 1, nilyr
       
       ! initialize new tracer
       trc2(k2) = c0
       
       ! calculate upper and lower boundary of new cell
       z2a = (k2 - 1) * hlyr2
       z2b = k2       * hlyr2

       ! calculate upper and lower boundary of added snow ice at top
       z1a = c0
       z1b = h2 - h1
       
       ! calculate overlap between added ice and new cell
       overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)
       
       ! aggregate added ice contribution to new cell
       trc2(k2) = trc2(k2) + overlap * trc0

       ! loop over old grid cells
       do k1 = 1, nilyr
          
          ! calculate upper and lower boundary of old cell
          z1a = (k1 - 1) * hlyr1 + h2 - h1
          z1b = k1       * hlyr1 + h2 - h1
          
          ! calculate overlap between old and new cell
          overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)
          
          ! aggregate old grid cell contribution to new cell
          trc2(k2) = trc2(k2) + overlap * trc(k1)
          
       enddo ! k1

       ! renormalize new grid cell
       trc2(k2) = trc2(k2) / hlyr2
       
    enddo ! k2
    
    ! update vertical tracer array with the adjusted tracer
    trc = trc2

  end subroutine update_vertical_tracers_ice

!=======================================================================
! Physical Quantities
!=======================================================================

  subroutine conductivity_mush_array(qin, Sin, km)

    ! detemine the conductivity of the mush from enthalpy and salinity
    
    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         qin , & ! ice layer enthalpy (J m-3) 
         Sin     ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(1:nilyr), intent(out) :: &
         km      ! ice layer conductivity (W m-1 K-1)
    
    integer(kind=int_kind) :: &
         k       ! ice layer index
    
    do k = 1, nilyr
      
       km(k) = &
            heat_conductivity(temperature_mush(qin(k), Sin(k)), Sin(k))
       
    enddo ! k

  end subroutine conductivity_mush_array

!=======================================================================

  function density_brine(Sbr) result(rho)
    
    ! density of brine from brine salinity

    real(kind=dbl_kind), intent(in) :: &
         Sbr ! brine salinity (ppt)

    real(kind=dbl_kind) :: &
         rho ! brine density (kg m-3)
    
    real(kind=dbl_kind), parameter :: &
         a = 1000.3_dbl_kind    , & ! zeroth empirical coefficient
         b = 0.78237_dbl_kind   , & ! linear empirical coefficient
         c = 2.8008e-4_dbl_kind     ! quadratic empirical coefficient
    
    rho = a + b * Sbr + c * Sbr**2
                
  end function density_brine

!=======================================================================
! Snow
!=======================================================================

  subroutine conductivity_snow_array(ks)

    ! heat conductivity of the snow

    real(kind=dbl_kind), dimension(1:nslyr), intent(out) :: &
         ks ! snow layer conductivity (W m-1 K-1)

    ks = ksno

  end subroutine conductivity_snow_array

!=======================================================================
  
  function enthalpy_snow(Tsn) result(qsn)
    
    ! enthalpy of snow from snow temperature

    real(kind=dbl_kind), intent(in) :: &
         Tsn ! snow layer temperature (C)

    real(kind=dbl_kind) :: &
         qsn ! snow layer enthalpy (J m-3) 
    
    qsn = -rhos * (-cp_ice * Tsn + Lfresh)
    
  end function enthalpy_snow

!=======================================================================
  
  function temperature_snow(qsn) result(Tsn)
    
    ! temperature of snow from the snow enthalpy

    real(kind=dbl_kind), intent(in) :: &
         qsn ! snow layer enthalpy (J m-3) 

    real(kind=dbl_kind) :: &
         Tsn ! snow layer temperature (C)
    
    real(kind=dbl_kind), parameter :: &
         A = c1 / (rhos * cp_ice) , &
         B = Lfresh / cp_ice
    
    Tsn = A * qsn + B

  end function temperature_snow

!=======================================================================
! Mushy Layer Formulation - Assur (1958) liquidus
!=======================================================================

  function liquidus_brine_salinity_mush(Tin) result(Sbr)

    ! liquidus relation: equilibrium brine salinity as function of temperature
    ! based on empirical data from Assur (1958)

    real(kind=dbl_kind), intent(in) :: &
         Tin          ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         Sbr          ! ice brine salinity (ppt)

    real(kind=dbl_kind) :: &
         t_high   , & ! mask for high temperature liquidus region
         lsubzero     ! mask for sub-zero temperatures

    t_high   = merge(c1, c0, (Tin > Tb_liq))
    lsubzero = merge(c1, c0, (Tin <= c0))

    Sbr = ((Tin + J1_liq) / (K1_liq * Tin + L1_liq)) * t_high + &
          ((Tin + J2_liq) / (K2_liq * Tin + L2_liq)) * (c1 - t_high)
!echmod    Sbr = -Tin/0.054 !echmod

    Sbr = Sbr * lsubzero

  end function liquidus_brine_salinity_mush

!=======================================================================

  function liquidus_temperature_mush(Sbr) result(Tin)

    ! liquidus relation: equilibrium temperature as function of brine salinity
    ! based on empirical data from Assur (1958)

    real(kind=dbl_kind), intent(in) :: &
         Sbr    ! ice brine salinity (ppt)

    real(kind=dbl_kind) :: &
         Tin    ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         t_high ! mask for high temperature liquidus region

    t_high = merge(c1, c0, (Sbr <= Sb_liq))

    Tin = ((Sbr / (M1_liq + N1_liq * Sbr)) + O1_liq) * t_high + &
          ((Sbr / (M2_liq + N2_liq * Sbr)) + O2_liq) * (c1 - t_high)
!echmod    Tin = -Sbr*0.054 !echmod

  end function liquidus_temperature_mush

!=======================================================================

  function enthalpy_mush(Tin, Sin) result(qin)

    ! enthalpy of mush from mush temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         Tin , & ! ice layer temperature (C)
         Sin     ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         qin     ! ice layer enthalpy (J m-3) 

    real(kind=dbl_kind) :: &
         phi     ! ice liquid fraction 

    phi = liquid_fraction(Tin, Sin)
    
    qin = phi * (cp_ocn * rhow - cp_ice * rhoi) * Tin + &
          rhoi * cp_ice * Tin - (c1 - phi) * rhoi * Lfresh

  end function enthalpy_mush

!=======================================================================

  function enthalpy_mush_liquid_fraction(Tin, phi) result(qin)

    ! enthalpy of mush from mush temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         Tin , & ! ice layer temperature (C)
         phi     ! liquid fraction

    real(kind=dbl_kind) :: &
         qin     ! ice layer enthalpy (J m-3) 

    qin = phi * (cp_ocn * rhow - cp_ice * rhoi) * Tin + &
          rhoi * cp_ice * Tin - (c1 - phi) * rhoi * Lfresh

  end function enthalpy_mush_liquid_fraction

!=======================================================================

  function enthalpy_of_melting(Sin) result(qm)

    ! enthalpy of melting of mush
    ! energy needed to fully melt mush (T < 0)

    real(kind=dbl_kind), intent(in) :: &
         Sin ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         qm  ! melting ice enthalpy (J m-3)

    qm = cp_ocn * rhow * liquidus_temperature_mush(Sin)

  end function enthalpy_of_melting

!=======================================================================

  function enthalpy_brine(Tin) result(qbr)

    ! enthalpy of brine (fully liquid)

    real(kind=dbl_kind), intent(in) :: &
         Tin ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         qbr ! brine enthalpy (J m-3)

    qbr = cp_ocn * rhow * Tin

  end function enthalpy_brine

!=======================================================================

  function temperature_mush(qin, Sin) result(Tin)

    ! temperature of mush from mush enthalpy

    real(kind=dbl_kind), intent(in) :: &
         qin    , & ! ice enthalpy (J m-3) 
         Sin        ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         Tin        ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         qb     , & ! liquidus break enthalpy
         q0     , & ! fully melted enthalpy
         A      , & ! quadratic equation A parameter
         B      , & ! quadratic equation B parameter
         C      , & ! quadratic equation C parameter
         S_low  , & ! mask for salinity less than the liquidus break salinity
         t_high , & ! mask for high temperature liquidus region
         t_low  , & ! mask for low temperature liquidus region
         q_melt     ! mask for all mush melted

    ! just melted enthalpy
    S_low = merge(c1, c0, (Sin < Sb_liq))
    q0 = ((F1_liq * Sin) / (G1_liq + Sin) + H1_liq) * S_low + &
         ((F2_liq * Sin) / (G2_liq + Sin) + H2_liq) * (c1 - S_low)
    q_melt = merge(c1, c0, (qin > q0))

    ! break enthalpy
    qb = D_liq * Sin + E_liq
    t_high = merge(c1, c0, (qin > qb))
    t_low = c1 - t_high

    ! quadratic values
    A = (AS1_liq * Sin                 + AC1_liq) * t_high + &
        (AS2_liq * Sin                 + AC2_liq) * t_low

    B = (BS1_liq * Sin + BQ1_liq * qin + BC1_liq) * t_high + &
        (BS2_liq * Sin + BQ2_liq * qin + BC2_liq) * t_low

    C = (CS1_liq * Sin + CQ1_liq * qin + CC1_liq) * t_high + &
        (CS2_liq * Sin + CQ2_liq * qin + CC2_liq) * t_low

    Tin = (-B + sqrt(max(B**2 - c4 * A * C,puny))) / (c2 * A)

    ! change T if all melted
    Tin = q_melt * qin * I_liq + (c1 - q_melt) * Tin

  end function temperature_mush

!=======================================================================

  function temperature_mush_liquid_fraction(qin, phi) result(Tin)

    ! temperature of mush from mush enthalpy

    real(kind=dbl_kind), intent(in) :: &
         qin    , & ! ice enthalpy (J m-3) 
         phi        ! liquid fraction

    real(kind=dbl_kind) :: &
         Tin        ! ice layer temperature (C)

    Tin = (qin + (c1 - phi) * rhoi * Lfresh) / &
          (phi * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice)

  end function temperature_mush_liquid_fraction

!=======================================================================

  function heat_conductivity(Tin, Sin) result(km)
    
    ! msuh heat conductivity from mush temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         Tin               , & ! ice layer temperature (C)
         Sin                   ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         km                    ! ice layer conductivity (W m-1 K-1)
    
    real(kind=dbl_kind), parameter :: &
         ki = 2.3_dbl_kind , & ! fresh ice conductivity (W m-1 K-1)
         kb = 0.3_dbl_kind     ! brine conductivity (W m-1 K-1)

    real(kind=dbl_kind) :: &
         phi                   ! liquid fraction

    phi = liquid_fraction(Tin, Sin)

    km = phi * (kb - ki) + ki

  end function heat_conductivity

  !=======================================================================

  function liquid_fraction(Tin, Sin) result(phi)

    ! liquid fraction of mush from mush temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         Tin , & ! ice layer temperature (C)
         Sin     ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         phi , & ! liquid fraction
         Sbr     ! brine salinity (ppt)

    Sbr = max(liquidus_brine_salinity_mush(Tin),puny)
    phi = Sin / max(Sbr, Sin)

  end function liquid_fraction

!=======================================================================
! debug
!=======================================================================

  subroutine temperature_changes_column_debug(istep1,   dt,       & 
                                              rhoa,     flw,      &
                                              potT,     Qa,       &
                                              shcoef,   lhcoef,   &
                                              fswsfc,   fswint,   &
                                              fswthrun,           &
                                              Sswabs,   Iswabs,   &
                                              hilyr,    hslyr,    &
                                              qin,      Tin,      &
                                              qsn,      Tsn,      &
                                              Sin,                &
                                              hpond,    apond,    & !!
                                              Tsf,      Tbot,     &
                                              sss,                & !!
                                              fsensn,   flatn,    &
                                              fswabsn,  flwoutn,  &
                                              fsurfn,             &
                                              fcondtopn,fcondbot, &
                                              fadvocn,  snoice,   & !!
                                              einit,    l_stop)

    ! this routine is useful for debugging situations when the Newton loop fails to converge
    ! run with this routine and with type="normal_run" - upon failure of the solver the state of
    ! that solver is saved to file. Rerun with type="load_saved_state" and cice will try to resolve
    ! the failed state the first time temperature_chnages is called. This vastly speeds debugging
    ! as you dont need to wait for the whole simulation until it fails again.

    use ice_fileunits, only: nu_jfnkdiag

    integer(kind=int_kind), intent(in) :: &
         istep1          ! time step index (diagnostic only)
    
    real(kind=dbl_kind), intent(in) :: &
         dt              ! time step
    
    real(kind=dbl_kind), intent(in) :: &
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         Tbot        , & ! ice bottom surface temperature (deg C)
         sss             ! sea surface salinity

    real(kind=dbl_kind), intent(inout) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         fswint      , & ! SW absorbed in ice interior below surface (W m-2)
         fswthrun        ! SW through ice to ocean         (W m-2)
    
    real(kind=dbl_kind), intent(inout) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr           ! snow layer thickness (m)

    real(kind=dbl_kind), intent(in) :: &
         einit           ! initial energy of melting (J m-2)
    
    real(kind=dbl_kind), dimension (nslyr), &
         intent(inout) :: &
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)
    
    real(kind=dbl_kind), dimension (nilyr), &
         intent(inout) :: &
         Iswabs          ! SW radiation absorbed in ice layers (W m-2)
    
    real(kind=dbl_kind), intent(inout):: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fswabsn     , & ! shortwave absorbed by ice (W m-2)
         flwoutn         ! upward LW at surface (W m-2)
    
    real(kind=dbl_kind), intent(out):: &
         fcondbot    , & ! downward cond flux at bottom surface (W m-2)
         fadvocn     , & ! flow of heat to ocean due to advection
         snoice          ! snow-ice formation

    real (kind=dbl_kind), intent(inout):: &
         Tsf         , & ! ice/snow surface temperature(C)
         hpond       , & ! melt pond depth (m)
         apond           ! melt pond area

    real(kind=dbl_kind), dimension (nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3)
         Tin         , & ! internal ice layer temperatures
         Sin             ! internal ice layer salinities
  
    real(kind=dbl_kind), dimension (nslyr), intent(inout) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         Tsn             ! internal snow layer temperatures
    
    logical(kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

    ! local versions
    integer(kind=int_kind) :: &
         istep1_in,    istep1_sv         ! time step index (diagnostic only)

    real(kind=dbl_kind) :: &
         dt_in,        dt_sv         , & ! time step
         rhoa_in,      rhoa_sv       , & ! air density (kg/m^3)
         flw_in,       flw_sv        , & ! incoming longwave radiation (W/m^2)
         potT_in,      potT_sv       , & ! air potential temperature  (K)
         Qa_in,        Qa_sv         , & ! specific humidity (kg/kg)
         shcoef_in,    shcoef_sv     , & ! transfer coefficient for sensible heat
         lhcoef_in,    lhcoef_sv     , & ! transfer coefficient for latent heat
         Tbot_in,      Tbot_sv       , & ! ice bottom surface temperature (deg C)
         sss_in,       sss_sv        , & ! sea surface salinity
         fswsfc_in,    fswsfc_sv     , & ! SW absorbed at ice/snow surface (W m-2)
         fswint_in,    fswint_sv     , & ! SW absorbed in ice interior below surface (W m-2)
         fswthrun_in,  fswthrun_sv   , & ! SW through ice to ocean         (W m-2)
         hilyr_in,     hilyr_sv      , & ! ice layer thickness (m)
         hslyr_in,     hslyr_sv      , & ! snow layer thickness (m)
         einit_in,     einit_sv      , & ! initial energy of melting (J m-2)
         fsurfn_in,    fsurfn_sv     , & ! net flux to top surface, excluding fcondtopn
         fcondtopn_in, fcondtopn_sv  , & ! downward cond flux at top surface (W m-2)
         fsensn_in,    fsensn_sv     , & ! surface downward sensible heat (W m-2)
         flatn_in,     flatn_sv      , & ! surface downward latent heat (W m-2)
         fswabsn_in,   fswabsn_sv    , & ! shortwave absorbed by ice (W m-2)
         flwoutn_in,   flwoutn_sv    , & ! upward LW at surface (W m-2)
         fcondbot_in,  fcondbot_sv   , & ! downward cond flux at bottom surface (W m-2)
         fadvocn_in,   fadvocn_sv    , & ! flow of heat to ocean due to advection
         snoice_in,    snoice_sv     , & ! snow-ice formation
         Tsf_in,       Tsf_sv        , & ! ice/snow surface temperature (C)
         hpond_in,     hpond_sv      , & ! melt pond thickness (m)
         apond_in,     apond_sv          ! melt pond area

    real(kind=dbl_kind), dimension (nslyr) :: &
         Sswabs_in, Sswabs_sv        , & ! SW radiation absorbed in snow layers (W m-2)
         qsn_in,    qsn_sv           , & ! snow layer enthalpy (J m-3)
         Tsn_in,    Tsn_sv               ! internal snow layer temperatures

    real(kind=dbl_kind), dimension (nilyr) :: &
         Iswabs_in, Iswabs_sv        , & ! SW radiation absorbed in ice layers (W m-2)
         qin_in,    qin_sv           , & ! ice layer enthalpy (J m-3)
         Tin_in,    Tin_sv           , & ! internal ice layer temperatures
         Sin_in,    Sin_sv               ! internal ice layer salinities

    logical(kind=log_kind) :: &
         l_stop_in, l_stop_sv            ! if true, print diagnostics and abort model

    logical(kind=log_kind) :: &
         lstop  ! solver failure flag 

    ! formatting type for save state
    character(len=200), parameter :: form_type = "unformatted"

    ! run type
#if load_save_state == 0
    character(len=200), parameter :: type = "normal_run"
#else
    character(len=200), parameter :: type = "load_saved_state"
#endif

    ! filename for file saving
    character(len=200) :: fileout

    ! always save the state
    logical(kind=log_kind), parameter :: l_save_state_always = .false.

    ! normal running conditions
    if (trim(type) == "normal_run") then

       ! make copies of input state
       istep1_in    = istep1    ; istep1_sv    = istep1
       dt_in        = dt        ; dt_sv        = dt
       rhoa_in      = rhoa      ; rhoa_sv      = rhoa
       flw_in       = flw       ; flw_sv       = flw
       potT_in      = potT      ; potT_sv      = potT
       Qa_in        = Qa        ; Qa_sv        = Qa
       shcoef_in    = shcoef    ; shcoef_sv    = shcoef
       lhcoef_in    = lhcoef    ; lhcoef_sv    = lhcoef
       Tbot_in      = Tbot      ; Tbot_sv      = Tbot
       sss_in       = sss       ; sss_sv       = sss
       fswsfc_in    = fswsfc    ; fswsfc_sv    = fswsfc
       fswint_in    = fswint    ; fswint_sv    = fswint
       fswthrun_in  = fswthrun  ; fswthrun_sv  = fswthrun
       hilyr_in     = hilyr     ; hilyr_sv     = hilyr
       hslyr_in     = hslyr     ; hslyr_sv     = hslyr
       einit_in     = einit     ; einit_sv     = einit
       Sswabs_in    = Sswabs    ; Sswabs_sv    = Sswabs
       Iswabs_in    = Iswabs    ; Iswabs_sv    = Iswabs
       hpond_in     = hpond     ; hpond_sv     = hpond
       apond_in     = apond     ; apond_sv     = apond
       fsurfn_in    = fsurfn    ; fsurfn_sv    = fsurfn
       fcondtopn_in = fcondtopn ; fcondtopn_sv = fcondtopn
       fsensn_in    = fsensn    ; fsensn_sv    = fsensn
       flatn_in     = flatn     ; flatn_sv     = flatn
       fswabsn_in   = fswabsn   ; fswabsn_sv   = fswabsn
       flwoutn_in   = flwoutn   ; flwoutn_sv   = flwoutn
       fcondbot_in  = fcondbot  ; fcondbot_sv  = fcondbot
       fadvocn_in   = fadvocn   ; fadvocn_sv   = fadvocn
       snoice_in    = snoice    ; snoice_sv    = snoice
       Tsf_in       = Tsf       ; Tsf_sv       = Tsf
       qin_in       = qin       ; qin_sv       = qin
       Tin_in       = Tin       ; Tin_sv       = Tin
       Sin_in       = Sin       ; Sin_sv       = Sin
       qsn_in       = qsn       ; qsn_sv       = qsn
       Tsn_in       = Tsn       ; Tsn_sv       = Tsn
       l_stop_in    = l_stop    ; l_stop_sv    = l_stop

    else
       
       ! start from saved state
       if (form_type == "formatted") then

          ! read in formatted saved state
          open(nu_jfnkdiag,file='save_state_form.dat',form="formatted",action='read')
          
          read(nu_jfnkdiag,*) istep1_in    ; istep1_sv    = istep1_in
          read(nu_jfnkdiag,*) dt_in        ; dt_sv        = dt_in
          read(nu_jfnkdiag,*) rhoa_in      ; rhoa_sv      = rhoa_in
          read(nu_jfnkdiag,*) flw_in       ; flw_sv       = flw_in
          read(nu_jfnkdiag,*) potT_in      ; potT_sv      = potT_in
          read(nu_jfnkdiag,*) Qa_in        ; Qa_sv        = Qa_in
          read(nu_jfnkdiag,*) shcoef_in    ; shcoef_sv    = shcoef_in
          read(nu_jfnkdiag,*) lhcoef_in    ; lhcoef_sv    = lhcoef_in
          read(nu_jfnkdiag,*) Tbot_in      ; Tbot_sv      = Tbot_in
          read(nu_jfnkdiag,*) sss_in       ; sss_sv       = sss_in
          read(nu_jfnkdiag,*) fswsfc_in    ; fswsfc_sv    = fswsfc_in
          read(nu_jfnkdiag,*) fswint_in    ; fswint_sv    = fswint_in
          read(nu_jfnkdiag,*) fswthrun_in  ; fswthrun_sv  = fswthrun_in
          read(nu_jfnkdiag,*) hilyr_in     ; hilyr_sv     = hilyr_in
          read(nu_jfnkdiag,*) hslyr_in     ; hslyr_sv     = hslyr_in
          read(nu_jfnkdiag,*) einit_in     ; einit_sv     = einit_in
          read(nu_jfnkdiag,*) Sswabs_in    ; Sswabs_sv    = Sswabs_in
          read(nu_jfnkdiag,*) Iswabs_in    ; Iswabs_sv    = Iswabs_in
          read(nu_jfnkdiag,*) hpond_in     ; hpond_sv     = hpond_in
          read(nu_jfnkdiag,*) apond_in     ; apond_sv     = apond_in
          read(nu_jfnkdiag,*) fsurfn_in    ; fsurfn_sv    = fsurfn_in
          read(nu_jfnkdiag,*) fcondtopn_in ; fcondtopn_sv = fcondtopn_in
          read(nu_jfnkdiag,*) fsensn_in    ; fsensn_sv    = fsensn_in
          read(nu_jfnkdiag,*) flatn_in     ; flatn_sv     = flatn_in
          read(nu_jfnkdiag,*) fswabsn_in   ; fswabsn_sv   = fswabsn_in
          read(nu_jfnkdiag,*) flwoutn_in   ; flwoutn_sv   = flwoutn_in
          read(nu_jfnkdiag,*) fcondbot_in  ; fcondbot_sv  = fcondbot_in
          read(nu_jfnkdiag,*) fadvocn_in   ; fadvocn_sv   = fadvocn
          read(nu_jfnkdiag,*) snoice_in    ; snoice_sv    = snoice
          read(nu_jfnkdiag,*) Tsf_in       ; Tsf_sv       = Tsf_in
          read(nu_jfnkdiag,*) qin_in       ; qin_sv       = qin_in
          read(nu_jfnkdiag,*) Tin_in       ; Tin_sv       = Tin_in
          read(nu_jfnkdiag,*) Sin_in       ; Sin_sv       = Sin_in
          read(nu_jfnkdiag,*) qsn_in       ; qsn_sv       = qsn_in
          read(nu_jfnkdiag,*) Tsn_in       ; Tsn_sv       = Tsn_in
          read(nu_jfnkdiag,*) l_stop_in    ; l_stop_sv    = l_stop_in        
          
          close(nu_jfnkdiag)
          
       else if (form_type == "unformatted") then

          ! read in binary saved state
          open(nu_jfnkdiag,file='save_state_unform.dat',form="unformatted",action='read')
          
          read(nu_jfnkdiag) istep1_in    ; istep1_sv    = istep1_in
          read(nu_jfnkdiag) dt_in        ; dt_sv        = dt_in
          read(nu_jfnkdiag) rhoa_in      ; rhoa_sv      = rhoa_in
          read(nu_jfnkdiag) flw_in       ; flw_sv       = flw_in
          read(nu_jfnkdiag) potT_in      ; potT_sv      = potT_in
          read(nu_jfnkdiag) Qa_in        ; Qa_sv        = Qa_in
          read(nu_jfnkdiag) shcoef_in    ; shcoef_sv    = shcoef_in
          read(nu_jfnkdiag) lhcoef_in    ; lhcoef_sv    = lhcoef_in
          read(nu_jfnkdiag) Tbot_in      ; Tbot_sv      = Tbot_in
          read(nu_jfnkdiag) sss_in       ; sss_sv       = sss_in
          read(nu_jfnkdiag) fswsfc_in    ; fswsfc_sv    = fswsfc_in
          read(nu_jfnkdiag) fswint_in    ; fswint_sv    = fswint_in
          read(nu_jfnkdiag) fswthrun_in  ; fswthrun_sv  = fswthrun_in
          read(nu_jfnkdiag) hilyr_in     ; hilyr_sv     = hilyr_in
          read(nu_jfnkdiag) hslyr_in     ; hslyr_sv     = hslyr_in
          read(nu_jfnkdiag) einit_in     ; einit_sv     = einit_in
          read(nu_jfnkdiag) Sswabs_in    ; Sswabs_sv    = Sswabs_in
          read(nu_jfnkdiag) Iswabs_in    ; Iswabs_sv    = Iswabs_in
          read(nu_jfnkdiag) hpond_in     ; hpond_sv     = hpond_in
          read(nu_jfnkdiag) apond_in     ; apond_sv     = apond_in
          read(nu_jfnkdiag) fsurfn_in    ; fsurfn_sv    = fsurfn_in
          read(nu_jfnkdiag) fcondtopn_in ; fcondtopn_sv = fcondtopn_in
          read(nu_jfnkdiag) fsensn_in    ; fsensn_sv    = fsensn_in
          read(nu_jfnkdiag) flatn_in     ; flatn_sv     = flatn_in
          read(nu_jfnkdiag) fswabsn_in   ; fswabsn_sv   = fswabsn_in
          read(nu_jfnkdiag) flwoutn_in   ; flwoutn_sv   = flwoutn_in
          read(nu_jfnkdiag) fcondbot_in  ; fcondbot_sv  = fcondbot_in
          read(nu_jfnkdiag) fadvocn_in   ; fadvocn_sv   = fadvocn
          read(nu_jfnkdiag) snoice_in    ; snoice_sv    = snoice
          read(nu_jfnkdiag) Tsf_in       ; Tsf_sv       = Tsf_in
          read(nu_jfnkdiag) qin_in       ; qin_sv       = qin_in
          read(nu_jfnkdiag) Tin_in       ; Tin_sv       = Tin_in
          read(nu_jfnkdiag) Sin_in       ; Sin_sv       = Sin_in
          read(nu_jfnkdiag) qsn_in       ; qsn_sv       = qsn_in
          read(nu_jfnkdiag) Tsn_in       ; Tsn_sv       = Tsn_in
          read(nu_jfnkdiag) l_stop_in    ; l_stop_sv    = l_stop_in        
          
          close(nu_jfnkdiag)

       endif

    endif

    ! perform the temperature changes calculation
    call temperature_changes_column(istep1_in,    dt_in,       & 
                                    rhoa_in,      flw_in,      &
                                    potT_in,      Qa_in,       &
                                    shcoef_in,    lhcoef_in,   &
                                    fswsfc_in,    fswint_in,   &
                                    fswthrun_in,               &
                                    Sswabs_in,    Iswabs_in,   &
                                    hilyr_in,     hslyr_in,    &
                                    qin_in,       Tin_in,      &
                                    qsn_in,       Tsn_in,      &
                                    Sin_in,                    &
                                    hpond_in,     apond_in,    &
                                    Tsf_in,       Tbot_in,     &
                                    sss_in,                    &
                                    fsensn_in,    flatn_in,    &
                                    fswabsn_in,   flwoutn_in,  &
                                    fsurfn_in,                 &
                                    fcondtopn_in, fcondbot_in, &
                                    fadvocn_in,   snoice_in,   &
                                    einit_in,     lstop)

    ! normal running conditions
    if (trim(type) == "normal_run") then
       if (.not. lstop .and. .not. l_save_state_always) then
       
          ! Solution found - copy out solution to output
          fswsfc    = fswsfc_in   
          fswint    = fswint_in   
          fswthrun  = fswthrun_in  
          Sswabs    = Sswabs_in   
          Iswabs    = Iswabs_in
          hpond     = hpond_in
          apond     = apond_in
          fsurfn    = fsurfn_in   
          fcondtopn = fcondtopn_in
          fsensn    = fsensn_in   
          flatn     = flatn_in    
          fswabsn   = fswabsn_in  
          flwoutn   = flwoutn_in  
          fcondbot  = fcondbot_in 
          fadvocn   = fadvocn_in
          snoice    = snoice_in
          Tsf       = Tsf_in      
          qin       = qin_in      
          Tin       = Tin_in      
          Sin       = Sin_in      
          qsn       = qsn_in      
          Tsn       = Tsn_in     
          hilyr     = hilyr_in
          hslyr     = hslyr_in
          l_stop    = l_stop_in       
          
       else if (lstop .or. l_save_state_always) then

          ! solution failed - save the input state
          if (form_type == "formatted") then

             ! save the input state with formatted output
             write(nu_diag,*) "save state formatted: istep, n, i, j: ", istep, g_n, g_i, g_j
             
             if (l_save_state_always) then
                write(fileout,fmt='(a,i5.5,a)') "save_state_unform_",istep1,".dat"
             else
                fileout = "save_state_unform.dat"
             endif

             open(nu_jfnkdiag,file=fileout,form="formatted",action='write')
             
             write(nu_jfnkdiag,*) istep1_sv   
             write(nu_jfnkdiag,*) dt_sv       
             write(nu_jfnkdiag,*) rhoa_sv     
             write(nu_jfnkdiag,*) flw_sv      
             write(nu_jfnkdiag,*) potT_sv     
             write(nu_jfnkdiag,*) Qa_sv       
             write(nu_jfnkdiag,*) shcoef_sv   
             write(nu_jfnkdiag,*) lhcoef_sv   
             write(nu_jfnkdiag,*) Tbot_sv
             write(nu_jfnkdiag,*) sss_sv     
             write(nu_jfnkdiag,*) fswsfc_sv   
             write(nu_jfnkdiag,*) fswint_sv   
             write(nu_jfnkdiag,*) fswthrun_sv 
             write(nu_jfnkdiag,*) hilyr_sv    
             write(nu_jfnkdiag,*) hslyr_sv    
             write(nu_jfnkdiag,*) einit_sv    
             write(nu_jfnkdiag,*) Sswabs_sv   
             write(nu_jfnkdiag,*) Iswabs_sv
             write(nu_jfnkdiag,*) hpond_sv
             write(nu_jfnkdiag,*) apond_sv
             write(nu_jfnkdiag,*) fsurfn_sv   
             write(nu_jfnkdiag,*) fcondtopn_sv
             write(nu_jfnkdiag,*) fsensn_sv   
             write(nu_jfnkdiag,*) flatn_sv    
             write(nu_jfnkdiag,*) fswabsn_sv  
             write(nu_jfnkdiag,*) flwoutn_sv  
             write(nu_jfnkdiag,*) fcondbot_sv 
             write(nu_jfnkdiag,*) fadvocn_sv
             write(nu_jfnkdiag,*) snoice_sv
             write(nu_jfnkdiag,*) Tsf_sv      
             write(nu_jfnkdiag,*) qin_sv      
             write(nu_jfnkdiag,*) Tin_sv      
             write(nu_jfnkdiag,*) Sin_sv 
             write(nu_jfnkdiag,*) qsn_sv      
             write(nu_jfnkdiag,*) Tsn_sv      
             write(nu_jfnkdiag,*) l_stop_sv   
             
             close(nu_jfnkdiag)
             
             if (.not. l_save_state_always) stop

          else if (form_type == "unformatted") then

             ! save the input state with binary output
             write(nu_diag,*) "save state formatted: istep, n, i, j: ", istep, g_n, g_i, g_j

             if (l_save_state_always) then
                write(fileout,fmt='(a,i5.5,a)') "save_state_unform_",istep1,".dat"
             else
                fileout = "save_state_unform.dat"
             endif             
             
             open(nu_jfnkdiag,file=fileout,form="unformatted",action='write')
             
             write(nu_jfnkdiag) istep1_sv   
             write(nu_jfnkdiag) dt_sv       
             write(nu_jfnkdiag) rhoa_sv     
             write(nu_jfnkdiag) flw_sv      
             write(nu_jfnkdiag) potT_sv     
             write(nu_jfnkdiag) Qa_sv       
             write(nu_jfnkdiag) shcoef_sv   
             write(nu_jfnkdiag) lhcoef_sv   
             write(nu_jfnkdiag) Tbot_sv
             write(nu_jfnkdiag) sss_sv     
             write(nu_jfnkdiag) fswsfc_sv   
             write(nu_jfnkdiag) fswint_sv   
             write(nu_jfnkdiag) fswthrun_sv 
             write(nu_jfnkdiag) hilyr_sv    
             write(nu_jfnkdiag) hslyr_sv    
             write(nu_jfnkdiag) einit_sv    
             write(nu_jfnkdiag) Sswabs_sv   
             write(nu_jfnkdiag) Iswabs_sv
             write(nu_jfnkdiag) hpond_sv
             write(nu_jfnkdiag) apond_sv
             write(nu_jfnkdiag) fsurfn_sv   
             write(nu_jfnkdiag) fcondtopn_sv
             write(nu_jfnkdiag) fsensn_sv   
             write(nu_jfnkdiag) flatn_sv    
             write(nu_jfnkdiag) fswabsn_sv  
             write(nu_jfnkdiag) flwoutn_sv  
             write(nu_jfnkdiag) fcondbot_sv 
             write(nu_jfnkdiag) fadvocn_sv
             write(nu_jfnkdiag) snoice_sv
             write(nu_jfnkdiag) Tsf_sv      
             write(nu_jfnkdiag) qin_sv      
             write(nu_jfnkdiag) Tin_sv      
             write(nu_jfnkdiag) Sin_sv     
             write(nu_jfnkdiag) qsn_sv      
             write(nu_jfnkdiag) Tsn_sv      
             write(nu_jfnkdiag) l_stop_sv   
             
             close(nu_jfnkdiag)
             
             if (.not. l_save_state_always) then

                call abort_ice("ice_therm_mushy: solver failure: state saved")

             endif

          endif
          
       endif
    endif

    ! after loading check to see result of solver
    if (trim(type) == "load_saved_state") then
       if (.not. l_stop_in) then
          
          if (my_task == master_task) then
             call abort_ice("Saved state worked")
          endif
          
       else
          
          if (my_task == master_task) then
             call abort_ice("Saved state failed")
          endif

       endif
    endif

  end subroutine temperature_changes_column_debug

!=======================================================================

end module ice_therm_mushy

!=======================================================================
