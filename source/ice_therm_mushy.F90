#define debug 0
#define diag 0
#define fully_converged 0
#define preconditioning 1
#define diag_nonconverge_always 0
#define gmres_accurate_solve 0

#define flushing 0
#define flushing_down 0
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
  public :: temperature_mush, temperature_snow, liquid_fraction, &
       liquidus_brine_salinity_mush, liquidus_temperature_mush, &
       temperature_changes_salinity, init_therm_mushy, enthalpy_mush, &
       enthalpy_of_melting, add_new_ice_mushy

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

  !real(kind=dbl_kind), parameter :: a_rapid_mode = 0.25e-3_dbl_kind
  !real(kind=dbl_kind), parameter :: a_rapid_mode = 0.5e-3_dbl_kind ! normal
  !real(kind=dbl_kind), parameter :: a_rapid_mode = 1.0e-3_dbl_kind
  
  !real(kind=dbl_kind), parameter :: Rac_rapid_mode = 10.0_dbl_kind ! normal
  !real(kind=dbl_kind), parameter :: Rac_rapid_mode = 1.0_dbl_kind
  !real(kind=dbl_kind), parameter :: Rac_rapid_mode = 20.0_dbl_kind
  
  !real(kind=dbl_kind), parameter :: dSdt_slow_mode = c0
  !real(kind=dbl_kind), parameter :: dSdt_slow_mode = -0.3e-7_dbl_kind
  !real(kind=dbl_kind), parameter :: dSdt_slow_mode = -1.5e-7_dbl_kind ! normal
  !real(kind=dbl_kind), parameter :: dSdt_slow_mode = -7.5e-7_dbl_kind
  
  !real(kind=dbl_kind), parameter :: phi_c_slow_mode = 0.025_dbl_kind ! for simulations
  !real(kind=dbl_kind), parameter :: phi_c_slow_mode = 0.05_dbl_kind ! normal
  !real(kind=dbl_kind), parameter :: phi_c_slow_mode = 0.1_dbl_kind

  !-----------------------------------------------------------------
  ! solver solution sizes
  !-----------------------------------------------------------------
  integer, parameter :: &
       nyn_nosnow_cold = 2 * nilyr + 1         , & ! solution size for cold, no snow conditions
       nyn_nosnow_melt = 2 * nilyr             , & ! solution size for melting, no snow conditions
       nyn_snow_cold   = 2 * nilyr + nslyr + 1 , & ! solution size for cold, snow conditions
       nyn_snow_melt   = 2 * nilyr + nslyr     , & ! solution size for melting, nsnow conditions
       nynmax          = nyn_snow_cold             ! maximum solution size

  !-----------------------------------------------------------------
  ! JFNK Solver global variables 
  !-----------------------------------------------------------------
  integer, parameter :: &
       nit_newton_max = 100                         , & ! maximum number of newton iterations
       nit_gmres_max  = 100                         , & ! maximum number of GMRES matrix-vector multiplies
       nsize_krylov   = nynmax+4                    , & ! size of Krylov subspace
       wsizemax       = (nynmax+3)*(nsize_krylov+2) + &
                        (nsize_krylov+1)*nsize_krylov/2 ! size of the GMRES workspace
 
  integer, dimension(16) :: &
       ipar ! GMRES integer input/output parameters

  real(kind=dbl_kind), dimension(16) :: &
       fpar ! GMRES real input/output parameters

  !-----------------------------------------------------------------
  ! physical parameters for passing to residual function through solver
  !-----------------------------------------------------------------
  real(kind=dbl_kind) :: &
       Tsf0     ! Surface temperature (C) at beginning of timestep

  real(kind=dbl_kind), dimension(1:nslyr) :: &
       qsn0 , & ! snow layer enthalpy (J m-3) at beginning of timestep
       ks0  , & ! snow layer conductivity (W m-1 K-1) at beginning of timestep
       Tsn0     ! snow layer temperature (C) at beginning of timestep

  real(kind=dbl_kind), dimension(1:nslyr+1) :: &
       ksstar0  ! snow interface conductivity (W m-1 K-1) at beginning of timestep

  real(kind=dbl_kind), dimension(1:nilyr) :: &
       qin0 , & ! ice layer enthalpy (J m-3) at beginning of timestep
       km0  , & ! ice layer conductivity (W m-1 K-1) at beginning of timestep
       Tin0 , & ! ice layer temperature (C) at beginning of timestep
       Sin0     ! ice layer bulk salinity (ppt) at beginning of timestep

  real(kind=dbl_kind), dimension(1:nilyr+1) :: &
       kmstar0  ! ice interface conductivity (W m-1 K-1) at beginning of timestep

  real(kind=dbl_kind) :: &
       g_dt      , & ! time step (s)
       g_hin     , & ! ice thickness (m)
       g_hsn     , & ! snow thickness (m)
       g_hilyr   , & ! ice layer thickness (m)
       g_hslyr   , & ! snow layer thickness (m)
       g_dti     , & ! inverse time step (s-1)
       g_hilyri  , & ! inverse ice layer thickness (m-1)
       g_hslyri  , & ! inverse snow layer thickness (m-1)
       g_hilyri2 , & ! square inverse ice layer thickness (m-2)
       g_hslyri2 , & ! square inverse snow layer thickness (m-2)
       g_fswsfc  , & ! SW absorbed at ice/snow surface (W m-2)
       g_rhoa    , & ! air density (kg/m^3)
       g_flw     , & ! incoming longwave radiation (W/m^2)
       g_potT    , & ! air potential temperature (K)
       g_Qa      , & ! specific humidity (kg/kg)
       g_shcoef  , & ! transfer coefficient for sensible heat
       g_lhcoef  , & ! transfer coefficient for latent heat
       g_Tbot    , & ! ice bottom surface temperature (deg C)
       g_hpond   , & ! melt pond depth 
       g_apond   , & ! melt pond area
       g_sss     , & ! sea surface salinity (ppt)
       g_qocn        ! sea water enthalpy (J m-3)

  real(kind=dbl_kind), dimension(1:nslyr) :: &
       g_Sswabs ! SW radiation absorbed in snow layers (W m-2)

  real(kind=dbl_kind), dimension(1:nilyr) :: &
       g_Iswabs , & ! SW radiation absorbed in ice layers (W m-2)
       g_dSdt       ! gravity drainage desalination rate for slow mode (ppt s-1)

  real(kind=dbl_kind), dimension(0:nilyr) :: &
       g_q ! upward interface vertical Darcy flow (m s-1)

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

  subroutine init_therm_mushy()

    ! initialize anything related to the mushy thermodynamics at the start of the simulation

    call init_jfnk_solver()

  end subroutine init_therm_mushy

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
! Solver with seperation of cold and melting phases
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
                                        fcondtopn,fcondbot, &
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
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
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
         km          , & ! ice conductivity (W m-1 K-1)
         phi         , & ! ice liquid fraction (-)
         dqin_dTin   , & ! derivative of ice enthalpy wrt ice temperature
         dTin_dSin       ! derivative of ice temperature wrt ice bulk salinity
    
    real(kind=dbl_kind), dimension(1:nilyr+1) :: &
         Sbr         , & ! brine salinity (ppt)
         qbr             ! brine enthalpy (J m-3)

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         ks          , & ! snow conductivity (W m-1 K-1)
         dqsn_dTsn       ! derivative of snow enthalpy wrt snow temperature

    real(kind=dbl_kind) :: &
         hin         , & ! ice thickness (m)
         hsn         , & ! snow thickness (m)
         hslyr_min   , & ! minimum snow layer thickness (m)
         w           , & ! vertical flushing Darcy velocity (m/s)
         einit       , & ! initial total energy (J)
         efinal      , & ! final total energy (J)
         sinit       , & ! initial total salt content (ppt m)
         sfinal      , & ! final total salt content
         fadvsalt        ! heat flux to ocean from brine advection (W m-2)

    integer(kind=int_kind) :: &
         k               ! ice/snow layer index

    logical(kind=log_kind) :: &
         lsnow           ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), dimension(nynmax) :: &
         yscale_cold , & ! solution scaling vector for cold conditions
         yscale_melt     ! solution scaling vector for melting conditions

    lstop = .false.
    fadvsalt = c0
    snoice = c0

#if defined oned
    if (hilyr <= c0) return
#endif

    call check_inputs_outputs(Tsf,qsn,qin,Sin,Tsn,Tin,hilyr,"start:")

    Tsf0 = Tsf
    qsn0 = qsn
    qin0 = qin
    Sin0 = Sin
    Tsn0 = Tsn
    Tin0 = Tin

    hslyr_min = hs_min / real(nslyr, dbl_kind)

    lsnow = (hslyr > hslyr_min)

    hin = hilyr * real(nilyr,dbl_kind)

    g_dt      = dt
    g_hin     = hin
    g_hilyr   = hilyr
    g_dti     = c1 / dt    
    g_hilyri  = c1 / hilyr
    g_hilyri2 = c1 / hilyr**2
    g_fswsfc  = fswsfc
    g_rhoa    = rhoa
    g_flw     = flw
    g_potT    = potT
    g_Qa      = Qa
    g_shcoef  = shcoef
    g_lhcoef  = lhcoef
    g_Tbot    = Tbot
    g_hpond   = hpond
    g_apond   = apond
    g_sss     = sss
    g_qocn    = enthalpy_brine(Tbot)

#if defined flushing_notz
    g_Iswabs = c0
#else
    g_Iswabs = Iswabs
#endif       

    if (lsnow) then
       hsn = hslyr * real(nslyr,dbl_kind)

       g_hslyr   = hslyr
       g_hsn     = hsn
       g_hslyri  = c1 / hslyr
       g_hslyri2 = c1 / hslyr**2
       g_Sswabs  = Sswabs
    else
       hsn = c0

       g_hslyr   = c0
       g_hsn     = c0
       g_hslyri  = c0
       g_hslyri2 = c0
       g_Sswabs  = c0
    endif

    if (lsnow) then
       ! case with snow
       
       ! calculate the total initial ice and snow energy content
       call sum_quantity_snow(qin, qsn, hilyr, hslyr, einit)
       call sum_quantity_nosnow(Sin, hilyr, sinit)

       ! calculate the conductivities
       call conductivity_mush_array(qin0, Sin, km0)
       call conductivity_snow_array(ks0)

       call interface_conductivities(km0, kmstar0, nilyr)
       call interface_conductivities(ks0, ksstar0, nslyr)

       ! calculate quantities related to drainage
       call explicit_flow_velocities(Sin,    qin,    &
                                     Tin,    Tsf,    &
                                     Tbot,   g_q,    &
                                     g_dSdt, Sbr,    &
                                     qbr,    dt,     &
                                     sss,    g_qocn, &
                                     hilyr,  hin)
       !g_q = c0 !; g_dSdt = c0 ; q_dSdt = c0

#if diag == 1
       call diagnose_salt_flux(g_q, g_dSdt)
#endif

       ! calculate drivatives needed for scaling
       do k = 1, nslyr
          dqsn_dTsn(k) = dqdT_snow()
       enddo ! k

       do k = 1, nilyr
          dqin_dTin(k) = dqdT(Tin(k), Sin(k))
          dTin_dSin(k) = -dqdS(Tin(k), Sin(k)) / dqin_dTin(k)
          if (dTin_dSin(k) == c0) then
             dTin_dSin(k) = c1
          endif
       enddo ! k

       ! set the scaling variables
       yscale_cold(1)                             = c1
       yscale_cold(2:nslyr+1)                     = dqsn_dTsn(1:nslyr)
       yscale_cold(nslyr+2:nilyr+nslyr+1)         = dqin_dTin(1:nilyr)
       yscale_cold(nilyr+nslyr+2:2*nilyr+nslyr+1) = &
                                               c1 / dTin_dSin(1:nilyr)

       yscale_melt(1:nslyr)                       = dqsn_dTsn(1:nslyr)
       yscale_melt(nslyr+1:nilyr+nslyr)           = dqin_dTin(1:nilyr)
       yscale_melt(nilyr+nslyr+1:2*nilyr+nslyr)   = &
                                               c1 / dTin_dSin(1:nilyr)

       ! run the two stage solver
       call two_stage_solver_snow(Tsf,         Tsf0,       &
                                  qsn,         qsn0,       &
                                  qin,         qin0,       &
                                  Sin,         Sin0,       &
                                  Tsn,         Tsn0,       &
                                  Tin,         Tin0,       &
                                  g_q,         fswsfc,     &
                                  rhoa,        flw,        &
                                  potT,        Qa,         &
                                  shcoef,      lhcoef,     &
                                  hilyr,       hslyr,      &
                                  lstop,                   &
                                  yscale_cold, yscale_melt)

       ! given the updated enthalpy and bulk salinity calculate other quantities
       do k = 1, nslyr
          Tsn(k) = temperature_snow(qsn(k))
       enddo ! k

       do k = 1, nilyr
          Tin(k) = temperature_mush(qin(k), Sin(k))
          Sbr(k) = liquidus_brine_salinity_mush(Tin(k)) 
          qbr(k) = enthalpy_brine(Tin(k))
          phi(k) = liquid_fraction(Tin(k), Sin(k))
       enddo ! k
          
       ! downward cond flux at top surface (W m-2)
       fcondtopn = (c2 * ks0(1) * (Tsf - Tsn(1))) / hslyr

       ! calculate the total final ice and snow energy content
       call sum_quantity_snow(qin, qsn, hilyr, hslyr, efinal)
       call sum_quantity_nosnow(Sin, hilyr, sfinal)
          
    else
       ! case without snow

       ! calculate the total initial ice energy content
       call sum_quantity_nosnow(qin, hilyr, einit)
       call sum_quantity_nosnow(Sin, hilyr, sinit)

       ! calculate the conductivities
       call conductivity_mush_array(qin0, Sin, km0)
       call interface_conductivities(km0, kmstar0, nilyr)

       ! calculate quantities related to drainage
       call explicit_flow_velocities(Sin,    qin,    &
                                     Tin,    Tsf,    &
                                     Tbot,   g_q,    &
                                     g_dSdt, Sbr,    &
                                     qbr,    dt,     &
                                     sss,    g_qocn, &
                                     hilyr,  hin)
       !g_q = c0 !; g_dSdt = c0 ; q_dSdt = c0

#if diag == 1
       call diagnose_salt_flux(g_q, g_dSdt)
#endif

       ! calculate drivatives needed for scaling
       do k = 1, nilyr
          dqin_dTin(k) = dqdT(Tin(k), Sin(k))
          dTin_dSin(k) = -dqdS(Tin(k), Sin(k)) / dqin_dTin(k)
          if (dTin_dSin(k) == c0) dTin_dSin(k) = c1
       enddo ! k

       ! set the scaling variables
       yscale_cold(1)                 = c1
       yscale_cold(2:nilyr+1)         = dqin_dTin(1:nilyr)
       yscale_cold(nilyr+2:2*nilyr+1) = c1 / dTin_dSin(1:nilyr)
    
       yscale_melt(1:nilyr)           = dqin_dTin(1:nilyr)
       yscale_melt(nilyr+1:2*nilyr)   = c1 / dTin_dSin(1:nilyr)

       ! run the two stage solver
       call two_stage_solver_nosnow(Tsf,         Tsf0,       &
                                    qin,         qin0,       &
                                    Sin,         Sin0,       &
                                    Tin,         Tin0,       &
                                    g_q,         fswsfc,     &
                                    rhoa,        flw,        &
                                    potT,        Qa,         &
                                    shcoef,      lhcoef,     &
                                    hilyr,       lstop,      &
                                    yscale_cold, yscale_melt)

       ! given the updated enthalpy and bulk salinity calculate other quantities
       do k = 1, nilyr
          Tin(k) = temperature_mush(qin(k), Sin(k))
          Sbr(k) = liquidus_brine_salinity_mush(Tin(k)) 
          qbr(k) = enthalpy_brine(Tin(k))
          phi(k) = liquid_fraction(Tin(k), Sin(k))
       enddo ! k

       ! downward cond flux at top surface (W m-2)
       fcondtopn = (c2 * km0(1) * (Tsf - Tin(1))) / hilyr
       
       ! calculate the total final ice and snow energy content
       call sum_quantity_nosnow(qin, hilyr, efinal)
       call sum_quantity_nosnow(Sin, hilyr, sfinal)

    endif

    if (lstop) then
       return
    end if

    ! surface fluxes
    call surface_heat_flux(Tsf,     fswsfc, &
                           rhoa,    flw,    &
                           potT,    Qa,     &
                           shcoef,  lhcoef, &
                           flwoutn, fsensn, &
                           flatn,   fsurfn)

    ! downward cond flux at bottom surface (W m-2)
    fcondbot = (c2 * km0(nilyr) * (Tin(nilyr) - Tbot)) / hilyr

#if flushing == 1
    ! calculate vertical bulk darcy flow
    call flushing_velocity(Tin,   Sin,   &
                           hin,   hsn,   &
                           hilyr,        &
                           hpond, apond, & 
                           dt,    w)

    open(55,file='history/w.txt',position='append')
    write(55,*) istep, abs(max(w,c0)), abs(min(w,c0))
    close(55)
#else
    w = c0
#endif

    call check_inputs_outputs(Tsf,qsn,qin,Sin,Tsn,Tin,hilyr,"end:")

    ! perform an energy conservation check on temperature changes part
    call drainage_fluxes(g_q,      g_dSdt,  &
                         qbr,      Sbr,     &
                         hilyr,    dt,      &
                         fadvheat, fadvsalt)

    call flushing_fluxes(w, qbr, Sbr, fadvheat, fadvsalt)

    call check_conservation_heat(einit,    efinal,    &
                                 fcondbot, fcondtopn, &
                                 fsurfn,   fadvheat,  &
                                 Iswabs,   Sswabs,    &
                                 dt)

    call check_conservation_salt(sinit, sfinal, fadvsalt, dt)

#if flushing == 1
    ! now perform final flushing related changes to the ice
    ! drain ponds from flushing 
    call flush_pond(w, hin, hpond, apond, dt)
       
    ! flood snow ice
    call flood_ice(w,     dt,    &
                   hsn,   hin,   &
                   hslyr, hilyr, & 
                   qsn,   qin,   &
                   Sin,   Sbr,   &
                   qbr,   snoice)
#endif

  end subroutine temperature_changes_column

!=======================================================================

  subroutine check_conservation_heat(einit,    efinal,   &
                                     fcondbot, fcondtop, &
                                     fsurfn,   fadvheat, &
                                     Iswabs,   Sswabs,   &
                                     dt)

    ! check the conservation of energy during the temperatue change solution
    
    real(kind=dbl_kind), intent(in) :: &
         einit    , & ! initial total energy (J m-2)
         efinal   , & ! final total energy (J m-2)
         fcondbot , & ! downward cond flux at bottom surface (W m-2)
         fcondtop , & ! downward cond flux at top surface (W m-2)
         fsurfn   , & ! net flux to top surface, excluding fcondtopn (W m-2)
         fadvheat , & ! heat flux to ocean from brine advection (W m-2)
         dt           ! time step (s)

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Iswabs ! SW radiation absorbed in ice layers (W m-2)
    
    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         Sswabs ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind) :: &
         internal_energy_change    , & ! change in internal energy during temperature_changes (J m-2)
         internal_energy_absorbed  , & ! SW absorbed energy during temperature_changes (J m-2)
         surface_energy_absorbed   , & ! energy absorbed at surfaces (J m-2)
         advective_energy_absorbed , & ! energy absorbed from brine advection (J m-2)
         en_uncon                      ! energy unconservation  (J m-2)

    ! change in internal energy of the mush/snow
    internal_energy_change = efinal - einit
    
    ! internal energy absorbed in the mush/snow
    internal_energy_absorbed = (sum(Iswabs) + sum(Sswabs)) * dt

    ! energy absorbed at the surfaces
    surface_energy_absorbed = fcondtop * dt - fcondbot * dt

    ! energy absorbed through the advection of brine
    advective_energy_absorbed = - fadvheat * dt

    ! unconservation of energy in the mush/snow
    en_uncon = internal_energy_change - &
               internal_energy_absorbed - &
               surface_energy_absorbed - &
               advective_energy_absorbed

    !open(55,file='history/econ.txt',position='append')
    !write(55,*) istep1, internal_energy_change, internal_energy_absorbed, surface_energy_absorbed, &
    !     advective_energy_absorbed, en_uncon
    !close(55)

  end subroutine check_conservation_heat

!=======================================================================

  subroutine check_conservation_salt(sinit, sfinal, fadvsalt, dt)

    ! check the conservation of salt during the temperatue change solution

    real(kind=dbl_kind), intent(in) :: &
         sinit    , & ! initial total salt content (ppt m)
         sfinal   , & ! final total salt content (ppt m)
         fadvsalt , & ! salt flux to ocean from brine advection (ppt m s-1)
         dt           ! time step (s)

    real(kind=dbl_kind) :: &
         internal_salt_change    , & ! change in internal salt content during temperature_changes (ppt m)
         advective_salt_absorbed , & ! salt content increase from brine advection (ppt m)
         sa_uncon                    ! salt unconservation unconservation  (ppt m)
    
    ! change in salt content of the mush
    internal_salt_change = sfinal - sinit

    ! salt retained through the advection of brine
    advective_salt_absorbed = - fadvsalt * dt

    ! unconservation of salt in the mush
    sa_uncon = internal_salt_change - advective_salt_absorbed

    !open(55,file='history/scon.txt',position='append')
    !write(55,*) istep1, internal_salt_change, advective_salt_absorbed, sa_uncon
    !close(55)

  end subroutine check_conservation_salt

!=======================================================================

  subroutine sum_quantity_snow(qin, qsn, hilyr, hslyr, quantity)

    ! perform a vertical thickness sum of an intensive quantity for case with snow

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qin ! intensive quantity in ice layers
    
    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         qsn ! intensive quantity in snow layers
    
    real(kind=dbl_kind), intent(in) :: &
         hilyr , & ! ice layer thickness (m)
         hslyr     ! snow layer thickness (m)

    real(kind=dbl_kind), intent(out) :: &
         quantity ! intensive quantity summed over thickness

    integer(kind=int_kind) :: &
         k ! ice/snow layer index

    quantity = c0

    ! sum intensive quantity over the snow
    do k = 1, nslyr
       quantity = quantity + qsn(k) * hslyr
    enddo ! k
 
    ! sum intensive quantity over the mush
    do k = 1, nilyr
       quantity = quantity + qin(k) * hilyr
    enddo ! k

  end subroutine sum_quantity_snow

!=======================================================================

  subroutine sum_quantity_nosnow(qin, hilyr, quantity)

    ! perform a vertical thickness sum of an intensive quantity for case with no snow

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qin ! intensive quantity in ice layers
    
    real(kind=dbl_kind), intent(in) :: &
         hilyr ! ice layer thickness (m)
    
    real(kind=dbl_kind), intent(out) :: &
         quantity ! intensive quantity summed over thickness

    integer(kind=int_kind) :: &
         k ! ice layer index

    quantity = c0
 
    ! sum intensive quantity over the mush
    do k = 1, nilyr
       quantity = quantity + qin(k) * hilyr
    enddo ! k

  end subroutine sum_quantity_nosnow

!=======================================================================

  subroutine check_inputs_outputs(Tsf,qsn,qin,Sin,Tsn,Tin,hilyr,message)

    ! check the physical reasonability of the input or output vertical thermo quantities
    
    real(kind=dbl_kind), intent(in) :: &
         Tsf    , & ! ice/snow surface temperature (C)
         hilyr      ! ice layer thickness (m)

    real(kind=dbl_kind), dimension(nslyr), intent(in) :: &
         qsn , & ! snow layer enthalpy (J m-3)
         Tsn     ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         qin , & ! ice layer enthalpy (J m-3)
         Sin , & ! ice layer bulk salinity (ppt)
         Tin     ! ice layer temperature (C)

    character(len=*), intent(in) :: &
         message ! message to print to abort message

    integer(kind=int_kind) :: &
         k  , &  ! ice/snow layer index
         kk      ! other ice/snow layer index

    ! check for NaNs
    if (isnan(Tsf))      call abort_ice("temperature_changes_column: "//trim(message)//" Tsf Nan")
    if (isnan(sum(qsn))) call abort_ice("temperature_changes_column: "//trim(message)//" qsn Nan")
    if (isnan(sum(qin))) call abort_ice("temperature_changes_column: "//trim(message)//" qin Nan")
    if (isnan(sum(Sin))) call abort_ice("temperature_changes_column: "//trim(message)//" Sin Nan")
    if (isnan(sum(Tsn))) call abort_ice("temperature_changes_column: "//trim(message)//" Tsn Nan")
    if (isnan(sum(Tin))) call abort_ice("temperature_changes_column: "//trim(message)//" Tin Nan")
    if (isnan(hilyr))    call abort_ice("temperature_changes_column: "//trim(message)//" hilyr Nan")

    ! check negative salinities
    do k = 1, nilyr

       if (Sin(k) < c0) then

          do kk = 1, nslyr
             write(*,*) kk, Sin(kk)
          enddo ! kk

          call abort_ice("temperature_changes_column: "//trim(message)//" Sin negative")

       endif

    enddo ! k

    ! check positive snow temperatures
    !do k = 1, nslyr

    !   if (Tsn(k) > c0) then

    !      write(*,*) k, Tsn(k)

    !      call abort_ice("temperature_changes_column: "//trim(message)//" Tsn positive")
          
    !   endif

    !enddo ! k

  end subroutine check_inputs_outputs

!=======================================================================

  subroutine two_stage_solver_snow(Tsf,         Tsf0,       &
                                   qsn,         qsn0,       &
                                   qin,         qin0,       &
                                   Sin,         Sin0,       &
                                   Tsn,         Tsn0,       &
                                   Tin,         Tin0,       &
                                   q,           fswsfc,     &
                                   rhoa,        flw,        &
                                   potT,        Qa,         &
                                   shcoef,      lhcoef,     &
                                   hilyr,       hslyr,      &
                                   lstop,                   &
                                   yscale_cold, yscale_melt)

    ! solve the vertical temperature and salt change for case with snow
    ! 1) determine what type of surface condition existed previously - cold or melting
    ! 2) solve the system assuming this condition persists
    ! 3) check the consistency of the surface condition of the solution
    ! 4) If the surface condition is inconsistent resolve for the other surface condition
    ! 5) If neither solution is consistent the resolve the inconsistency

    real(kind=dbl_kind), intent(inout) :: &
         Tsf         , & ! snow surface temperature (C)
         Tsf0            ! snow surface temperature (C) at beginning of timestep

    real(kind=dbl_kind), dimension(1:nslyr), intent(inout) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         qsn0        , & ! snow layer enthalpy (J m-3) at beginning of timestep
         Tsn         , & ! snow layer temperature (C)
         Tsn0            ! snow layer temperature (C) at beginning of timestep

    real(kind=dbl_kind), dimension(1:nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         qin0        , & ! ice layer enthalpy (J m-3) at beginning of timestep
         Sin         , & ! ice layer bulk salinity (ppt)
         Sin0        , & ! ice layer bulk salinity (ppt) at beginning of timestep
         Tin         , & ! ice layer temperature (C)
         Tin0            ! ice layer temperature (C) at beginning of timestep

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

    logical(kind=log_kind), intent(inout) :: &
         lstop           ! solver failure flag

    real(kind=dbl_kind), dimension(nynmax), intent(in) :: &
         yscale_cold , & ! solution scaling vector for cold conditions
         yscale_melt     ! solution scaling vector for melting conditions

    real(kind=dbl_kind) :: &
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         Tsn1            ! upper snow layer temperature (C)

    real(kind=dbl_kind), dimension(nynmax) :: &
         yn_cold     , & ! solution vector for cold conditions
         y0_cold     , & ! initial solution vector for cold conditions
         yn_melt     , & ! solution vector for melting conditions
         y0_melt     , & ! initial solution vector for melting conditions
         fscale_cold , & ! residual scaling vector for cold conditions
         fscale_melt , & ! residual scaling vector for melting conditions
         Jac_a_cold  , & ! Jacobian lower off-diagonal for cold conditions
         Jac_b_cold  , & ! Jacobian diagonal for cold conditions
         Jac_c_cold  , & ! Jacobian upper off-diagonal for cold conditions
         Jac_a_melt  , & ! Jacobian lower off-diagonal for melting conditions
         Jac_b_melt  , & ! Jacobian diagonal for melting conditions
         Jac_c_melt      ! Jacobian upper off-diagonal for melting conditions

    ! determine if surface is initially cold or melting
    if (Tsf < c0) then

       ! initially cold

#if debug == 1
       write(*,*) "Two stage snow cold: 1"
#endif

       ! solve the system for cold and snow
       call solver_snow_cold(yn_cold,     y0_cold,     &
                             Tsf,         Sin,         &
                             qsn,         qin,         &
                             Tsn,         Tin,         &
                             q,           fswsfc,      &
                             rhoa,        flw,         &
                             potT,        Qa,          &
                             shcoef,      lhcoef,      &
                             hilyr,       hslyr,       &
                             yscale_cold, fscale_cold, &
                             Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                             lstop)

       ! halt if solver failed
       if (lstop) return

#if debug == 1
       write(*,*) "Two stage snow: 1", Tsf, Tsf < c0
#endif

       ! check if solution is consistent - surface should still be cold
       if (Tsf < dTsf_errcon) then

          ! solution is consistent - have solution so finish
          return

       else
          
          ! solution is inconsistent - surface is warmer than melting 
          ! resolve assuming surface is melting

#if debug == 1
          write(*,*) "Two stage snow melt: 2"
#endif

          ! reset the solution to initial values
          Tsf = Tsf0
          qsn = qsn0
          qin = qin0
          Sin = Sin0

          ! solve the system for melting and snow
          call solver_snow_melt(yn_melt,     y0_melt,     &
                                Sin,                      &
                                qsn,         qin,         &
                                Tsn,         Tin,         &
                                q,           fswsfc,      &
                                rhoa,        flw,         &
                                potT,        Qa,          &
                                shcoef,      lhcoef,      &
                                hilyr,       hslyr,       &
                                yscale_melt, fscale_melt, &
                                Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                lstop)

          ! halt if solver failed
          if (lstop) return

          ! calculate the surface condition
          Tsf = c0
          Tsn1 = temperature_snow(qsn(1))
          fcondtopn = (c2 * ks0(1) * (Tsf - Tsn1)) / hslyr

          call surface_heat_flux(Tsf,     fswsfc, &
                                 rhoa,    flw,    &
                                 potT,    Qa,     &
                                 shcoef,  lhcoef, &
                                 flwoutn, fsensn, &
                                 flatn,   fsurfn)

#if debug == 1
          write(*,*) "Two stage snow: 2", fsurfn, fcondtopn, fsurfn > fcondtopn
#endif

          ! check if solution is consistent - surface conductive heat flux should be less than incoming surface heat flux
          if (fcondtopn - fsurfn < ferrcon) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent
             ! failed to find a solution so need to refine solutions until consistency found

             call resolve_inconsistency_snow(Tsf,         Sin,                   &
                                             qsn,         qin,                   &
                                             yn_cold,     y0_cold,               &
                                             yn_melt,     y0_melt,               &
                                             yscale_cold, yscale_melt,           &
                                             fscale_cold, fscale_melt,           &
                                             Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                             Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                             fswsfc,                             &
                                             rhoa,       flw,                    &
                                             potT,       Qa,                     &
                                             shcoef,     lhcoef,                 &
                                             hilyr,      hslyr,                  &
                                             lstop)

          endif

       endif

    else

       ! initially melting

#if debug == 1
       write(*,*) "Two stage snow melt: 3"
#endif

       ! solve the system for melting and snow
       call solver_snow_melt(yn_melt,     y0_melt,     &
                             Sin,                      &
                             qsn,         qin,         &
                             Tsn,         Tin,         &
                             q,           fswsfc,      &
                             rhoa,        flw,         &
                             potT,        Qa,          &
                             shcoef,      lhcoef,      &
                             hilyr,       hslyr,       &
                             yscale_melt, fscale_melt, &
                             Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                             lstop)

       ! halt if solver failed
       if (lstop) return

       ! calculate the surface condition
       Tsf = c0
       Tsn1 = temperature_snow(qsn(1))
       fcondtopn = (c2 * ks0(1) * (Tsf - Tsn1)) / hslyr

       call surface_heat_flux(Tsf,     fswsfc, &
                              rhoa,    flw,    &
                              potT,    Qa,     &
                              shcoef,  lhcoef, &
                              flwoutn, fsensn, &
                              flatn,   fsurfn)
       
#if debug == 1
       write(*,*) "Two stage snow: 3", fsurfn, fcondtopn, fsurfn > fcondtopn
#endif
       
       ! check if solution is consistent - surface conductive heat flux should be less than incoming surface heat flux
       if (fcondtopn - fsurfn < ferrcon) then

          ! solution is consistent - have solution so finish
          return

       else

          ! solution is inconsistent - resolve assuming other surface condition
          ! assume surface is cold

#if debug == 1
          write(*,*) "Two stage snow cold: 4"
#endif

          ! reset the solution to initial values
          Tsf = Tsf0
          qsn = qsn0
          qin = qin0
          Sin = Sin0

          ! solve the system for cold and snow
          call solver_snow_cold(yn_cold,     y0_cold,     &
                                Tsf,         Sin,         &
                                qsn,         qin,         &
                                Tsn,         Tin,         &
                                q,           fswsfc,      &
                                rhoa,        flw,         &
                                potT,        Qa,          &
                                shcoef,      lhcoef,      &
                                hilyr,       hslyr,       &
                                yscale_cold, fscale_cold, &
                                Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                lstop)

          ! halt if solver failed
          if (lstop) return

#if debug == 1
          write(*,*) "Two stage snow: 4", Tsf, Tsf < c0
#endif

          ! check if solution is consistent - surface should be cold
          if (Tsf < dTsf_errcon) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent
             ! failed to find a solution so need to refine solutions until consistency found

             call resolve_inconsistency_snow(Tsf,         Sin,                   &
                                             qsn,         qin,                   &
                                             yn_cold,     y0_cold,               &
                                             yn_melt,     y0_melt,               &
                                             yscale_cold, yscale_melt,           &
                                             fscale_cold, fscale_melt,           &
                                             Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                             Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                             fswsfc,                             &
                                             rhoa,       flw,                    &
                                             potT,       Qa,                     &
                                             shcoef,     lhcoef,                 &
                                             hilyr,      hslyr,                  &
                                             lstop)

          endif
          
       endif

    endif

  end subroutine two_stage_solver_snow

!=======================================================================

  subroutine solver_snow_cold(yn_cold,     y0_cold,     &
                              Tsf,         Sin,         & 
                              qsn,         qin,         &
                              Tsn,         Tin,         &
                              q,           fswsfc,      &
                              rhoa,        flw,         &
                              potT,        Qa,          &
                              shcoef,      lhcoef,      &
                              hilyr,       hslyr,       &
                              yscale_cold, fscale_cold, &
                              Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                              lstop)
  
    ! solve the heat and temperature equations simultaneously for case with a cold snow covered surface

    real(kind=dbl_kind), intent(inout) :: &
         Tsf             ! snow surface temperature (C)

    real(kind=dbl_kind), dimension(1:nslyr), intent(inout) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         Tsn             ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(1:nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         Sin         , & ! ice layer bulk salinity (ppt)
         Tin             ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

    real(kind=dbl_kind), dimension(nynmax), intent(in) :: &
         yscale_cold     ! solution scaling vector

    real(kind=dbl_kind), dimension(nynmax), intent(out) :: &
         yn_cold     , & ! solution vector
         y0_cold     , & ! initial solution vector
         fscale_cold , & ! residual scaling vector
         Jac_a_cold  , & ! Jacobian lower off-diagonal
         Jac_b_cold  , & ! Jacobian diagonal
         Jac_c_cold      ! Jacobian upper off-diagonal

     logical(kind=log_kind), intent(inout) :: &
         lstop           ! solver failure flag

    integer(kind=int_kind) :: &
         l, &            ! solution vector index
         solve_type      ! 1: nosnow, cold condition
                         ! 2: nosnow, melt condition
                         ! 3: snow, cold condition
                         ! 4: snow, melt condition

    solve_type = 3

    ! calculte the jacobian 
    call calculate_jacobian_snow_cold(Tsf,    Tsn,    &
                                      Tin,    Sin,    &
                                      q,              &
                                      hilyr,  hslyr,  &
                                      dt,     fswsfc, &
                                      rhoa,   flw,    &
                                      potT,   Qa,     &
                                      shcoef, lhcoef, &
                                      Jac_a_cold, Jac_b_cold, Jac_c_cold)

    ! set the residual scaling vector by the jacobian main diagonal
    do l = 1, nyn_snow_cold
       fscale_cold(l) = c1 / (Jac_b_cold(l) * yscale_cold(l))
    enddo ! l
    
    ! rescale the jacobian given the values of the solution and residual scaling vectors
    call rescale_jacobian_main_diagonal(Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                        yscale_cold, fscale_cold, nyn_snow_cold)
    
    ! calculate the initial solution vector
    call physical_to_state_snow_cold(Tsf, qsn, qin, Sin, &
                                     yn_cold, y0_cold, yscale_cold)
    
    ! solve the system of equations to get a new solution vector
    call newton_loop(yn_cold(1:nyn_snow_cold), &
                      y0_cold(1:nyn_snow_cold), &
                      yscale_cold(1:nyn_snow_cold), &
                      fscale_cold(1:nyn_snow_cold), &
                      nyn_snow_cold, solve_type, lstop, &
                      Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                      residual_function_snow_cold, &
                      converged_snow_cold)

    ! convert to physical values from the final solution vector
    call state_to_physical_snow_cold(yn_cold, yscale_cold, &
                                     Tsf, qsn, qin, Sin)

  end subroutine solver_snow_cold

!=======================================================================

  subroutine solver_snow_melt(yn_melt,     y0_melt,     &
                              Sin,                      &
                              qsn,         qin,         &
                              Tsn,         Tin,         &
                              q,           fswsfc,      &
                              rhoa,        flw,         &
                              potT,        Qa,          &
                              shcoef,      lhcoef,      &
                              hilyr,       hslyr,       &
                              yscale_melt, fscale_melt, &
                              Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                              lstop)

    ! solve the heat and temperature equations simultaneously for case with a melting snow covered surface

    real(kind=dbl_kind), dimension(1:nslyr), intent(inout) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         Tsn             ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(1:nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         Sin         , & ! ice layer bulk salinity (ppt)
         Tin             ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

    real(kind=dbl_kind), dimension(nynmax), intent(in) :: &
         yscale_melt     ! solution scaling vector for melting conditions

    real(kind=dbl_kind), dimension(nynmax), intent(out) :: &
         yn_melt     , & ! solution vector
         y0_melt     , & ! initial solution vector
         fscale_melt , & ! residual scaling vector
         Jac_a_melt  , & ! Jacobian lower off-diagonal
         Jac_b_melt  , & ! Jacobian diagonal
         Jac_c_melt      ! Jacobian upper off-diagonal

    logical(kind=log_kind), intent(inout) :: &
         lstop           ! solver failure flag

    integer(kind=int_kind) :: &
         l, &            ! solution vector index
         solve_type      ! 1: nosnow, cold condition
                         ! 2: nosnow, melt condition
                         ! 3: snow, cold condition
                         ! 4: snow, melt condition

    solve_type = 4

    ! calculte the jacobian 
    call calculate_jacobian_snow_melt(Tsn,   Tin,   &
                                      Sin,   q,     &
                                      hilyr, hslyr, &
                                      dt,           &
                                      Jac_a_melt, Jac_b_melt, Jac_c_melt)

    ! set the residual scaling vector by the jacobian main diagonal
    do l = 1, nyn_snow_melt
       fscale_melt(l) = c1 / (Jac_b_melt(l) * yscale_melt(l))
    enddo ! l

    ! rescale the jacobian given the values of the solution and residual scaling vectors
    call rescale_jacobian_main_diagonal(Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                        yscale_melt, fscale_melt, nyn_snow_melt)

    ! calculate the initial solution vector
    call physical_to_state_snow_melt(qsn, qin, Sin, &
                                     yn_melt, y0_melt, yscale_melt)

    ! solve the system of equations to get a new solution vector
    call newton_loop(yn_melt(1:nyn_snow_melt), &
                     y0_melt(1:nyn_snow_melt), &
                     yscale_melt(1:nyn_snow_melt), &
                     fscale_melt(1:nyn_snow_melt), &
                     nyn_snow_melt, solve_type, lstop, &
                     Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                     residual_function_snow_melt, &
                     converged_snow_melt)

    ! convert to physical values from the final solution vector
    call state_to_physical_snow_melt(yn_melt, yscale_melt, &
                                     qsn, qin, Sin)

  end subroutine solver_snow_melt

!=======================================================================

  subroutine two_stage_solver_nosnow(Tsf,         Tsf0,       &
                                     qin,         qin0,       &
                                     Sin,         Sin0,       &
                                     Tin,         Tin0,       &
                                     q,           fswsfc,     &
                                     rhoa,        flw,        &
                                     potT,        Qa,         &
                                     shcoef,      lhcoef,     &
                                     hilyr,       lstop,      &
                                     yscale_cold, yscale_melt)
    
    ! solve the vertical temperature and salt change for case with no snow
    ! 1) determine what type of surface condition existed previously - cold or melting
    ! 2) solve the system assuming this condition persists
    ! 3) check the consistency of the surface condition of the solution
    ! 4) If the surface condition is inconsistent resolve for the other surface condition
    ! 5) If neither solution is consistent the resolve the inconsistency

    real(kind=dbl_kind), intent(inout) :: &
         Tsf         , & ! ice surface temperature (C)
         Tsf0            ! ice surface temperature (C) at beginning of timestep
    
    real(kind=dbl_kind), dimension(1:nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         qin0        , & ! ice layer enthalpy (J m-3) at beginning of timestep
         Sin         , & ! ice layer bulk salinity (ppt)
         Sin0        , & ! ice layer bulk salinity (ppt) at beginning of timestep
         Tin         , & ! ice layer temperature (C)
         Tin0            ! ice layer temperature (C) at beginning of timestep

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         hilyr           ! ice layer thickness (m)

    logical, intent(inout) :: &
         lstop           ! solver failure flag

    real(kind=dbl_kind), dimension(nynmax), intent(in) :: &
         yscale_cold , & ! solution scaling vector for cold conditions
         yscale_melt     ! solution scaling vector for melting conditions 

    real(kind=dbl_kind) :: &
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         Tmlt        , & ! upper ice layer melting temperature (C)
         Tin1            ! upper ice layer temperature (C)
    
    real(kind=dbl_kind), dimension(nynmax) :: &
         yn_cold     , & ! solution vector for cold conditions
         y0_cold     , & ! initial solution vector for cold conditions
         yn_melt     , & ! solution vector for melting conditions
         y0_melt     , & ! initial solution vector for melting conditions
         fscale_cold , & ! residual scaling vector for cold conditions
         fscale_melt , & ! residual scaling vector for melting conditions
         Jac_a_cold  , & ! Jacobian lower off-diagonal for cold conditions
         Jac_b_cold  , & ! Jacobian diagonal for cold conditions
         Jac_c_cold  , & ! Jacobian upper off-diagonal for cold conditions
         Jac_a_melt  , & ! Jacobian lower off-diagonal for melting conditions
         Jac_b_melt  , & ! Jacobian diagonal for melting conditions
         Jac_c_melt      ! Jacobian upper off-diagonal for melting conditions

    ! initial surface melting temperature
    Tmlt = liquidus_temperature_mush(Sin(1))

#if defined flushing_notz
    if (.true.) then
#else
    ! determine if surface is initially cold or melting
    if (Tsf < Tmlt) then
#endif
       ! initially cold

#if debug == 1
       write(*,*) "Two stage nosnow cold: 1"
#endif

       ! solve the system for cold and no snow
       call solver_nosnow_cold(yn_cold,     y0_cold,     &
                               Tsf,         Sin,         &
                               qin,         Tin,         &
                               q,           fswsfc,      &
                               rhoa,        flw,         &
                               potT,        Qa,          &
                               shcoef,      lhcoef,      &
                               hilyr,                    &
                               yscale_cold, fscale_cold, & 
                               Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                               lstop)

       ! halt if solver failed
       if (lstop) return

       ! surface melting temperature
       Tmlt = liquidus_temperature_mush(Sin(1))

#if debug == 1
       write(*,*) "Two stage nosnow: 1", Tsf, Tmlt, Tsf < Tmlt, Tsf - Tmlt
#endif

#if defined flushing_notz
       return
#endif

       ! check if solution is consistent - surface should still be cold
       if (Tsf < Tmlt + dTsf_errcon) then

          ! solution is consistent - have solution so finish
          return

       else

          ! solution is inconsistent - surface is warmer than melting 
          ! resolve assuming surface is melting

#if debug == 1
          write(*,*) "Two stage nosnow melt: 2: "
#endif

          ! reset the solution to initial values
          Tsf = Tsf0
          qin = qin0
          Sin = Sin0

          ! solve the system for melt and no snow
          call solver_nosnow_melt(yn_melt,     y0_melt,     &
                                  Sin,                      &
                                  qin,         Tin,         &
                                  q,           fswsfc,      &
                                  rhoa,        flw,         &
                                  potT,        Qa,          &
                                  shcoef,      lhcoef,      &
                                  hilyr,                    &
                                  yscale_melt, fscale_melt, & 
                                  Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                  lstop)

          ! halt if solver failed
          if (lstop) return

          ! calculate the surface condition
          Tsf = liquidus_temperature_mush(Sin(1))
          Tin1 = temperature_mush(qin(1), Sin(1))
          fcondtopn = (c2 * km0(1) * (Tsf - Tin1)) / hilyr

          call surface_heat_flux(Tsf,     fswsfc, &
                                 rhoa,    flw,    &
                                 potT,    Qa,     &
                                 shcoef,  lhcoef, &
                                 flwoutn, fsensn, &
                                 flatn,   fsurfn)

#if debug == 1
          write(*,*) "Two stage nosnow: 2", fsurfn, fcondtopn, fsurfn > fcondtopn, fsurfn - fcondtopn
#endif

          ! check if solution is consistent - surface conductive heat flux should be less than incoming surface heat flux
          if (fcondtopn - fsurfn < ferrcon) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent
             ! failed to find a solution so need to refine solutions until consistency found

             call resolve_inconsistency_nosnow(Tsf,         Sin,                   &
                                               qin,                                &
                                               yn_cold,     y0_cold,               &
                                               yn_melt,     y0_melt,               &
                                               yscale_cold, yscale_melt,           &
                                               fscale_cold, fscale_melt,           &
                                               Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                               Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                               fswsfc,                             &
                                               rhoa,        flw,                   &
                                               potT,        Qa,                    &
                                               shcoef,      lhcoef,                &
                                               hilyr,       lstop)


          endif

       endif

    else
       ! initially melting

#if debug == 1
       write(*,*) "Two stage nosnow melt: 3"
#endif

       ! solve the system for melt and no snow
       call solver_nosnow_melt(yn_melt,     y0_melt,     &
                               Sin,                      &
                               qin,         Tin,         &
                               q,           fswsfc,      &
                               rhoa,        flw,         & 
                               potT,        Qa,          &
                               shcoef,      lhcoef,      &
                               hilyr,                    &
                               yscale_melt, fscale_melt, &
                               Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                               lstop)

       ! halt if solver failed
       if (lstop) return

       ! calculate the surface condition
       Tsf = liquidus_temperature_mush(Sin(1))
       Tin1 = temperature_mush(qin(1), Sin(1))
       fcondtopn = (c2 * km0(1) * (Tsf - Tin1)) / hilyr

       call surface_heat_flux(Tsf,     fswsfc, &
                              rhoa,    flw,    &
                              potT,    Qa,     &
                              shcoef,  lhcoef, &
                              flwoutn, fsensn, &
                              flatn,   fsurfn)

#if debug == 1
       write(*,*) "Two stage nosnow: 3", fsurfn, fcondtopn, fsurfn > fcondtopn
#endif

       ! check if solution is consistent - surface conductive heat flux should be less than incoming surface heat flux
       if (fcondtopn - fsurfn < ferrcon) then

          ! solution is consistent - have solution so finish
          return

       else

          ! solution is inconsistent - resolve assuming other surface condition
          ! assume surface is cold

#if debug == 1
          write(*,*) "Two stage nosnow cold: 4"
#endif

          ! reset the solution to initial values
          Tsf = Tsf0
          qin = qin0
          Sin = Sin0

          ! solve the system for cold and no snow
          call solver_nosnow_cold(yn_cold,     y0_cold,     &
                                  Tsf,         Sin,         &
                                  qin,         Tin,         &
                                  q,           fswsfc,      &
                                  rhoa,        flw,         &
                                  potT,        Qa,          &
                                  shcoef,      lhcoef,      &
                                  hilyr,                    & 
                                  yscale_cold, fscale_cold, &
                                  Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                  lstop)

          ! halt if solver failed
          if (lstop) return

#if debug == 1
          write(*,*) "Two stage nosnow: 4", Tsf, Tsf < c0
#endif

          ! surface melting temperature
          Tmlt = liquidus_temperature_mush(Sin(1))

          ! check if solution is consistent - surface should be cold
          if (Tsf < Tmlt + dTsf_errcon) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent
             ! failed to find a solution so need to refine solutions until consistency found

             call resolve_inconsistency_nosnow(Tsf,         Sin,                   &
                                               qin,                                &
                                               yn_cold,     y0_cold,               &
                                               yn_melt,     y0_melt,               &
                                               yscale_cold, yscale_melt,           &
                                               fscale_cold, fscale_melt,           &
                                               Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                               Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                               fswsfc,                             &
                                               rhoa,        flw,                   &
                                               potT,        Qa,                    &
                                               shcoef,      lhcoef,                &
                                               hilyr,       lstop)
             
          endif
          
       endif

    endif

  end subroutine two_stage_solver_nosnow

!=======================================================================

  subroutine solver_nosnow_cold(yn_cold,     y0_cold,     &
                                Tsf,         Sin,         &
                                qin,         Tin,         &
                                q,           fswsfc,      &
                                rhoa,        flw,         &
                                potT,        Qa,          &
                                shcoef,      lhcoef,      &
                                hilyr,                    &
                                yscale_cold, fscale_cold, & 
                                Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                lstop)

    ! solve the heat and temperature equations simultaneously for case with a cold no snow surface

    real(kind=dbl_kind), intent(inout) :: &
         Tsf             ! ice surface temperature (C)

    real(kind=dbl_kind), dimension(1:nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         Sin         , & ! ice layer bulk salinity (ppt)
         Tin             ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

    real(kind=dbl_kind), dimension(nynmax), intent(in) :: &
         yscale_cold     ! solution scaling vector for cold conditions

    real(kind=dbl_kind), dimension(nynmax), intent(out) :: &
         yn_cold     , & ! solution vector
         y0_cold     , & ! initial solution vector
         fscale_cold , & ! residual scaling vector
         Jac_a_cold  , & ! Jacobian lower off-diagonal
         Jac_b_cold  , & ! Jacobian diagonal
         Jac_c_cold      ! Jacobian upper off-diagonal

    logical(kind=log_kind), intent(inout) :: &
         lstop           ! solver failure flag

    integer(kind=int_kind) :: &
         l           , & ! solution vector index
         solve_type      ! 1: nosnow, cold condition
                         ! 2: nosnow, melt condition
                         ! 3: snow, cold condition
                         ! 4: snow, melt condition

    solve_type = 1

    ! calculte the jacobian 
    call calculate_jacobian_nosnow_cold(Tsf,    Tin,    &
                                        Sin,            &
                                        q,      hilyr,  &
                                        dt,     fswsfc, &
                                        rhoa,   flw,    &
                                        potT,   Qa,     &
                                        shcoef, lhcoef, &
                                        Jac_a_cold, Jac_b_cold, Jac_c_cold)

    ! set the residual scaling vector by the jacobian main diagonal
    do l = 1, nyn_nosnow_cold
       fscale_cold(l) = c1 / (Jac_b_cold(l) * yscale_cold(l))
    enddo ! l
    
    ! rescale the jacobian given the values of the solution and residual scaling vectors
    call rescale_jacobian_main_diagonal(Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                        yscale_cold, fscale_cold, nyn_nosnow_cold)

    ! calculate the initial solution vector
    call physical_to_state_nosnow_cold(Tsf, qin, Sin, &
                                       yn_cold, y0_cold, yscale_cold)

    ! solve the system of equations to get a new solution vector
    call newton_loop(yn_cold(1:nyn_nosnow_cold), &
                     y0_cold(1:nyn_nosnow_cold), &
                     yscale_cold(1:nyn_nosnow_cold), &
                     fscale_cold(1:nyn_nosnow_cold), &
                     nyn_nosnow_cold, solve_type, lstop, &
                     Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                     residual_function_nosnow_cold, &
                     converged_nosnow_cold)
    
    ! convert to physical values from the final solution vector
    call state_to_physical_nosnow_cold(yn_cold, yscale_cold, &
                                       Tsf, qin, Sin)

  end subroutine solver_nosnow_cold

!=======================================================================

  subroutine solver_nosnow_melt(yn_melt,     y0_melt,     &
                                Sin,                      &
                                qin,         Tin,         &
                                q,           fswsfc,      &
                                rhoa,        flw,         &
                                potT,        Qa,          &
                                shcoef,      lhcoef,      &
                                hilyr,                    &
                                yscale_melt, fscale_melt, &
                                Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                lstop)

    ! solve the heat and temperature equations simultaneously for case with a melting no snow surface

    real(kind=dbl_kind), dimension(1:nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         Sin         , & ! ice layer bulk salinity (ppt)
         Tin             ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr       , & ! ice layer thickness (m)
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

    real(kind=dbl_kind), dimension(nynmax), intent(in) :: &
         yscale_melt     ! solution scaling vector for melting conditions

    real(kind=dbl_kind), dimension(nynmax), intent(out) :: &
         yn_melt     , & ! solution vector
         y0_melt     , & ! initial solution vector
         fscale_melt , & ! residual scaling vector
         Jac_a_melt  , & ! Jacobian lower off-diagonal
         Jac_b_melt  , & ! Jacobian diagonal
         Jac_c_melt      ! Jacobian upper off-diagonal

    logical(kind=log_kind), intent(inout) :: &
         lstop           ! solver failure flag

    integer(kind=int_kind) :: &
         l           , & ! solution vector index
         solve_type      ! 1: nosnow, cold condition
                         ! 2: nosnow, melt condition
                         ! 3: snow, cold condition
                         ! 4: snow, melt condition

    solve_type = 2

    ! calculte the jacobian
    call calculate_jacobian_nosnow_melt(Tin, Sin,   &
                                        q,   hilyr, &
                                        dt,         &
                                        Jac_a_melt, Jac_b_melt, Jac_c_melt)
    
    ! set the residual scaling vector by the jacobian main diagonal
    do l = 1, nyn_nosnow_melt
       fscale_melt(l) = c1 / (Jac_b_melt(l) * yscale_melt(l))
    enddo ! l
          
    ! rescale the jacobian given the values of the solution and residual scaling vectors
    call rescale_jacobian_main_diagonal(Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                        yscale_melt, fscale_melt, nyn_nosnow_melt)

    ! calculate the initial solution vector
    call physical_to_state_nosnow_melt(qin, Sin, &
                                       yn_melt, y0_melt, yscale_melt)
          
    ! solve the system of equations to get a new solution vector
    call newton_loop(yn_melt(1:nyn_nosnow_melt), &
                     y0_melt(1:nyn_nosnow_melt), &
                     yscale_melt(1:nyn_nosnow_melt), &
                     fscale_melt(1:nyn_nosnow_melt), &
                     nyn_nosnow_melt, solve_type, lstop, &
                     Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                     residual_function_nosnow_melt, &
                     converged_nosnow_melt)

    ! convert to physical values from the final solution vector
    call state_to_physical_nosnow_melt(yn_melt, yscale_melt, &
                                       qin, Sin)

  end subroutine solver_nosnow_melt

!=======================================================================

  subroutine resolve_inconsistency_snow(Tsf,         Sin,                   &
                                        qsn,         qin,                   &
                                        yn_cold,     y0_cold,               &
                                        yn_melt,     y0_melt,               &
                                        yscale_cold, yscale_melt,           &
                                        fscale_cold, fscale_melt,           &
                                        Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                        Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                        fswsfc,                             &
                                        rhoa,        flw,                   &
                                        potT,        Qa,                    &
                                        shcoef,      lhcoef,                &
                                        hilyr,       hslyr, &
                                        lstop)

    ! if neither surface condition for the two stage solver wuth snow is consistent then resolve that inconsistency.
    ! take both the solutions from the two stage solver and continue performing newton iterations on them
    ! both at the same time until one solution becomes consistent - use this solution

    real(kind=dbl_kind), intent(out) :: &
         Tsf                 ! snow surface temperature (C)

    real(kind=dbl_kind), dimension(nslyr), intent(out) :: &
         qsn                 ! snow layer enthalpy (J m-3) 

    real(kind=dbl_kind), dimension(nilyr), intent(out) :: &
         qin             , & ! ice layer enthalpy (J m-3) 
         Sin                 ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         yn_cold         , & ! solution vector for cold conditions
         yn_melt             ! solution vector for melting conditions

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y0_cold         , & ! initial solution vector for cold conditions
         y0_melt         , & ! initial solution vector for melting conditions
         yscale_cold     , & ! solution scaling vector for cold conditions
         yscale_melt     , & ! solution scaling vector for melting conditions
         fscale_cold     , & ! residual scaling vector for cold conditions
         fscale_melt     , & ! residual scaling vector for melting conditions
         Jac_a_cold      , & ! Jacobian lower off-diagonal for cold conditions
         Jac_b_cold      , & ! Jacobian diagonal for cold conditions
         Jac_c_cold      , & ! Jacobian upper off-diagonal for cold conditions
         Jac_a_melt      , & ! Jacobian lower off-diagonal for melting conditions
         Jac_b_melt      , & ! Jacobian diagonal for melting conditions
         Jac_c_melt          ! Jacobian upper off-diagonal for melting conditions

    real(kind=dbl_kind), intent(in) :: &
         hilyr           , & ! ice layer thickness (m)
         hslyr           , & ! snow layer thickness (m)
         fswsfc          , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa            , & ! air density (kg/m^3)
         flw             , & ! incoming longwave radiation (W/m^2)
         potT            , & ! air potential temperature (K)
         Qa              , & ! specific humidity (kg/kg)
         shcoef          , & ! transfer coefficient for sensible heat
         lhcoef              ! transfer coefficient for latent heat

    logical(kind=log_kind), intent(inout) :: &
         lstop               ! solver failure flag

    real(kind=dbl_kind) :: &
         Tsf_melt        , & ! snow surface temperature (C) for cold conditions
         Tsf_cold            ! snow surface temperature (C) for melting conditions

    real(kind=dbl_kind), dimension(nslyr) :: &
         qsn_cold        , & ! snow layer enthalpy (J m-3) for cold conditions
         qsn_melt            ! snow layer enthalpy (J m-3) for melting conditions

    real(kind=dbl_kind), dimension(nilyr) :: &
         qin_cold        , & ! ice layer enthalpy (J m-3) for cold conditions
         qin_melt        , & ! ice layer enthalpy (J m-3) for melting conditions
         Sin_cold        , & ! ice layer bulk salinity (ppt) for cold conditions
         Sin_melt            ! ice layer bulk salinity (ppt) for melting conditions

    real(kind=dbl_kind), dimension(nynmax) :: &
         delta_yn_cold   , & ! solution vector change from Newton iteration for cold conditions
         delta_yn_melt   , & ! solution vector change from Newton iteration for melting conditions
         rhs_cold        , & ! residual vector for cold conditions
         rhs_melt            ! residual vector for melting conditions

    real(kind=dbl_kind) :: &
         Tsn1            , & ! upper snow layer temperature (C)
         fcondtopn       , & ! downward cond flux at top surface (W m-2)
         flwoutn         , & ! upward LW at surface (W m-2)
         fsensn          , & ! surface downward sensible heat (W m-2)                        
         flatn           , & ! surface downward latent heat (W m-2)
         fsurfn              ! net flux to top surface, excluding fcondtopn

    integer(kind=int_kind) :: &
         nit             , & ! Newton iteration counter
         solve_type_cold , & ! solution type for snow, cold conditions
         solve_type_melt     ! solution type for snow, melting conditions

    solve_type_cold = 3
    solve_type_melt = 4
    
    ! iterate Newton solver until consistent solution found
    do nit = 1, nit_newton_max
       
       !----------------------------------------------------
       ! advance the cold solution first
       !----------------------------------------------------
 
       ! calculate the rhs to solve in GMRES
       call residual_function_snow_cold(yn_cold,         y0_cold,     &
                                        yscale_cold,     fscale_cold, & 
                                        rhs_cold,        nyn_snow_cold)
    

       ! solve GMRES for cold
       call gmres_loop(yn_cold,      y0_cold,             &
                      -rhs_cold,     delta_yn_cold,       &
                      yscale_cold,   fscale_cold,         &
                      nyn_snow_cold, solve_type_cold,     &
                      Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                      residual_function_snow_cold,        &
                      lstop)

       ! update the cold solution
       yn_cold = yn_cold + delta_yn_cold(1:nyn_snow_cold)

       ! calculate the physical state for the cold solution
       call state_to_physical_snow_cold(yn_cold, yscale_cold, &
                                        Tsf_cold, qsn_cold, qin_cold, Sin_cold)

       ! check for consistency of cold solution
       if (Tsf_cold < dTsf_errcon) then

          ! solution is consistent - record solution and return
          Tsf = Tsf_cold
          qsn = qsn_cold
          qin = qin_cold
          Sin = Sin_cold
          return

       endif

       !----------------------------------------------------
       ! advance the melt solution second
       !----------------------------------------------------

       ! calculate the rhs to solve in GMRES
       call residual_function_snow_melt(yn_melt,         y0_melt,     &
                                        yscale_melt,     fscale_melt, & 
                                        rhs_melt,        nyn_snow_melt)
       
       ! solve GMRES for melt
       call gmres_loop(yn_melt,      y0_melt,             &
                      -rhs_melt,     delta_yn_melt,       &
                      yscale_melt,   fscale_melt,         &
                      nyn_snow_melt, solve_type_melt,     &
                      Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                      residual_function_snow_melt,        &
                      lstop)

       ! update the melt solution
       yn_melt = yn_melt + delta_yn_melt(1:nyn_snow_melt)

       ! calculate the physical state for the melt solution
       call state_to_physical_snow_melt(yn_melt, yscale_melt, &
                                        qsn_melt, qin_melt, Sin_melt)

       ! check for consistency of melt solution
       Tsf_melt = c0
       Tsn1 = temperature_snow(qsn_melt(1))
       fcondtopn = (c2 * ks0(1) * (Tsf_melt - Tsn1)) / hslyr
       
       call surface_heat_flux(Tsf_melt, fswsfc, &
                              rhoa,     flw,    &
                              potT,     Qa,     &
                              shcoef,   lhcoef, &
                              flwoutn,  fsensn, &
                              flatn,    fsurfn)
       
       if (fcondtopn - fsurfn < ferrcon) then      
          
          ! solution is consistent - record solution and return
          Tsf = Tsf_melt
          qsn = qsn_melt
          qin = qin_melt
          Sin = Sin_melt
          return
          
       endif

    enddo ! nit

    write(nu_diag,*) "ice_therm_mushy: inconsistency failure: snow"
    lstop = .true.

  end subroutine resolve_inconsistency_snow

!=======================================================================

  subroutine resolve_inconsistency_nosnow(Tsf,         Sin,                   &
                                          qin,                                &
                                          yn_cold,     y0_cold,               &
                                          yn_melt,     y0_melt,               &
                                          yscale_cold, yscale_melt,           &
                                          fscale_cold, fscale_melt,           &
                                          Jac_a_cold, Jac_b_cold, Jac_c_cold, &
                                          Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                                          fswsfc,                             &
                                          rhoa,        flw,                   &
                                          potT,        Qa,                    &
                                          shcoef,      lhcoef,                &
                                          hilyr,       lstop)

    ! if neither surface condition for the two stage solver wuth snow is consistent then resolve that inconsistency.
    ! take both the solutions from the two stage solver and continue performing newton iterations on them
    ! both at the same time until one solution becomes consistent - use this solution

    real(kind=dbl_kind), intent(out) :: &
         Tsf                 ! ice surface temperature (C)

    real(kind=dbl_kind), dimension(nilyr), intent(out) :: &
         qin             , & ! ice layer enthalpy (J m-3) 
         Sin                 ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         yn_cold         , & ! solution vector for cold conditions
         yn_melt             ! solution vector for melting conditions

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y0_cold         , & ! initial solution vector for cold conditions
         y0_melt         , & ! initial solution vector for melting conditions
         yscale_cold     , & ! solution scaling vector for cold conditions
         yscale_melt     , & ! solution scaling vector for melting conditions                     
         fscale_cold     , & ! residual scaling vector for cold conditions
         fscale_melt     , & ! residual scaling vector for melting conditions
         Jac_a_cold      , & ! Jacobian lower off-diagonal for cold conditions
         Jac_b_cold      , & ! Jacobian diagonal for cold conditions
         Jac_c_cold      , & ! Jacobian upper off-diagonal for cold conditions
         Jac_a_melt      , & ! Jacobian lower off-diagonal for melting conditions
         Jac_b_melt      , & ! Jacobian diagonal for melting conditions
         Jac_c_melt          ! Jacobian upper off-diagonal for melting conditions

    real(kind=dbl_kind), intent(in) :: &
         hilyr           , & ! ice layer thickness (m)
         fswsfc          , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa            , & ! air density (kg/m^3)
         flw             , & ! incoming longwave radiation (W/m^2)
         potT            , & ! air potential temperature (K)
         Qa              , & ! specific humidity (kg/kg)
         shcoef          , & ! transfer coefficient for sensible heat
         lhcoef              ! transfer coefficient for latent heat

    logical, intent(inout) :: &
         lstop               ! solver failure flag

    real(kind=dbl_kind) :: &
         Tsf_melt        , & ! ice surface temperature (C) for cold conditions
         Tsf_cold            ! ice surface temperature (C) for melting conditions

    real(kind=dbl_kind), dimension(nilyr) :: &
         qin_cold        , & ! ice layer enthalpy (J m-3) for cold conditions
         qin_melt        , & ! ice layer enthalpy (J m-3) for melting conditions
         Sin_cold        , & ! ice layer bulk salinity (ppt) for cold conditions
         Sin_melt            ! ice layer bulk salinity (ppt) for melting conditions

    real(kind=dbl_kind), dimension(nynmax) :: &
         delta_yn_cold   , & ! solution vector change from Newton iteration for cold conditions
         delta_yn_melt   , & ! solution vector change from Newton iteration for melting conditions
         rhs_cold        , & ! residual vector for cold conditions
         rhs_melt            ! residual vector for melting conditions

    real(kind=dbl_kind) :: &
         Tmlt            , & ! upper ice layer melting temperature (C)
         Tin1            , & ! upper ice layer temperature (C)
         fcondtopn       , & ! downward cond flux at top surface (W m-2)
         flwoutn         , & ! upward LW at surface (W m-2)
         fsensn          , & ! surface downward sensible heat (W m-2)                          
         flatn           , & ! surface downward latent heat (W m-2)
         fsurfn              ! net flux to top surface, excluding fcondtopn

    integer(kind=int_kind) :: &
         nit , & ! Newton iteration counter
         solve_type_cold , & ! solver type for no snow, cold conditions
         solve_type_melt     ! solver type for no snow, melting conditions

    solve_type_cold = 1
    solve_type_melt = 2
    
    ! iterate Newton solver until consistent solution found
    do nit = 1, nit_newton_max
       
       !----------------------------------------------------
       ! advance the cold solution first
       !----------------------------------------------------
 
       ! calculate the rhs to solve in GMRES
       call residual_function_nosnow_cold(yn_cold,         y0_cold,     &
                                          yscale_cold,     fscale_cold, & 
                                          rhs_cold,        nyn_nosnow_cold)
    

       ! solve GMRES for cold
       call gmres_loop(yn_cold,        y0_cold,             &
                      -rhs_cold,       delta_yn_cold,       &
                      yscale_cold,     fscale_cold,         &
                      nyn_nosnow_cold, solve_type_cold,     &
                      Jac_a_cold, Jac_b_cold, Jac_c_cold,   &
                      residual_function_nosnow_cold,        &
                      lstop)

       ! update the cold solution
       yn_cold = yn_cold + delta_yn_cold(1:nyn_nosnow_cold)

       ! calculate the physical state for the cold solution
       call state_to_physical_nosnow_cold(yn_cold, yscale_cold, &
                                          Tsf_cold, qin_cold, Sin_cold)

       ! check for consistency of cold solution
       Tmlt = liquidus_temperature_mush(Sin_cold(1))

       if (Tsf_cold < Tmlt + dTsf_errcon) then
  
          ! solution is consistent - record solution and return
          Tsf = Tsf_cold
          qin = qin_cold
          Sin = Sin_cold
          return

       endif

       !----------------------------------------------------
       ! advance the melt solution second
       !----------------------------------------------------

       ! calculate the rhs to solve in GMRES
       call residual_function_nosnow_melt(yn_melt,         y0_melt,     &
                                          yscale_melt,     fscale_melt, & 
                                          rhs_melt,        nyn_nosnow_melt)
       
       ! solve GMRES for melt
       call gmres_loop(yn_melt,        y0_melt,           &
                      -rhs_melt,       delta_yn_melt,     &
                      yscale_melt,     fscale_melt,       &
                      nyn_nosnow_melt, solve_type_melt,   &
                      Jac_a_melt, Jac_b_melt, Jac_c_melt, &
                      residual_function_nosnow_melt,      &
                      lstop)

       ! update the melt solution
       yn_melt = yn_melt + delta_yn_melt(1:nyn_nosnow_melt)

       ! calculate the physical state for the melt solution
       call state_to_physical_nosnow_melt(yn_melt, yscale_melt, &
                                          qin_melt, Sin_melt)

       ! check for consistency of melt solution
       Tsf_melt = liquidus_temperature_mush(Sin_melt(1))
       Tin1 = temperature_mush(qin_melt(1), Sin_melt(1))
       fcondtopn = (c2 * km0(1) * (Tsf_melt - Tin1)) / hilyr
   
       call surface_heat_flux(Tsf_melt, fswsfc, &
                              rhoa,     flw,    &
                              potT,     Qa,     &
                              shcoef,   lhcoef, &
                              flwoutn,  fsensn, &
                              flatn,    fsurfn)
       
       if (fcondtopn - fsurfn < ferrcon) then
 
          ! solution is consistent - record solution and return
          Tsf = Tsf_melt
          qin = qin_melt
          Sin = Sin_melt
          return
          
       endif

    enddo ! nit

    write(nu_diag,*) "ice_therm_mushy: inconsistency failure: snow"
    lstop = .true.

  end subroutine resolve_inconsistency_nosnow

!=======================================================================
! Conversion between state and physical variables
!=======================================================================

  subroutine physical_to_state_snow_cold(Tsf, qsn, qin, Sin, &
                                         yn, y0, yscale)

    ! convert from physical to state variables for case with snow and cold surface

    real(kind=dbl_kind), intent(in) :: &
         Tsf     ! snow surface temperature (C)
    
    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         qsn     ! snow layer enthalpy (J m-3) 

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         qin , & ! ice layer enthalpy (J m-3) 
         Sin     ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         yn  , & ! solution vector
         y0      ! initial solution vector

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         yscale  ! solution scaling vector

    integer(kind=int_kind) :: &
         k   , & ! ice/snow layer index
         l1  , & ! solution vector index for temperature/enthalpy
         l2      ! solution vector index for bulk salinity

    l1 = 1
    yn(l1) = Tsf / yscale(l1)
    y0(l1) = Tsf / yscale(l1)

    do k = 1, nslyr
       l1 = k + 1
       yn(l1) = qsn(k) / yscale(l1)
       y0(l1) = qsn(k) / yscale(l1)
    enddo ! k
    
    do k = 1, nilyr

       l1 = k + nslyr + 1
       yn(l1) = qin(k) / yscale(l1)
       y0(l1) = qin(k) / yscale(l1)

       l2 = k + nilyr + nslyr + 1
       yn(l2) = Sin(k) / yscale(l2)
       y0(l2) = Sin(k) / yscale(l2)

    enddo ! k

  end subroutine physical_to_state_snow_cold

!=======================================================================

  subroutine physical_to_state_snow_melt(qsn, qin, Sin, &
                                         yn, y0, yscale)

    ! convert from physical to state variables for case with snow and melting surface

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         qsn     ! snow layer enthalpy (J m-3) 

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         qin , & ! ice layer enthalpy (J m-3) 
         Sin     ! ice layer bulk salinity (ppt)
    
    real(kind=dbl_kind), dimension(:), intent(out) :: &
         yn  , & ! solution vector
         y0      ! initial solution vector

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         yscale  ! solution scaling vector

    integer(kind=int_kind) :: &
         k   , & ! ice/snow layer index
         l1  , & ! solution vector index for enthalpy
         l2      ! solution vector index for bulk salinity

    do k = 1, nslyr
       l1 = k
       yn(l1) = qsn(k) / yscale(l1)
       y0(l1) = qsn(k) / yscale(l1)
    enddo ! k
    
    do k = 1, nilyr

       l1 = k + nslyr
       yn(l1) = qin(k) / yscale(l1)
       y0(l1) = qin(k) / yscale(l1)

       l2 = k + nilyr + nslyr
       yn(l2) = Sin(k) / yscale(l2)
       y0(l2) = Sin(k) / yscale(l2)

    enddo ! k

  end subroutine physical_to_state_snow_melt

!=======================================================================

  subroutine physical_to_state_nosnow_cold(Tsf, qin, Sin, &
                                           yn, y0, yscale)

    ! convert from physical to state variables for case without snow and cold surface

    real(kind=dbl_kind), intent(in) :: &
         Tsf     ! ice surface temperature (C)

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         qin , & ! ice layer enthalpy (J m-3) 
         Sin     ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         yn  , & ! solution vector
         y0      ! initial solution vector

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         yscale  ! solution scaling vector

    integer :: &
         k   , & ! ice layer index
         l1  , & ! solution vector index for temperature/enthalpy
         l2      ! solution vector index for bulk salinity

    l1 = 1
    yn(l1) = Tsf / yscale(l1)
    y0(l1) = Tsf / yscale(l1)   

    do k = 1, nilyr

       l1 = k + 1
       yn(l1) = qin(k) / yscale(l1)
       y0(l1) = qin(k) / yscale(l1)

       l2 = k + nilyr + 1
       yn(l2) = Sin(k) / yscale(l2)
       y0(l2) = Sin(k) / yscale(l2)

    enddo ! k
    
  end subroutine physical_to_state_nosnow_cold

!=======================================================================

  subroutine physical_to_state_nosnow_melt(qin, Sin, &
                                           yn, y0, yscale)

    ! convert from physical to state variables for case without snow and melting surface

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         qin , & ! ice layer enthalpy (J m-3) 
         Sin     ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         yn  , & ! solution vector
         y0      ! initial solution vector

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         yscale  ! solution scaling vector

    integer :: &
         k   , & ! ice layer index
         l1  , & ! solution vector index for enthalpy
         l2      ! solution vector index for bulk salinity

    do k = 1, nilyr

       l1 = k
       yn(l1) = qin(k) / yscale(l1)
       y0(l1) = qin(k) / yscale(l1)

       l2 = k + nilyr
       yn(l2) = Sin(k) / yscale(l2)
       y0(l2) = Sin(k) / yscale(l2)

    enddo ! k

  end subroutine physical_to_state_nosnow_melt

!=======================================================================

  subroutine state_to_physical_snow_cold(y, yscale, &
                                         Tsf, qsn, qin, Sin)

    ! convert from state to physical variables for case with snow and cold surface

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y      , & ! solution vector
         yscale     ! solution scaling vector

    real(kind=dbl_kind), intent(out) :: &
         Tsf        ! snow surface temperature (C)

    real(kind=dbl_kind), dimension(1:nslyr), intent(out) :: &
         qsn        ! snow layer enthalpy (J m-3) 

    real(kind=dbl_kind), dimension(1:nilyr), intent(out) :: &
         qin    , & ! ice layer enthalpy (J m-3) 
         Sin        ! ice layer bulk salinity (ppt)

    integer(kind=int_kind) :: &
         k      , & ! ice/snow layer index
         l1     , & ! solution vector index for temperature/enthalpy
         l2         ! solution vector index for bulk salinity

    l1 = 1
    Tsf = y(l1) * yscale(l1)

    do k = 1, nslyr
       l1 = k + 1
       qsn(k) = y(l1) * yscale(l1)
    enddo ! k

    do k = 1, nilyr

       l1 = k + nslyr + 1
       qin(k) = y(l1) * yscale(l1)

       l2 = k + nilyr + nslyr + 1
       Sin(k) = y(l2) * yscale(l2)

    enddo ! k

  end subroutine state_to_physical_snow_cold

!=======================================================================

  subroutine state_to_physical_snow_melt(y, yscale, &
                                         qsn, qin, Sin)

    ! convert from state to physical variables for case with snow and melting surface

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y      , & ! solution vector
         yscale     ! solution scaling vector

    real(kind=dbl_kind), dimension(1:nslyr), intent(out) :: &
         qsn        ! snow layer enthalpy (J m-3) 

    real(kind=dbl_kind), dimension(1:nilyr), intent(out) :: &
         qin    , & ! ice layer enthalpy (J m-3) 
         Sin        ! ice layer bulk salinity (ppt)

    integer(kind=int_kind) :: &
         k      , & ! ice/snow layer index
         l1     , & ! solution vector index for enthalpy
         l2         ! solution vector index for bulk salinity

    do k = 1, nslyr
       l1 = k
       qsn(k) = y(l1) * yscale(l1)
    enddo ! k

    do k = 1, nilyr

       l1 = k + nslyr
       qin(k) = y(l1) * yscale(l1)

       l2 = k + nilyr + nslyr
       Sin(k) = y(l2) * yscale(l2)

    enddo ! k

  end subroutine state_to_physical_snow_melt

!=======================================================================

  subroutine state_to_physical_nosnow_cold(y, yscale, &
                                           Tsf, qin, Sin)

    ! convert from state to physical variables for case without snow and cold surface

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y      , & ! solution vector
         yscale     ! solution scaling vector

    real(kind=dbl_kind), intent(out) :: &
         Tsf        ! ice surface temperature (C)

    real(kind=dbl_kind), dimension(1:nilyr), intent(out) :: &
         qin    , & ! ice layer enthalpy (J m-3) 
         Sin        ! ice layer bulk salinity (ppt)

    integer(kind=int_kind) :: &
         k      , & ! ice layer index
         l1     , & ! solution vector index for temperature/enthalpy
         l2         ! solution vector index for bulk salinity

    l1 = 1
    Tsf = y(l1) * yscale(l1)

    do k = 1, nilyr

       l1 = k + 1
       qin(k) = y(l1) * yscale(l1)

       l2 = k + nilyr + 1
       Sin(k) = y(l2) * yscale(l2)

    enddo ! k

  end subroutine state_to_physical_nosnow_cold

!=======================================================================

  subroutine state_to_physical_nosnow_melt(y, yscale, &
                                           qin, Sin)

    ! convert from state to physical variables for case with snow and melting surface

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y      , & ! solution vector
         yscale     ! solution scaling vector

    real(kind=dbl_kind), dimension(1:nilyr), intent(out) :: &
         qin    , & ! ice layer enthalpy (J m-3) 
         Sin        ! ice layer bulk salinity (ppt)

    integer(kind=int_kind) :: &
         k      , & ! ice layer index
         l1     , & ! solution vector index for enthalpy
         l2         ! solution vector index for bulk salinity

    do k = 1, nilyr

       l1 = k
       qin(k) = y(l1) * yscale(l1)

       l2 = k + nilyr
       Sin(k) = y(l2) * yscale(l2)

    enddo ! k

  end subroutine state_to_physical_nosnow_melt

!=======================================================================
! Residual functions
!=======================================================================

  subroutine residual_function_snow_cold(y,           y0,          &
                                         yscale_cold, fscale_cold, &
                                         f,           nyn)
    
    ! calculate the full residual vector for cold snow covered conditions
    
    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y           , & ! solution vector
         y0          , & ! initial solution vector
         yscale_cold , & ! solution scaling vector
         fscale_cold     ! residual scaling vector

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         f               ! residual vector

    integer(kind=int_kind), intent(in) :: &
         nyn             ! solution vector size

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         Sin         , & ! ice layer bulk salinity (ppt)
         phi         , & ! ice layer liquid fraction
         Sbr         , & ! ice layer brine salinity (ppt)
         km              ! ice layer conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         qsn         , & ! snow layer enthalpy (J m-3) 
         ks              ! snow layer conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(0:nilyr+1) :: &
         Tin             ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(0:nslyr+1) :: &
         Tsn             ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(1:nilyr+1) :: &
         kmstar          ! ice interface conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(1:nslyr+1) :: &
         ksstar          ! snow interface conductivity (W m-1 K-1)

    real(kind=dbl_kind) :: &
         Tsf         , & ! snow surface temperature (C)
         Tis         , & ! ice/snow interface temperature (C)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn          ! net flux to top surface, excluding fcondtopn

    integer(kind=int_kind) :: &
         l_Tsf       , & ! position in solution vector of surface condition
         l_qsn1      , & ! start position of snow enthalpy in solution vector
         l_qsnn      , & ! end position of snow enthalpy in solution vector
         l_Sin1      , & ! start position of bulk salinity in solution vector
         l_Sinn      , & ! end position of bulk salinity in solution vector
         l_qin1      , & ! start position of ice enthalpy in solution vector
         l_qinn          ! end position of ice enthalpy in solution vector

    ! unpack primary variables
    call state_to_physical_snow_cold(y, yscale_cold, &
                                     Tsf, qsn, qin, Sin)

    ! calculate secondary variables
    call mush_derived_quantities(qin, Sin, Tin(1:nilyr), Sbr, phi)
    call snow_derived_quantities(qsn, Tsn(1:nslyr))

    ! calculate snow/ice interface temperature
    call ice_snow_interface_temperature(Tis,                  &
                                        Tin(1),   Tsn(nslyr), &
                                        km0(1),   ks0(nslyr), &
                                        g_hilyri, g_hslyri)
    
    ! temperature boundary values
    Tsn(0)       = Tsf
    Tsn(nslyr+1) = Tis
    Tin(0)       = Tis
    Tin(nilyr+1) = g_Tbot

    ! surface flux
    l_Tsf = 1
    call surface_heat_flux(Tsf,      g_fswsfc, &
                           g_rhoa,   g_flw,    &
                           g_potT,   g_Qa,     &
                           g_shcoef, g_lhcoef, &
                           flwoutn,  fsensn,   &
                           flatn,    fsurfn)

    ! surface residual
    call residual_surface(Tsf,      fscale_cold(l_Tsf), &
                          f(l_Tsf), fsurfn,             &
                          ks0(1),   Tsn(1),             &
                          c0,       g_hslyri)

    ! heat diffusion in snow
    l_qsn1 = 2 ; l_qsnn = nslyr + 1
    call residual_heat_diffusion(qsn,              qsn0,                       &
                                 Tsn,              ksstar0,                    &
                                 g_Sswabs,         fscale_cold(l_qsn1:l_qsnn), &
                                 f(l_qsn1:l_qsnn), nslyr,                      &
                                 g_hslyri,         g_hslyri2,                  &
                                 g_dti)

    ! heat diffusion in mush
    l_qin1 = nslyr + 2 ; l_qinn = nslyr + nilyr + 1
    call residual_heat_diffusion(qin,              qin0,                       &
                                 Tin,              kmstar0,                    &
                                 g_Iswabs,         fscale_cold(l_qin1:l_qinn), &
                                 f(l_qin1:l_qinn), nilyr,                      &
                                 g_hilyri,         g_hilyri2,                  &
                                 g_dti)

    ! bulk salinity part of residual function
    l_Sin1 = nslyr + nilyr + 2 ; l_Sinn = nslyr + nilyr*2 + 1
    ! gravity drainage
    call salinity_residual(Sin,                        Sin0,                      &
                           Tin(1:nilyr),               g_dt,                      &
                           g_q,                        g_dSdt,                    &
                           f(l_Sin1:l_Sinn),           f(l_qin1:l_qinn),          &
                           fscale_cold(l_Sin1:l_Sinn), fscale_cold(l_qin1:l_qinn))

#if flushing == 1
    ! flushing
    call flushing_advection(Tin(1:nilyr),               Sin(1:nilyr),              &
                            g_hin,                      g_hsn,                     &
                            g_hilyr,                    g_hilyri,                  &
                            g_hpond,                    g_apond,                   &
                            g_sss,                      g_qocn,                    &
                            g_dt,                                                  &
                            f(l_Sin1:l_Sinn),           f(l_qin1:l_qinn),          &
                            fscale_cold(l_Sin1:l_Sinn), fscale_cold(l_qin1:l_qinn))
#endif
    !f(l_Sin1:l_Sinn) = c0

  end subroutine residual_function_snow_cold

!=======================================================================

  subroutine residual_function_snow_melt(y,           y0,          &
                                         yscale_melt, fscale_melt, &
                                         f,           nyn)
    
    ! calculate the full residual vector for melting snow covered conditions

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y           , & ! solution vector
         y0          , & ! initial solution vector
         yscale_melt , & ! solution scaling vector
         fscale_melt     ! residual scaling vector

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         f               ! residual vector

    integer(kind=int_kind), intent(in) :: &
         nyn             ! solution vector size

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         Sin         , & ! ice layer bulk salinity (ppt)
         phi         , & ! ice layer liquid fraction
         Sbr         , & ! ice layer brine salinity (ppt)
         km              ! ice layer conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         ks              ! snow layer conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(0:nilyr+1) :: &
         Tin             ! ice layer temperature (C)
    
    real(kind=dbl_kind), dimension(0:nslyr+1) :: &
         Tsn             ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(1:nilyr+1) :: &
         kmstar          ! ice interface conductivity (W m-1 K-1)
    
    real(kind=dbl_kind), dimension(1:nslyr+1) :: &
         ksstar          ! snow interface conductivity (W m-1 K-1) 

    real(kind=dbl_kind) :: &
         Tsf         , & ! snow surface temperature (C)
         Tis         , & ! ice/snow interface temperature (C)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn          ! net flux to top surface, excluding fcondtopn

    integer(kind=int_kind) :: &
         l_qsn1      , & ! start position of snow enthalpy in solution vector
         l_qsnn      , & ! end position of snow enthalpy in solution vector
         l_Sin1      , & ! start position of bulk salinity in solution vector 
         l_Sinn      , & ! end position of bulk salinity in solution vector 
         l_qin1      , & ! start position of ice enthalpy in solution vector 
         l_qinn          ! end position of ice enthalpy in solution vector 

    ! unpack primary variables
    call state_to_physical_snow_melt(y, yscale_melt, &
                                     qsn, qin, Sin)

    ! calculate secondary variables
    call mush_derived_quantities(qin, Sin, Tin(1:nilyr), Sbr, phi)
    call snow_derived_quantities(qsn, Tsn(1:nslyr))

    ! calculate snow/ice interface temperature
    call ice_snow_interface_temperature(Tis,                  &
                                        Tin(1),   Tsn(nslyr), &
                                        km0(1),   ks0(nslyr), &
                                        g_hilyri, g_hslyri)
    
    ! temperature boundary values
    Tsn(0)       = c0
    Tsn(nslyr+1) = Tis
    Tin(0)       = Tis
    Tin(nilyr+1) = g_Tbot

    ! heat diffusion in snow
    l_qsn1 = 1 ; l_qsnn = nslyr
    call residual_heat_diffusion(qsn,              qsn0,                       &
                                 Tsn,              ksstar0,                    &
                                 g_Sswabs,         fscale_melt(l_qsn1:l_qsnn), &
                                 f(l_qsn1:l_qsnn), nslyr,                      &
                                 g_hslyri,         g_hslyri2,                  &
                                 g_dti)

    ! heat diffusion in mush
    l_qin1 = nslyr + 1 ; l_qinn = nslyr + nilyr
    call residual_heat_diffusion(qin,              qin0,                       &
                                 Tin,              kmstar0,                    &
                                 g_Iswabs,         fscale_melt(l_qin1:l_qinn), &
                                 f(l_qin1:l_qinn), nilyr,                      &
                                 g_hilyri,         g_hilyri2,                  &
                                 g_dti)

    ! bulk salinity part of residual function
    l_Sin1 = nslyr + nilyr + 1 ; l_Sinn = nslyr + nilyr*2
    ! bulk salinity
    call salinity_residual(Sin,                        Sin0,                      &
                           Tin(1:nilyr),               g_dt,                      &
                           g_q,                        g_dSdt,                    &
                           f(l_Sin1:l_Sinn),           f(l_qin1:l_qinn),          &
                           fscale_melt(l_Sin1:l_Sinn), fscale_melt(l_qin1:l_qinn))

#if flushing == 1
    ! flushing
    call flushing_advection(Tin(1:nilyr),               Sin(1:nilyr),              &
                            g_hin,                      g_hsn,                     &
                            g_hilyr,                    g_hilyri,                  &
                            g_hpond,                    g_apond,                   &
                            g_sss,                      g_qocn,                    &
                            g_dt,                                                  &
                            f(l_Sin1:l_Sinn),           f(l_qin1:l_qinn),          &
                            fscale_melt(l_Sin1:l_Sinn), fscale_melt(l_qin1:l_qinn))
#endif
    !f(l_Sin1:l_Sinn) = c0

  end subroutine residual_function_snow_melt

!=======================================================================

  subroutine residual_function_nosnow_cold(y,           y0,          &
                                           yscale_cold, fscale_cold, &
                                           f,           nyn)
    
    ! calculate the full residual vector for cold no snow conditions

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y           , & ! solution vector
         y0          , & ! initial solution vector
         yscale_cold , & ! solution scaling vector
         fscale_cold     ! residual scaling vector

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         f               ! residual vector

    integer(kind=int_kind), intent(in) :: &
         nyn             ! solution vector size

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         Sin         , & ! ice layer bulk salinity (ppt)
         Sbr         , & ! ice layer brine salinity (ppt) 
         phi         , & ! ice layer liquid fraction
         km              ! ice layer conductivity (W m-1 K-1)
    
    real(kind=dbl_kind), dimension(0:nilyr+1) :: &
         Tin             ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(1:nilyr+1) :: &
         kmstar          ! ice interface conductivity (W m-1 K-1)

    real(kind=dbl_kind) :: &
         Tsf         , & ! ice surface temperature (C)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn          ! net flux to top surface, excluding fcondtopn

    integer(kind=int_kind) :: & 
         l_Tsf       , & ! position in solution vector of surface condition
         l_Sin1      , & ! start position of bulk salinity in solution vector 
         l_Sinn      , & ! end position of bulk salinity in solution vector 
         l_qin1      , & ! start position of ice enthalpy in solution vector 
         l_qinn          ! end position of ice enthalpy in solution vector 

    ! unpack primary variables
    call state_to_physical_nosnow_cold(y, yscale_cold, &
                                       Tsf, qin, Sin)

    ! calculate secondary variables
    call mush_derived_quantities(qin, Sin, Tin(1:nilyr), Sbr, phi)

    ! temperature boundary values
    Tin(0)       = Tsf
    Tin(nilyr+1) = g_Tbot

    ! surface flux
    l_Tsf = 1
    call surface_heat_flux(Tsf,      g_fswsfc, &
                           g_rhoa,   g_flw,    &
                           g_potT,   g_Qa,     &
                           g_shcoef, g_lhcoef, &
                           flwoutn,  fsensn,   &
                           flatn,    fsurfn)
    
    ! surface residual
    call residual_surface(Tsf,      fscale_cold(l_Tsf), &
                          f(l_Tsf), fsurfn,             &
                          km0(1),   Tin(1),             &
                          Sin(1),   g_hilyri)

    ! heat diffusion in mush
    l_qin1 = 2 ; l_qinn = nilyr + 1
    call residual_heat_diffusion(qin,              qin0,                       &
                                 Tin,              kmstar0,                    &
                                 g_Iswabs,         fscale_cold(l_qin1:l_qinn), &
                                 f(l_qin1:l_qinn), nilyr,                      &
                                 g_hilyri,         g_hilyri2,                  &
                                 g_dti)

    ! bulk salinity part of residual function
    l_Sin1 = nilyr + 2 ; l_Sinn = nilyr*2 + 1
    ! gravity drainage
    call salinity_residual(Sin,                  Sin0,                &
                           Tin(1:nilyr),         g_dt,                &
                           g_q,                  g_dSdt,              &
                           f(l_Sin1:l_Sinn),           f(l_qin1:l_qinn),          &
                           fscale_cold(l_Sin1:l_Sinn), fscale_cold(l_qin1:l_qinn))

#if flushing == 1
    ! flushing
    call flushing_advection(Tin(1:nilyr),               Sin(1:nilyr),              &
                            g_hin,                      g_hsn,                     &
                            g_hilyr,                    g_hilyri,                  &
                            g_hpond,                    g_apond,                   & 
                            g_sss,                      g_qocn,                    &
                            g_dt,                                                  &
                            f(l_Sin1:l_Sinn),           f(l_qin1:l_qinn),          &
                            fscale_cold(l_Sin1:l_Sinn), fscale_cold(l_qin1:l_qinn))
#endif
    !f(l_Sin1:l_Sinn) = c0

  end subroutine residual_function_nosnow_cold

!=======================================================================

  subroutine residual_function_nosnow_melt(y,           y0,          &
                                           yscale_melt, fscale_melt, &
                                           f,           nyn)

    ! calculate the full residual vector for melting no snow conditions
    
    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y           , & ! solution vector
         y0          , & ! initial solution vector
         yscale_melt , & ! solution scaling vector
         fscale_melt     ! residual scaling vector

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         f               ! residual vector

    integer(kind=int_kind), intent(in) :: &
         nyn             ! solution vector size

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         qin         , & ! ice layer enthalpy (J m-3) 
         Sin         , & ! ice layer bulk salinity (ppt)
         Sbr         , & ! ice layer brine salinity (ppt)
         phi         , & ! ice layer liquid fraction
         km              ! ice layer conductivity (W m-1 K-1)
    
    real(kind=dbl_kind), dimension(0:nilyr+1) :: &
         Tin             ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(1:nilyr+1) :: &
         kmstar          ! ice interface conductivity (W m-1 K-1)

    real(kind=dbl_kind) :: &
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn          ! net flux to top surface, excluding fcondtopn

    integer(kind=int_kind) :: &
         l_Sin1      , & ! start position of bulk salinity in solution vector 
         l_Sinn      , & ! end position of bulk salinity in solution vector 
         l_qin1      , & ! start position of ice enthalpy in solution vector 
         l_qinn          ! end position of ice enthalpy in solution vector 

    ! unpack primary variables
    call state_to_physical_nosnow_melt(y, yscale_melt, &
                                       qin, Sin)

    ! calculate secondary variables
    call mush_derived_quantities(qin, Sin, Tin(1:nilyr), Sbr, phi)

    ! temperature boundary values
    Tin(0)       = liquidus_temperature_mush(Sin(1))
    Tin(nilyr+1) = g_Tbot

    ! heat diffusion in mush
    l_qin1 = 1 ; l_qinn = nilyr
    call residual_heat_diffusion(qin,              qin0,                       &
                                 Tin,              kmstar0,                    &
                                 g_Iswabs,         fscale_melt(l_qin1:l_qinn), &
                                 f(l_qin1:l_qinn), nilyr,                      &
                                 g_hilyri,         g_hilyri2,                  &
                                 g_dti)

    ! bulk salinity part of residual function 
    l_Sin1 = nilyr + 1 ; l_Sinn = nilyr*2
    call salinity_residual(Sin,                        Sin0,                      &
                           Tin(1:nilyr),               g_dt,                      &
                           g_q,                        g_dSdt,                    &
                           f(l_Sin1:l_Sinn),           f(l_qin1:l_qinn),          &
                           fscale_melt(l_Sin1:l_Sinn), fscale_melt(l_qin1:l_qinn))

#if flushing == 1
    ! flushing
    call flushing_advection(Tin(1:nilyr),               Sin(1:nilyr),              &
                            g_hin,                      g_hsn,                     &
                            g_hilyr,                    g_hilyri,                  &
                            g_hpond,                    g_apond,                   &
                            g_sss,                      g_qocn,                    &
                            g_dt,                                                  &
                            f(l_Sin1:l_Sinn),           f(l_qin1:l_qinn),          &
                            fscale_melt(l_Sin1:l_Sinn), fscale_melt(l_qin1:l_qinn))
#endif
    !f(l_Sin1:l_Sinn) = c0

  end subroutine residual_function_nosnow_melt

!=======================================================================
! Full Jacobian
!=======================================================================

  subroutine calculate_jacobian_snow_cold(Tsf,    Tsn,    &
                                          Tin,    Sin,    &
                                          q,              &
                                          hilyr,  hslyr,  &
                                          dt,     fswsfc, &
                                          rhoa,   flw,    &
                                          potT,   Qa,     &
                                          shcoef, lhcoef, &
                                          Jac_a_cold, Jac_b_cold, Jac_c_cold)

    ! calculate the jacobian for cold snow covered conditions
    
    real(kind=dbl_kind), intent(in) :: &
         Tsf        , & ! snow surface temperature (C)
         hilyr      , & ! ice layer thickness (m)
         hslyr      , & ! snow layer thickness (m)
         dt         , & ! time step (s)
         fswsfc     , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa       , & ! air density (kg/m^3)
         flw        , & ! incoming longwave radiation (W/m^2)
         potT       , & ! air potential temperature (K)
         Qa         , & ! specific humidity (kg/kg)
         shcoef     , & ! transfer coefficient for sensible heat
         lhcoef         ! transfer coefficient for latent heat

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         Jac_a_cold , & ! Jacobian lower off-diagonal
         Jac_b_cold , & ! Jacobian diagonal
         Jac_c_cold     ! Jacobian upper off-diagonal

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         Tsn            ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         Tin        , & ! ice layer temperature (C)
         Sin            ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q              ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         dTsn_dqsn      ! derivative of snow temperature wrt snow enthalpy (K m3 J-1)
 
    real(kind=dbl_kind), dimension(1:nilyr) :: &
         dTin_dqin  , & ! derivative of ice temperature wrt ice enthalpy (K m3 J-1)
         dqin_dSin  , & ! derivative of ice enthalpy wrt bulk salinity (J m-3 ppt-1)
         dTin_dSin  , & ! derivative of ice temperature wrt bulk salinity (K ppt-1)
         dqbr_dqin  , & ! derivative of brine enthalpy wrt ice enthalpy (-)
         dqbr_dSin  , & ! derivative of brine enthalpy wrt bulk salinity (J m-3 ppt-1)
         dSbr_dqin  , & ! derivative of brine salinity wrt ice enthalpy (ppt m3 J-1)
         dSbr_dSin      ! derivative of brine salinity wrt bulk salinity (-)
 
    integer(kind=int_kind) :: &
         l_Tsf      , & ! position in solution vector of surface condition
         l_qsn1     , & ! start position of snow enthalpy in solution vector
         l_qsnn     , & ! end position of snow enthalpy in solution vector
         l_qin1     , & ! start position of ice enthalpy in solution vector
         l_qinn     , & ! end position of ice enthalpy in solution vector
         l_Sin1     , & ! start position of bulk salinity in solution vector
         l_Sinn         ! end position of bulk salinity in solution vector

    l_Tsf = 1
    l_qsn1 = 2
    l_qsnn = 1 + nslyr
    l_qin1 = 2 + nslyr
    l_qinn = 1 + nslyr + nilyr
    l_Sin1 = 2 + nslyr + nilyr
    l_Sinn = 1 + nslyr + nilyr*2

    ! calculate derivatives
    call calculate_derivatives_for_jacobian_snow(dTsn_dqsn)
    call calculate_derivatives_for_jacobian_ice(Tin,       Sin,       &
                                                dTin_dqin, dqin_dSin, &
                                                dTin_dSin, dqbr_dqin, &
                                                dqbr_dSin, dSbr_dqin, & 
                                                dSbr_dSin)

    ! initialize jacobian vectors
    call jacobian_initialize(Jac_a_cold, Jac_b_cold, Jac_c_cold)

    ! surface
    call jacobian_Tsf_snow(Tsf,    fswsfc,    &
                           rhoa,   flw,       &
                           potT,   Qa,        &
                           shcoef, lhcoef,    &
                           hslyr,  dTsn_dqsn, &
                           l_Tsf,  l_qsn1,    &
                           Jac_a_cold, Jac_b_cold, Jac_c_cold)
    
    ! heat conduction
    call jacobian_dfqsn_dqsn(dTsn_dqsn, hslyr,     &
                             dt,        l_qsn1,    &
                             Jac_a_cold, Jac_b_cold, Jac_c_cold)

    call jacobian_dfqin_dqin(dTin_dqin, hilyr,     &
                             dt,        l_qin1,    &
                             Jac_a_cold, Jac_b_cold, Jac_c_cold)

    call jacobian_qin_qsn_couple(Tsn,       Tin,        &
                                 dTsn_dqsn, dTin_dqin,  &
                                 dTin_dSin,             &
                                 hilyr,     hslyr,      &
                                 l_qsnn,    l_qin1,     & 
                                 Jac_a_cold, Jac_b_cold, Jac_c_cold)

    ! main salt diagonal
    call jacobian_dfSin_dSin(dt, l_Sin1, Jac_b_cold)

    ! salt advection
    call jacobian_salt_advection(dqbr_dqin, dqbr_dSin, &
                                 dSbr_dqin, dSbr_dSin, &
                                 hilyr,     q,         &
                                 l_qin1,    l_Sin1,    &
                                 Jac_a_cold, Jac_b_cold, Jac_c_cold)

  end subroutine calculate_jacobian_snow_cold

!=======================================================================

  subroutine calculate_jacobian_snow_melt(Tsn,   Tin,   &
                                          Sin,   q,     &
                                          hilyr, hslyr, &
                                          dt,           &
                                          Jac_a_melt, Jac_b_melt, Jac_c_melt)

    ! calculate the jacobian for melting snow covered conditions

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         Tsn            ! snow layer temperature (C)
         
    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         Tin        , & ! ice layer temperature (C)
         Sin            ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q              ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr      , & ! ice layer thickness (m)
         hslyr      , & ! snow layer thickness (m)
         dt             ! time step (s)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         Jac_a_melt , & ! Jacobian lower off-diagonal
         Jac_b_melt , & ! Jacobian diagonal
         Jac_c_melt     ! Jacobian upper off-diagonal

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         dTsn_dqsn      ! derivative of snow temperature wrt snow enthalpy (K m3 J-1)

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         dTin_dqin  , & ! derivative of ice temperature wrt ice enthalpy (K m3 J-1)
         dqin_dSin  , & ! derivative of ice enthalpy wrt bulk salinity (J m-3 ppt-1)
         dTin_dSin  , & ! derivative of ice temperature wrt bulk salinity (K ppt-1)
         dqbr_dqin  , & ! derivative of brine enthalpy wrt ice enthalpy (-)
         dqbr_dSin  , & ! derivative of brine enthalpy wrt bulk salinity (J m-3 ppt-1)
         dSbr_dqin  , & ! derivative of brine salinity wrt ice enthalpy (ppt m3 J-1)
         dSbr_dSin      ! derivative of brine salinity wrt bulk salinity (-)
 
    integer(kind=int_kind) :: &
         l_qsn1     , & ! start position of snow enthalpy in solution vector
         l_qsnn     , & ! end position of snow enthalpy in solution vector
         l_qin1     , & ! start position of ice enthalpy in solution vector
         l_qinn     , & ! end position of ice enthalpy in solution vector
         l_Sin1     , & ! start position of bulk salinity in solution vector
         l_Sinn         ! end position of bulk salinity in solution vector

    l_qsn1 = 1
    l_qsnn = nslyr
    l_qin1 = 1 + nslyr
    l_qinn = nslyr + nilyr
    l_Sin1 = 1 + nslyr + nilyr
    l_Sinn = nslyr + nilyr*2

    ! calculate derivatives
    call calculate_derivatives_for_jacobian_snow(dTsn_dqsn)
    call calculate_derivatives_for_jacobian_ice(Tin,       Sin,       &
                                                dTin_dqin, dqin_dSin, &
                                                dTin_dSin, dqbr_dqin, &
                                                dqbr_dSin, dSbr_dqin, & 
                                                dSbr_dSin)

    ! initialize jacobian vectors
    call jacobian_initialize(Jac_a_melt, Jac_b_melt, Jac_c_melt)
    
    ! heat conduction
    call jacobian_dfqsn_dqsn(dTsn_dqsn, hslyr,     &
                             dt,        l_qsn1,    &
                             Jac_a_melt, Jac_b_melt, Jac_c_melt)

    call jacobian_dfqin_dqin(dTin_dqin, hilyr,     &
                             dt,        l_qin1,    &
                             Jac_a_melt, Jac_b_melt, Jac_c_melt)

    call jacobian_qin_qsn_couple(Tsn,       Tin,        &
                                 dTsn_dqsn, dTin_dqin,  &
                                 dTin_dSin,             &
                                 hilyr,     hslyr,      &
                                 l_qsnn,    l_qin1,     &
                                 Jac_a_melt, Jac_b_melt, Jac_c_melt)

    ! main salt diagonal
    call jacobian_dfSin_dSin(dt, l_Sin1, Jac_b_melt)

    ! salt advection
    call jacobian_salt_advection(dqbr_dqin, dqbr_dSin, &
                                 dSbr_dqin, dSbr_dSin, &
                                 hilyr,     q,         &
                                 l_qin1,    l_Sin1,    &
                                 Jac_a_melt, Jac_b_melt, Jac_c_melt)

  end subroutine calculate_jacobian_snow_melt

!=======================================================================

  subroutine calculate_jacobian_nosnow_cold(Tsf,    Tin,    &
                                            Sin,            &
                                            q,      hilyr,  &
                                            dt,     fswsfc, &
                                            rhoa,   flw,    & 
                                            potT,   Qa,     &
                                            shcoef, lhcoef, &
                                            Jac_a_cold, Jac_b_cold, Jac_c_cold)

    ! calculate the jacobian for cold no snow conditions

    real(kind=dbl_kind), intent(in) :: &
         Tsf        , & ! ice surface temperature (C)
         hilyr      , & ! ice layer thickness (m)
         dt         , & ! time step (s)
         fswsfc     , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa       , & ! air density (kg/m^3)
         flw        , & ! incoming longwave radiation (W/m^2)
         potT       , & ! air potential temperature (K)
         Qa         , & ! specific humidity (kg/kg)
         shcoef     , & ! transfer coefficient for sensible heat
         lhcoef         ! transfer coefficient for latent heat

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         Tin        , & ! ice layer temperature (C)
         Sin            ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q              ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(:), intent(out) :: &  
         Jac_a_cold , & ! Jacobian lower off-diagonal
         Jac_b_cold , & ! Jacobian diagonal
         Jac_c_cold     ! Jacobian upper off-diagonal

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         dTin_dqin  , & ! derivative of ice temperature wrt ice enthalpy (K m3 J-1)
         dqin_dSin  , & ! derivative of ice enthalpy wrt bulk salinity (J m-3 ppt-1)
         dTin_dSin  , & ! derivative of ice temperature wrt bulk salinity (K ppt-1)
         dqbr_dqin  , & ! derivative of brine enthalpy wrt ice enthalpy (-)
         dqbr_dSin  , & ! derivative of brine enthalpy wrt bulk salinity (J m-3 ppt-1)
         dSbr_dqin  , & ! derivative of brine salinity wrt ice enthalpy (ppt m3 J-1)
         dSbr_dSin      ! derivative of brine salinity wrt bulk salinity (-)
 
    integer(kind=int_kind) :: &
         l_Tsf      , & ! position in solution vector of surface condition
         l_qin1     , & ! start position of ice enthalpy in solution vector
         l_qinn     , & ! end position of ice enthalpy in solution vector
         l_Sin1     , & ! start position of bulk salinity in solution vector
         l_Sinn         ! end position of bulk salinity in solution vector

    l_Tsf = 1
    l_qin1 = 2
    l_qinn = 1 + nilyr
    l_Sin1 = 2 + nilyr
    l_Sinn = 1 + nilyr*2

    ! calculate derivatives
    call calculate_derivatives_for_jacobian_ice(Tin,       Sin,       &
                                                dTin_dqin, dqin_dSin, &
                                                dTin_dSin, dqbr_dqin, &
                                                dqbr_dSin, dSbr_dqin, & 
                                                dSbr_dSin)

    ! initialize jacobian vectors
    call jacobian_initialize(Jac_a_cold, Jac_b_cold, Jac_c_cold)

    ! surface
    call jacobian_Tsf_nosnow(Tsf,                  &
                             hilyr,     fswsfc,    &
                             rhoa,      flw,       &
                             potT,      Qa,        &
                             shcoef,    lhcoef,    &
                             dTin_dqin, dTin_dSin, & 
                             l_Tsf,     l_qin1,    &
                             Jac_a_cold, Jac_b_cold, Jac_c_cold)
    
    ! heat conduction
    call jacobian_dfqin_dqin(dTin_dqin, hilyr,     &
                             dt,        l_qin1,    &
                             Jac_a_cold, Jac_b_cold, Jac_c_cold)

    ! main salt diagonal
    call jacobian_dfSin_dSin(dt, l_Sin1, Jac_b_cold)

    ! salt advection
    call jacobian_salt_advection(dqbr_dqin, dqbr_dSin, &
                                 dSbr_dqin, dSbr_dSin, &
                                 hilyr,     q,         &
                                 l_qin1,    l_Sin1,    &
                                 Jac_a_cold, Jac_b_cold, Jac_c_cold)

  end subroutine calculate_jacobian_nosnow_cold

!=======================================================================

  subroutine calculate_jacobian_nosnow_melt(Tin, Sin,   &
                                            q,   hilyr, &
                                            dt,         &
                                            Jac_a_melt, Jac_b_melt, Jac_c_melt)

    ! calculate the jacobian for melting no snow conditions

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         Tin        , & ! ice layer temperature (C)
         Sin            ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q              ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr      , & ! ice layer thickness (m)
         dt             ! time step (s)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         Jac_a_melt , & ! Jacobian lower off-diagonal
         Jac_b_melt , & ! Jacobian diagonal
         Jac_c_melt     ! Jacobian upper off-diagonal

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         dTin_dqin  , & ! derivative of ice temperature wrt ice enthalpy (K m3 J-1)
         dqin_dSin  , & ! derivative of ice enthalpy wrt bulk salinity (J m-3 ppt-1)
         dTin_dSin  , & ! derivative of ice temperature wrt bulk salinity (K ppt-1)
         dqbr_dqin  , & ! derivative of brine enthalpy wrt ice enthalpy (-)
         dqbr_dSin  , & ! derivative of brine enthalpy wrt bulk salinity (J m-3 ppt-1)
         dSbr_dqin  , & ! derivative of brine salinity wrt ice enthalpy (ppt m3 J-1)
         dSbr_dSin      ! derivative of brine salinity wrt bulk salinity (-)
 
    integer(kind=int_kind) :: &
         l_qin1     , & ! start position of ice enthalpy in solution vector
         l_qinn     , & ! end position of ice enthalpy in solution vector
         l_Sin1     , & ! start position of bulk salinity in solution vector
         l_Sinn         ! end position of bulk salinity in solution vector

    l_qin1 = 1
    l_qinn = nilyr
    l_Sin1 = 1 + nilyr
    l_Sinn = nilyr*2

    ! calculate derivatives
    call calculate_derivatives_for_jacobian_ice(Tin,       Sin,       &
                                                dTin_dqin, dqin_dSin, &
                                                dTin_dSin, dqbr_dqin, &
                                                dqbr_dSin, dSbr_dqin, &
                                                dSbr_dSin)

    ! initialize jacobian vectors
    call jacobian_initialize(Jac_a_melt, Jac_b_melt, Jac_c_melt)
    
    ! heat conduction
    call jacobian_dfqin_dqin(dTin_dqin, hilyr,     &
                             dt,        l_qin1,    &
                             Jac_a_melt, Jac_b_melt, Jac_c_melt)

    ! main salt diagonal
    call jacobian_dfSin_dSin(dt, l_Sin1, Jac_b_melt)

    ! salt advection
    call jacobian_salt_advection(dqbr_dqin, dqbr_dSin, &
                                 dSbr_dqin, dSbr_dSin, &
                                 hilyr,     q,         &
                                 l_qin1,    l_Sin1,    &
                                 Jac_a_melt, Jac_b_melt, Jac_c_melt)

  end subroutine calculate_jacobian_nosnow_melt

!=======================================================================

  subroutine calculate_derivatives_for_jacobian_snow(dTsn_dqsn)

    ! calculate snow derivatives needed for calculating the Jacobian

    real(kind=dbl_kind), dimension(nslyr), intent(out) :: &
         dTsn_dqsn ! derivative of snow temperature wrt snow enthalpy (K m3 J-1)

    integer(kind=int_kind) :: &
         k         ! snow layer index

    ! snow derivaties
    do k = 1, nslyr
       dTsn_dqsn(k) = c1 / dqdT_snow()
    enddo ! k

  end subroutine calculate_derivatives_for_jacobian_snow

!=======================================================================

  subroutine calculate_derivatives_for_jacobian_ice &
                                    (Tin,       Sin, &
                                     dTin_dqin, dqin_dSin, &
                                     dTin_dSin, dqbr_dqin, &
                                     dqbr_dSin, dSbr_dqin, &
                                     dSbr_dSin)

    ! calculate ice derivatives needed for calculating the Jacobian

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Tin       , & ! ice layer temperature (C)
         Sin           ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(nilyr), intent(out) :: &
         dTin_dqin , & ! derivative of ice temperature wrt ice enthalpy (K m3 J-1)
         dqin_dSin , & ! derivative of ice enthalpy wrt bulk salinity (J m-3 ppt-1)
         dTin_dSin , & ! derivative of ice temperature wrt bulk salinity (K ppt-1)
         dqbr_dqin , & ! derivative of brine enthalpy wrt ice enthalpy (-)
         dqbr_dSin , & ! derivative of brine enthalpy wrt bulk salinity (J m-3 ppt-1)
         dSbr_dqin , & ! derivative of brine salinity wrt ice enthalpy (ppt m3 J-1)
         dSbr_dSin     ! derivative of brine salinity wrt bulk salinity (-)

    ! derivatives needed for main terms
    ! dTsn_dqsn, dTin_dqin
    
    ! dTsn_dqsn = 1 / dqsn_dTsn
    ! dTin_dqin = 1 / dqin_dTin


    ! derivatives needed for other terms:
    ! dTin_dSin, dqbr_dqin, dqbr_dSin, dSbr_dqin, dSbr_dSin

    ! dTin_dSin = dTin_dqin * dqin_dSin
    ! dqbr_dqin = dqbr_dTin * dTin_dqin
    ! dqbr_dSin = dqbr_dTin * dTin_dSin = dqbr_dTin * dTin_dqin * dqin_dSin
    ! dSbr_dqin = dSbr_dTin * dTin_dqin
    ! dSbr_dSin = dSbr_dTin * dTin_dSin = dSbr_dTin * dTin_dqin * dqin_dSin

    real(kind=dbl_kind) :: &
         dqbr_dTin , & ! derivative of brine enthalpy wrt ice temperature (J m-3 K-1)
         dSbr_dTin     ! derivative of brine salinity wrt ice temperature (ppt K-1)

    integer(kind=int_kind) :: &
         k ! ice layer index

    ! ice derivatives
    do k = 1, nilyr

       dTin_dqin(k) = c1 / dqdT(Tin(k), Sin(k))

       dqin_dSin(k) = dqdS(Tin(k), Sin(k))

       dTin_dSin(k) = -dTin_dqin(k) * dqin_dSin(k)

       dqbr_dTin = denthalpy_brine_dT()

       dqbr_dqin(k) = dqbr_dTin * dTin_dqin(k)

       dqbr_dSin(k) = dqbr_dTin * dTin_dSin(k)

       dSbr_dTin = dliquidus_brine_salinity_mush_dT(Tin(k))

       dSbr_dqin(k) = dSbr_dTin * dTin_dqin(k)

       dSbr_dSin(k) = dSbr_dTin * dTin_dSin(k)

    enddo ! k

  end subroutine calculate_derivatives_for_jacobian_ice

!=======================================================================

  subroutine jacobian_initialize(Jac_a, Jac_b, Jac_c)

    ! initialize the elements of the Jacobian

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         Jac_a , & ! Jacobian lower off-diagonal
         Jac_b , & ! Jacobian diagonal
         Jac_c     ! Jacobian upper off-diagonal
    
    integer(kind=int_kind) :: &
         l         ! solution vector index

    do l = 1, nynmax

       Jac_a(l) = c0
       Jac_b(l) = c0
       Jac_c(l) = c0

    enddo ! l

  end subroutine jacobian_initialize

!=======================================================================

  subroutine jacobian_Tsf_snow(Tsf,    fswsfc,     &
                               rhoa,   flw,        &
                               potT,   Qa,         &
                               shcoef, lhcoef,     & 
                               hslyr,  dTsn_dqsn,  &
                               l_Tsf,  l_qsn1,     & 
                               Jac_a, Jac_b, Jac_c)

    ! calculate the elements of the Jacobian related to the surface condition for case with snow

    real(kind=dbl_kind), intent(in) :: &
         Tsf      , & ! snow surface temperature (C)
         hslyr    , & ! snow layer thickness (m)
         fswsfc   , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa     , & ! air density (kg/m^3)
         flw      , & ! incoming longwave radiation (W/m^2)
         potT     , & ! air potential temperature (K)
         Qa       , & ! specific humidity (kg/kg)
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         dTsn_dqsn    ! derivative of snow temperature wrt snow enthalpy (K m3 J-1)

    integer(kind=int_kind), intent(in) :: &
         l_Tsf    , & ! position in solution vector of surface condition
         l_qsn1       ! start position of snow enthalpy in solution vector

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         Jac_a    , & ! Jacobian lower off-diagonal
         Jac_b    , & ! Jacobian diagonal
         Jac_c        ! Jacobian upper off-diagonal

    real(kind=dbl_kind) :: &
         dfsurfn_dTsf ! derivative of surface flux wrt surface temperature (W m-2 K-1)

    ! dfsurfn_dTsf
    call dsurface_heat_flux_dTsf(Tsf,          fswsfc, &
                                 rhoa,         flw,    &
                                 potT,         Qa,     &
                                 shcoef,       lhcoef, &
                                 dfsurfn_dTsf)

    ! df(Tsf)_dTsf
    Jac_b(l_Tsf)  = Jac_b(l_Tsf)  + dfsurfn_dTsf - (c2 * ks0(1)) / hslyr

    ! df(Tsf)_dqsn(1)
    Jac_c(l_Tsf)  = Jac_c(l_Tsf)  + (c2 * ks0(1) * dTsn_dqsn(1)) / hslyr

    ! df(qsn(1))_dTsf
    Jac_a(l_qsn1) = Jac_a(l_qsn1) - (c2 * ks0(1)) / hslyr**2

  end subroutine jacobian_Tsf_snow

!=======================================================================

  subroutine jacobian_Tsf_nosnow(Tsf,                  &
                                 hilyr,     fswsfc,    &
                                 rhoa,      flw,       & 
                                 potT,      Qa,        &
                                 shcoef,    lhcoef,    &
                                 dTin_dqin, dTin_dSin, &
                                 l_Tsf,     l_qin1,    &
                                 Jac_a, Jac_b, Jac_c)

    ! calculate the elements of the Jacobian related to the surface condition for case without snow

    real(kind=dbl_kind), intent(in) :: &
         Tsf       , & ! ice surface temperature (C)
         hilyr     , & ! ice layer thickness (m)
         fswsfc    , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa      , & ! air density (kg/m^3)
         flw       , & ! incoming longwave radiation (W/m^2)
         potT      , & ! air potential temperature (K)
         Qa        , & ! specific humidity (kg/kg)
         shcoef    , & ! transfer coefficient for sensible heat
         lhcoef        ! transfer coefficient for latent heat

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         dTin_dqin , & ! derivative of ice temperature wrt ice enthalpy (K m3 J-1)
         dTin_dSin     ! derivative of ice temperature wrt bulk salinity (K ppt-1)

    integer(kind=int_kind), intent(in) :: &
         l_Tsf     , & ! position in solution vector of surface condition
         l_qin1        ! start position of ice enthalpy in solution vector

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         Jac_a     , & ! Jacobian lower off-diagonal
         Jac_b     , & ! Jacobian diagonal
         Jac_c         ! Jacobian upper off-diagonal

    real(kind=dbl_kind) :: &
         dfsurfn_dTsf  ! derivative of surface flux wrt surface temperature (W m-2 K-1)

    ! dfsurfn_dTsf
    call dsurface_heat_flux_dTsf(Tsf,          fswsfc, &
                                 rhoa,         flw,    &
                                 potT,         Qa,     &
                                 shcoef,       lhcoef, &
                                 dfsurfn_dTsf)

    ! df(Tsf)_dTsf
    Jac_b(l_Tsf)  = Jac_b(l_Tsf)  + dfsurfn_dTsf - (c2 * km0(1)) / hilyr

    ! df(Tsf)_dqin(1)
    Jac_c(l_Tsf)  = Jac_c(l_Tsf)  + (c2 * km0(1) * dTin_dqin(1)) / hilyr

    ! df(qin(1))_dTsf
    Jac_a(l_qin1) = Jac_a(l_qin1) - (c2 * km0(1)) / hilyr**2

  end subroutine jacobian_Tsf_nosnow

!=======================================================================

  subroutine jacobian_dfqsn_dqsn(dTsn_dqsn, hslyr,     &
                                 dt,        l_qsn1,    &
                                 Jac_a, Jac_b, Jac_c)

    ! calulate the d(f(qsn))/d(qsn) elements of the Jacobian 

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         dTsn_dqsn ! derivative of snow temperature wrt snow enthalpy (K m3 J-1)

    real(kind=dbl_kind), intent(in) :: &
         hslyr , & ! snow layer thickness (m)
         dt        ! time step (s)

    integer(kind=int_kind), intent(in) :: &
         l_qsn1    ! start position of snow enthalpy in solution vector

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         Jac_a , & ! Jacobian lower off-diagonal
         Jac_b , & ! Jacobian diagonal
         Jac_c     ! Jacobian upper off-diagonal

    integer(kind=int_kind) :: &
         k     , & ! snow layer index
         l         ! solution vector index
    
    ! diagonal
    do k = 1, nslyr

       l = l_qsn1 + k - 1
       Jac_b(l) = Jac_b(l) + (c1 / dt) + &
               (dTsn_dqsn(k) * (ksstar0(k+1) + ksstar0(k))) / hslyr**2

    enddo ! k

    ! lower off-diagonal
    do k = 2, nslyr

       l = l_qsn1 + k - 1
       Jac_a(l) = Jac_a(l) - (dTsn_dqsn(k-1) * ksstar0(k)) / hslyr**2

    enddo ! k

    ! upper-off diagonal
    do k = 1, nslyr-1

       l = l_qsn1 + k - 1
       Jac_c(l) = Jac_c(l) - (dTsn_dqsn(k+1) * ksstar0(k+1)) / hslyr**2

    enddo ! k

  end subroutine jacobian_dfqsn_dqsn

!=======================================================================

  subroutine jacobian_qin_qsn_couple(Tsn,       Tin,       & 
                                     dTsn_dqsn, dTin_dqin, &
                                     dTin_dSin,            &
                                     hilyr,     hslyr,     &
                                     l_qsnn,    l_qin1,    &
                                     Jac_a, Jac_b, Jac_c)

    ! calulate the d(f(qin))/d(qsn) and d(f(qsn))/d(qin) elements of the Jacobian 

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         Tsn       , & ! snow layer temperature (C)
         dTsn_dqsn     ! derivative of snow temperature wrt snow enthalpy (K m3 J-1)

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         Tin       , & ! ice layer temperature (C)
         dTin_dqin , & ! derivative of ice temperature wrt ice enthalpy (K m3 J-1)
         dTin_dSin     ! derivative of ice temperature wrt bulk salinity (K ppt-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr     , & ! ice layer thickness (m)
         hslyr         ! snow layer thickness (m)

    integer(kind=int_kind), intent(in) :: &
         l_qsnn    , & ! end position of snow enthalpy in solution vector
         l_qin1        ! start position of ice enthalpy in solution vector

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         Jac_a     , & ! Jacobian lower off-diagonal
         Jac_b     , & ! Jacobian diagonal
         Jac_c         ! Jacobian upper off-diagonal

    real(kind=dbl_kind) :: &
         dTisdTi   , & ! derivative of ice/snow interface temperature wrt upper ice layer temperature (-)
         dTisdTs       ! derivative of ice/snow interface temperature wrt lower snow layer temperature (-)

    call dTisdTiTs_fixedk(dTisdTi, dTisdTs,    &
                          Tin(1),  Tsn(nslyr), &
                          km0(1),  ks0(nslyr), &
                          hilyr,   hslyr)

    ! df(qsn(n))_dqsn(n)
    Jac_b(l_qsnn) = Jac_b(l_qsnn) - (ksstar0(nslyr+1) &
                                  * dTisdTs * dTsn_dqsn(nslyr)) / hslyr**2

    ! df(qin(1))_dqin(1)
    Jac_b(l_qin1) = Jac_b(l_qin1) - (kmstar0(1) &
                                  * dTisdTi * dTin_dqin(1)) / hilyr**2

    ! df(qsn(n))_dqin(1)
    Jac_c(l_qsnn) = Jac_c(l_qsnn) - (ksstar0(nslyr+1) &
                                  * dTisdTi * dTin_dqin(1)) / hslyr**2

    ! df(qin(1))_dqsn(n)
    Jac_a(l_qin1) = Jac_a(l_qin1) - (kmstar0(1) &
                                  * dTisdTs * dTsn_dqsn(nslyr)) / hilyr**2

  end subroutine jacobian_qin_qsn_couple

!=======================================================================

  subroutine jacobian_dfqin_dqin(dTin_dqin, hilyr,     &
                                 dt,        l_qin1,    &
                                 Jac_a, Jac_b, Jac_c)

    ! calulate the d(f(qin))/d(qin) elements of the Jacobian 

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         dTin_dqin ! derivative of ice temperature wrt ice enthalpy (K m3 J-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr , & ! ice layer thickness (m)
         dt        ! time step (s)

    integer(kind=int_kind), intent(in) :: &
         l_qin1    ! start position of ice enthalpy in solution vector

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         Jac_a , & ! Jacobian lower off-diagonal
         Jac_b , & ! Jacobian diagonal
         Jac_c     ! Jacobian upper off-diagonal

    integer(kind=int_kind) :: &
         k     , & ! ice layer index
         l         ! solution vector index
    
    ! diagonal
    do k = 1, nilyr

       l = l_qin1 + k - 1
       Jac_b(l) = Jac_b(l) + (c1 / dt) + &
               (dTin_dqin(k) * (kmstar0(k+1) + kmstar0(k))) / hilyr**2

    enddo ! k

    ! lower off-diagonal
    do k = 2, nilyr

       l = l_qin1 + k - 1
       Jac_a(l) = Jac_a(l) - (dTin_dqin(k-1) * kmstar0(k)) / hilyr**2

    enddo ! k

    ! upper-off diagonal
    do k = 1, nilyr-1

       l = l_qin1 + k - 1
       Jac_c(l) = Jac_c(l) - (dTin_dqin(k+1) * kmstar0(k+1)) / hilyr**2

    enddo ! k

  end subroutine jacobian_dfqin_dqin

!=======================================================================

  subroutine jacobian_dfSin_dSin(dt, l_Sin1, Jac_b)

    ! calulate the d(f(Sin))/d(Sin) elements of the Jacobian 

    real(kind=dbl_kind), intent(in) :: &
         dt     ! time step (s)

    integer(kind=int_kind), intent(in) :: &
         l_Sin1 ! start position of bulk salinity in solution vector

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         Jac_b  ! Jacobian diagonal

    integer(kind=int_kind) :: &
         k  , & ! ice layer index
         l      ! solution vector index

    do k = 1, nilyr

       l = l_Sin1 + k - 1
       Jac_b(l) = Jac_b(l) + (c1 / dt)

    enddo ! k

  end subroutine jacobian_dfSin_dSin

!=======================================================================

  subroutine jacobian_salt_advection(dqbr_dqin, dqbr_dSin, & 
                                     dSbr_dqin, dSbr_dSin, &
                                     hilyr,     q,         &
                                     l_qin1,    l_Sin1,    &
                                     Jac_a, Jac_b, Jac_c)

    ! calulate elements of the Jacobian related to salt advection
    
    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         dqbr_dqin , & ! derivative of brine enthalpy wrt ice enthalpy (-)
         dqbr_dSin , & ! derivative of brine enthalpy wrt bulk salinity (J m-3 ppt-1)
         dSbr_dqin , & ! derivative of brine salinity wrt ice enthalpy (ppt m3 J-1)
         dSbr_dSin     ! derivative of brine salinity wrt bulk salinity (-)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q             ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr         ! ice layer thickness (m)

    integer(kind=int_kind), intent(in) :: &
         l_qin1    , & ! start position of ice enthalpy in solution vector
         l_Sin1        ! start position of bulk salinity in solution vector

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         Jac_a     , & ! Jacobian lower off-diagonal
         Jac_b     , & ! Jacobian diagonal
         Jac_c         ! Jacobian upper off-diagonal
    
    integer(kind=int_kind) :: &
         k         , & ! ice layer index
         l             ! solution vector index

    ! df(qin)_dqin - through qbr in adv_q
    ! diagonal
    do k = 1, nilyr
       l = l_qin1 + k - 1
       Jac_b(l) = Jac_b(l) + (q(k) * dqbr_dqin(k)) / hilyr
    enddo ! k

    ! upper off diagonal (no lower off diagonal since upwind)
    do k = 1, nilyr-1
       l = l_qin1 + k - 1
       Jac_c(l) = Jac_c(l) - (q(k) * dqbr_dqin(k+1)) / hilyr
    enddo ! k


    ! df(Sin)_dSin - through Sbr in adv_s
    ! diagonal
    do k = 1, nilyr
       l = l_Sin1 + k - 1
       Jac_b(l) = Jac_b(l) + (q(k) * dSbr_dSin(k)) / hilyr
    enddo ! k

    ! upper off diagonal (no lower off diagonal since upwind)
    do k = 1, nilyr-1
       l = l_Sin1 + k - 1
       Jac_c(l) = Jac_c(l) - (q(k) * dSbr_dSin(k+1)) / hilyr
    enddo ! k

  end subroutine jacobian_salt_advection

!=======================================================================

  subroutine rescale_jacobian_main_diagonal(Jac_a, Jac_b, Jac_c, &
                                            yscale, fscale, nyn)
    
    ! rescale the elements of the Jacobian by the residual and state scaling vectors.

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         Jac_a  , & ! Jacobian lower off-diagonal
         Jac_b  , & ! Jacobian diagonal
         Jac_c      ! Jacobian upper off-diagonal

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         yscale , & ! solution scaling vector
         fscale     ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn        ! solution vector size

    integer(kind=int_kind) :: &
         l          ! solution vector index

    l = 1
    Jac_b(l) = Jac_b(l) * fscale(l) * yscale(l)
    Jac_c(l) = Jac_c(l) * fscale(l) * yscale(l+1)
    
    do l = 2, nyn-1
       Jac_a(l) = Jac_a(l) * fscale(l) * yscale(l-1)
       Jac_b(l) = Jac_b(l) * fscale(l) * yscale(l)
       Jac_c(l) = Jac_c(l) * fscale(l) * yscale(l+1)
    enddo ! l
    
    l = nyn
    Jac_a(l) = Jac_a(l) * fscale(l) * yscale(l-1)
    Jac_b(l) = Jac_b(l) * fscale(l) * yscale(l)
    
  end subroutine rescale_jacobian_main_diagonal

!=======================================================================
! Convergence tests
!=======================================================================

  function converged_snow_cold(fk,          f0,          &
                               yn,          y0,          &
                               yscale_cold, fscale_cold, &
                               nyn,         nit_newton)  &
                               result(lconverged)

    ! check for convergence of the Newton solver for cold snow covered conditions

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         fk              , & ! residual vector
         f0              , & ! initial residual vector
         yn              , & ! solution vector
         y0              , & ! initial solution vector
         yscale_cold     , & ! solution scaling vector
         fscale_cold         ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn             , & ! solution vector size
         nit_newton          ! Newton iteration counter
    
    logical(kind=log_kind) :: &
         lconverged          ! solver convergence test

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         qsn                 ! snow layer enthalpy (J m-3)
    
    real(kind=dbl_kind), dimension(1:nilyr) :: &
         qin             , & ! ice layer enthalpy (J m-3) 
         Sin                 ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         Tsf             , & ! snow surface temperature (C)
         en_uncon        , & ! solution energy unconservation (J m-2)
         flux_error      , & ! solution flux error (W m-2)
         sa_uncon        , & ! solution salt unconservation (ppt m)
         salt_flux_error , & ! solution salt flux unconservation (ppt m s-1)
         dTsf            , & ! estimated surface temperature error (C)
         surf_flux_error     ! solution surface flux error (W m-2)

    ! calculate the physical state from state vector
    call state_to_physical_snow_cold(yn, yscale_cold, &
                                     Tsf, qsn, qin, Sin)

    ! internal heat conservation error
    call energy_unconservation_cold(fk,       fscale_cold, &
                                    nyn,                   &
                                    g_hslyr,  g_hilyr,     &
                                    g_dt,     g_dti,       &
                                    .true.,                &
                                    en_uncon, flux_error)

    ! salt conservation error
    call salt_unconservation_cold(fk,       fscale_cold,    &
                                  nyn,                      &
                                  g_hilyr,                  &
                                  g_dt,     g_dti,          &
                                  .true.,                   &
                                  sa_uncon, salt_flux_error)

    ! estimated surface temperature error
    call surface_temperature_error(fk(1),                    &
                                   Tsf,      ks0(1),         & 
                                   nslyr,    g_hslyri,       &
                                   g_shcoef, g_lhcoef,       &
                                   g_rhoa,   fscale_cold(1), &
                                   dTsf)

    ! surface heat flux error
    call surface_flux_error(fk(1), fscale_cold(1), surf_flux_error)

    ! convergence test
    lconverged = (flux_error + surf_flux_error < ferrmax .and. &
                  salt_flux_error              < serrmax .and. &
                  dTsf                         < dTsf_errmax)

    !write(*,*) istep1, 1, nit_newton, lconverged, flux_error + surf_flux_error, salt_flux_error, dTsf

  end function converged_snow_cold

!=======================================================================

  function converged_snow_melt(fk,          f0,          &
                               yn,          y0,          &
                               yscale_melt, fscale_melt, &
                               nyn,         nit_newton)  &
                               result(lconverged)

    ! check for convergence of the Newton solver for melting snow covered conditions

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         fk              , & ! residual vector
         f0              , & ! initial residual vector
         yn              , & ! solution vector
         y0              , & ! initial solution vector
         yscale_melt     , & ! solution scaling vector
         fscale_melt         ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn             , & ! solution vector size
         nit_newton          ! Newton iteration counter
    
    logical(kind=log_kind) :: &
         lconverged          ! solver convergence test

    real(kind=dbl_kind) :: &
         en_uncon        , & ! solution energy unconservation (J m-2)
         flux_error      , & ! solution flux error (W m-2)
         sa_uncon        , & ! solution salt unconservation (ppt m)
         salt_flux_error     ! solution salt flux unconservation (ppt m s-1)

    ! internal heat conservation error
    call energy_unconservation_melt(fk,       fscale_melt, &
                                    nyn,                   &
                                    g_hslyr,  g_hilyr,     &
                                    g_dt,     g_dti,       &
                                    .true.,                &
                                    en_uncon, flux_error)

    ! salt conservation error
    call salt_unconservation_melt(fk,       fscale_melt,    & 
                                  nyn,      g_hilyr,        &
                                  g_dt,     g_dti,          &
                                  .true.,                   &
                                  sa_uncon, salt_flux_error)

    ! convergence test
    lconverged = (flux_error      < ferrmax .and. &
                  salt_flux_error < serrmax)

    !write(*,*) istep1, 2, nit_newton, lconverged, flux_error, salt_flux_error

  end function converged_snow_melt

!=======================================================================

  function converged_nosnow_cold(fk,          f0,          &
                                 yn,          y0,          &
                                 yscale_cold, fscale_cold, &
                                 nyn,         nit_newton)  &
                                 result(lconverged)
    
    ! check for convergence of the Newton solver for cold no snow conditions

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         fk              , & ! residual vector
         f0              , & ! initial residual vector
         yn              , & ! solution vector
         y0              , & ! initial solution vector
         yscale_cold     , & ! solution scaling vector
         fscale_cold         ! residual scaling vector
    
    integer(kind=int_kind), intent(in) :: &
         nyn             , & ! solution vector size
         nit_newton          ! Newton iteration counter
    
    logical(kind=log_kind) :: &
         lconverged          ! solver convergence test

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         qin             , & ! ice layer enthalpy (J m-3) 
         Sin                 ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         Tsf             , & ! ice surface temperature (C)
         en_uncon        , & ! solution energy unconservation (J m-2)
         flux_error      , & ! solution flux error (W m-2)
         sa_uncon        , & ! solution salt unconservation (ppt m)
         salt_flux_error , & ! solution salt flux unconservation (ppt m s-1)
         dTsf            , & ! estimated surface temperature error (C)
         surf_flux_error     ! solution surface flux error (W m-2)

    ! calculate the physical state from state vector
    call state_to_physical_nosnow_cold(yn, yscale_cold, &
                                       Tsf, qin, Sin)

    ! internal heat conservation error
    call energy_unconservation_cold(fk,       fscale_cold , &
                                    nyn,                    &
                                    g_hslyr,  g_hilyr,      &
                                    g_dt,     g_dti,        &
                                    .false.,                &
                                    en_uncon, flux_error)
    
    ! salt conservation error
    call salt_unconservation_cold(fk,       fscale_cold,    &
                                  nyn,      g_hilyr,        &
                                  g_dt,     g_dti,          &
                                  .false.,                  &
                                  sa_uncon, salt_flux_error)

    ! estimated surface temperature error
    call surface_temperature_error(fk(1),                    &
                                   Tsf,      km0(1),         &
                                   nilyr,    g_hilyri,       &
                                   g_shcoef, g_lhcoef,       &
                                   g_rhoa,   fscale_cold(1), &
                                   dTsf)

    ! surface heat flux error
    call surface_flux_error(fk(1), fscale_cold(1), surf_flux_error)

    ! convergence test
    lconverged = (flux_error + surf_flux_error < ferrmax .and. &
                  salt_flux_error              < serrmax .and.&
                  dTsf                         < dTsf_errmax)

    !write(*,*) istep1, 3, nit_newton, lconverged, flux_error + surf_flux_error, salt_flux_error, dTsf

  end function converged_nosnow_cold

!=======================================================================
  
  function converged_nosnow_melt(fk,          f0,          &
                                 yn,          y0,          &
                                 yscale_melt, fscale_melt, &
                                 nyn,         nit_newton)  &
                                 result(lconverged)
    
    ! check for convergence of the Newton solver for melting no snow conditions

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         fk              , & ! residual vector
         f0              , & ! initial residual vector
         yn              , & ! solution vector
         y0              , & ! initial solution vector
         yscale_melt     , & ! solution scaling vector
         fscale_melt         ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn             , & ! solution vector size
         nit_newton          ! Newton iteration counter
    
    logical(kind=log_kind) :: &
         lconverged          ! solver convergence test

    real(kind=dbl_kind) :: &
         en_uncon        , & ! solution energy unconservation (J m-2)
         flux_error      , & ! solution flux error (W m-2)
         sa_uncon        , & ! solution salt unconservation (ppt m)
         salt_flux_error     ! solution salt flux unconservation (ppt m s-1)

    ! internal heat conservation error
    call energy_unconservation_melt(fk,       fscale_melt , &
                                    nyn,                    &
                                    g_hslyr,  g_hilyr,      &
                                    g_dt,     g_dti,        &
                                    .false.,                &
                                    en_uncon, flux_error)

    ! salt conservation error
    call salt_unconservation_melt(fk,       fscale_melt ,   &
                                  nyn,      g_hilyr,        &
                                  g_dt,     g_dti,          &
                                  .false.,                  &
                                  sa_uncon, salt_flux_error)

    ! convergence test
    lconverged = (flux_error      < ferrmax .and. &
                  salt_flux_error < serrmax)

    !write(*,*) istep1, 4, nit_newton, lconverged, flux_error, salt_flux_error

  end function converged_nosnow_melt

!=======================================================================
  
  subroutine energy_unconservation_cold(f,        fscale_cold, &
                                        nyn,                   & 
                                        hslyr,    hilyr,       & 
                                        dt,       dti,         & 
                                        lsnow,                 & 
                                        en_uncon, flux_error)

    ! determine the energy unconservation of the snow/mush interior for cold conditions

    real(kind=dbl_kind), dimension(nyn), intent(in) :: & 
         f           , & ! residual vector
         fscale_cold     ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn             ! solution vector size

    real(kind=dbl_kind), intent(in) :: &
         hslyr       , & ! ice layer thickness (m)
         hilyr       , & ! snow layer thickness (m)
         dt          , & ! time step (s)
         dti             ! inverse time step (s-1)

    logical(kind=log_kind), intent(in) :: &
         lsnow           ! snow presence

    real(kind=dbl_kind), intent(out) :: &
         en_uncon    , & ! solution energy unconservation (J m-2)
         flux_error      ! solution flux error (W m-2)

    integer(kind=int_kind) :: &
         l           , & ! solution vector index
         k               ! ice/snow layer index

    en_uncon = c0

    if (lsnow) then

       ! internal heat conservation error in snow
       do k = 1, nslyr
          
          l = k + 1
          en_uncon = en_uncon + abs((f(l) * dt * hslyr) / fscale_cold(l))

       enddo ! k 

       ! internal heat conservation error in ice
       do k = 1, nilyr
          
          l = k + nslyr + 1
          en_uncon = en_uncon + abs((f(l) * dt * hilyr) / fscale_cold(l))

       enddo ! k 

    else

       ! internal heat conservation error in ice
       do k = 1, nilyr
          
          l = k + 1
          en_uncon = en_uncon + abs((f(l) * dt * hilyr) / fscale_cold(l))

       enddo ! k 

    endif

    flux_error = en_uncon * dti

    !open(55,file='./history/conservation.txt',position="append")
    !write(55,*) istep, 1, en_uncon
    !close(55)

  end subroutine energy_unconservation_cold

!=======================================================================
  
  subroutine energy_unconservation_melt(f,        fscale_melt, &
                                        nyn,                   &
                                        hslyr,    hilyr,       & 
                                        dt,       dti,         &
                                        lsnow,                 &
                                        en_uncon, flux_error)

    ! determine the energy unconservation of the snow/mush interior for melting conditions

    real(kind=dbl_kind), dimension(nyn), intent(in) :: &
         f           , & ! residual vector
         fscale_melt     ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn             ! solution vector size

    real(kind=dbl_kind), intent(in) :: &
         hslyr       , & ! ice layer thickness (m)
         hilyr       , & ! snow layer thickness (m)
         dt          , & ! time step (s)
         dti             ! inverse time step (s-1)

    logical(kind=log_kind), intent(in) :: &
         lsnow           ! snow presence

    real(kind=dbl_kind), intent(out) :: &
         en_uncon    , & ! solution energy unconservation (J m-2)
         flux_error      ! solution flux error (W m-2)

    integer(kind=int_kind) :: &
         l           , & ! solution vector index
         k               ! ice/snow layer index

    en_uncon = c0

    if (lsnow) then

       ! internal heat conservation error in snow
       do k = 1, nslyr
          
          l = k
          en_uncon = en_uncon + abs((f(l) * dt * hslyr) / fscale_melt(l))

       enddo ! k 

       ! internal heat conservation error in ice
       do k = 1, nilyr
          
          l = k + nslyr
          en_uncon = en_uncon + abs((f(l) * dt * hilyr) / fscale_melt(l))

       enddo ! k 

    else

       ! internal heat conservation error in ice
       do k = 1, nilyr
          
          l = k
          en_uncon = en_uncon + abs((f(l) * dt * hilyr) / fscale_melt(l))

       enddo ! k 

    endif

    flux_error = en_uncon * dti

    !open(55,file='./history/conservation.txt',position="append")
    !write(55,*) istep, en_uncon
    !close(55)

  end subroutine energy_unconservation_melt

  !=======================================================================

  subroutine salt_unconservation_cold(f,        fscale_cold,    &
                                      nyn,      hilyr,          &
                                      dt,       dti,            & 
                                      lsnow,                    &
                                      sa_uncon, salt_flux_error)

    ! determine the salt unconservation of the mush interior for cold conditions

    real(kind=dbl_kind), dimension(nyn), intent(in) :: &
         f               , & ! residual vector
         fscale_cold         ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn                 ! solution vector size

    real(kind=dbl_kind), intent(in) :: &
         hilyr           , & ! ice layer thickness (m)
         dt              , & ! time step (s)
         dti                 ! inverse time step (s-1)

    logical(kind=log_kind), intent(in) :: &
         lsnow               ! snow presence

    real(kind=dbl_kind), intent(out) :: &
         sa_uncon        , & ! solution salt unconservation (ppt m)
         salt_flux_error     ! solution salt flux unconservation (ppt m s-1)

    integer(kind=int_kind) :: &
         k               , & ! ice layer index
         l                   ! solution vector index

    sa_uncon = c0

    if (lsnow) then

       ! salt conservation error in ice
       do k = 1, nilyr

          l = k + nilyr + nslyr + 1
          sa_uncon = sa_uncon + abs((f(l) * dt * hilyr) / fscale_cold(l))

       enddo ! k

    else

       ! salt conservation error in ice
       do k = 1, nilyr

          l = k + nilyr + 1
          sa_uncon = sa_uncon + abs((f(l) * dt * hilyr) / fscale_cold(l))

       enddo ! k

    endif

    salt_flux_error = sa_uncon * dti

  end subroutine salt_unconservation_cold

  !=======================================================================

  subroutine salt_unconservation_melt(f,        fscale_melt,    &
                                      nyn,      hilyr,          &
                                      dt,       dti,            &
                                      lsnow,                    &
                                      sa_uncon, salt_flux_error)

    ! determine the salt unconservation of the mush interior for cold conditions

    real(kind=dbl_kind), dimension(nyn), intent(in) :: &
         f               , & ! residual vector
         fscale_melt         ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn                 ! solution vector size

    real(kind=dbl_kind), intent(in) :: &
         hilyr           , & ! ice layer thickness (m)
         dt              , & ! time step (s)
         dti                 ! inverse time step (s-1)

    logical(kind=log_kind), intent(in) :: &
         lsnow               ! snow presence

    real(kind=dbl_kind), intent(out) :: &
         sa_uncon        , & ! solution salt unconservation (ppt m)
         salt_flux_error     ! solution salt flux unconservation (ppt m s-1)

    integer(kind=int_kind) :: &
         k               , & ! ice layer index
         l                   ! solution vector index

    sa_uncon = c0

    if (lsnow) then

       ! salt conservation error in ice
       do k = 1, nilyr

          l = k + nilyr + nslyr
          sa_uncon = sa_uncon + abs((f(l) * dt * hilyr) / fscale_melt(l))

       enddo ! k

    else

       ! salt conservation error in ice
       do k = 1, nilyr

          l = k + nilyr
          sa_uncon = sa_uncon + abs((f(l) * dt * hilyr) / fscale_melt(l))

       enddo ! k

    endif

    salt_flux_error = sa_uncon * dti

  end subroutine salt_unconservation_melt

  !=======================================================================

  subroutine surface_temperature_error(fsf,               &
                                       Tsf,    km,        & 
                                       nlyr,   hlyri,     &
                                       shcoef, lhcoef,    &
                                       rhoa,   fsf_scale, &
                                       dTsf)

    ! determine the approximate error in the surface temperature for cold conditions
    ! this takes into account only the error from f(Tsf) not f(Tin(1)/f(Tsn(1)))

    real(kind=dbl_kind), intent(in) :: &
         fsf       , & ! surface residual
         Tsf       , & ! ice/snow surface temperature (C)
         km        , & ! upper ice/snow layer conductivity (W m-1 K-1)
         hlyri     , & ! ice/snow layer inverse thickness (m-1)
         fsf_scale , & ! surface residual scaling
         shcoef    , & ! transfer coefficient for sensible heat
         lhcoef    , & ! transfer coefficient for latent heat
         rhoa          ! air density (kg/m^3)

    integer(kind=int_kind), intent(in) :: &
         nlyr          ! number of ice/snow layers

    real(kind=dbl_kind), intent(out) :: &
         dTsf          ! estimated surface temperature error (C)

    real(kind=dbl_kind) :: &
         A         , & ! factor of surface temperature error
         B         , & ! factor of surface temperature error
         D         , & ! factor of surface temperature error
         E         , & ! factor of surface temperature error
         F         , & ! factor of surface temperature error
         Tsfk      , & ! surface absolute temperature (K)
         denom         ! denominator of surface temperature error

    ! surface temperature in Kelvin
    Tsfk = max(Tsf + Tffresh, puny)

    A = fsf_scale

    B = fsf_scale * c2 * km * hlyri

    D = -emissivity * stefan_boltzmann

    E = -shcoef

    F = (-lhcoef * qqqice) / rhoa

    denom = c4 * A * D + A * E + &
            A * F * (TTTice / Tsfk**2) * exp(-TTTice / Tsfk) - B

    dTsf = abs(fsf / denom)

  end subroutine surface_temperature_error

  !=======================================================================

  subroutine surface_flux_error(fsf, fsf_fscale, surf_flux_error)

    ! determine the error in the surface energy flux

    real(kind=dbl_kind), intent(in) :: &
         fsf         , & ! surface residual
         fsf_fscale      ! surface residual scaling

    real(kind=dbl_kind), intent(out) :: &
         surf_flux_error ! solution surface flux error (W m-2)

    surf_flux_error = abs(fsf / fsf_fscale)

  end subroutine surface_flux_error

  !=======================================================================
  ! Residual componets
  !=======================================================================

  subroutine residual_heat_diffusion(q,     q0,     &
                                     T,     kcstar, &
                                     swabs, fscale, &
                                     f,     nlyr,   &
                                     hlyri, hlyri2, &
                                     dti)

    ! residual related to internal heat conservation of either snow or mush

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         q      , & ! layer enthalpy (J m-3)
         q0     , & ! layer enthalpy (J m-3) at beginning of timestep 
         kcstar , & ! interface conductivity (W m-1 K-1)
         swabs  , & ! SW radiation absorbed in layers (W m-2)
         fscale     ! residual scaling vector

    real(kind=dbl_kind), dimension(0:), intent(in) :: &
         T          ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         f          ! heat diffusion residual vector

    integer(kind=int_kind), intent(in) :: &
         nlyr       ! number of layers

    real(kind=dbl_kind), intent(in) :: &
         hlyri  , & ! inverse layer thickness (m-1)
         hlyri2 , & ! inverse squared layer thickness (m-2)
         dti        ! inverse time step (s-1)
    
    integer(kind=int_kind) :: &
         k          ! ice/snow layer index
    
    do k = 1, nlyr
       
       f(k) = ((q(k) - q0(k)) * dti &
           
            ! short wave absorption
            - swabs(k) * hlyri &
 
            ! heat diffusion
            - (kcstar(k+1) * (T(k+1) - T(k)  ) - &
               kcstar(k)   * (T(k)   - T(k-1))) * hlyri2 &

            ! rescale residual
            ) * fscale(k)
            
    enddo ! k
    
  end subroutine residual_heat_diffusion

!=======================================================================
  
  subroutine interface_conductivities(kc, kcstar, nlyr)
    
    ! calculate the interface heat conductivities from cell heat conductivities
    ! harmonic mean of layer conductivities

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         kc     ! layer conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         kcstar ! interface conductivity (W m-1 K-1)

    integer(kind=int_kind), intent(in) :: &
         nlyr   ! number of layers
    
    integer(kind=int_kind) :: &
         k      ! ice/snow layer index
    
    kcstar(1) = c2 * kc(1)
    
    do k = 2, nlyr
       
       kcstar(k) = (c2 * kc(k) * kc(k-1)) / (kc(k) + kc(k-1))
       
    enddo ! k
    
    kcstar(nlyr+1) = c2 * kc(nlyr)
    
  end subroutine interface_conductivities

  !=======================================================================
     
  subroutine residual_surface(Tsf, fscale, &
                              f,   fsurfn, &
                              km,  Tin,    &
                              Sin, hlyri)

    ! calculate the resiual related to the surface temperature

    real(kind=dbl_kind), intent(in) :: &
         Tsf    , & ! ice/snow surface temperature (C)
         fscale , & ! surface residual scaling
         fsurfn , & ! net flux to top surface, excluding fcondtopn
         km     , & ! upper layer conductivity (W m-1 K-1)
         Tin    , & ! upper ice layer temperature (C)
         Sin    , & ! upper ice layer bulk salinity (ppt)
         hlyri      ! inverse layer thickness (m-1)

    real(kind=dbl_kind), intent(out) :: &
         f          ! surface residual

#if defined notz_experiment

    real(kind=dbl_kind), parameter :: &
         Tsurf_exp = -10.0_dbl_kind

    f = (Tsf - Tsurf_exp) * fscale

#elif defined flushing_notz

    real(kind=dbl_kind) :: &
         Tsurf_exp 

    Tsurf_exp = liquidus_temperature_mush(Sin)

    f = (Tsf - Tsurf_exp) * fscale

#else
    f = (fsurfn - km * c2 * (Tsf - Tin) * hlyri) * fscale
#endif

  end subroutine residual_surface

!=======================================================================

  subroutine surface_heat_flux(Tsf,     fswsfc, &
                               rhoa,    flw,    &
                               potT,    Qa,     &
                               shcoef,  lhcoef, &
                               flwoutn, fsensn, &
                               flatn,   fsurfn)       
    
    ! heat flux into ice
    
    ! input surface temperature
    real(kind=dbl_kind), intent(in) :: &
         Tsf             ! ice/snow surface temperature (C)
    
    ! input variables
    real(kind=dbl_kind), intent(in) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat
    
    ! output
    real(kind=dbl_kind), intent(out) :: &
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsurfn          ! net flux to top surface, excluding fcondtopn
    
    ! local variables
    real(kind=dbl_kind) :: &
         TsfK        , & ! ice/snow surface temperature (K)
         Qsfc        , & ! saturated surface specific humidity (kg/kg)
         qsat        , & ! the saturation humidity of air (kg/m^3)
         flwdabs     , & ! downward longwave absorbed heat flx (W/m^2)
         tmpvar          ! 1/TsfK
    
    
    ! ice surface temperature in Kelvin
    TsfK = max(Tsf + Tffresh, c1)
    tmpvar = c1/TsfK
    
    ! saturation humidity
    qsat    = qqqice * exp(-TTTice*tmpvar)
    Qsfc    = qsat / rhoa
    
    ! longwave radiative flux
    flwdabs =  emissivity * flw
    flwoutn = -emissivity * stefan_boltzmann * TsfK**4
    
    ! downward latent and sensible heat fluxes
#if !defined notz_experiment
    fsensn = shcoef * (potT - TsfK)
#elif defined  notz_experiment
    fsensn = c0
#else
    fsensn = heat_coeff_air_global * (potT - TsfK)
#endif
    flatn  = lhcoef * (Qa - Qsfc)
    
    ! combine fluxes
#if !defined notz_experiment
    fsurfn = fswsfc + flwdabs + flwoutn + fsensn + flatn
#elif defined notz_experiment
    fsurfn = c0
    flwdabs = c0
    flatn = c0
    flwoutn = c0
#else
    flatn = c0
    fsurfn = fsensn + flwoutn * longwave_fudge_global + flwdabs
#endif

  end subroutine surface_heat_flux

  !=======================================================================
  
  subroutine dsurface_heat_flux_dTsf(Tsf,     fswsfc, &
                                     rhoa,    flw,    &
                                     potT,    Qa,     &
                                     shcoef,  lhcoef, &
                                     dfsurfn_dTsf)
    
    ! input surface temperature
    real(kind=dbl_kind), intent(in) :: &
         Tsf               ! ice/snow surface temperature (C)
    
    ! input variables
    real(kind=dbl_kind), intent(in) :: &
         fswsfc        , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa          , & ! air density (kg/m^3)
         flw           , & ! incoming longwave radiation (W/m^2)
         potT          , & ! air potential temperature  (K)
         Qa            , & ! specific humidity (kg/kg)
         shcoef        , & ! transfer coefficient for sensible heat
         lhcoef            ! transfer coefficient for latent heat
    
    ! output
    real(kind=dbl_kind), intent(out) :: &
         dfsurfn_dTsf      ! net flux to top surface, excluding fcondtopn
    
    ! local variables
    real(kind=dbl_kind) :: &
         TsfK          , & ! ice/snow surface temperature (K)
         dQsfc_dTsf    , & ! saturated surface specific humidity (kg/kg)
         qsat          , & ! the saturation humidity of air (kg/m^3)
         flwdabs       , & ! downward longwave absorbed heat flx (W/m^2)
         tmpvar            ! 1/TsfK
    
    real(kind=dbl_kind) :: &
         dflwoutn_dTsf , & ! derivative of longwave flux wrt surface temperature
         dfsensn_dTsf  , & ! derivative of sensible heat flux wrt surface temperature
         dflatn_dTsf       ! derivative of latent heat flux wrt surface temperature
    
    ! ice surface temperature in Kelvin
    TsfK = max(Tsf + Tffresh, c1)
    tmpvar = c1/TsfK
    
    ! saturation humidity
    qsat          = qqqice * exp(-TTTice*tmpvar)
    dQsfc_dTsf    = (qsat / rhoa) * TTTice * tmpvar * tmpvar
    
    ! longwave radiative flux
    dflwoutn_dTsf = -c4 * emissivity * stefan_boltzmann * TsfK**3
    
    ! downward latent and sensible heat fluxes
    dfsensn_dTsf = -shcoef
    dflatn_dTsf  = - lhcoef * dQsfc_dTsf
    
    ! combine fluxes
    dfsurfn_dTsf = dflwoutn_dTsf + dfsensn_dTsf + dflatn_dTsf
    
  end subroutine dsurface_heat_flux_dTsf

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
         Sin_min       = 0.1_dbl_kind                          ! minimum bulk salinity (ppt)

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
       dS_guess = ((q(k) * (Sbr(k+1) - Sbr(k))) / hilyr + dSdt(k)) * dt
       
       alpha = (Sin_min - Sin(k)) / dS_guess

       if (alpha < c0 .or. alpha > c1) alpha = c1

       q(k)    = q(k)    * alpha
       dSdt(k) = dSdt(k) * alpha

    enddo ! k

  end subroutine explicit_flow_velocities

!=======================================================================

  subroutine salinity_residual(Sin,     Sin0,   &
                               Tin,     dt,     &
                               q,       dSdt,   &
                               fs,      fq,     &
                               fsscale, fqscale)

    ! given strength of gravity drainage determine contribution to residual 
    ! in heat and bulk salinity conservation equations

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         Sin     , & ! ice layer bulk salinity (ppt)
         Sin0    , & ! ice layer bulk salinity at beginning of timestep (ppt)
         Tin     , & ! ice layer temperature (C)
         dSdt    , & ! slow mode drainage rate (ppt s-1)
         fsscale , & ! residual scaling vector for bulk salinity
         fqscale     ! residual scaling vector for ice enthalpy

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q           ! rapid mode upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(1:nilyr), intent(inout) :: &
         fs      , & ! residual vector for bulk salinity
         fq          ! residual vector for ice enthalpy

    real(kind=dbl_kind), intent(in) :: &
         dt          ! time step (s)

    real(kind=dbl_kind) :: &
         adv_s   , & ! advective contribution for bulk salinity 
         adv_q       ! advective contribution for ice enthalpy

    real(kind=dbl_kind), dimension(1:nilyr+1) :: &
         Sbr     , & ! ice layer brine salinity (ppt)
         qbr         ! ice layer brine enthalpy (J m-3)
    
    integer(kind=int_kind) :: &
         k           ! ice layer index

    ! ocean properties
    Sbr(nilyr+1) = g_sss
    qbr(nilyr+1) = g_qocn

    ! calculate residual from base upwards
    do k = nilyr, 1, -1

       ! brine properties
       Sbr(k) = max(liquidus_brine_salinity_mush(Tin(k)), Sin(k))

       qbr(k) = enthalpy_brine(Tin(k))

       ! advection
       adv_s = (q(k) * g_hilyri) * (Sbr(k+1) - Sbr(k))

       adv_q = (q(k) * g_hilyri) * (qbr(k+1) - qbr(k))

       ! residuals 
       fs(k) = ((Sin(k) - Sin0(k)) / dt - adv_s &
                                        - dSdt(k)) * fsscale(k)

       fq(k) = fq(k) - adv_q * fqscale(k)

    enddo ! k

  end subroutine salinity_residual

!=======================================================================

  subroutine drainage_fluxes(q,        dSdt,    &
                             qbr,      Sbr,     &
                             hilyr,    dt,      &
                             fadvheat, fadvsalt)

    ! calculate heat and salt flux associated with gravity drainage

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q            ! rapid mode upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         dSdt         ! slow mode drainage rate (ppt s-1)

    real(kind=dbl_kind), dimension(1:nilyr+1), intent(in) :: &
         qbr      , & ! ice layer brine enthalpy (J m-3)
         Sbr          ! ice layer brine salinity (ppt)

    real(kind=dbl_kind), intent(in) :: &
         hilyr    , & ! ice layer thickness (m)
         dt           ! time step (s)

    real(kind=dbl_kind), intent(inout) :: &
         fadvheat , & ! heat flux to ocean from brine advection (W m-2)
         fadvsalt     ! salt flux to ocean from brine advection 
        
    integer(kind=int_kind) :: &
         k            ! ice layer index

    ! calculate fluxes from base upwards
    do k = nilyr, 1, -1

       fadvsalt = fadvsalt - q(k) * (Sbr(k+1) - Sbr(k)) &
                           - dSdt(k) * hilyr

       fadvheat = fadvheat - q(k) * (qbr(k+1) - qbr(k))

    enddo ! k

    !open(55,file='history/q_conservation.txt',position="append")
    !write(55,*) istep, fadvocn*g_dt
    !close(55)

  end subroutine drainage_fluxes

!=======================================================================
! Flushing
!=======================================================================

  subroutine flushing_velocity(Tin,   Sin,   &
                               hin,   hsn,   &
                               hilyr,        &
                               hpond, apond, &
                               dt,    w)
   
    ! calculate the vertical flushing Darcy velocity
    ! negative - downward flushing
    ! positive - upward flushing

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Tin       , & ! ice layer temperature (C)
         Sin           ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), intent(in) :: &
         hilyr     , & ! ice layer thickness (m)
         hpond     , & ! melt pond thickness (m)
         apond     , & ! melt pond area (-)
         hsn       , & ! snow thickness (m)
         hin       , & ! ice thickness (m)
         dt            ! time step (s)

    real(kind=dbl_kind), intent(out) :: &
         w             ! vertical flushing Darcy flow rate (m s-1)

    real(kind=dbl_kind), parameter :: &
         viscosity_dyn = 1.79e-3_dbl_kind ! dynamic viscosity

    real(kind=dbl_kind) :: &
         phi       , & ! ice layer liquid fraction
         perm      , & ! ice layer permeability (m2)
         Mice      , & ! mass of ice (kg m-2)
         perm_harm , & ! harmonic mean of ice permeability (m2)
         hocn      , & ! ocean surface height above ice base (m)
         hbrine    , & ! brine surface height above ice base (m)
         dhhead    , & ! hydraulic head (m)
         w_up_max  , & ! maximum upward flushing Darcy flow rate (m s-1) 
         w_down_max    ! maximum downward flushing Darcy flow rate (m s-1) 

    integer(kind=int_kind) :: &
         k             ! ice layer index

    ! initialize
    Mice = c0
    perm_harm = c0

    do k = 1, nilyr

       ! liquid fraction
       phi = liquid_fraction(Tin(k), Sin(k))

       ! permeability
       perm = permeability(phi)
       !perm = 3.0e-8 * phi**3
       !perm = 3.0e-8 * (0.05_dbl_kind)**3

       ! ice mass
       Mice = Mice + hilyr * (phi * rhow + (c1 - phi) * rhoi)

       ! permeability harmonic mean
       perm_harm = perm_harm + c1 / (perm + 1e-30_dbl_kind)

    enddo ! k

    perm_harm = real(nilyr,dbl_kind) / perm_harm 

    ! calculate ocean surface height above bottom of ice
    hocn = (Mice + hpond * apond * rhow + hsn * rhos) / rhow

    ! calculate brine height above bottom of ice
    hbrine = hin + hpond

    ! hydraulic head
    dhhead = hbrine - hocn
    if (hbrine == hin) then
       dhhead = min(dhhead,c0)
    endif

    ! darcy flow through ice
    w = -(perm_harm * rhow * gravit * (dhhead / hin)) / viscosity_dyn

    ! maximum down flow to drain pond
    w_down_max = -(hpond * apond) / dt
    
    ! maximum up flow to flood to sea level
    w_up_max = max(((Mice + hsn * rhos - hin * rhow) * (rhoi - rhos)) / (dt * rhos * rhow), c0)

    !open(44,file='./history/wup.txt',position='append')
    !write(44,*) istep, w, w_up_max, min(max(w,w_down_max),w_up_max), hocn, hin
    !close(44)

    !open(44,file='./history/wdn.txt',position='append')
    !write(44,*) istep, w, abs(min(w,c0)), w_down_max, hpond, apond, &
    !     min(max(w,w_down_max),w_up_max), perm_harm, hbrine, hocn, dhhead, hbrine - hocn
    !close(44)

    ! limit flow
    w = min(max(w,w_down_max),w_up_max)
    !w = c0!

#if flushing_up == 0
    ! no flooding
    w = min(w,c0)
#endif

#if flushing_down == 0
    ! no flushing
    w = max(w,c0)
#endif

  end subroutine flushing_velocity

!=======================================================================

  subroutine flushing_advection(Tin,     Sin,    &
                                hin,     hsn,    &
                                hilyr,   hilyri, &
                                hpond,   apond,  &
                                sss,     qocn,   &
                                dt,              &
                                fs,      fq,     &
                                fsscale, fqscale)

    ! given strength of flushing determine contribution to residual 
    ! in bulk salinity and enthalpy conservation equations

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Tin     , & ! ice layer temperature (C)
         Sin     , & ! ice layer bulk salinity (ppt)
         fsscale , & ! residual scaling vector for bulk salinity
         fqscale     ! residual scaling vector for ice enthalpy

    real(kind=dbl_kind), intent(in) :: &
         hin     , & ! ice thickness (m)
         hsn     , & ! snow thickness (m)
         hilyr   , & ! ice layer thickness (m)
         hilyri  , & ! snow layer thickness (m)
         hpond   , & ! melt pond thickness (m)
         apond   , & ! melt pond area (-)
         sss     , & ! sea surface salinity (ppt)
         qocn    , & ! ocean enthalpy (J m-3)
         dt          ! time step (s)

    real(kind=dbl_kind), dimension(nilyr), intent(inout) :: &
         fs      , & ! residual vector for bulk salinity
         fq          ! residual vector for ice enthalpy

    integer(kind=int_kind) :: &
         k           ! ice layer index

    real(kind=dbl_kind), dimension(nilyr) :: &
         Sbr     , & ! ice layer brine salinity (ppt)
         qbr         ! ice layer brine enthalpy (J m-2)

    real(kind=dbl_kind) :: &
         Spond   , & ! melt pond salinity (ppt)
         qpond   , & ! melt pond enthalpy (J m-2)
         w       , & ! vertical flushing Darcy flow rate (m s-1)
         wdown   , & ! downward flushing Darcy flow rate (m s-1) 
         wup         ! upward flushing Darcy flow rate (m s-1) 

    ! brine properties
    do k = 1, nilyr

       Sbr(k) = max(liquidus_brine_salinity_mush(Tin(k)), Sin(k))
       qbr(k) = enthalpy_brine(Tin(k))

    enddo ! k

    ! melt pond properties
    Spond = 0.0_dbl_kind
    qpond = enthalpy_brine(c0)

    ! calculate the flushing velocity
    call flushing_velocity(Tin,   Sin,   &
                           hin,   hsn,   &
                           hilyr,        &
                           hpond, apond, &
                           dt,    w)
    
    wdown = -min(w,c0) * hilyri
    wup   =  max(w,c0) * hilyri

    ! calculate residuals
    k = 1
    fs(k) = fs(k) - (wdown * (Spond    - Sbr(k)) + &
                     wup   * (Sbr(k+1) - Sbr(k))) * fsscale(k)
    fq(k) = fq(k) - (wdown * (qpond    - qbr(k)) + &
                     wup   * (qbr(k+1) - qbr(k))) * fqscale(k)

    do k = 2, nilyr-1
       
       fs(k) = fs(k) - (wdown * (Sbr(k-1) - Sbr(k)) + &
                        wup   * (Sbr(k+1) - Sbr(k))) * fsscale(k)
       fq(k) = fq(k) - (wdown * (qbr(k-1) - qbr(k)) + &
                        wup   * (qbr(k+1) - qbr(k))) * fqscale(k)

    enddo ! k

    k = nilyr
    fs(k) = fs(k) - (wdown * (Sbr(k-1) - Sbr(k)) + &
                     wup   * (sss      - Sbr(k))) * fsscale(k)
    fq(k) = fq(k) - (wdown * (qbr(k-1) - qbr(k)) + &
                     wup   * (qocn - qbr(k))) * fqscale(k)

  end subroutine flushing_advection

!=======================================================================

  subroutine flushing_fluxes(w, qbr, Sbr, fadvheat, fadvsalt)

    ! determine heat and salt flux to ocean from flushing
    
    real(kind=dbl_kind), intent(in) :: &
         w            ! vertical flushing Darcy flow rate (m s-1)
    
    real(kind=dbl_kind), dimension(1:nilyr+1), intent(in) :: &
         qbr      , & ! ice layer brine enthalpy (J m-2)
         Sbr          ! ice layer brine salinity (ppt)

    real(kind=dbl_kind), intent(inout) :: &
         fadvheat , & ! heat flux to ocean from brine advection (W m-2)
         fadvsalt     ! salt flux to ocean from brine advection

    integer(kind=int_kind) :: &
         k            ! ice layer index

    real(kind=dbl_kind) :: &
         Spond    , & ! melt pond salinity (ppt)
         qpond    , & ! melt pond enthalpy (J m-2)
         wdown    , & ! downward flushing Darcy flow rate (m s-1) 
         wup          ! upward flushing Darcy flow rate (m s-1) 

    ! pond properties
    Spond = 0.0_dbl_kind
    qpond = enthalpy_brine(c0)

    wdown = -min(w,c0)
    wup   =  max(w,c0)
    
    ! fluxes
    k = 1
    fadvsalt = fadvsalt - (wdown * (Spond    - Sbr(k)) + &
                           wup   * (Sbr(k+1) - Sbr(k)))
    fadvheat = fadvheat - (wdown * (qpond    - qbr(k)) + &
                           wup   * (qbr(k+1) - qbr(k)))

    do k = 2, nilyr
       
       fadvsalt = fadvsalt - (wdown * (Sbr(k-1) - Sbr(k)) + &
                              wup   * (Sbr(k+1) - Sbr(k)))
       fadvheat = fadvheat - (wdown * (qbr(k-1) - qbr(k)) + &
                              wup   * (qbr(k+1) - qbr(k)))

    enddo ! k

    !open(55,file='history/flush_conservation.txt',position="append")
    !write(55,*) istep, flux_heat, flux_heat*g_dt
    !close(55)

  end subroutine flushing_fluxes

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

    if (apond > c0 .and. hpond > c0 .and. w < c0) then

       hpond = hpond + (w * dt) / apond
       
       hpond = max(hpond, c0)

    endif

    open(55,file="history/flush_pond.txt",position='append')
    write(55,*) istep, w, abs(min(w,c0)), (min(w,c0) * dt), (min(w,c0) * dt) / apond, hpond
    close(55)

  end subroutine flush_pond

 !=======================================================================

  subroutine flood_ice(w,     dt,    &
                       hsn,   hin,   &
                       hslyr, hilyr, & 
                       qsn,   qin,   &
                       Sin,   Sbr,   &
                       qbr,   snoice)

    ! given upwards flushing brine flow calculate amount of snow ice and
    ! convert snow to ice with appropriate properties

    real(kind=dbl_kind), intent(in) :: &
         w           , & ! vertical flushing Darcy flow rate (m s-1)
         dt          , & ! time step (s)
         hsn         , & ! snow thickness (m)
         hin             ! ice thickness (m)

    real(kind=dbl_kind), dimension(nslyr), intent(inout) :: &
         qsn             ! snow layer enthalpy (J m-2)

    real(kind=dbl_kind), dimension(nilyr), intent(inout) :: &
         qin         , & ! ice layer enthalpy (J m-2)
         Sin             ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(nilyr), intent(in) :: &
         Sbr         , & ! ice layer brine salinity (ppt)
         qbr             ! ice layer brine enthalpy (J m-2)

    real(kind=dbl_kind), intent(out) :: &
         hslyr       , & ! ice layer thickness (m)
         hilyr       , & ! snow layer thickness (m)
         snoice          ! snow ice formation

    real(kind=dbl_kind) :: &
         hin2        , & ! new ice thickness (m)
         hsn2        , & ! new snow thickness (m)
         hilyr2      , & ! new ice layer thickness (m)
         hslyr2      , & ! new snow layer thickness (m)
         dh          , & ! thickness of snowice formation (m)
         phi_snowice , & ! liquid fraction of new snow ice
         Sin_snowice , & ! bulk salinity of new snowice (ppt)
         qin_snowice , & ! ice enthalpy of new snowice (J m-2)
         qsn_snowice     ! snow enthalpy of snow thats becoming snowice (J m-2)

    if (w > c0) then 

       !integer :: k
       !real(kind=dbl_kind) :: einit, sinit
       !real(kind=dbl_kind) :: efinal, sfinal
       !real(kind=dbl_kind) :: eadded, sadded
       !real(kind=dbl_kind) :: euncon, suncon
       
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
       
       ! liquid fraction of newly formed snow ice
       phi_snowice = (c1 - rhos / rhoi)
       
       ! calculate thickness of new ice added
       dh = (w * dt) / phi_snowice
       
       ! enthalpy of snow thats becoming snowice
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
       
       ! change thicknesses
       hilyr = hilyr2
       hslyr = hslyr2
       snoice = dh
       
       ! initial energy and salt
       !efinal = c0
       !sfinal = c0
       
       !do k = 1, nslyr
       !   efinal = efinal + qsn(k) * hslyr
       !enddo ! k
       
       !do k = 1, nilyr
       !   efinal = efinal + qin(k) * hilyr
       !   sfinal = sfinal + Sin(k) * hilyr
       !enddo ! k
       
       !eadded = w * dt * qbr(1)
       !sadded = w * dt * Sbr(1)
       
       !euncon = efinal - einit - eadded
       !suncon = sfinal - sinit - sadded
       
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

  subroutine ice_snow_interface_temperature(Tis,           &
                                            Tin,    Tsn,   &
                                            ki,     ks,    &
                                            hilyri, hslyri)
    
    ! determine the mush/snow interface temperature from heat flux conservation

    real(kind=dbl_kind), intent(out) :: &
         Tis        ! ice/snow interface temperature (C)

    real(kind=dbl_kind), intent(in) :: &
         Tin    , & ! top ice layer temperature (C)
         Tsn    , & ! bottom snow layer temperature (C)
         ki     , & ! top ice layer conductivity (W m-1 K-1)
         ks     , & ! bottom snow layer conductivity (W m-1 K-1)
         hilyri , & ! top ice layer inverse thickness (m-1)
         hslyri     ! bottom snow layer inverse thickness (m-1)
    
    real(kind=dbl_kind) :: &
         A      , & ! ice conductivity divided by ice layer thickness (W m-2 K-1)
         B          ! snow conductivity divided by snow layer thickness (W m-2 K-1)
    
    A = ki * hilyri
    B = ks * hslyri
    
    Tis = (Tin * A + Tsn * B) / (A + B)
    
  end subroutine ice_snow_interface_temperature

  !=======================================================================
  
  subroutine dTisdTiTs(dTisdTi, dTisdTs, &
                       Tin,     Tsn,     &
                       ki,      ks,      &
                       hilyr,   hslyr,   &
                       dkidTi,  dksdTs)
    
    ! derivative of the mush/snow interface temperature

    real(kind=dbl_kind), intent(out) :: &
         dTisdTi , & ! derivative of ice/snow interface temperature wrt ice temperature (-)
         dTisdTs     ! derivative of ice/snow interface temperature wrt snow temperature (-)

    real(kind=dbl_kind), intent(in) :: &
         Tin     , & ! top ice layer temperature (C)
         Tsn     , & ! bottom ice layer temperature (C)
         ki      , & ! top ice layer conductivity (W m-1 K-1)
         ks      , & ! bottom snow layer conductivity (W m-1 K-1)
         hilyr   , & ! ice layer thickness (m)
         hslyr   , & ! snow layer thickness (m)
         dkidTi  , & ! derivative of ice conductivity wrt ice temperature (W m-1 K-2)
         dksdTs      ! derivative of snow conductivity wrt snow temperature (W m-1 K-2)
    
    real(kind=dbl_kind) :: &
         A       , & ! ice conductivity divided by ice layer thickness (W m-2 K-1)
         B           ! snow conductivity divided by snow layer thickness (W m-2 K-1)
    
    A = ki / hilyr
    B = ks / hslyr
    
    dTisdTi = (A + ((Tin * dkidTi) / hilyr)) / (A + B) - &
              (Tin * A + Tsn * B) * ((dkidTi / hilyr) / (A + B)**2)
    
    dTisdTs = (B + ((Tsn * dksdTs) / hslyr)) / (A + B) - &
              (Tin * A + Tsn * B) * ((dksdTs / hslyr) / (A + B)**2)
    
  end subroutine dTisdTiTs

!=======================================================================

  subroutine dTisdTiTs_fixedk(dTisdTi, dTisdTs, &
                              Tin,     Tsn,     &
                              ki,      ks,      &
                              hilyr,   hslyr)
    
    ! derivative of the mush/snow interface temperature assuming fixed heat conductivities

    real(kind=dbl_kind), intent(out) :: &
         dTisdTi , & ! derivative of ice/snow interface temperature wrt upper ice layer temperature (-) 
         dTisdTs     ! derivative of ice/snow interface temperature wrt lower snow layer temperature (-) 

    real(kind=dbl_kind), intent(in) :: &
         Tin     , & ! top ice layer temperature (C)
         Tsn     , & ! bottom ice layer temperature (C)
         ki      , & ! top ice layer conductivity (W m-1 K-1)
         ks      , & ! bottom snow layer conductivity (W m-1 K-1)
         hilyr   , & ! ice layer thickness (m)
         hslyr       ! snow layer thickness (m)
    
    real(kind=dbl_kind) :: &
         A       , & ! ice conductivity divided by ice layer thickness (W m-2 K-1)
         B           ! snow conductivity divided by snow layer thickness (W m-2 K-1)
    
    A = ki / hilyr
    B = ks / hslyr
    
    dTisdTi = A / (A + B)
    
    dTisdTs = B / (A + B)
    
  end subroutine dTisdTiTs_fixedk

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

  subroutine snow_derived_quantities(qsn, Tsn)

    ! temperature of the snow layers from the snow layer enthalpies

    real(kind=dbl_kind), dimension(1:nslyr), intent(in) :: &
         qsn ! snow layer enthalpy (J m-3) 

    real(kind=dbl_kind), dimension(1:nslyr), intent(out) :: &
         Tsn ! snow layer temperature (C)
    
    integer(kind=int_kind) :: &
         k   ! snow layer index
    
    do k = 1, nslyr
       
       Tsn(k) = temperature_snow(qsn(k))
           
    enddo ! k
    
  end subroutine snow_derived_quantities

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
    
    !if (Tsn <= c0) then
       
       qsn = -rhos * (-cp_ice * Tsn + Lfresh)
       
    !else
       
    !   qsn = rhos * cp_ice * Tsn
       
    !endif

    
  end function enthalpy_snow

!=======================================================================
  
  function dqdT_snow() result(dq_dT)
    
    ! derivative of snow enthalpy wrt snow temperature

    real(kind=dbl_kind) :: &
         dq_dT ! derivative of snow enthalpy wrt to snow temperature (J m-3 K-1)
    
    dq_dT = rhos * cp_ice
    
  end function dqdT_snow

!=======================================================================
  
  function dqdT_ice() result(dq_dT)
    
    ! derivative of pure ice enthalpy wrt ice temperature

    real(kind=dbl_kind) :: &
         dq_dT ! derivative of ice enthalpy wrt to ice temperature (J m-3 K-1)
    
    dq_dT = rhoi * cp_ice
    
  end function dqdT_ice
  
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

    !real(kind=dbl_kind) :: &
    !     qm
    
    !qm = -rhos * Lfresh
    
    !if (qsn < qm) then
       
    !   Tsn = (qsn / (-rhos) - Lfresh) / (-cp_ice)

       Tsn = A * qsn + B
       
    !else if (qsn > c0) then
       
    !   Tsn = qsn / (rhos * cp_ice)
       
    !else
       
    !   Tsn = c0
       
    !endif

  end function temperature_snow

!=======================================================================
! Mushy Layer Formulation - Assur (1958) liquidus
!=======================================================================

 subroutine mush_derived_quantities(qin, Sin, Tin, Sbr, phi)

   ! physical mush quantities drived from the mush enthalpy and bulk salinity

    real(kind=dbl_kind), dimension(1:nilyr), intent(in) :: &
         qin , & ! ice layer enthalpy (J m-3) 
         Sin     ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(1:nilyr), intent(out) :: &
         Tin , & ! ice layer temperature (C)
         Sbr , & ! ice layer brine salinity (ppt)
         phi     ! ice layer liquid fraction
    
    integer(kind=int_kind) :: &
         k       ! ice layer index
    
    do k = 1, nilyr
       
       ! derived state variables
       Tin(k) = temperature_mush(qin(k), Sin(k))
       Sbr(k) = liquidus_brine_salinity_mush(Tin(k))
       phi(k) = liquid_fraction(Tin(k), Sin(k))

    enddo ! k
    
  end subroutine mush_derived_quantities

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

    Sbr = Sbr * lsubzero

  end function liquidus_brine_salinity_mush

!=======================================================================

  function dliquidus_brine_salinity_mush_dT(Tin) result(dSbr_dTin)

    ! derivative of liquidus salinity wrt temperature
    ! based on empirical data from Assur (1958)

    real(kind=dbl_kind), intent(in) :: &
         Tin          ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         dSbr_dTin    ! derivative of brine salinity wrt to ice temperature (ppt K-1)

    real(kind=dbl_kind) :: &
         t_high   , & ! mask for high temperature liquidus region
         lsubzero     ! mask for sub-zero temperatures

    t_high   = merge(c1, c0, (Tin > Tb_liq))
    lsubzero = merge(c1, c0, (Tin <= c0))

    dSbr_dTin = ((L1_liq - J1_liq * K1_liq) / &
                (K1_liq * Tin + L1_liq)**2) * t_high + &
                ((L2_liq - J2_liq * K2_liq) / &
                (K2_liq * Tin + L2_liq)**2) * (c1 - t_high)

    dSbr_dTin = dSbr_dTin * lsubzero

  end function dliquidus_brine_salinity_mush_dT

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

  function denthalpy_brine_dT() result(dqbr_dTin)

    ! derivative of brine enthalpy wrt brine temperature

    real(kind=dbl_kind) :: &
         dqbr_dTin ! derivative of brine enthalpy wrt ice temperature (J m-3 K-1)

    dqbr_dTin = cp_ocn * rhow

  end function denthalpy_brine_dT

!=======================================================================

  function dqdT(Tin, Sin) result(dqin_dTin)
    
    ! derivative of mush enthalpy wrt mush temperature

    real(kind=dbl_kind), intent(in) :: &
         Tin       , & ! ice layer temperature (C)
         Sin           ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         dqin_dTin     ! derivative of ice enthalpy wrt ice temperature (J m-3 K-1)

    real(kind=dbl_kind) :: &
         phi       , & ! liquid fraction
         dphi_dTin , & ! derivative of liquid fraction wrt temperature (C-1)
         l_mush        ! mask for presence of mush
 
    phi = liquid_fraction(Tin, Sin)

    dphi_dTin = dliquid_fraction_dT(Tin, Sin)

    l_mush = merge(c1, c0, (phi < c1))

    dqin_dTin = (cp_ocn * rhow - cp_ice * rhoi) * (phi + Tin * dphi_dTin) + &
                 rhoi * cp_ice + rhoi * Lfresh * dphi_dTin

    dqin_dTin = dqin_dTin * l_mush + cp_ocn * rhow * (c1 - l_mush)

  end function dqdT

!=======================================================================

  function dqdS(Tin, Sin) result(dqin_dSin)

    ! derivative of mush enthalpy wrt mush bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         Tin       , & ! ice layer temperature (C)
         Sin           ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         dqin_dSin     ! derivative of ice enthalpy wrt bulk salinity (J m-3 ppt-1)

    real(kind=dbl_kind) :: &
         phi       , & ! liquid fraction
         dphi_dSin , & ! derivative of liquid fraction wrt temperature (C-1)
         l_mush        ! mask for presence of mush

    phi = liquid_fraction(Tin, Sin)

    dphi_dSin = dliquid_fraction_dS(Tin, Sin)

    l_mush = merge(c1, c0, (phi < c1))

    dqin_dSin = ((cp_ocn * rhow - cp_ice * rhoi) * Tin + &
                 rhoi * Lfresh) * dphi_dSin

    dqin_dSin = dqin_dSin * l_mush

  end function dqdS

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

    Tin = (-B + sqrt(B**2 - c4 * A * C)) / (c2 * A)

    ! change T if all melted
    Tin = q_melt * qin * I_liq + (c1 - q_melt) * Tin

  end function temperature_mush

!=======================================================================

  function heat_conductivity(Tin, Sin) result(km)
    
    ! msuh heat conductivity from mush temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         Tin               , & ! ice layer temperature (C)
         Sin                   ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         km                    ! ice layer conductivity (W m-1 K-1)
    
    real(kind=dbl_kind), parameter :: &
         ki = 2.3_dbl_kind , & ! ice conductivity (W m-1 K-1)
         kb = 0.3_dbl_kind     ! brine conductivity (W m-1 K-1)

    real(kind=dbl_kind) :: &
         phi                   ! liquid fraction

    phi = liquid_fraction(Tin, Sin)

    km = phi * (kb - ki) + ki

  end function heat_conductivity

!=======================================================================
  
  function dheat_conductivity_dT(Tin, Sin) result(dkm_dTin)
    
    ! derivative of mush heat conductivity wrt temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         Tin               , & ! ice layer temperature (C)
         Sin                   ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         dkm_dTin              ! derivative of ice conductivity wrt temperature (W m-1 K-2)
    
    real(kind=dbl_kind), parameter :: &
         ki = 2.3_dbl_kind , & ! ice conductivity (W m-1 K-1)
         kb = 0.3_dbl_kind     ! brine conductivity (W m-1 K-1)

    real(kind=dbl_kind) :: &
         phi               , & ! liquid fraction
         dphi_dTin         , & ! derivative of liquid fraction wrt temperature (C-1)
         l_mush                ! mask for presence of mush

    phi = liquid_fraction(Tin, Sin)

    dphi_dTin = dliquid_fraction_dT(Tin, Sin)

    l_mush = merge(c1, c0, (phi < c1))

    dkm_dTin = (kb - ki) * dphi_dTin * l_mush

  end function dheat_conductivity_dT

  !=======================================================================

  function liquid_fraction(Tin, Sin) result(phi)

    ! liquid fraction of mush from mush temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         Tin , & ! ice layer temperature (C)
         Sin     ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         phi , & ! liquid fraction
         Sbr     ! brine salinity (ppt)

    Sbr = liquidus_brine_salinity_mush(Tin)
    phi = Sin / max(Sbr, Sin)

  end function liquid_fraction

!=======================================================================

  function dliquid_fraction_dT(Tin, Sin) result(dphi_dTin)

    ! derivative of liquid fraction wrt mush temperature

    real(kind=dbl_kind), intent(in) :: &
         Tin       , & ! ice layer temperature (C)
         Sin           ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         dphi_dTin , & ! derivative of liquid fraction wrt temperature (C-1)
         Sbr       , & ! brine salinity (ppt)
         dSbr_dTin , & ! derivative of brine salinity wrt temperature (ppt C-1)
         phi       , & ! liquid fraction
         lmush         ! mask for presence of mush

    ! original
    Sbr       = liquidus_brine_salinity_mush(Tin)
    dSbr_dTin = dliquidus_brine_salinity_mush_dT(Tin)
    phi = Sin / max(Sbr, Sin)
    lmush = merge(c1, c0, (phi < c1))

    dphi_dTin = lmush * (-phi**2 * dSbr_dTin) / Sin

 end function dliquid_fraction_dT

!=======================================================================

  function dliquid_fraction_dS(Tin, Sin) result(dphi_dSin)

    ! derivative of liquid fraction wrt mush bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         Tin       , & ! ice layer temperature (C)
         Sin           ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         dphi_dSin , & ! derivative of liquid fraction wrt bulk salinity (ppt-1)
         Sbr       , & ! brine salinity (ppt)
         dSbr_dTin , & ! derivative of brine salinity wrt temperature (ppt C-1)
         phi       , & ! liquid fraction
         lmush         ! mask for presence of mush

    Sbr   = liquidus_brine_salinity_mush(Tin)
    phi   = Sin / max(Sbr, Sin)
    lmush = merge(c1, c0, (phi < c1))

    dphi_dSin = lmush / max(Sbr, Sin)

  end function dliquid_fraction_dS

!=======================================================================
! Jacobian Free Newton Krylov Solver
!=======================================================================

  subroutine init_jfnk_solver()

    ! initialize parameters used by the GMRES solver within the JFNK solver
    
    ! permanant GMRES parameters
    
    ! set up ipar parameters
    ipar(1) = 0 ! starting condition
#if preconditioning == 1
    ipar(2) = 2 ! status of the preconditioning
#else
    ipar(2) = 0 ! status of the preconditioning
#endif
    !       0 == no preconditioning
    !       1 == left preconditioning only
    !       2 == right preconditioning only
    !       3 == both left and right preconditioning
    
    ! stopping criteria 
    ipar(3) = 1 
    !       1 == || residual || <= rtol * || initial residual || + atol
    !ipar(3) = 999
#if gmres_accurate_solve == 0
    fpar(1) = 0.001_dbl_kind ! the relative tolerance - default 0.001
#else
    fpar(1) = 1.0e-8_dbl_kind
#endif
    fpar(2) = c0   ! the absolute tolerance
    
    ipar(5) = nsize_krylov ! size of the Krylov subspace - 15 is the default
    ipar(6) = ipar(5)*nit_gmres_max! maximum number of matrix-vector multiplies
    
  end subroutine init_jfnk_solver
  
!=======================================================================

  subroutine newton_loop(yn, y0, yscale, fscale, &
                         nyn, solve_type, lstop, &
                         Jac_a, Jac_b, Jac_c, &
                         residual_function, converged_function)

    ! Jacobian Free Newton Krylov (JFNK) solver - this level is the outer Newton loop

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         yn             ! solution vector

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         y0         , & ! initial solution vector
         yscale     , & ! solution scaling vector
         fscale     , & ! residual scaling vector
         Jac_a      , & ! Jacobian lower off-diagonal
         Jac_b      , & ! Jacobian diagonal
         Jac_c          ! Jacobian upper off-diagonal
    
    integer(kind=int_kind), intent(in) :: &
         nyn        , & ! solution vector size
         solve_type     ! 1: nosnow, cold condition
                        ! 2: nosnow, melt condition
                        ! 3: snow, cold condition
                        ! 4: snow, melt condition
    
    logical(kind=log_kind), intent(inout) :: &
         lstop          ! solver failure flag 

    ! function pointer to the non-linear residual function
    interface 
       subroutine residual_function(yn,     y0,     &
                                    yscale, fscale, &
                                    f,      nyn)
         use ice_kinds_mod

         real(kind=dbl_kind), dimension(:), intent(in) :: &
              yn     , & ! solution vector
              y0     , & ! initial solution vector
              yscale , & ! solution scaling vector
              fscale     ! residual scaling vector

         real(kind=dbl_kind), dimension(:), intent(out) :: &
              f          ! residual vector

         integer(kind=int_kind), intent(in) :: &
              nyn        ! solution vector size

       end subroutine residual_function
    end interface
    
    ! function pointer to the non-linear residual function
    interface 
       function converged_function(fk,     f0,         &
                                   yn,     y0,         &
                                   yscale, fscale,     &
                                   nyn,    nit_newton) &
                                   result(lconverged)
         use ice_kinds_mod

         real(kind=dbl_kind), dimension(:), intent(in) :: &
              fk          , & ! residual vector
              f0          , & ! initial residual vector
              yn          , & ! solution vector
              y0          , & ! initial solution vector
              yscale      , & ! solution scaling vector
              fscale          ! residual scaling vector

         integer(kind=int_kind), intent(in) :: &
              nyn        , &  ! solution vector size
              nit_newton      ! Newton iteration counter
         
         logical(kind=log_kind) :: &
              lconverged      ! solver convergence test

       end function converged_function
    end interface

    logical(kind=log_kind) :: &
         lconv_newton ! solver convergence test

    integer(kind=int_kind) :: &
         nit_newton   ! Newton iteration counter
    
    real(kind=dbl_kind), dimension(1:nynmax) :: &
         rhs      , & ! residual vector
         rhs0     , & ! initial residual vector
         delta_yn     ! solution vector change from Newton iteration
    
    nit_newton = 0
    lconv_newton = .false.

    ! initial call to residual function
    call residual_function(yn,         y0,     &
                           yscale,     fscale, & 
                           rhs(1:nyn), nyn)

    rhs0(1:nyn) = rhs(1:nyn)

    ! perform Newton iterations
    newton_do: do while (.not. lconv_newton .and. .not. lstop)
       
       ! run the gmres loop
       nit_newton = nit_newton + 1

       call gmres_loop(yn,         y0,              &
                      -rhs(1:nyn), delta_yn(1:nyn), &
                      yscale,      fscale,          &
                      nyn,         solve_type,      &
                      Jac_a, Jac_b, Jac_c,          &
                      residual_function,            &
                      lstop)

       ! update the solution with the result of the gmres loop
       yn = yn + delta_yn(1:nyn)
       
       ! test convergence
       call residual_function(yn,         y0,     &
                              yscale,     fscale, & 
                              rhs(1:nyn), nyn)
       
       call newton_stopping_condition(rhs(1:nyn),    rhs0(1:nyn),   &
                                      yn(1:nyn),     y0(1:nyn),     &
                                      yscale(1:nyn), fscale(1:nyn), &
                                      nyn,           nit_newton,    &
                                      lconv_newton,  lstop,         &
                                      converged_function)
       
    enddo newton_do

    ! if not converged write out diagnostics
#if diag_nonconverge_always == 0
    if (.not. lconv_newton) then
#endif
       call non_converged_newton_loop(yn,     y0,         &
                                      yscale, fscale,     &
                                      nyn,    solve_type, &
                                      residual_function)
       
#if diag_nonconverge_always == 0
    end if
#endif
    
  end subroutine newton_loop

!=======================================================================

  subroutine newton_stopping_condition(fk,        f0,         & 
                                       yn,        y0,         &
                                       yscale,    fscale,     &
                                       nyn,       nit_newton, &
                                       lconverge, lstop,      &
                                       converged_function)

    ! test stopping condition of the Newton loop in the JFNK solver

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         fk         , & ! residual vector
         f0         , & ! initial residual vector
         yn         , & ! solution vector
         y0         , & ! initial solution vector
         yscale     , & ! solution scaling vector
         fscale         ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn        , & ! solution vector size
         nit_newton     ! Newton iteration counter

    logical(kind=log_kind), intent(inout) :: &
         lconverge  , & ! solver convergence test
         lstop          ! solver failure flag 

    ! function pointer to the non-linear residual function
    interface 
       function converged_function(fk,     f0,         &
                                   yn,     y0,         &
                                   yscale, fscale,     &
                                   nyn,    nit_newton) &
                                   result(lconverged)
         use ice_kinds_mod
 
         real(kind=dbl_kind), dimension(:), intent(in) :: &
              fk         , & ! residual vector
              f0         , & ! initial residual vector
              yn         , & ! solution vector
              y0         , & ! initial solution vector
              yscale     , & ! solution scaling vector
              fscale         ! residual scaling vector

         integer(kind=int_kind), intent(in) :: &
              nyn        , & ! solution vector size
              nit_newton     ! Newton iteration counter
         
         logical(kind=log_kind) :: &
              lconverged     ! solver convergence test

       end function converged_function
    end interface
   
    ! get convergence test
    lconverge = converged_function(fk,     f0,        &
                                   yn,     y0,        &
                                   yscale, fscale,    &
                                   nyn,    nit_newton)
    
    ! never allow newton loop to stop
#if fully_converged == 1
    lconverge = .false.
    if (nit_newton >= nit_newton_max) then
       lconverge = .true.
    endif
    return
#endif
    
    ! return if converged
    if (lconverge) then
       return
    end if
    
    ! check if have exceeded Newton iteration limit
    if (nit_newton >= nit_newton_max) then
       lstop = .true.
       return
    end if
    
  end subroutine newton_stopping_condition

!=======================================================================

  subroutine gmres_loop(yn,     y0,          &
                        rhs,    sol,         &
                        yscale, fscale,      &
                        nyn,    solve_type,  &
                        Jac_a, Jac_b, Jac_c, &
                        residual_function,   &
                        lstop)

    ! perform a solution with GMRES of the linear problem in JFNK

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         rhs        , & ! residual vector to form right hand side for GMRES solution
         yn         , & ! solution vector
         y0         , & ! initial solution vector
         yscale     , & ! solution scaling vector
         fscale     , & ! residual scaling vector
         Jac_a      , & ! Jacobian lower off-diagonal
         Jac_b      , & ! Jacobian diagonal
         Jac_c          ! Jacobian upper off-diagonal

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         sol            ! solution of the inner GMRES iteration

    integer(kind=int_kind), intent(in) :: &
         nyn        , & ! solution vector size
         solve_type     ! 1: nosnow, cold condition
                        ! 2: nosnow, melt condition
                        ! 3: snow, cold condition
                        ! 4: snow, melt condition

    logical(kind=log_kind), intent(inout) :: &
         lstop          ! solver failure flag 
    
    
    ! function pointer to the non-linear residual function
    interface 
       subroutine residual_function(yn,     y0,     &
                                    yscale, fscale, &
                                    f,      nyn)
         use ice_kinds_mod

         real(kind=dbl_kind), dimension(:), intent(in) :: &
              yn     , & ! solution vector
              y0     , & ! initial solution vector
              yscale , & ! solution scaling vector
              fscale     ! residual scaling vector

         real(kind=dbl_kind), dimension(:), intent(out) :: &
              f          ! residual vector

         integer(kind=int_kind), intent(in) :: &
              nyn        ! solution vector size

       end subroutine residual_function
    end interface
   
    integer(kind=int_kind) :: &
         wsize ! size of the GMRES work space vector

    real(kind=dbl_kind), dimension(1:wsizemax) :: &
         w     ! work space for GMRES
    
    ! required storage space
    wsize = (nyn+3)*(nsize_krylov+2) + (nsize_krylov+1)*nsize_krylov/2 ! size of the workspace
    
    ! workspace size
    ipar(4) = wsize ! number of elements in the array 'w'
    
    ! set up fpar parameters
    fpar(11) = c0 ! number of floating-point operations
    
    sol = c0
    
    ipar(1) = 0

    ! loop until GMRES is finished - reverse communication protocol
    gmres_do: do

       !call gmres(nyn, rhs, sol, ipar, fpar, w(1:wsize))
       call gmres_fast(nyn, rhs, sol, ipar, fpar, w(1:wsize))
       
       select case (ipar(1))
       case (0)

          ! gmres finished
          return

       case (-3)

          ! return due to anticipated break-down
          if (ipar(12) == -3) then
             ! accidental exact solution
             return
          else
             ! error in GMRES
             exit
          endif

       case (1)
          
          ! perform an approximated matrix multiplication
          call gmres_matvec_multiply(w(ipar(8):ipar(8)+nyn-1), &
                                     w(ipar(9):ipar(9)+nyn-1), &
                                     yn,          y0,          &
                                     yscale,      fscale,      &
                                     -rhs(1:nyn), nyn,         &
                                     residual_function)
          
       case (5)
          
          ! perform a right preconditioning multiplication
          call gmres_right_precond_multiply(w(ipar(8):ipar(8)+nyn-1), &
                                            w(ipar(9):ipar(9)+nyn-1), &
                                            Jac_a, Jac_b, Jac_c, &
                                            nyn)
          
       case default
          
          ! error condition
          exit
          
       end select
       
    enddo gmres_do

    ! failure of gmres
    lstop = .true. 
    call gmres_fail(ipar, fpar)
    
  end subroutine gmres_loop

!=======================================================================

  subroutine gmres_matvec_multiply(u,      v,      &
                                   yn,     y0,     &
                                   yscale, fscale, &
                                   f2,     nyn,    &
                                   residual_function)

    ! perform an approximated matrix-vector multiply for GMRES

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         u      , & ! input vector of matrix multiplication
         yn     , & ! solution vector
         y0     , & ! initial solution vector
         f2     , & ! first residual vector for residual difference
         yscale , & ! solution scaling vector
         fscale     ! residual scaling vector

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         v          ! output vector of matrix multiplication

    integer(kind=int_kind), intent(in) :: &
         nyn        ! solution vector size
    
    ! function pointer to the non-linear residual function
    interface 
       subroutine residual_function(yn,     y0,     &
                                    yscale, fscale, &
                                    f, nyn)
         use ice_kinds_mod

         real(kind=dbl_kind), dimension(:), intent(in) :: &
              yn     , & ! solution vector
              y0     , & ! initial solution vector
              yscale , & ! solution scaling vector
              fscale     ! residual scaling vector

         real(kind=dbl_kind), dimension(:), intent(out) :: &
              f          ! residual vector

         integer(kind=int_kind), intent(in) :: &
              nyn        ! solution vector size

       end subroutine residual_function
    end interface

    real(kind=dbl_kind) :: &
         epsilon , & ! matrix multiplication approximation perturbation
         u_norm      ! L2 norm of the input vector for matrix multiplication

    real(kind=dbl_kind), dimension(1:nynmax) :: &
         y       , & ! solution vector for second residual vector
         f1          ! second residual vector for residual difference
    
    u_norm = L2_norm(u, nyn)
    
    if (u_norm > c0) then ! check v is not zero
       
       epsilon = calc_epsilon(u_norm, yn, nyn)
       
       y(1:nyn) = yn + epsilon * u
       
       call residual_function(y(1:nyn),      y0,            &
                              yscale(1:nyn), fscale(1:nyn), &
                              f1(1:nyn),     nyn)                 
       
       v = (f1(1:nyn) - f2(1:nyn)) / epsilon
       
    else
       
       v = c0
       
    endif
    
  end subroutine gmres_matvec_multiply

!=======================================================================

  function calc_epsilon(u_norm, yn, nyn) result(epsilon)

    ! calculate the matrix-vector multiplication approximation perturbation

    real(kind=dbl_kind), intent(in) :: &
         u_norm  ! L2 norm of the input vector for matrix multiplication

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         yn      ! solution vector

    integer(kind=int_kind), intent(in) :: &
         nyn     ! solution vector size

    real(kind=dbl_kind) :: &
         epsilon ! matrix multiplication approximation perturbation

    integer(kind=int_kind) :: &
         l       ! solution vector index

    ! few orders of magnitude bigger than square root of machine precision 
    real(kind=dbl_kind), parameter :: &
         sqrtmpb = sqrt(c2**(-53)) * 100.0_dbl_kind
    
    epsilon = c0
    
    do l = 1, nyn
       epsilon = epsilon + (abs(yn(l)) + c1) * sqrtmpb
    enddo
    
    epsilon = epsilon / (real(nyn, dbl_kind) * u_norm)
    
  end function calc_epsilon

!=======================================================================

  function L2_norm(v, n) result(L2)

    ! calculate the L2 norm of a vector

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         v  ! input vector

    integer(kind=int_kind), intent(in) :: &
         n  ! input vector size

    real(kind=dbl_kind) :: &
         L2 ! L2 norm

    integer(kind=int_kind) :: &
         l  ! solution vector index
    
    L2 = c0
    do l = 1, n
       L2 = L2 + v(l)**2
    enddo
    L2 = sqrt(L2)
    
  end function L2_norm

!=======================================================================

  subroutine gmres_fail(ipar, fpar)

    ! error messages upon failure of GMRES

    integer(kind=int_kind), dimension(16), intent(in) :: &
         ipar        ! GMRES integer input/output parameters

    real(kind=dbl_kind), dimension(16), intent(in) :: &
         fpar        ! GMRES real input/output parameters

    character(len=200) :: &
         gmres_error ! GMRES error string

    write(nu_diag,*) "GMRES failure"
    write(nu_diag,*) "-------------"
    write(nu_diag,*) "ipar(1) : status of the call/return       :", ipar(1)
    write(nu_diag,*) "ipar(2) : status of preconditioning       :", ipar(2)
    write(nu_diag,*) "ipar(3) : stopping criteria               :", ipar(3)
    write(nu_diag,*) "ipar(4) : w array size                    :", ipar(4)
    write(nu_diag,*) "ipar(5) : size of the Krylov subspace     :", ipar(5)
    write(nu_diag,*) "ipar(6) : maximum no. of matrix multiplies:", ipar(6)
    write(nu_diag,*) "ipar(7) : current no. of matrix multiplies:", ipar(7)
    write(nu_diag,*) "ipar(12): error code of MGSRO             :", ipar(12)
    write(nu_diag,*) "ipar(13): number of initializations       :", ipar(13)
    write(nu_diag,*) ""
    write(nu_diag,*) "fpar(1) : the relative tolerance             :", fpar(1)
    write(nu_diag,*) "fpar(2) : the absolute tolerance             :", fpar(2)
    write(nu_diag,*) "fpar(3) : initial residual/error norm        :", fpar(3)
    write(nu_diag,*) "fpar(4) : target residual/error norm         :", fpar(4)
    write(nu_diag,*) "fpar(5) : current residual norm              :", fpar(5)
    write(nu_diag,*) "fpar(6) : current residual/error norm        :", fpar(6)
    write(nu_diag,*) "fpar(7) : convergence rate                   :", fpar(7)
    write(nu_diag,*) "fpar(11): number of floating-point operations:", fpar(11)
    write(nu_diag,*) ""

    ! write out gmres error message
    call gmres_error_message(ipar(1), gmres_error)
    write(nu_diag,*) trim(gmres_error)
    if (ipar(1) == -3) then
       call mgsro_error_message(ipar(12), gmres_error)
       write(nu_diag,*) trim(gmres_error)
    endif
    write(nu_diag,*) ""

  end subroutine gmres_fail

!=======================================================================

  subroutine gmres_error_message(ipar1, gmres_error)

    ! GMRES error messages from error code

    integer(kind=int_kind), intent(in) :: &
         ipar1       ! first GMRES integer input/output parameter (ipar(1))

    character(len=200), intent(out) :: &
         gmres_error ! GMRES error string

    character(len=5) :: &
         strerror    ! ipar(1) as a string
    
    select case (ipar1)
    case (2) ! request a matvec with A^T
       gmres_error = "gmres: 2: gmres should not have requested a matvec with A^T"
    case (3) ! request a left preconditioner solve (Ml^{-1})
       gmres_error = "gmres: 3: gmres should not have requested a left preconditioner solve (Ml^{-1})"
    case (4) ! request a left preconditioner transposed solve (Ml^{-T})
       gmres_error = "gmres: 4: gmres should not have requested a left preconditioner transposed solve (Ml^{-T})"
    case (6) ! request a right preconditioner transposed solve (Mr^{-T})
       gmres_error = "gmres: 6: gmres should not have requested a right preconditioner transposed solve (Mr^{-T})"
    case (10) ! request the caller to perform stopping test
       gmres_error = "gmres: 10: gmres should not have requested the caller to perform stopping test"
    case (-1) ! termination because iteration number is greater than the preset limit       
       gmres_error = "gmres: -1: termination because iteration number is greater than the preset limit"
    case (-2) ! return due to insufficient work space
       gmres_error = "gmres: -2: return due to insufficient work space"
    case (-3) ! return due to anticipated break-down / divide by zero, in the case where Arnoldi procedure is used
       gmres_error = "gmres: -3: return due to anticipated break-down: ipar(12):"
    case (-4) ! the values of fpar(1) and fpar(2) are both <= 0, the valid ranges are 0 <= fpar(1) < 1, 0 <= fpar(2), and they cannot be zero at the same time
       gmres_error = "gmres: -4: error with values of fpar(1) and fpar(2)"
    case (-9) ! while trying to detect a break-down, an abnormal number is detected.
       gmres_error = "gmres: -9: while trying to detect a break-down, an abnormal number is detected."
    case (-10) ! return due to some non-numerical reasons, e.g. invalid floating-point numbers etc.
       gmres_error = "gmres: -10: return due to some non-numerical reasons, e.g. invalid floating-point numbers etc."
    case default
       write(strerror,fmt='(i3)') ipar1
       gmres_error = "gmres: unknown error state: "//trim(strerror)
    end select
    
  end subroutine gmres_error_message

!=======================================================================

  subroutine mgsro_error_message(ipar12, mgsro_error)

    ! MGSRO error messages from code

    integer(kind=int_kind), intent(in) :: &
         ipar12      ! twelth GMRES integer input/output parameter (ipar(12))

    character(len=200), intent(out) :: &
         mgsro_error ! MGSRO error string

    character(len=5) :: &
         strerror    ! ipar(12) as a string

    select case (ipar12)
    case(-1)
       mgsro_error = "-1: zero input vector"
    case(-2)
       mgsro_error = "-2: input vector contains abnormal numbers"
    case(-3)
       mgsro_error = "-3: input vector is a linear combination of others"
    case(-4)
       mgsro_error = "-4: trianguler system in GMRES/FOM/etc. has nul rank"
    case default
       write(strerror,fmt='(i3)') ipar12
       mgsro_error = "msgro: unknown error state: "//trim(strerror)
    end select

  end subroutine mgsro_error_message

!=======================================================================
! preconditioning stuff
!=======================================================================

  subroutine gmres_right_precond_multiply(u, v, Jac_a, Jac_b, Jac_c, nyn)

    ! perform a matrix multiply with right preconditioner

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         u     , & ! input vector
         Jac_a , & ! Jacobian lower off-diagonal
         Jac_b , & ! Jacobian diagonal
         Jac_c     ! Jacobian upper off-diagonal

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         v         ! output vector

    integer(kind=int_kind), intent(in) :: &
         nyn       ! solution vector size

    call tdma_solve_sparse(Jac_a(1:nyn), Jac_b(1:nyn), Jac_c(1:nyn), &
                           u(1:nyn),     v(1:nyn),     nyn)
    
  end subroutine gmres_right_precond_multiply
  
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
    
    real(kind=dbl_kind), dimension(nynmax) :: &
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
! Functionality for non-converging situations
!=======================================================================

  subroutine non_converged_newton_loop(yn,     y0,         &
                                       yscale, fscale,     &
                                       nyn,    solve_type, &
                                       residual_function)
    
    ! diagnostics if JFNK fails to find a solution

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         yn         , & ! solution vector
         y0         , & ! initial solution vector
         yscale     , & ! solution scaling vector
         fscale         ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn        , & ! solution vector size
         solve_type     ! 1: nosnow, cold condition
                        ! 2: nosnow, melt condition
                        ! 3: snow, cold condition
                        ! 4: snow, melt condition

    ! function pointer to the non-linear residual function
    interface 
       subroutine residual_function(yn,     y0,     &
                                    yscale, fscale, &
                                    f,      nyn)
         use ice_kinds_mod

         real(kind=dbl_kind), dimension(:), intent(in) :: &
              yn     , & ! solution vector
              y0     , & ! initial solution vector
              yscale , & ! solution scaling vector
              fscale     ! residual scaling vector

         real(kind=dbl_kind), dimension(:), intent(out) :: &
              f          ! residual vector

         integer(kind=int_kind), intent(in) :: &
              nyn        ! solution vector size

       end subroutine residual_function
    end interface
   
    real(kind=dbl_kind), dimension(1:nynmax) :: &
         rhs ! residual vector

    integer(kind=int_kind) :: &
         l   ! solution vector index

    ! basic stopping information
    write(nu_diag,*) "Newton loop not converged!"
    write(nu_diag,*) "istep1:  ", istep1
    write(nu_diag,*) "my_task: ", my_task
    write(nu_diag,*) "g_n:     ", g_n
    write(nu_diag,*) "g_i:     ", g_i
    write(nu_diag,*) "g_j:     ", g_j
    write(nu_diag,*) "nyn:     ", nyn
    write(nu_diag,*)
    
    ! write out global physical prameters passed to the residual function
    call diag_global_variables()
    
    ! write out the physical state from the current and initial solution vectors
    call writeout_physical_state(yn, y0, yscale, fscale, &
                                 nyn, solve_type)

    ! write out the final residual
    write(nu_diag,*) "Final residual: istep1: ", istep1
    call residual_function(yn,         y0,     &
                           yscale,     fscale, &
                           rhs(1:nyn), nyn)
    do l = 1, nyn
       write(nu_diag,*) l, rhs(l), yn(l), y0(l)
    enddo ! l
    write(nu_diag,*)
    
    ! write out the initial residual
    write(nu_diag,*) "Initial residual:"
    call residual_function(y0,         y0,     &
                           yscale,     fscale, &
                           rhs(1:nyn), nyn)
    do l = 1, nyn
       write(nu_diag,*) l, rhs(l), y0(l), y0(l)
    enddo ! l
    write(nu_diag,*)

  end subroutine non_converged_newton_loop

!=======================================================================

  subroutine writeout_physical_state(yn, y0, yscale, fscale, &
                                     nyn, solve_type)

    ! given the state variables write out the physical state of the system

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         yn          , & ! solution vector
         y0          , & ! initial solution vector
         yscale      , & ! solution scaling vector
         fscale          ! residual scaling vector

    integer(kind=int_kind), intent(in) :: &
         nyn        , &  ! solution vector size
         solve_type      ! 1: nosnow, cold condition
                         ! 2: nosnow, melt condition
                         ! 3: snow, cold condition
                         ! 4: snow, melt condition

    real(kind=dbl_kind) :: &
         Tsf         , & ! ice/snow surface temperature (C)
         Tsf0            ! ice/snow surface temperature (C) at beginning of time step

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         qsn         , & ! snow layer enthalpy (J m-3)
         qsn0        , & ! snow layer enthalpy (J m-3) at beginning of time step
         Tsn         , & ! snow layer temperature (C)
         Tsn0            ! snow layer temperature (C) at beginning of time step

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         qin         , & ! ice layer enthalpy (J m-3)
         qin0        , & ! ice layer enthalpy (J m-3) at beginning of time step
         Tin         , & ! ice layer temperature (C)
         Tin0        , & ! ice layer temperature (C) at beginning of time step
         Sin         , & ! ice layer bulk salinity (ppt)
         Sin0        , & ! ice layer bulk salinity (ppt) at beginning of time step
         Tmlt        , & ! ice layer melting temperature (C)
         Tmlt0           ! ice layer melting temperature (C) at beginning of time step

    integer(kind=int_kind) :: &
         k               ! ice/snow layer index

    logical(kind=log_kind) :: &
         lsnow           ! flag for snow presence
       
    if (solve_type == 1) then

       ! calculate physical state for cold conditions and no snow
       call state_to_physical_nosnow_cold(yn, yscale, &
                                          Tsf,  qin,  Sin)
       call state_to_physical_nosnow_cold(y0, yscale, &
                                          Tsf0, qin0, Sin0)

       do k = 1, nilyr
          Tin(k)   = temperature_mush(qin(k), Sin(k))
          Tin0(k)  = temperature_mush(qin0(k), Sin0(k))
          Tmlt(k)  = liquidus_temperature_mush(Sin(k))
          Tmlt0(k) = liquidus_temperature_mush(Sin0(k))
       enddo ! k

       lsnow = .false.

    else if (solve_type == 2) then

       ! calculate physical state for melting conditions and no snow
       call state_to_physical_nosnow_melt(yn, yscale, &
                                          qin,  Sin)
       call state_to_physical_nosnow_melt(y0, yscale, &
                                          qin0, Sin0)

       Tsf = -999.9_dbl_kind ; Tsf0 = -999.9_dbl_kind

       do k = 1, nilyr
          Tin(k)   = temperature_mush(qin(k), Sin(k))
          Tin0(k)  = temperature_mush(qin0(k), Sin0(k))
          Tmlt(k)  = liquidus_temperature_mush(Sin(k))
          Tmlt0(k) = liquidus_temperature_mush(Sin0(k))
       enddo ! k

       lsnow = .false.

    else if (solve_type == 3) then

       ! calculate physical state for cold conditions and snow present
       call state_to_physical_snow_cold(yn, yscale, &
                                        Tsf,  qsn,  qin,  Sin)
       call state_to_physical_snow_cold(y0, yscale, &
                                        Tsf0, qsn0, qin0, Sin0)

       do k = 1, nslyr
          Tsn(k)  = temperature_snow(qsn(k))
          Tsn0(k) = temperature_snow(qsn0(k))
       enddo ! k

       do k = 1, nilyr
          Tin(k)   = temperature_mush(qin(k), Sin(k))
          Tin0(k)  = temperature_mush(qin0(k), Sin0(k))
          Tmlt(k)  = liquidus_temperature_mush(Sin(k))
          Tmlt0(k) = liquidus_temperature_mush(Sin0(k))
       enddo ! k

       lsnow = .true.

    else if (solve_type == 4) then

       ! calculate physical state for melting conditions and snow present
       call state_to_physical_snow_melt(yn, yscale, &
                                        qsn,  qin,  Sin)
       call state_to_physical_snow_melt(y0, yscale, &
                                        qsn0, qin0, Sin0)
       Tsf = -999.9_dbl_kind   ; Tsf0 = -999.9_dbl_kind

       do k = 1, nslyr
          Tsn(k)  = temperature_snow(qsn(k))
          Tsn0(k) = temperature_snow(qsn0(k))
       enddo ! k

       do k = 1, nilyr
          Tin(k)   = temperature_mush(qin(k), Sin(k))
          Tin0(k)  = temperature_mush(qin0(k), Sin0(k))
          Tmlt(k)  = liquidus_temperature_mush(Sin(k))
          Tmlt0(k) = liquidus_temperature_mush(Sin0(k))
       enddo ! k

       lsnow = .true.

    endif
   
    ! write out physical state
    write(nu_diag,*) "-----------------------------------"
    write(nu_diag,*) "Physical state: ", solve_type
    write(nu_diag,*) "-----------------------------------"
    write(nu_diag,*)
    write(nu_diag,*) "Tsf  : ", 0, Tsf, Tsf0
    write(nu_diag,*)
    if (lsnow) then
       do k = 1, nslyr
          write(nu_diag,*) "Snow : ", k, Tsn(k), Tsn0(k), &
                                         qsn(k), qsn0(k)
       enddo ! k
       write(nu_diag,*)
    endif ! lsnow
    do k = 1, nilyr
       write(nu_diag,*) "Ice  : ", k, Tin(k),   Tin0(k),  &
                                      qin(k),   qin0(k),  &
                                      Sin(k),   Sin0(k),  &
                                      Tmlt(k),  Tmlt0(k), &
                                      Tin(k)  - Tmlt(k),  &
                                      Tin0(k) - Tmlt0(k)
    enddo ! k
    write(nu_diag,*)

    ! write out the scaling factors used
    write(nu_diag,*) "-----------------------------------"
    write(nu_diag,*) "Scalings: k, yscale, fscale"
    write(nu_diag,*) "-----------------------------------"
    write(nu_diag,*)

    if (solve_type == 1 .or. solve_type == 3) then
       
       do k = 1, nyn
          write(nu_diag,*) k, yscale(k), fscale(k)
       enddo ! k

    else if (solve_type == 2 .or. solve_type == 4) then

       do k = 1, nyn
          write(nu_diag,*) k, yscale(k), fscale(k)
       enddo ! k

    endif

    write(nu_diag,*)

  end subroutine writeout_physical_state

!=======================================================================

  subroutine diag_global_variables()
        
    ! write out the global variables used by the residual functions

    write(nu_diag,*) "Global variables: "
    write(nu_diag,*) "g_dt:      ", g_dt     
    write(nu_diag,*) "g_hin:     ", g_hin    
    write(nu_diag,*) "g_hilyr:   ", g_hilyr  
    write(nu_diag,*) "g_hslyr:   ", g_hslyr  
    write(nu_diag,*) "g_dti:     ", g_dti    
    write(nu_diag,*) "g_hilyri:  ", g_hilyri 
    write(nu_diag,*) "g_hslyri:  ", g_hslyri 
    write(nu_diag,*) "g_hilyri2: ", g_hilyri2
    write(nu_diag,*) "g_hslyri2: ", g_hslyri2
    write(nu_diag,*) "g_fswsfc:  ", g_fswsfc 
    write(nu_diag,*) "g_rhoa:    ", g_rhoa   
    write(nu_diag,*) "g_flw:     ", g_flw    
    write(nu_diag,*) "g_potT:    ", g_potT   
    write(nu_diag,*) "g_Qa:      ", g_Qa     
    write(nu_diag,*) "g_shcoef:  ", g_shcoef 
    write(nu_diag,*) "g_lhcoef:  ", g_lhcoef 
    write(nu_diag,*) "g_Tbot:    ", g_Tbot   
    write(nu_diag,*) "g_hpond:   ", g_hpond  
    write(nu_diag,*) "g_apond:   ", g_apond  
    write(nu_diag,*) "g_sss:     ", g_sss    
    write(nu_diag,*) "g_qocn:    ", g_qocn   
    write(nu_diag,*)
    
  end subroutine diag_global_variables

!=======================================================================

!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!         Basic Iterative Solvers with Reverse Communication           c
!----------------------------------------------------------------------c
!     This file currently has several basic iterative linear system    c
!     solvers. They are:                                               c
!     CG       -- Conjugate Gradient Method                            c
!     CGNR     -- Conjugate Gradient Method (Normal Residual equation) c
!     BCG      -- Bi-Conjugate Gradient Method                         c
!     DBCG     -- BCG with partial pivoting                            c
!     BCGSTAB  -- BCG stabilized                                       c
!     TFQMR    -- Transpose-Free Quasi-Minimum Residual method         c
!     FOM      -- Full Orthogonalization Method                        c
!     GMRES    -- Generalized Minimum RESidual method                  c
!     FGMRES   -- Flexible version of Generalized Minimum              c
!                 RESidual method                                      c
!     DQGMRES  -- Direct versions of Quasi Generalize Minimum          c
!                 Residual method                                      c
!----------------------------------------------------------------------c
!     They all have the following calling sequence:
!      subroutine solver(n, rhs, sol, ipar, fpar, w)
!      integer n, ipar(16)
!      real*8 rhs(n), sol(n), fpar(16), w(*)
!     Where
!     (1) 'n' is the size of the linear system,
!     (2) 'rhs' is the right-hand side of the linear system,
!     (3) 'sol' is the solution to the linear system,
!     (4) 'ipar' is an integer parameter array for the reverse
!     communication protocol,
!     (5) 'fpar' is an floating-point parameter array storing
!     information to and from the iterative solvers.
!     (6) 'w' is the work space (size is specified in ipar)
!
!     They are preconditioned iterative solvers with reverse
!     communication. The preconditioners can be applied from either
!     from left or right or both (specified by ipar(2), see below).
!
!     Author: Kesheng John Wu (kewu@mail.cs.umn.edu) 1993
!
!     NOTES:
!
!     (1) Work space required by each of the iterative solver
!     routines is as follows:
!       CG      == 5 * n
!       CGNR    == 5 * n
!       BCG     == 7 * n
!       DBCG    == 11 * n
!       BCGSTAB == 8 * n
!       TFQMR   == 11 * n
!       FOM     == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
!       GMRES   == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
!       FGMRES  == 2*n*(m+1) + (m+1)*m/2 + 3*m + 2 (m = ipar(5),
!                  default m=15)
!       DQGMRES == n + lb * (2*n+4) (lb=ipar(5)+1, default lb = 16)
!
!     (2) ALL iterative solvers require a user-supplied DOT-product
!     routine named DISTDOT. The prototype of DISTDOT is
!
!     real(kind=dbl_kind) function distdot(n,x,ix,y,iy)
!     integer n, ix, iy
!     real(kind=dbl_kind) x(1+(n-1)*ix), y(1+(n-1)*iy)
!
!     This interface of DISTDOT is exactly the same as that of
!     DDOT (or SDOT if real == real(kind=dbl_kind)) from BLAS-1. It should have
!     same functionality as DDOT on a single processor machine. On a
!     parallel/distributed environment, each processor can perform
!     DDOT on the data it has, then perform a summation on all the
!     partial results.
!
!     (3) To use this set of routines under SPMD/MIMD program paradigm,
!     several things are to be noted: (a) 'n' should be the number of
!     vector elements of 'rhs' that is present on the local processor.
!     (b) if RHS(i) is on processor j, it is expected that SOL(i)
!     will be on the same processor, i.e. the vectors are distributed
!     to each processor in the same way. (c) the preconditioning and
!     stopping criteria specifications have to be the same on all
!     processor involved, ipar and fpar have to be the same on each
!     processor. (d) DISTDOT should be replaced by a distributed
!     dot-product function.
!
!     ..................................................................
!     Reverse Communication Protocols
!
!     When a reverse-communication routine returns, it could be either
!     that the routine has terminated or it simply requires the caller
!     to perform one matrix-vector multiplication. The possible matrices
!     that involve in the matrix-vector multiplications are:
!     A       (the matrix of the linear system),
!     A^T     (A transposed),
!     Ml^{-1} (inverse of the left preconditioner),
!     Ml^{-T} (inverse of the left preconditioner transposed),
!     Mr^{-1} (inverse of the right preconditioner),
!     Mr^{-T} (inverse of the right preconditioner transposed).
!     For all the matrix vector multiplication, v = A u. The input and
!     output vectors are supposed to be part of the work space 'w', and
!     the starting positions of them are stored in ipar(8:9), see below.
!
!     The array 'ipar' is used to store the information about the solver.
!     Here is the list of what each element represents:
!
!     ipar(1) -- status of the call/return.
!     A call to the solver with ipar(1) == 0 will initialize the
!     iterative solver. On return from the iterative solver, ipar(1)
!     carries the status flag which indicates the condition of the
!     return. The status information is divided into two categories,
!     (1) a positive value indicates the solver requires a matrix-vector
!     multiplication,
!     (2) a non-positive value indicates termination of the solver.
!     Here is the current definition:
!       1 == request a matvec with A,
!       2 == request a matvec with A^T,
!       3 == request a left preconditioner solve (Ml^{-1}),
!       4 == request a left preconditioner transposed solve (Ml^{-T}),
!       5 == request a right preconditioner solve (Mr^{-1}),
!       6 == request a right preconditioner transposed solve (Mr^{-T}),
!      10 == request the caller to perform stopping test,
!       0 == normal termination of the solver, satisfied the stopping
!            criteria,
!      -1 == termination because iteration number is greater than the
!            preset limit,
!      -2 == return due to insufficient work space,
!      -3 == return due to anticipated break-down / divide by zero,
!            in the case where Arnoldi procedure is used, additional
!            error code can be found in ipar(12), where ipar(12) is
!            the error code of orthogonalization procedure MGSRO:
!               -1: zero input vector
!               -2: input vector contains abnormal numbers
!               -3: input vector is a linear combination of others
!               -4: trianguler system in GMRES/FOM/etc. has nul rank
!      -4 == the values of fpar(1) and fpar(2) are both <= 0, the valid
!            ranges are 0 <= fpar(1) < 1, 0 <= fpar(2), and they can
!            not be zero at the same time
!      -9 == while trying to detect a break-down, an abnormal number is
!            detected.
!     -10 == return due to some non-numerical reasons, e.g. invalid
!            floating-point numbers etc.
!
!     ipar(2) -- status of the preconditioning:
!       0 == no preconditioning
!       1 == left preconditioning only
!       2 == right preconditioning only
!       3 == both left and right preconditioning
!
!     ipar(3) -- stopping criteria (details of this will be
!     discussed later).
!
!     ipar(4) -- number of elements in the array 'w'. if this is less
!     than the desired size, it will be over-written with the minimum
!     requirement. In which case the status flag ipar(1) = -2.
!
!     ipar(5) -- size of the Krylov subspace (used by GMRES and its
!     variants), e.g. GMRES(ipar(5)), FGMRES(ipar(5)),
!     DQGMRES(ipar(5)).
!
!     ipar(6) -- maximum number of matrix-vector multiplies, if not a
!     positive number the iterative solver will run till convergence
!     test is satisfied.
!
!     ipar(7) -- current number of matrix-vector multiplies. It is
!     incremented after each matrix-vector multiplication. If there
!     is preconditioning, the counter is incremented after the
!     preconditioning associated with each matrix-vector multiplication.
!
!     ipar(8) -- pointer to the input vector to the requested matrix-
!     vector multiplication.
!
!     ipar(9) -- pointer to the output vector of the requested matrix-
!     vector multiplication.
!
!     To perform v = A * u, it is assumed that u is w(ipar(8):ipar(8)+n-1)
!     and v is stored as w(ipar(9):ipar(9)+n-1).
!
!     ipar(10) -- the return address (used to determine where to go to
!     inside the iterative solvers after the caller has performed the
!     requested services).
!
!     ipar(11) -- the result of the external convergence test
!     On final return from the iterative solvers, this value
!     will be reflected by ipar(1) = 0 (details discussed later)
!
!     ipar(12) -- error code of MGSRO, it is
!                  1 if the input vector to MGSRO is linear combination
!                    of others,
!                  0 if MGSRO was successful,
!                 -1 if the input vector to MGSRO is zero,
!                 -2 if the input vector contains invalid number.
!
!     ipar(13) -- number of initializations. During each initilization
!                 residual norm is computed directly from M_l(b - A x).
!
!     ipar(14) to ipar(16) are NOT defined, they are NOT USED by
!     any iterative solver at this time.
!
!     Information about the error and tolerance are stored in the array
!     FPAR. So are some internal variables that need to be saved from
!     one iteration to the next one. Since the internal variables are
!     not the same for each routine, we only define the common ones.
!
!     The first two are input parameters:
!     fpar(1) -- the relative tolerance,
!     fpar(2) -- the absolute tolerance (details discussed later),
!
!     When the iterative solver terminates,
!     fpar(3) -- initial residual/error norm,
!     fpar(4) -- target residual/error norm,
!     fpar(5) -- current residual norm (if available),
!     fpar(6) -- current residual/error norm,
!     fpar(7) -- convergence rate,
!
!     fpar(8:10) are used by some of the iterative solvers to save some
!     internal information.
!
!     fpar(11) -- number of floating-point operations. The iterative
!     solvers will add the number of FLOPS they used to this variable,
!     but they do NOT initialize it, nor add the number of FLOPS due to
!     matrix-vector multiplications (since matvec is outside of the
!     iterative solvers). To insure the correct FLOPS count, the
!     caller should set fpar(11) = 0 before invoking the iterative
!     solvers and account for the number of FLOPS from matrix-vector
!     multiplications and preconditioners.
!
!     fpar(12:16) are not used in current implementation.
!
!     Whether the content of fpar(3), fpar(4) and fpar(6) are residual
!     norms or error norms depends on ipar(3). If the requested
!     convergence test is based on the residual norm, they will be
!     residual norms. If the caller want to test convergence based the
!     error norms (estimated by the norm of the modifications applied
!     to the approximate solution), they will be error norms.
!     Convergence rate is defined by (Fortran 77 statement)
!     fpar(7) = log10(fpar(3) / fpar(6)) / (ipar(7)-ipar(13))
!     If fpar(7) = 0.5, it means that approximately every 2 (= 1/0.5)
!     steps the residual/error norm decrease by a factor of 10.
!
!     ..................................................................
!     Stopping criteria,
!
!     An iterative solver may be terminated due to (1) satisfying
!     convergence test; (2) exceeding iteration limit; (3) insufficient
!     work space; (4) break-down. Checking of the work space is
!     only done in the initialization stage, i.e. when it is called with
!     ipar(1) == 0. A complete convergence test is done after each
!     update of the solutions. Other conditions are monitored
!     continuously.
!
!     With regard to the number of iteration, when ipar(6) is positive,
!     the current iteration number will be checked against it. If
!     current iteration number is greater the ipar(6) than the solver
!     will return with status -1. If ipar(6) is not positive, the
!     iteration will continue until convergence test is satisfied.
!
!     Two things may be used in the convergence tests, one is the
!     residual 2-norm, the other one is 2-norm of the change in the
!     approximate solution. The residual and the change in approximate
!     solution are from the preconditioned system (if preconditioning
!     is applied). The DQGMRES and TFQMR use two estimates for the
!     residual norms. The estimates are not accurate, but they are
!     acceptable in most of the cases. Generally speaking, the error
!     of the TFQMRs estimate is less accurate.
!
!     The convergence test type is indicated by ipar(3). There are four
!     type convergence tests: (1) tests based on the residual norm;
!     (2) tests based on change in approximate solution; (3) caller
!     does not care, the solver choose one from above two on its own;
!     (4) caller will perform the test, the solver should simply continue.
!     Here is the complete definition:
!      -2 == || dx(i) || <= rtol * || rhs || + atol
!      -1 == || dx(i) || <= rtol * || dx(1) || + atol
!       0 == solver will choose test 1 (next)
!       1 == || residual || <= rtol * || initial residual || + atol
!       2 == || residual || <= rtol * || rhs || + atol
!     999 == caller will perform the test
!     where dx(i) denote the change in the solution at the ith update.
!     ||.|| denotes 2-norm. rtol = fpar(1) and atol = fpar(2).
!
!     If the caller is to perform the convergence test, the outcome
!     should be stored in ipar(11).
!     ipar(11) = 0 -- failed the convergence test, iterative solver
!     should continue
!     ipar(11) = 1 -- satisfied convergence test, iterative solver
!     should perform the clean up job and stop.
!
!     Upon return with ipar(1) = 10,
!     ipar(8)  points to the starting position of the change in
!              solution Sx, where the actual solution of the step is
!              x_j = x_0 + M_r^{-1} Sx.
!              Exception: ipar(8) < 0, Sx = 0. It is mostly used by
!              GMRES and variants to indicate (1) Sx was not necessary,
!              (2) intermediate result of Sx is not computed.
!     ipar(9)  points to the starting position of a work vector that
!              can be used by the caller.
!
!     NOTE: the caller should allow the iterative solver to perform
!     clean up job after the external convergence test is satisfied,
!     since some of the iterative solvers do not directly
!     update the 'sol' array. A typical clean-up stage includes
!     performing the final update of the approximate solution and
!     computing the convergence information (e.g. values of fpar(3:7)).
!
!     NOTE: fpar(4) and fpar(6) are not set by the accelerators (the
!     routines implemented here) if ipar(3) = 999.
!
!     ..................................................................
!     Usage:
!
!     To start solving a linear system, the user needs to specify
!     first 6 elements of the ipar, and first 2 elements of fpar.
!     The user may optionally set fpar(11) = 0 if one wants to count
!     the number of floating-point operations. (Note: the iterative
!     solvers will only add the floating-point operations inside
!     themselves, the caller will have to add the FLOPS from the
!     matrix-vector multiplication routines and the preconditioning
!     routines in order to account for all the arithmetic operations.)
!
!     Here is an example:
!     ipar(1) = 0	! always 0 to start an iterative solver
!     ipar(2) = 2	! right preconditioning
!     ipar(3) = 1	! use convergence test scheme 1
!     ipar(4) = 10000	! the 'w' has 10,000 elements
!     ipar(5) = 10	! use *GMRES(10) (e.g. FGMRES(10))
!     ipar(6) = 100	! use at most 100 matvecs
!     fpar(1) = 1.0E-6	! relative tolerance 1.0E-6
!     fpar(2) = 1.0E-10 ! absolute tolerance 1.0E-10
!     fpar(11) = 0.0	! clearing the FLOPS counter
!
!     After the above specifications, one can start to call an iterative
!     solver, say BCG. Here is a piece of pseudo-code showing how it can
!     be done,
!
! 10   call bcg(n,rhs,sol,ipar,fpar,w)
!      if (ipar(1).eq.1) then
!         call amux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
!         goto 10
!      else if (ipar(1).eq.2) then
!         call atmux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
!         goto 10
!      else if (ipar(1).eq.3) then
!         left preconditioner solver
!         goto 10
!      else if (ipar(1).eq.4) then
!         left preconditioner transposed solve
!         goto 10
!      else if (ipar(1).eq.5) then
!         right preconditioner solve
!         goto 10
!      else if (ipar(1).eq.6) then
!         right preconditioner transposed solve
!         goto 10
!      else if (ipar(1).eq.10) then
!         call my own stopping test routine
!         goto 10
!      else if (ipar(1).gt.0) then
!         ipar(1) is an unspecified code
!      else
!         the iterative solver terminated with code = ipar(1)
!      endif
!
!     This segment of pseudo-code assumes the matrix is in CSR format,
!     AMUX and ATMUX are two routines from the SPARSKIT MATVEC module.
!     They perform matrix-vector multiplications for CSR matrices,
!     where w(ipar(8)) is the first element of the input vectors to the
!     two routines, and w(ipar(9)) is the first element of the output
!     vectors from them. For simplicity, we did not show the name of
!     the routine that performs the preconditioning operations or the
!     convergence tests.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine gmres(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real(kind=dbl_kind) :: rhs(n), sol(n), fpar(16), w(*)
!-----------------------------------------------------------------------
!     This a version of GMRES implemented with reverse communication.
!     It is a simple restart version of the GMRES algorithm.
!
!     ipar(5) == the dimension of the Krylov subspace
!     after every ipar(5) iterations, the GMRES will restart with
!     the updated solution and recomputed residual vector.
!
!     the space of the 'w' is used as follows:
!     (1) the basis for the Krylov subspace, size n*(m+1);
!     (2) the Hessenberg matrix, only the upper triangular
!     portion of the matrix is stored, size (m+1)*m/2 + 1
!     (3) three vectors, all are of size m, they are
!     the cosine and sine of the Givens rotations, the third one holds
!     the residuals, it is of size m+1.
!
!     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
!     Note: m == ipar(5). The default value for this is 15 if
!     ipar(5) <= 1.
!-----------------------------------------------------------------------
!
!     local variables, ptr and p2 are temporary pointers,
!     hess points to the Hessenberg matrix,
!     vc, vs point to the cosines and sines of the Givens rotations
!     vrn points to the vectors of residual norms, more precisely
!     the right hand side of the least square problem solved.
!
      integer i,ii,idx,k,m,ptr,p2,hess,vc,vs,vrn
      real(kind=dbl_kind) :: alpha, c, s
      logical lp, rp
      save
!
!     check the status of the call
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60, 70) ipar(10)
!
!     initialization
!
      if (ipar(5).le.1) then
         m = 15
      else
         m = ipar(5)
      endif
      idx = n * (m+1)
      hess = idx + n
      vc = hess + (m+1) * m / 2 + 1
      vs = vc + m
      vrn = vs + m
      i = vrn + m + 1
      call bisinit(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
!
!     request for matrix vector multiplication A*x in the initialization
!
 100  ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = 1
      ipar(10) = 1
      k = 0
      do i = 1, n
         w(n+i) = sol(i)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      if (lp) then
         do i = 1, n
            w(n+i) = rhs(i) - w(i)
         enddo
         ipar(1) = 3
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i) = rhs(i) - w(i)
         enddo
      endif
      fpar(11) = fpar(11) + n
!
 20   alpha = sqrt(distdot(n,w,1,w,1))
      !write(*,*) alpha
      fpar(11) = fpar(11) + 2*n
      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
         if (abs(ipar(3)).eq.2) then
            fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
            fpar(11) = fpar(11) + 2*n
         else
            fpar(4) = fpar(1) * alpha + fpar(2)
         endif
         fpar(3) = alpha
      endif
      fpar(5) = alpha
      w(vrn+1) = alpha
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 300
      endif
      alpha = c1 / alpha
      do ii = 1, n
         w(ii) = alpha * w(ii)
      enddo
      fpar(11) = fpar(11) + n
!
!     request for (1) right preconditioning
!     (2) matrix vector multiplication
!     (3) left preconditioning
!
 110  k = k + 1
      if (rp) then
         ipar(1) = 5
         ipar(8) = k*n - n + 1
         if (lp) then
            ipar(9) = k*n + 1
         else
            ipar(9) = idx + 1
         endif
         ipar(10) = 3
         return
      endif
!
 30   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (k-1)*n + 1
      endif
      if (lp) then
         ipar(9) = idx + 1
      else
         ipar(9) = 1 + k*n
      endif
      ipar(10) = 4
      return
!
 40   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = k*n + 1
         ipar(10) = 5
         return
      endif
!
!     Modified Gram-Schmidt orthogonalization procedure
!     temporary pointer 'ptr' is pointing to the current column of the
!     Hessenberg matrix. 'p2' points to the new basis vector
!
 50   ipar(7) = ipar(7) + 1
      ptr = k * (k - 1) / 2 + hess
      p2 = ipar(9)
      call mgsro(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1), &
           ipar(12))
      if (ipar(12).lt.0) goto 200
!
!     apply previous Givens rotations and generate a new one to eliminate
!     the subdiagonal element.
!
      p2 = ptr + 1
      do i = 1, k-1
         ptr = p2
         p2 = p2 + 1
         alpha = w(ptr)
         c = w(vc+i)
         s = w(vs+i)
         w(ptr) = c * alpha + s * w(p2)
         w(p2) = c * w(p2) - s * alpha
      enddo
      call givens(w(p2), w(p2+1), c, s)
      w(vc+k) = c
      w(vs+k) = s
      p2 = vrn + k
      alpha = - s * w(p2)
      w(p2) = c * w(p2)
      w(p2+1) = alpha
!
!     end of one Arnoldi iteration, alpha will store the estimated
!     residual norm at current stage
!
      fpar(11) = fpar(11) + 6*k + 2
      alpha = abs(alpha)
      fpar(5) = alpha
      if (k.lt.m .and. .not.(ipar(3).ge.0 .and. alpha.le.fpar(4)) &
           .and. (ipar(6).le.0 .or. ipar(7).lt.ipar(6))) goto 110
!
!     update the approximate solution, first solve the upper triangular
!     system, temporary pointer ptr points to the Hessenberg matrix,
!     p2 points to the right-hand-side (also the solution) of the system.
!
 200  ptr = hess + k * (k + 1) / 2
      p2 = vrn + k
      if (w(ptr).eq.c0) then
!
!     if the diagonal elements of the last column is zero, reduce k by 1
!     so that a smaller trianguler system is solved [It should only
!     happen when the matrix is singular, and at most once!]
!
         k = k - 1
         if (k.gt.0) then
            goto 200
         else
            ipar(1) = -3
            ipar(12) = -4
            goto 300
         endif
      endif
      w(p2) = w(p2) / w(ptr)
      do i = k-1, 1, -1
         ptr = ptr - i - 1
         do ii = 1, i
            w(vrn+ii) = w(vrn+ii) - w(p2) * w(ptr+ii)
         enddo
         p2 = p2 - 1
         w(p2) = w(p2) / w(ptr)
      enddo
!
      do ii = 1, n
         w(ii) = w(ii) * w(p2)
      enddo
      do i = 1, k-1
         ptr = i*n
         p2 = p2 + 1
         do ii = 1, n
            w(ii) = w(ii) + w(p2) * w(ptr+ii)
         enddo
      enddo
      fpar(11) = fpar(11) + 2*k*n - n + k*(k+1)
!
      if (rp) then
         ipar(1) = 5
         ipar(8) = 1
         ipar(9) = idx + 1
         ipar(10) = 6
         return
      endif
!
60    if (rp) then
         do i = 1, n
            sol(i) = sol(i) + w(idx+i)
         enddo
      else
         do i = 1, n
            sol(i) = sol(i) + w(i)
         enddo
      endif
      fpar(11) = fpar(11) + n
!
!     process the complete stopping criteria
!
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = -1
         ipar(9) = idx + 1
         ipar(10) = 7
         return
      else if (ipar(3).lt.0) then
         if (ipar(7).le.m+1) then
            fpar(3) = abs(w(vrn+1))
            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
         endif
         fpar(6) = abs(w(vrn+k))
      else
         fpar(6) = fpar(5)
      endif
!
!     do we need to restart ?
!
 70   if (ipar(12).ne.0) then
         ipar(1) = -3
         goto 300
      endif
      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0) .and. &
           ((ipar(3).eq.999.and.ipar(11).eq.0) .or. &
           (ipar(3).ne.999.and.fpar(6).gt.fpar(4)))) goto 100
!
!     termination, set error code, compute convergence rate
!
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
            ipar(1) = 0
         else if (ipar(3).ne.999 .and. fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
 300  if (fpar(3).ne.c0 .and. fpar(6).ne.c0 .and. &
           ipar(7).gt.ipar(13)) then
         fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
      else
         fpar(7) = c0
      endif
      return
      end subroutine gmres
!-----end-of-gmres
!-----------------------------------------------------------------------
      subroutine gmres_fast(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real(kind=dbl_kind) :: rhs(n), sol(n), fpar(16), w(*)
!
      integer i,ii,idx,k,m,ptr,p2,hess,vc,vs,vrn
      real(kind=dbl_kind) :: alpha, c, s
      logical lp, rp
      save
!
!     check the status of the call
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60, 70) ipar(10)
!
!     initialization
!
      !!if (ipar(5).le.1) then
      !!   m = 15
      !!else
         m = ipar(5)
      !!endif
      idx = n * (m+1)
      hess = idx + n
      vc = hess + (m+1) * m / 2 + 1
      vs = vc + m
      vrn = vs + m
      i = vrn + m + 1
      call bisinit_fast(ipar,fpar,i,1,lp,rp,w)
      !!if (ipar(1).lt.0) return
!
!     request for matrix vector multiplication A*x in the initialization
!
 100  ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = 1
      ipar(10) = 1
      k = 0
      do i = 1, n
         w(n+i) = sol(i)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      !!if (lp) then
      !!   do i = 1, n
      !!      w(n+i) = rhs(i) - w(i)
      !!   enddo
      !!   ipar(1) = 3
      !!   ipar(10) = 2
      !!   return
      !!else
         do i = 1, n
            w(i) = rhs(i) - w(i)
         enddo
      !!endif
      !!fpar(11) = fpar(11) + n
!
 20   alpha = sqrt(distdot(n,w,1,w,1))
      !!fpar(11) = fpar(11) + 2*n
      !!if (ipar(7).eq.1 .and. ipar(3).ne.999) then
      if (ipar(7).eq.1) then
         !!if (abs(ipar(3)).eq.2) then
         !!   fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         !!   fpar(11) = fpar(11) + 2*n
         !!else
            fpar(4) = fpar(1) * alpha + fpar(2)
         !!endif
         fpar(3) = alpha
      endif
      fpar(5) = alpha
      w(vrn+1) = alpha
      !!if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
      if (alpha.le.fpar(4)) then
         ipar(1) = 0
         fpar(6) = alpha
         !!goto 300
         return
      endif
      alpha = c1 / alpha
      do ii = 1, n
         w(ii) = alpha * w(ii)
      enddo
      !!fpar(11) = fpar(11) + n
!
!     request for (1) right preconditioning
!     (2) matrix vector multiplication
!     (3) left preconditioning
!
 110  k = k + 1
      !!if (rp) then
         ipar(1) = 5
         ipar(8) = k*n - n + 1
         !!if (lp) then
         !!   ipar(9) = k*n + 1
         !!else
            ipar(9) = idx + 1
         !!endif
         ipar(10) = 3
         return
      !!endif
!
 30   ipar(1) = 1
      !!if (rp) then
         ipar(8) = ipar(9)
      !!else
      !!   ipar(8) = (k-1)*n + 1
      !!endif
      !!if (lp) then
      !!   ipar(9) = idx + 1
      !!else
         ipar(9) = 1 + k*n
      !!endif
      ipar(10) = 4
      return
!
 40   continue
      !!if (lp) then
      !!   ipar(1) = 3
      !!   ipar(8) = ipar(9)
      !!   ipar(9) = k*n + 1
      !!   ipar(10) = 5
      !!   return
      !!endif
!
!     Modified Gram-Schmidt orthogonalization procedure
!     temporary pointer 'ptr' is pointing to the current column of the
!     Hessenberg matrix. 'p2' points to the new basis vector
!
 50   ipar(7) = ipar(7) + 1
      ptr = k * (k - 1) / 2 + hess
      p2 = ipar(9)
      call mgsro_fast(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1), &
           ipar(12))
      if (ipar(12).lt.0) goto 200
!
!     apply previous Givens rotations and generate a new one to eliminate
!     the subdiagonal element.
!
      p2 = ptr + 1
      do i = 1, k-1
         ptr = p2
         p2 = p2 + 1
         alpha = w(ptr)
         c = w(vc+i)
         s = w(vs+i)
         w(ptr) = c * alpha + s * w(p2)
         w(p2) = c * w(p2) - s * alpha
      enddo
      call givens(w(p2), w(p2+1), c, s)
      w(vc+k) = c
      w(vs+k) = s
      p2 = vrn + k
      alpha = - s * w(p2)
      w(p2) = c * w(p2)
      w(p2+1) = alpha
!
!     end of one Arnoldi iteration, alpha will store the estimated
!     residual norm at current stage
!
      !!fpar(11) = fpar(11) + 6*k + 2
      alpha = abs(alpha)
      fpar(5) = alpha
      !!if (k.lt.m .and. .not.(ipar(3).ge.0 .and. alpha.le.fpar(4)) &
      !!     .and. (ipar(6).le.0 .or. ipar(7).lt.ipar(6))) goto 110
      if (k.lt.m .and. .not.(alpha.le.fpar(4)) &
           .and. (ipar(6).le.0 .or. ipar(7).lt.ipar(6))) goto 110
!
!     update the approximate solution, first solve the upper triangular
!     system, temporary pointer ptr points to the Hessenberg matrix,
!     p2 points to the right-hand-side (also the solution) of the system.
!
 200  ptr = hess + k * (k + 1) / 2
      p2 = vrn + k
      if (w(ptr).eq.c0) then
!
!     if the diagonal elements of the last column is zero, reduce k by 1
!     so that a smaller trianguler system is solved [It should only
!     happen when the matrix is singular, and at most once!]
!
         k = k - 1
         if (k.gt.0) then
            goto 200
         else
            ipar(1) = -3
            ipar(12) = -4
            !!goto 300
            return
         endif
      endif
      w(p2) = w(p2) / w(ptr)
      do i = k-1, 1, -1
         ptr = ptr - i - 1
         do ii = 1, i
            w(vrn+ii) = w(vrn+ii) - w(p2) * w(ptr+ii)
         enddo
         p2 = p2 - 1
         w(p2) = w(p2) / w(ptr)
      enddo
!
      do ii = 1, n
         w(ii) = w(ii) * w(p2)
      enddo
      do i = 1, k-1
         ptr = i*n
         p2 = p2 + 1
         do ii = 1, n
            w(ii) = w(ii) + w(p2) * w(ptr+ii)
         enddo
      enddo
      !!fpar(11) = fpar(11) + 2*k*n - n + k*(k+1)
!
      !!if (rp) then
         ipar(1) = 5
         ipar(8) = 1
         ipar(9) = idx + 1
         ipar(10) = 6
         return
      !!endif
!
60    continue
      !!if (rp) then
         do i = 1, n
            sol(i) = sol(i) + w(idx+i)
         enddo
      !!else
      !!   do i = 1, n
      !!      sol(i) = sol(i) + w(i)
      !!   enddo
      !!endif
      !!fpar(11) = fpar(11) + n
!
!     process the complete stopping criteria
!
      !!if (ipar(3).eq.999) then
      !!   ipar(1) = 10
      !!   ipar(8) = -1
      !!   ipar(9) = idx + 1
      !!   ipar(10) = 7
      !!   return
      !!else if (ipar(3).lt.0) then
      !!   if (ipar(7).le.m+1) then
      !!      fpar(3) = abs(w(vrn+1))
      !!      if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
      !!   endif
      !!   fpar(6) = abs(w(vrn+k))
      !!else
         fpar(6) = fpar(5)
      !!endif
!
!     do we need to restart ?
!
 70   if (ipar(12).ne.0) then
         ipar(1) = -3
         !!goto 300
         return
      endif
      !!if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0) .and. &
      !!     ((ipar(3).eq.999.and.ipar(11).eq.0) .or. &
      !!     (ipar(3).ne.999.and.fpar(6).gt.fpar(4)))) goto 100
      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0) .and. fpar(6).gt.fpar(4)) goto 100
!
!     termination, set error code, compute convergence rate
!
      if (ipar(1).gt.0) then
         !!if (ipar(3).eq.999 .and. ipar(11).eq.1) then
         !!   ipar(1) = 0
         !!else if (ipar(3).ne.999 .and. fpar(6).le.fpar(4)) then
         if (fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
!! 300  continue
      !!if (fpar(3).ne.c0 .and. fpar(6).ne.c0 .and. &
      !!     ipar(7).gt.ipar(13)) then
      !!   fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
      !!else
      !!   fpar(7) = c0
      !!endif
      return
      end subroutine gmres_fast
!-----end-of-gmres
!-----------------------------------------------------------------------
      subroutine givens(x,y,c,s)
      real(kind=dbl_kind) :: x,y,c,s
!-----------------------------------------------------------------------
!     Given x and y, this subroutine generates a Givens rotation c, s.
!     And apply the rotation on (x,y) ==> (sqrt(x**2 + y**2), 0).
!     (See P 202 of "matrix computation" by Golub and van Loan.)
!-----------------------------------------------------------------------
      real(kind=dbl_kind) :: t
!
      if (x.eq.c0 .and. y.eq.c0) then
         c = c1
         s = c0
      else if (abs(y).gt.abs(x)) then
         t = x / y
         x = sqrt(c1+t*t)
         s = sign(c1 / x, y)
         c = t*s
      else if (abs(y).le.abs(x)) then
         t = y / x
         y = sqrt(c1+t*t)
         c = sign(c1 / y, x)
         s = t*c
      else
!
!     X or Y must be an invalid floating-point number, set both to zero
!
         x = c0
         y = c0
         c = c1
         s = c0
      endif
      x = abs(x*y)
!
!     end of givens
!
      return
      end subroutine givens
!-----end-of-givens
!-----------------------------------------------------------------------
      subroutine bisinit(ipar,fpar,wksize,dsc,lp,rp,wk)
      implicit none
      integer i,ipar(16),wksize,dsc
      logical lp,rp
      real(kind=dbl_kind) :: fpar(16),wk(*)
!-----------------------------------------------------------------------
!     some common initializations for the iterative solvers
!-----------------------------------------------------------------------
!
!     ipar(1) = -2 inidcate that there are not enough space in the work
!     array
!
      if (ipar(4).lt.wksize) then
         ipar(1) = -2
         ipar(4) = wksize
         return
      endif
!
      if (ipar(2).gt.2) then
         lp = .true.
         rp = .true.
      else if (ipar(2).eq.2) then
         lp = .false.
         rp = .true.
      else if (ipar(2).eq.1) then
         lp = .true.
         rp = .false.
      else
         lp = .false.
         rp = .false.
      endif
      if (ipar(3).eq.0) ipar(3) = dsc
!     .. clear the ipar elements used
      ipar(7) = 0
      ipar(8) = 0
      ipar(9) = 0
      ipar(10) = 0
      ipar(11) = 0
      ipar(12) = 0
      ipar(13) = 0
!
!     fpar(1) must be between (0, 1), fpar(2) must be positive,
!     fpar(1) and fpar(2) can NOT both be zero
!     Normally return ipar(1) = -4 to indicate any of above error
!
      if (fpar(1).lt.c0 .or. fpar(1).ge.c1 .or. fpar(2).lt.c0 .or. &
           (fpar(1).eq.c0 .and. fpar(2).eq.c0)) then
         if (ipar(1).eq.0) then
            ipar(1) = -4
            return
         else
            fpar(1) = 1.0D-6
            fpar(2) = 1.0D-16
         endif
      endif
!     .. clear the fpar elements
      do i = 3, 10
         fpar(i) = c0
      enddo
      if (fpar(11).lt.c0) fpar(11) = c0
!     .. clear the used portion of the work array to zero
      do i = 1, wksize
         wk(i) = c0
      enddo
!
      return
!-----end-of-bisinit
      end subroutine bisinit
!-----------------------------------------------------------------------
      subroutine bisinit_fast(ipar,fpar,wksize,dsc,lp,rp,wk)
      implicit none
      integer i,ipar(16),wksize,dsc
      logical lp,rp
      real(kind=dbl_kind) :: fpar(16),wk(*)
!-----------------------------------------------------------------------
!     some common initializations for the iterative solvers
!-----------------------------------------------------------------------
!
!     ipar(1) = -2 inidcate that there are not enough space in the work
!     array
!
      !!if (ipar(4).lt.wksize) then
      !!   ipar(1) = -2
      !!   ipar(4) = wksize
      !!   return
      !!endif
!
      !!if (ipar(2).gt.2) then
      !!   lp = .true.
      !!   rp = .true.
      !!else if (ipar(2).eq.2) then
         lp = .false.
         rp = .true.
      !!else if (ipar(2).eq.1) then
      !!   lp = .true.
      !!   rp = .false.
      !!else
      !!   lp = .false.
      !!   rp = .false.
      !!endif
      !!if (ipar(3).eq.0) ipar(3) = dsc
!     .. clear the ipar elements used
      ipar(7) = 0
      ipar(8) = 0
      ipar(9) = 0
      ipar(10) = 0
      ipar(11) = 0
      ipar(12) = 0
      ipar(13) = 0
!
!     fpar(1) must be between (0, 1), fpar(2) must be positive,
!     fpar(1) and fpar(2) can NOT both be zero
!     Normally return ipar(1) = -4 to indicate any of above error
!
      !!if (fpar(1).lt.c0 .or. fpar(1).ge.c1 .or. fpar(2).lt.c0 .or. &
      !!     (fpar(1).eq.c0 .and. fpar(2).eq.c0)) then
      !!   if (ipar(1).eq.0) then
      !!      ipar(1) = -4
      !!      return
      !!   else
      !!      fpar(1) = 1.0D-6
      !!      fpar(2) = 1.0D-16
      !!   endif
      !!endif
!     .. clear the fpar elements
      do i = 3, 10
         fpar(i) = c0
      enddo
      !!if (fpar(11).lt.c0) fpar(11) = c0
!     .. clear the used portion of the work array to zero
      do i = 1, wksize
         wk(i) = c0
      enddo
!
      return
!-----end-of-bisinit
      end subroutine bisinit_fast
!-----------------------------------------------------------------------
      subroutine mgsro(full,lda,n,m,ind,ops,vec,hh,ierr)
      implicit none
      logical full
      integer lda,m,n,ind,ierr
      real(kind=dbl_kind) :: ops,hh(m),vec(lda,m)
!-----------------------------------------------------------------------
!     MGSRO  -- Modified Gram-Schmidt procedure with Selective Re-
!               Orthogonalization
!     The indth vector of VEC is orthogonalized against the rest of
!     the vectors.
!
!     The test for performing re-orthogonalization is performed for
!     each indivadual vectors. If the cosine between the two vectors
!     is greater than 0.99 (REORTH = 0.99**2), re-orthogonalization is
!     performed. The norm of the 'new' vector is kept in variable NRM0,
!     and updated after operating with each vector.
!
!     full   -- .ture. if it is necessary to orthogonalize the indth
!               against all the vectors vec(:,1:ind-1), vec(:,ind+2:m)
!               .false. only orthogonalize againt vec(:,1:ind-1)
!     lda    -- the leading dimension of VEC
!     n      -- length of the vector in VEC
!     m      -- number of vectors can be stored in VEC
!     ind    -- index to the vector to be changed
!     ops    -- operation counts
!     vec    -- vector of LDA X M storing the vectors
!     hh     -- coefficient of the orthogonalization
!     ierr   -- error code
!               0 : successful return
!               -1: zero input vector
!               -2: input vector contains abnormal numbers
!               -3: input vector is a linear combination of others
!
!     External routines used: real(kind=dbl_kind) distdot
!-----------------------------------------------------------------------
      integer i,k
      real(kind=dbl_kind) :: nrm0, nrm1, fct, thr
      real(kind=dbl_kind), parameter :: reorth = 0.98_dbl_kind
!
!     compute the norm of the input vector
!
      nrm0 = distdot(n,vec(1,ind),1,vec(1,ind),1)
      ops = ops + n + n
      thr = nrm0 * reorth
      if (nrm0.le.c0) then
         ierr = - 1
         return
      else if (nrm0.gt.c0 .and. c1/nrm0.gt.c0) then
         ierr = 0
      else
         ierr = -2
         return
      endif
!
!     Modified Gram-Schmidt loop
!
      if (full) then
         do 40 i = ind+1, m
            fct = distdot(n,vec(1,ind),1,vec(1,i),1)
            hh(i) = fct
            do 20 k = 1, n
               vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 20         continue
            ops = ops + 4 * n + 2
            if (fct*fct.gt.thr) then
               fct = distdot(n,vec(1,ind),1,vec(1,i),1)
               hh(i) = hh(i) + fct
               do 30 k = 1, n
                  vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 30            continue
               ops = ops + 4*n + 1
            endif
            nrm0 = nrm0 - hh(i) * hh(i)
            if (nrm0.lt.c0) nrm0 = c0
            thr = nrm0 * reorth
 40      continue
      endif
!
      do 70 i = 1, ind-1
         fct = distdot(n,vec(1,ind),1,vec(1,i),1)
         hh(i) = fct
         do 50 k = 1, n
            vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 50      continue
         ops = ops + 4 * n + 2
         if (fct*fct.gt.thr) then
            fct = distdot(n,vec(1,ind),1,vec(1,i),1)
            hh(i) = hh(i) + fct
            do 60 k = 1, n
               vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 60         continue
            ops = ops + 4*n + 1
         endif
         nrm0 = nrm0 - hh(i) * hh(i)
         if (nrm0.lt.c0) nrm0 = c0
         thr = nrm0 * reorth
 70   continue
!
!     test the resulting vector
!
      nrm1 = sqrt(distdot(n,vec(1,ind),1,vec(1,ind),1))
      ops = ops + n + n
 75   hh(ind) = nrm1
      if (nrm1.le.c0) then
         ierr = -3
         return
      endif
!
!     scale the resulting vector
!
      fct = c1 / nrm1
      do 80 k = 1, n
         vec(k,ind) = vec(k,ind) * fct
 80   continue
      ops = ops + n + 1
!
!     normal return
!
      ierr = 0
      return
!     end surbotine mgsro
      end subroutine mgsro
!-----------------------------------------------------------------------
      subroutine mgsro_fast(full,lda,n,m,ind,ops,vec,hh,ierr)
      implicit none
      logical full
      integer lda,m,n,ind,ierr
      real(kind=dbl_kind) :: ops,hh(m),vec(lda,m)
!-----------------------------------------------------------------------
!     MGSRO  -- Modified Gram-Schmidt procedure with Selective Re-
!               Orthogonalization
!     The indth vector of VEC is orthogonalized against the rest of
!     the vectors.
!
!     The test for performing re-orthogonalization is performed for
!     each indivadual vectors. If the cosine between the two vectors
!     is greater than 0.99 (REORTH = 0.99**2), re-orthogonalization is
!     performed. The norm of the 'new' vector is kept in variable NRM0,
!     and updated after operating with each vector.
!
!     full   -- .ture. if it is necessary to orthogonalize the indth
!               against all the vectors vec(:,1:ind-1), vec(:,ind+2:m)
!               .false. only orthogonalize againt vec(:,1:ind-1)
!     lda    -- the leading dimension of VEC
!     n      -- length of the vector in VEC
!     m      -- number of vectors can be stored in VEC
!     ind    -- index to the vector to be changed
!     ops    -- operation counts
!     vec    -- vector of LDA X M storing the vectors
!     hh     -- coefficient of the orthogonalization
!     ierr   -- error code
!               0 : successful return
!               -1: zero input vector
!               -2: input vector contains abnormal numbers
!               -3: input vector is a linear combination of others
!
!     External routines used: real(kind=dbl_kind) distdot
!-----------------------------------------------------------------------
      integer i,k
      real(kind=dbl_kind) :: nrm0, nrm1, fct, thr
      real(kind=dbl_kind), parameter :: reorth = 0.98_dbl_kind
!
!     compute the norm of the input vector
!
      nrm0 = distdot(n,vec(1,ind),1,vec(1,ind),1)
      !!ops = ops + n + n
      thr = nrm0 * reorth
      if (nrm0.le.c0) then
         ierr = - 1
         return
      else if (nrm0.gt.c0 .and. c1/nrm0.gt.c0) then
         ierr = 0
      else
         ierr = -2
         return
      endif
!
!     Modified Gram-Schmidt loop
!
      if (full) then
         do 40 i = ind+1, m
            fct = distdot(n,vec(1,ind),1,vec(1,i),1)
            hh(i) = fct
            do 20 k = 1, n
               vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 20         continue
            !!ops = ops + 4 * n + 2
            if (fct*fct.gt.thr) then
               fct = distdot(n,vec(1,ind),1,vec(1,i),1)
               hh(i) = hh(i) + fct
               do 30 k = 1, n
                  vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 30            continue
               !!ops = ops + 4*n + 1
            endif
            nrm0 = nrm0 - hh(i) * hh(i)
            if (nrm0.lt.c0) nrm0 = c0
            thr = nrm0 * reorth
 40      continue
      endif
!
      do 70 i = 1, ind-1
         fct = distdot(n,vec(1,ind),1,vec(1,i),1)
         hh(i) = fct
         do 50 k = 1, n
            vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 50      continue
         !!ops = ops + 4 * n + 2
         if (fct*fct.gt.thr) then
            fct = distdot(n,vec(1,ind),1,vec(1,i),1)
            hh(i) = hh(i) + fct
            do 60 k = 1, n
               vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 60         continue
            !!ops = ops + 4*n + 1
         endif
         nrm0 = nrm0 - hh(i) * hh(i)
         if (nrm0.lt.c0) nrm0 = c0
         thr = nrm0 * reorth
 70   continue
!
!     test the resulting vector
!
      nrm1 = sqrt(distdot(n,vec(1,ind),1,vec(1,ind),1))
      !!ops = ops + n + n
 75   hh(ind) = nrm1
      if (nrm1.le.c0) then
         ierr = -3
         return
      endif
!
!     scale the resulting vector
!
      fct = c1 / nrm1
      do 80 k = 1, n
         vec(k,ind) = vec(k,ind) * fct
 80   continue
      !!ops = ops + n + 1
!
!     normal return
!
      ierr = 0
      return
!     end surbotine mgsro
      end subroutine mgsro_fast

!-----------------------------------------------------------------------

      function distdot(N,DX,INCX,DY,INCY) result(ddot)

        integer :: N
        real(kind=dbl_kind) :: DX(*)
        integer :: INCX
        real(kind=dbl_kind) :: DY(*)
        integer :: INCY

        real(kind=dbl_kind) :: ddot

        integer :: i

        ddot = 0.0_dbl_kind

        do i = 1, N

           ddot = ddot + DX(i) * DY(i)

        enddo ! i

      end function distdot

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
