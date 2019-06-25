c $Id: ice.h,v 1.5.2.1 1998/04/20 20:36:13 eclare Exp $
c.. Common blocks for almost everything that the sea ice model passes around.
c-----------------------------------------------------------------------
c.. Grid parameters
      integer imt     ! number of gridpts in x direction, not incl. ghost pts. 
      integer jmt     ! number of gridpts in y direction, not incl. ghost pts. 
      integer nkmax   ! maximum number of ice layers

      parameter (imt=192, jmt=128, nkmax=2) !POP <4/3>
c      parameter (imt=384, jmt=288, nkmax=2) !POP augmented <2/3>
c      parameter (imt=107, jmt=49, nkmax=2)  !arctic
c      parameter (imt=192, jmt=30, nkmax=2)  !antarctic
c      parameter (imt=94, jmt=50, nkmax=2)   !NRL

c.. Flags and bounds
      integer iplot   ! how often output is written (10 = once per 10 dt)
      integer kread   ! 0  standard start
                      ! 1  restart from file dump.d
      integer kstress ! type of input data used - dynamics
                      ! 0  flux coupler
      integer kth     ! 0  thin ice
                      ! 1  thick ice
      integer ktherm  ! type of input data used - thermodynamics
                      ! 0  flux coupler 
      real max_comp   ! max(compact) = 0.999
      real tiny       ! a small positive number = 10^(-11)

      common /odd/ tiny, max_comp, iplot, ktherm, kstress, kread, kth

c.. Constants and parameters associated with time
      integer istep   ! step counter for time loop
      integer istep0  ! step counter, dumped value
      integer istep1  ! step counter, after dump = istep + istep0
      integer month   ! 1 \le month \le 12
      integer monthl  ! value of month at the last time step
      integer year    ! initial year
      integer idate   ! date (yyyymmdd)
      integer sec     ! elapsed seconds into date
      integer ndte    ! number of subcycles:  ndte=dt/dte
      integer npt     ! total number of time steps (dt) 

      real dt         ! thermo/transport timestep (s)
      real dte        ! subcycling timestep for EVP dynamics (s)
      real dtei       ! 1/dte
      real secyr      ! number of seconds in a year  (3.1536e7 = 365-day year)
      real time       ! time (s)
      real yday       ! day of the year

      common /times/ dt, dte, dtei, secyr, time, yday,
     & istep, istep0, istep1, month, monthl, year, idate, sec, 
     & ndte, npt

c.. Grid variables
      integer nkij(imt,jmt,0:1)   ! number of ice layers
                                  ! nk= 0  1 ice layer, < 25 cm thick
                                  ! nk=-1  no ice  (open ocean)

      real dxt(0:imt+1,0:jmt+1)   ! width of T-cell through the middle (cm)
      real dxtr(imt,jmt)          ! 1/dxt
      real dyt(0:imt+1,0:jmt+1)   ! height of T-cell through the middle (cm)
      real dytr(imt,jmt)          ! 1/dyt
      real dxu(0:imt+1,0:jmt+1)   ! width of U-cell through the middle (cm)
      real dyu(0:imt+1,0:jmt+1)   ! height of U-cell through the middle (cm)
      real HTE(0:imt+1,0:jmt+1)   ! length of eastern edge of T-cell (cm)
      real HTN(0:imt+1,0:jmt+1)   ! length of northern edge of T-cell (cm)
      real tarea(0:imt+1,0:jmt+1) ! area of T-cell (cm^2)
      real ULONG(imt,jmt)         ! longitude of velocity gridpoints (radians)
      real ULAT(imt,jmt)          ! latitude of velocity gridpoints  (radians)
      real ANGLE(imt,jmt)         ! for conversions betwn POP grid and lat/lng
      real xymin                  ! min(dxt,dyt)^2 (cm^2)

      common /grids/ dxt, dxtr, dyt, dytr, dxu, dyu,
     & HTN, HTE, ULONG, ULAT, ANGLE, tarea, xymin, nkij

c.. Masks
      real hm(0:imt+1,0:jmt+1)       ! land/boundary mask, thickness (T-cell)
      real uvm(0:imt,0:jmt)          ! land/boundary mask, velocity (U-cell)
      logical icehm(0:imt+1,0:jmt+1) ! ice extent mask (T-cell)
      logical iceuvm(0:imt,0:jmt)    ! ice extent mask (U-cell)
      logical iceuvmp(imt,jmt)       ! ice extent mask (U-cell), w/ N,E points

      common /masks/ hm, uvm, icehm, iceuvm, iceuvmp

c.. Transport component
      real cpct       ! initial value for compact
      real compact(0:imt+1,0:jmt+1)  ! fractional area of a T-cell 
                                     ! occupied by thick ice
      real fice(imt,jmt)             ! d(hi)/dt from thermodynamics (cm/s)
      real fsnow(imt,jmt)            ! d(hs)/dt from thermodynamics (cm/s)
      real gamma      ! snow aging timescale, ~55 days (2.e-7 /s)
      real hzero      ! classification thickness between thick/thin ice (cm)
      real qdrown     ! drowning rates for snow on thin and thick ice (cm/s)
      real rhoiw      ! rhow - rhoi            (g/cm^3)
      real rhosiw     ! rhos/(rhoi*rhow*dt)    (cm^3/g s)

      common /transp/ 
     & cpct, compact, fice, fsnow, 
     & gamma, hzero, qdrown, rhoiw, rhosiw

      real hice                      ! initial value for hithick  (cm)
      real hiold(0:imt+1,0:jmt+1)    ! temp variable: previous ice thickness
      real hitherm(0:imt+1,0:jmt+1)  ! ice thickness predicted by thermodynam.
      real hithick(0:imt+1,0:jmt+1)  ! thickness of thick ice (T-cell)
      real hithin(0:imt+1,0:jmt+1)   ! thickness of thin ice (T-cell)
      real hsnow                     ! initial value for hsthick 
      real hsold(0:imt+1,0:jmt+1)    ! temp variable: previous snow depth
      real hstherm(0:imt+1,0:jmt+1)  ! snow thickness predicted by thermodynam.
      real hsthick(0:imt+1,0:jmt+1)  ! thickness of snow on thick ice (T-cell)
      real hsthin(0:imt+1,0:jmt+1)   ! thickness of snow on thin ice (T-cell)

      common /thickness/ hice, hsnow, hiold, hsold,
     & hitherm, hstherm, hithick, hithin, hsthick, hsthin

c.. Dynamics component
      real u(0:imt,0:jmt)    ! x-component of velocity (cm/s)
      real v(0:imt,0:jmt)    ! y-component of velocity (cm/s)
      real divu(0:imt,0:jmt) ! velocity divergence (diagnostic) (1/s)
      real fcor(imt,jmt)     ! Coriolis parameter (1/s)
      real fm(imt,jmt)       ! Coriolis parameter * mass in U-cell (g/s)
      real gwatx(imt,jmt)    ! ocean current, x direction (cm/s)
      real gwaty(imt,jmt)    ! ocean current, x direction (cm/s)
      real strairx(imt,jmt)  ! stress on ice by air, x-direction (g/cm s^2)
      real strairy(imt,jmt)  ! stress on ice by air, y-direction (g/cm s^2)
      real strocnx(imt,jmt)  ! ice-ocean stress, x-direction (U-cell)
      real strocny(imt,jmt)  ! ice-ocean stress, y-direction (U-cell)
      real tiltx(imt,jmt)    ! sea surface slope, x-direction
      real tilty(imt,jmt)    ! sea surface slope, y-direction
      real waterx(imt,jmt)   ! temp variable for ice-ocean stress (cm/s)
      real watery(imt,jmt)   ! temp variable for ice-ocean stress (cm/s)
      real dragw             ! drag coefficient for water on ice * rhow 
      real ecc2              ! 1/e^2 for the elliptical yield curve 
      real ecc2m             ! 2.*(1.-ecc2)
      real ecc2p             ! 1.+ecc2
      real cw                ! cosine of the turning angle in water (phiwater)
      real sw                ! sine of the turning angle in water (phiwater)
      real rhoi              ! density of ice (g/cm^3)
      real rhos              ! density of snow (g/cm^3)
      real u0                ! constant coefficient for initial u field (cm/s)
      real v0                ! constant coefficient for initial v field (cm/s)

      common /dyn1/ u, v, divu, fcor, fm, gwatx, gwaty,
     & strairx, strairy, strocnx, strocny,
     & tiltx, tilty, waterx, watery,
     & dragw, ecc2, ecc2m, ecc2p, cw, sw,
     & rhoi, rhos, u0, v0

      real eyc               ! coefficient for calculating the parameter E
      real cst               ! c*, parameter in the defn of prss
      real pst               ! P*, parameter in the defn of prss (dyn/cm^2)
      real zetamin           ! minimum value of zetan(s,e,w) (g/s)
      real etan(imt,jmt)     ! shear viscosity in north triangle of T-cell 
      real etas(imt,jmt)     ! shear viscosity in south triangle of T-cell 
      real etae(imt,jmt)     ! shear viscosity in east  triangle of T-cell 
      real etaw(imt,jmt)     ! shear viscosity in west  triangle of T-cell 
      real prss(imt,jmt)     ! pressure P (centered in T-cell) (dyn/cm)
      real umass(0:imt,0:jmt)! mass of U-cell (g)
      real sig11n(imt,jmt)   ! internal ice stress tensor, sigma_11, north
      real sig11e(imt,jmt)   !   sigma_11, east       (g/s^2)
      real sig11s(imt,jmt)   !   sigma_11, south
      real sig11w(imt,jmt)   !   sigma_11, west
      real sig12n(imt,jmt)   !   sigma_12, north
      real sig12e(imt,jmt)   !   sigma_12, east
      real sig12s(imt,jmt)   !   sigma_12, south
      real sig12w(imt,jmt)   !   sigma_12, west
      real sig22n(imt,jmt)   !   sigma_22, north
      real sig22e(imt,jmt)   !   sigma_22, east
      real sig22s(imt,jmt)   !   sigma_22, south
      real sig22w(imt,jmt)   !   sigma_22, west
      real zetan(imt,jmt)    ! bulk viscosity in north triangle of T-cell (g/s)
      real zetas(imt,jmt)    ! bulk viscosity in south triangle of T-cell (g/s)
      real zetae(imt,jmt)    ! bulk viscosity in east  triangle of T-cell (g/s)
      real zetaw(imt,jmt)    ! bulk viscosity in west  triangle of T-cell (g/s)

      common /dyn2/ eyc, cst, pst, zetamin, 
     & etan, etas, etae, etaw, prss, umass, 
     & sig11n, sig11e, sig11s, sig11w, 
     & sig12n, sig12e, sig12s, sig12w, 
     & sig22n, sig22e, sig22s, sig22w,
     & zetan, zetas, zetae, zetaw

                             ! temporary variables for dynamics
      real
     & a2na(imt,jmt), a2sa(imt,jmt), a2ea(imt,jmt), a2wa(imt,jmt), 
     & b2n(imt,jmt), b2s(imt,jmt), b2e(imt,jmt), b2w(imt,jmt),
     & dxt8(imt,jmt), dyt8(imt,jmt), 
     & edy(imt,jmt), edx(imt,jmt), eHN(imt,jmt), eHE(imt,jmt),
     & eHNm(imt,jmt), eHEm(imt,jmt), HTN4(imt,jmt), HTE4(imt,jmt),
     & h2n(imt,jmt), h2s(imt,jmt), h2e(imt,jmt), h2w(imt,jmt),
     & prssn(imt,jmt), prsss(imt,jmt), prsse(imt,jmt), prssw(imt,jmt)

      common /dyn3/ 
     & a2na, a2sa, a2ea, a2wa, b2n, b2s, b2e, b2w,
     & dxt8, dyt8, edy, edx, eHN, eHE, eHNm, eHEm, HTN4, HTE4,
     & h2n, h2s, h2e, h2w, prssn, prsss, prsse, prssw

c.. Thermodynamics component
      real albice     ! ice albedo
      real albsnow    ! snow albedo
      real albw       ! ocean albedo
      real hmin       ! minimum allowable ice layer thickness (cm)
      real qs         ! heat of fusion of snow (cal/cm^3)
      real qi         ! heat of fusion of surface ice (cal/cm^3)
      real qb         ! heat of fusion of bottom ice (cal/cm^3)
      real rhocsi     ! 1/(density)*(specific heat) of snow (cm^3 C/cal)
      real rhoci      ! (density)*(specific heat) of ice (cal/cm^3 C)
      real rhocii     ! 1/rhoci
      real salinity   ! salinity of ice 
      real sigma      ! Stephan-Boltzmann constant (cal/s cm^2 K^4)
      real smin       ! minimum allowable snow layer thickness (cm)
      real solarIo    ! fraction of penetrating solar radiation
      real sk         ! conductivity of snow (cal/cm s C)
      real yk         ! conductivity of ice (cal/cm s C)

      common /therm1/ albice, albsnow, albw,
     & hmin, qs, qi, qb, rhocsi, rhoci, rhocii, 
     & salinity, sigma, smin, solarIo, sk, yk

      real qstorij(0:imt+1,0:jmt+1,0:1) ! heat stored in brine (cal/cm^2)
      real tij(imt,jmt,0:1,0:nkmax)     ! temperature of ice layers (C)
      real toceanij(imt,jmt,0:1)        ! ocean surface temperature (C)
      real tsfcij(imt,jmt,0:1)          ! surface temperature of ice/snow (C)

      common /therm2/ qstorij, tij, toceanij, tsfcij

      real flw(imt,jmt)    ! incoming longwave radiation (cal/s cm^2)
      real fsw(imt,jmt)    ! incoming shortwave radiation (cal/s cm^2)
      real fw(imt,jmt)     ! upward oceanic heat flux (cal/s cm^2)
      real Qa(imt,jmt)     ! specific humidity at 10 m 
      real snow(imt,jmt)   ! snowfall rate (cm/s)
      real sss(imt,jmt)    ! sea surface salinity (ppt)
      real sst(imt,jmt)    ! sea surface temperature (C)
      real tair(imt,jmt)   ! air temperature at 10 m (C)
      real tf(imt,jmt)     ! freezing temperature (C)
      real wind(imt,jmt)   ! wind speed (cm/s)

      common /fluxin/ flw, fsw, fw, Qa, snow, sss, sst,
     & tair, tf, wind

      real fsensible(imt,jmt,0:1)  ! sensible heat flux (cal/s cm^2)
      real flatent(imt,jmt,0:1)    ! latent heat flux (cal/s cm^2)
      real flwout(imt,jmt,0:1)     ! outgoing longwave radiation  (cal/s cm^2)
      real albij(imt,jmt,0:1)      ! surface albedo 
      real fresh(imt,jmt,0:1)      ! fresh water flux (g/cm^2 s)
      real fhnet(imt,jmt,0:1)      ! net heat flux to ocean

      common /fluxout/ fsensible, flatent, flwout, albij,
     & fresh, fhnet

      real salbedo(12)             ! monthly snow albedo
      real sfalltmp(12)            ! monthly snowfall rate (cm/s)

      common /Semtner/ salbedo, sfalltmp

      real sstd(imt,jmt,12)        ! sea surface temperature (C)
      real gwatxd(imt,jmt,12)      ! surface ocean current, x-direction (cm/s)
      real gwatyd(imt,jmt,12)      ! surface ocean current, y-direction (cm/s)
      real strxd(imt,jmt,12)       ! wind stress, x-direction (g/cm s^2)
      real stryd(imt,jmt,12)       ! wind stress, y-direction (g/cm s^2)
      real c0intp                  ! interpolation coefficients
      real c1intp
      real c2intp
      real c3intp

      common /POP/ sstd, gwatxd, gwatyd, strxd, stryd,
     & c0intp, c1intp, c2intp, c3intp

c.. File names
      character (80) :: 
     &       mask_file,         ! output file for mask info
     &       compact_file,      ! output file for compactness
     &       thick_file,        ! output file for thickness
     &       vel_file,          ! output file for velocities
     &       temps_file,        ! output file for temps
     &       fluxes_file,       ! output file for fluxes
     &       grid_file,         ! output file for POP grid info
     &       dump_file,         ! output file for restart dump
     &       restrt_file        ! input file for restarting
 
      common /filenames/ mask_file, compact_file, thick_file,
     &       vel_file, temps_file, fluxes_file, grid_file,
     &       dump_file, restrt_file
