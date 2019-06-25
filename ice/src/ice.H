c $Id: ice.H,v 2.2 2000/01/14 23:07:07 eclare Exp $
c----------------------------------------------------------------------------
c.. Common blocks for almost everything that the sea ice model passes around.
c----------------------------------------------------------------------------
c.. Global grid parameters 
      integer nlyr      ! maximum number of ice layers
      real rnlyr          ! real value of nlyr
      parameter (nlyr = 2)
      parameter (rnlyr = 2.)

      integer ncat     ! number of ice categories  !!! ncat > 1 !!!
      parameter (ncat = 2)

      integer imt_global ! number of gridpts in x direction, global domain,
                         !   not incl. ghost cells
      integer jmt_global ! number of gridpts in y direction, global domain,
                         !   not incl. ghost cells

      parameter (imt_global = 192, jmt_global = 128) !POP <4/3>

c.. Grid parameters for grid decomposition (local domains)
      integer nproc_x    ! number of processors assigned in x direction
      integer nproc_y    ! number of processors assigned in y direction
      integer nproc_s    ! total number of processors assigned, both directions
      integer nghost     ! number of ghost cell rows along local domain edges
      integer imt_local  ! number of gridpts in x direction, local domain,
                         !   including ghost cells
      integer jmt_local  ! number of gridpts in y direction, local domain,
                         !   including ghost cells
      integer ilo, ihi   ! first/last i indices of physical domain (local)
      integer jlo, jhi   ! first/last j indices of physical domain (local)

      parameter (nproc_x = NPROC_X, nproc_y = NPROC_Y)  ! set in Makefile
      parameter (nproc_s = nproc_x*nproc_y)
      parameter (nghost = 1)

      parameter (imt_local = (imt_global-1)/nproc_x+1+2*nghost)     
      parameter (jmt_local = (jmt_global-1)/nproc_y+1+2*nghost)     

      parameter (ilo = nghost+1, jlo = nghost+1)
      parameter (ihi = imt_local-nghost, jhi = jmt_local-nghost)

c.. Flags and bounds
      real tiny       ! a small positive number 
      parameter(tiny = 1.e-11)

      real c_max      ! maximum ice compactness
      real c_min      ! minimum ice compactness
      parameter(c_max = 0.99)    
      parameter(c_min = 1. - c_max)

      integer kread   ! 0  standard start
                      ! 1  restart from file dump.d
      integer kdata   ! type of input data used 
                      ! 0  flux coupler
      integer ktherm  ! type of thermodynamics (future implementation)

      common /odd/  ktherm, kdata, kread

c.. Constants and parameters associated with time
      real secyr      ! number of seconds in a year  
      parameter (secyr = 3.1536e7)     ! 365-day year
      integer daymo(12)                ! number of days in each month
      integer daycal(13)               ! day number at end of month
      data daymo/31,28,31,30,31,30,31,31,30,31,30,31/
      data daycal/0,31,59,90,120,151,181,212,243,273,304,334,365/

      integer istep   ! local step counter for time loop
      integer istep0  ! step counter, number of steps taken in previous run
      integer istep1  ! step counter, number of steps at current timestep
      integer mday    ! day of the month
      integer week    ! week of the year
      integer month   ! 1 \le month \le 12
      integer monthp  ! last month
      integer year    ! initial year
      integer nyr     ! year number
      integer idate   ! date (yyyymmdd)
      integer sec     ! elapsed seconds into date
      integer ndte    ! number of subcycles:  ndte=dt/dte
      integer npt     ! total number of time steps (dt) 

      logical new_year         ! new year = .true.
      logical new_month        ! new month = .true.
      logical new_week         ! new week = .true.
      logical new_day          ! new day = .true.
      logical write_history    ! write history now
      logical write_restart    ! write restart now
      logical write_ic         ! write initial condition now
      logical update_forcing   ! update forcing now
      character*1 outfreq      !  output frequency, 'y','m','d','1'
      character*1 dumpfreq     ! restart frequency, 'y','m','d'
      integer diagfreq         ! diagnostic output frequency
                               !    (10 = once per 10 dt)

      real dt         ! thermo/transport timestep (s)
      real dtei       ! 1/dte, where dte is the EVP timestep (1/s)
      real time       ! total elapsed time (s)
      real time_forc  ! time of last forcing update (s)
      real yday       ! day of the year

      common /times/ dt, dtei, time, time_forc, yday, monthp,
     & istep, istep0, istep1, mday, week, month, year, nyr, idate, sec,
     & ndte, npt, diagfreq, new_year, new_month, new_week, new_day, 
     & write_history, write_restart, update_forcing, outfreq, dumpfreq,
     & write_ic

c.. Constants and parameters associated with ice, snow and seawater
      real rhos                      ! density of snow (kg/m^3)
      real rhoi                      ! density of ice (kg/m^3)
      real rhoml                     ! density of seawater (kg/m^3)
      real albsnodry                 ! albedo of dry snow (Tsfc < 0)
      real albsnowet                 ! albedo of melting snow (Tsfc = 0)
      real albicemax                 ! max albedo of bare ice
      real albicemin                 ! minimum ice albedo (thin melting ice)
      real halb                      ! melt rate parameter for albedo

      parameter (rhos = 330.)
      parameter (rhoi = 900.)
      parameter (rhoml = 1026.)
      parameter (albsnodry = 0.85)
      parameter (albsnowet = 0.75)
      parameter (albicemax = 0.65)
      parameter (albicemin = 0.20)
      parameter (halb = .5)

      real rhoiw                     ! ocean-ice density difference (kg/m^3)
      parameter(rhoiw = rhoml - rhoi)

      real salinity                  ! salinity of ice (o/oo)
      parameter (salinity = 4.) 

c.. Grid variables
      real dxt(imt_local,jmt_local)  ! width of T-cell through the middle (m)
      real dxtr(ilo:ihi,jlo:jhi)     ! 0.5/dxt
      real dxtr4(ilo:ihi,jlo:jhi)    ! 0.25/dxt
      real dyt(imt_local,jmt_local)  ! height of T-cell through the middle (m)
      real dytr(ilo:ihi,jlo:jhi)     ! 0.5/dyt
      real dytr4(ilo:ihi,jlo:jhi)    ! 0.25/dyt
      real dxu(imt_local,jmt_local)  ! width of U-cell through the middle (m)
      real dyu(imt_local,jmt_local)  ! height of U-cell through the middle (m)
      real HTE(imt_local,jmt_local)  ! length of eastern edge of T-cell (m)
      real HTN(imt_local,jmt_local)  ! length of northern edge of T-cell (m)
      real tarea(imt_local,jmt_local)! area of T-cell (m^2)
      real uarea(imt_local,jmt_local)! area of U-cell (m^2)
      real ULONG(imt_local,jmt_local)! longitude of velocity pts (radians)
      real ULAT(imt_local,jmt_local) ! latitude of velocity pts (radians)
      real ANGLE(ilo:ihi,jlo:jhi)    ! for conversions between
                                     !    POP grid and lat/lon

      common /grids/ dxt, dxtr, dyt, dytr, dxu, dyu, dxtr4, dytr4,
     & HTN, HTE, ULONG, ULAT, ANGLE, tarea, uarea

c.. Masks
      real hm(imt_local,jmt_local)   ! land/boundary mask, thickness (T-cell)
      real uvm(imt_local,jmt_local)  ! land/boundary mask, velocity (U-cell)
      logical icehm(imt_local,jmt_local)   ! ice extent mask (T-cell)
      logical iceuvm(ilo:ihi,jlo:jhi)      ! ice extent mask (U-cell)
      real mask_n(imt_local,jmt_local)   ! northern hemisphere
      real mask_s(imt_local,jmt_local)   ! southern hemisphere

      common /masks/ hm, uvm, icehm, iceuvm, mask_n, mask_s

c.. Transport component
      real lateral                   ! fraction of new ice volume used
      parameter (lateral=0.5)        !   for lateral ice growth

      real c0(ncat)                  ! initial values for compact
      real compact(imt_local,jmt_local,0:ncat) ! fractional area of a T-cell 
                                     ! occupied by ice of the given category
      real ifrc(ilo:ihi,jlo:jhi)     ! total ice fraction
      real hice0(ncat)               ! initial values for ice thickness (m)
      real hsnow0(ncat)              ! initial values for snow thickness (m)
      real hice(imt_local,jmt_local,0:ncat)  ! ice thickness (T-cell) (m)
      real hsnow(imt_local,jmt_local,0:ncat) ! snow depth (m)
      real hzero(ncat-1)     ! ice category demarcation thicknesses (m)

      common /transp/  c0, compact, hice0, hsnow0, hice, hsnow, hzero,
     &  ifrc

c.. Dynamics component
      real grav              ! gravitational acceleration (m/s^2)
      real dragw             ! drag coefficient for water on ice * rhoml
      real cosw              ! cosine of the turning angle in water (phiwater)
      real sinw              ! sine of the turning angle in water (phiwater)

      parameter(grav = 9.80616)      
      parameter(dragw = 0.0055*rhoml)! (kg/m^3)
      parameter(cosw = 1.)           ! phiwater = 0 deg
      parameter(sinw = 0.)
c      phiwater = 25.*acos(-1.)/180. ! 25 deg converted to radians
c      sinw = sin(phiwater)
c      cosw = cos(phiwater)

      real u(imt_local,jmt_local)    ! x-component of velocity (m/s)
      real v(imt_local,jmt_local)    ! y-component of velocity (m/s)
      real divu(ilo:ihi,jlo:jhi)     ! velocity divergence (diagnostic) (1/s)
      real umassdtei(ilo:ihi,jlo:jhi)! mass of U-cell/dte (kg/m^2 s)
      real fcor(ilo:ihi,jlo:jhi)     ! Coriolis parameter (1/s)
      real fm(ilo:ihi,jlo:jhi)       ! Coriolis param. * mass in U-cell (kg/s)
      real prss(ilo:ihi,jlo:jhi)     ! pressure P (centered in T-cell) (kg/s)
      real cst               ! c*, parameter in the defn of prss
      real pst               ! P*, parameter in the defn of prss (kg/m s)
      real u0                ! constant coefficient for initial u field (m/s)
      real v0                ! constant coefficient for initial v field (m/s)
      real ecci              ! 4/e^2 for the elliptical yield curve 
      real eccm              ! 2.-0.5*ecci
      real eccp              ! 1.+ecci*0.25
      real eyc               ! coefficient for calculating the parameter E
      real acoef             ! cryptic combinations of stress equation
      real ccoef             !   coefficients
      real bcoef             ! 
      real Tdtei             ! 2(wave damping timescale)/dte

      common /dyn1/ u, v, divu, umassdtei, fcor, fm, prss, cst, pst,
     & u0, v0, ecci, eccm, eccp, eyc, acoef, bcoef, ccoef, Tdtei

      real sig11ne(ilo:ihi,jlo:jhi) ! ice stress tensor: sigma_11 north
      real sig11sw(ilo:ihi,jlo:jhi) !   sigma_11, east       (kg/s^2)
      real sig11nw(ilo:ihi,jlo:jhi) !   sigma_11, south
      real sig11se(ilo:ihi,jlo:jhi) !   sigma_11, west
      real sig12ne(ilo:ihi,jlo:jhi) !   sigma_12, north
      real sig12sw(ilo:ihi,jlo:jhi) !   sigma_12, east
      real sig12nw(ilo:ihi,jlo:jhi) !   sigma_12, south
      real sig12se(ilo:ihi,jlo:jhi) !   sigma_12, west
      real sig22ne(ilo:ihi,jlo:jhi) !   sigma_22, north
      real sig22sw(ilo:ihi,jlo:jhi) !   sigma_22, east
      real sig22nw(ilo:ihi,jlo:jhi) !   sigma_22, south
      real sig22se(ilo:ihi,jlo:jhi) !   sigma_22, west
      real prs_sig(ilo:ihi,jlo:jhi) ! replacement pressure, for stress calc
      real sig1   (ilo:ihi,jlo:jhi) ! principal stress component (diagnostic)
      real sig2   (ilo:ihi,jlo:jhi) ! principal stress component (diagnostic)

      common /dyn_stress/  prs_sig, sig1, sig2,
     & sig11ne, sig11se, sig11sw, sig11nw, 
     & sig12ne, sig12se, sig12sw, sig12nw, 
     & sig22ne, sig22se, sig22sw, sig22nw

      real gwatx(ilo:ihi,jlo:jhi)    ! ocean current, x-direction (m/s)
      real gwaty(ilo:ihi,jlo:jhi)    ! ocean current, y-direction (m/s)
      real strairx(ilo:ihi,jlo:jhi)  ! stress on ice by air, x-direction
      real strairy(ilo:ihi,jlo:jhi)  ! stress on ice by air, y-direction
                                          ! (kg/m s^2)
      real strocnx(ilo:ihi,jlo:jhi)  ! ice-ocean stress, x-dir. (U-cell)
      real strocny(ilo:ihi,jlo:jhi)  ! ice-ocean stress, y-dir. (U-cell)
      real tiltx(ilo:ihi,jlo:jhi)    ! sea surface slope, x-direction
      real tilty(ilo:ihi,jlo:jhi)    ! sea surface slope, y-direction
      real waterx(ilo:ihi,jlo:jhi)   ! temp variable for ice-ocean stress 
      real watery(ilo:ihi,jlo:jhi)   ! temp variable for ice-ocean stress 
                                          ! (m/s)

      common /dyn_flux/ gwatx, gwaty, strairx, strairy, 
     & strocnx, strocny, tiltx, tilty, waterx, watery

c.. Thermodynamics component
      real mu_Tf                          ! Tf:brine salinity ratio (C/ppt)
      parameter (mu_Tf = 0.054)

      real iceruf(ncat)                   ! ice surface roughness (m)
      real Tmlt (nlyr)                    ! melting temp; depends on S (C)
      real Tsfc(imt_local,jmt_local,0:ncat)! ice/snow surface temperature (C)
      real qice (imt_local,jmt_local,0:ncat,0:nlyr)
                                          !enthalpy of ice layers (J/kg)
      real shcoef(ilo:ihi,jlo:jhi,ncat) ! transfer coefficient for sensible ht
      real lhcoef(ilo:ihi,jlo:jhi,ncat) ! transfer coefficient for latent heat
      real fw   (ilo:ihi,jlo:jhi)  ! upward oceanic heat flux (W/m^2)

      common /therm_save/ Tmlt, Tsfc, qice, shcoef, lhcoef, iceruf, fw

      real flw   (ilo:ihi,jlo:jhi)  ! incoming longwave radiation (W/m^2)
      real fsw   (ilo:ihi,jlo:jhi)  ! incoming shortwave radiation (W/m^2)
      real frzmlt(ilo:ihi,jlo:jhi)  ! freezing/melting potential (W/m^2)
      real rain  (ilo:ihi,jlo:jhi)  ! rainfall rate (kg/m^2 s)
      real snow  (ilo:ihi,jlo:jhi)  ! snowfall rate (m/s)
      real sss   (ilo:ihi,jlo:jhi)  ! sea surface salinity (ppt)
      real sst   (ilo:ihi,jlo:jhi)  ! sea surface temperature (C)
      real Tf    (ilo:ihi,jlo:jhi)  ! freezing temperature (C)
      real wind  (ilo:ihi,jlo:jhi)  ! wind speed (m/s)

      common /therm_in/ flw, fsw, frzmlt, rain, snow, sss, sst,
     & Tf, wind

      real fsensible(ilo:ihi,jlo:jhi,ncat)  ! sensible heat flux (W/m^2)
      real flatent  (ilo:ihi,jlo:jhi,ncat)  ! latent heat flux (W/m^2)
      real flwout   (ilo:ihi,jlo:jhi,ncat)  ! outgoing longwave radiation  
                                            !    (W/m^2)
      real albedo    (ilo:ihi,jlo:jhi,ncat) ! surface albedo 
      real fresh    (ilo:ihi,jlo:jhi,ncat)  ! fresh water flux (kg/m^2 s)
      real fhnet    (ilo:ihi,jlo:jhi,ncat)  ! net heat flux to ocean
      real fswthru  (ilo:ihi,jlo:jhi,ncat)  ! shortwave thru ice to ocean
      real fsalt    (ilo:ihi,jlo:jhi,ncat)  ! net salt flux to ocean

      common /therm_out/ fsensible, flatent, flwout, albedo,
     & fresh, fhnet, fswthru, fsalt

c.. Atmospheric state/fluxes from/to coupler
      real zlvl (ilo:ihi,jlo:jhi)  ! atm level height (m)
      real uatm (ilo:ihi,jlo:jhi)  ! wind velocity    (m/s)
      real vatm (ilo:ihi,jlo:jhi)
      real potT (ilo:ihi,jlo:jhi)  ! air potential temperature  (K)
      real Tair (ilo:ihi,jlo:jhi)  ! air temperature  (K)
      real Qa   (ilo:ihi,jlo:jhi)  ! specific humidity (kg/kg)
      real rhoa (ilo:ihi,jlo:jhi)  ! air density (kg/m^3)
      real swvdr(ilo:ihi,jlo:jhi)  ! sw down, visible, direct  (W/m^2)
      real swvdf(ilo:ihi,jlo:jhi)  ! sw down, visible, diffuse (W/m^2)
      real swidr(ilo:ihi,jlo:jhi)  ! sw down, near IR, direct  (W/m^2)
      real swidf(ilo:ihi,jlo:jhi)  ! sw down, near IR, diffuse (W/m^2)
      real Tref (ilo:ihi,jlo:jhi)  ! 2m atm reference temperature (K)

      common /atmo/ zlvl, uatm, vatm, potT, Tair, Qa, rhoa,
     & swvdr, swvdf, swidr, swidf, Tref

      real c0intp                  ! interpolation coefficients
      real c1intp
      real c2intp
      real c3intp

      common /interp/
     & c0intp, c1intp, c2intp, c3intp

c.. File names
      character (80) :: 
     &       grid_file,         !  input file for POP grid info
     &       kmt_file,          !  input file for POP grid info
     &       dump_file,         ! output file for restart dump
     &       restrt_file,       !  input file for restarting
     &       history_file       ! output file for history "tape"
 
      common /filenames/ grid_file, kmt_file,
     &       dump_file, restrt_file, history_file


      character (72) ::
     &       sst_file,          ! filenames, sea surface temperature
     &       sss_file,          ! filenames, sea surface salinity
     &       fsw_file,          ! filenames, shortwave
     &       flw_file,          ! filenames, longwave
     &       rain_file,         ! filenames, rainfall rate
     &       tair_file,         ! filenames, sfc air temperature
     &       strx_file,         ! filenames, wind stress, x
     &       stry_file,         ! filenames, wind stress, y 
     &       wind_file,         ! filenames, wind speed
     &       humid_file,        ! filenames, specific humidity
     &       height_file,       ! filenames, atmo level height
     &       uwind_file,        ! filenames, u wind component
     &       vwind_file,        ! filenames, v wind component
     &       potT_file,         ! filenames, air potential temperature
     &       rhoa_file          ! filenames, air density
 
      common /Weddell/ fsw_file, flw_file, rain_file, 
     &       tair_file, strx_file, stry_file, wind_file, humid_file, 
     &       sst_file, sss_file, height_file, uwind_file, vwind_file,
     &       potT_file, rhoa_file


c.. Parallelization stuff
      real ::
     &   ice_global_real_minval,! global minimum function 
     &   ice_global_real_maxval,! global maximum function
     &   ice_global_real_sum    ! global sum
      integer ::
     &   my_task,               ! task ID for current processor
     &   master_task,           ! task ID for controlling processor
     &   nb_tasks,              ! total number of tasks
     &   ierr                   ! general-use error flag
      integer ::
     &   MPI_COMM_ICE           ! communicator for internal ice comms

      common / setup /
     &   my_task, nb_tasks, master_task,
     &   MPI_COMM_ICE




