#define diag 1

module ice_therm_oned

  use ice_kinds_mod
  use ice_constants
  use ice_domain_size
  use ice_therm_mushy

  implicit none

  real(kind=dbl_kind) :: z_interface_snow ! snow surface
  real(kind=dbl_kind) :: z_interface_snic ! snow ice surface
  real(kind=dbl_kind) :: z_interface_ices ! ice surface
  real(kind=dbl_kind) :: z_interface_iceb ! ice base

  real(kind=dbl_kind) :: z_interface_snow_float ! snow surface relative to ocean surface
  real(kind=dbl_kind) :: z_interface_snic_float ! snow ice surface relative to ocean surface
  real(kind=dbl_kind) :: z_interface_ices_float ! ice surface relative to ocean surface
  real(kind=dbl_kind) :: z_interface_iceb_float ! ice base relative to ocean surface

  integer :: istep1_prev = 0

#ifdef notz_fieldwork
  !-----------------------------------------------------------------
  ! Notz fieldwork
  !-----------------------------------------------------------------

  ! adventdalen weather station data
  integer, parameter :: ndata_notz = 4582
  real(kind=dbl_kind), dimension(ndata_notz) :: time_notz
  real(kind=dbl_kind), dimension(ndata_notz) :: tair_notz
  real(kind=dbl_kind), dimension(ndata_notz) :: relh_notz
  real(kind=dbl_kind), dimension(ndata_notz) :: wind_notz

  ! ERA data
  integer, parameter :: ndata_era = 249
  real(kind=dbl_kind), dimension(ndata_era) :: time_era
  real(kind=dbl_kind), dimension(ndata_era) :: tcc_era
  real(kind=dbl_kind), dimension(ndata_era) :: ssrd_era
  real(kind=dbl_kind), dimension(ndata_era) :: strd_era

  !real(kind=dbl_kind), parameter :: time0_notz = 5313600.0_dbl_kind ! simulation time of notz t=0
  !real(kind=dbl_kind), parameter :: time0_notz = c0 ! simulation time of notz t=0
  !real(kind=dbl_kind), parameter :: time1_notz = 471900.0_dbl_kind  ! first restart in notz time
  !real(kind=dbl_kind), parameter :: time2_notz = 1000500.0_dbl_kind ! second restart in notz time
  
  ! in secs since start of the year
  real(kind=dbl_kind), parameter :: time1_notz = 5741712.0_dbl_kind  ! first restart in notz time
  real(kind=dbl_kind), parameter :: time2_notz = 6270048.0_dbl_kind ! second restart in notz time

  integer :: istep_notz1
  integer :: istep_notz2

  logical :: lrestart_notz1 = .false.
  logical :: lrestart_notz2 = .false.

  real(kind=dbl_kind), parameter :: lat_advent = 78.241077
  real(kind=dbl_kind), parameter :: lon_advent = 15.595265
  real(kind=dbl_kind), parameter :: lat_advent_rad = lat_advent / rad_to_deg
  real(kind=dbl_kind), parameter :: lon_advent_rad = lon_advent / rad_to_deg

#endif

#if defined era-interim

  ! ERA-Interim data
  integer, parameter :: ntimes1_era = 1460
  integer, parameter :: ntimes2_era = 2920
  real(kind=dbl_kind), dimension(ntimes1_era) :: times1_era
  real(kind=dbl_kind), dimension(ntimes1_era) :: t2m_era
  real(kind=dbl_kind), dimension(ntimes1_era) :: d2m_era
  real(kind=dbl_kind), dimension(ntimes1_era) :: u10_era
  real(kind=dbl_kind), dimension(ntimes1_era) :: v10_era
  real(kind=dbl_kind), dimension(ntimes1_era) :: tcc_era

  real(kind=dbl_kind), dimension(ntimes2_era) :: times2_era
  real(kind=dbl_kind), dimension(ntimes2_era) :: tp_era
  real(kind=dbl_kind), dimension(ntimes2_era) :: sf_era

  ! ocean data
  integer, parameter :: ntimes_ocean = 14
  real(kind=dbl_kind), dimension(ntimes_ocean) :: times_ocean
  real(kind=dbl_kind), dimension(ntimes_ocean) :: S_ocean
  real(kind=dbl_kind), dimension(ntimes_ocean) :: T_ocean
  real(kind=dbl_kind), dimension(ntimes_ocean) :: U_ocean
  real(kind=dbl_kind), dimension(ntimes_ocean) :: V_ocean
  real(kind=dbl_kind), dimension(ntimes_ocean) :: dhdx_ocean
  real(kind=dbl_kind), dimension(ntimes_ocean) :: dhdy_ocean
  real(kind=dbl_kind), dimension(ntimes_ocean) :: hblt_ocean
  real(kind=dbl_kind), dimension(ntimes_ocean) :: qdp_ocean

#if defined snowice_maksym

  real(kind=dbl_kind), parameter :: time_snowice_maksym = 121.0_dbl_kind * 24.0_dbl_kind * 3600.0_dbl_kind

  logical :: lrestart_snowice_maksym = .false.

#endif
#endif

contains

!=======================================================================

  subroutine init_therm_oned()

#ifdef notz_fieldwork
    call load_notz_fieldwork_data()
#endif
    
#ifdef notz_experiment
    call calculate_notz_fbot_coeff(0.01_dbl_kind)
#endif
    
#if defined era-interim
    call load_era_data()
    call load_ocean_data()
#endif

  end subroutine init_therm_oned

#ifdef notz_fieldwork
!=======================================================================
! Notz fieldwork forcing
!=======================================================================

  subroutine restart_simulation_notz_fieldwork_times(time, nx_block, ny_block, vicen, vsnon, trcrn, phi_init)

    real(kind=dbl_kind), intent(in) :: time
    integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block  ! block dimensions
    real(kind=dbl_kind), dimension(nx_block,ny_block,ncat), intent(inout) :: &
         vicen, vsnon
    real(kind=dbl_kind), dimension(nx_block,ny_block,max_ntrcr,ncat), intent(inout) :: &
         trcrn
    real(kind=dbl_kind), intent(in) :: phi_init
   
    !write(*,*) time
    
    ! first experiment
    if (time >= time1_notz .and. .not. lrestart_notz1) then
       
       ! sort out diagnostics
       !call notz_diagnostics()

       ! restart the simulation
       call restart_simulation_notz_fieldwork(time, nx_block, ny_block, vicen, vsnon, trcrn, phi_init)

       ! record the restarting timestep
       istep_notz1 = istep1
       
       ! lets not do it again!
       lrestart_notz1 = .true.
       
    endif

    ! second experiment
    if (time >= time2_notz .and. .not. lrestart_notz2) then
       
       ! restart the simulation
       call restart_simulation_notz_fieldwork(time, nx_block, ny_block, vicen, vsnon, trcrn, phi_init)

       ! record the restarting timestep
       istep_notz2 = istep1
       
       ! lets not do it again!
       lrestart_notz2 = .true.
       
    endif

! plot 'hin_gnuplot.txt' using ($1*900):2 with lines, 'hin_gnuplot.txt' using ($1*900):3 with lines, '/Users/akt/Work/Convection/NotzThesis/fieldwork/fieldwork1/h.txt' using ($1+5785500):2 with lines, '/Users/akt/Work/Convection/NotzThesis/fieldwork/fieldwork2/h.txt' using ($1+6314100):2 with lines

! plot 'Tin_cice_1.txt' using ($1/3600):2 with lines, 'Tin_cice_1.txt' using ($1/3600):3 with lines, 'Tin_cice_1.txt' using ($1/3600):4 with lines, 'Tin_cice_1.txt' using ($1/3600):5 with lines, 'Tin_cice_1.txt' using ($1/3600):6 with lines, 'Tin_cice_1.txt' using ($1/3600):7 with lines, 'Tin_cice_1.txt' using ($1/3600):8 with lines, 'Tin_cice_1.txt' using ($1/3600):9 with lines, 'Tin_cice_1.txt' using ($1/3600):10 with lines, 'Tin_cice_1.txt' using ($1/3600):11 with lines, 'Tin_cice_1.txt' using ($1/3600):12 with lines, 'Tin_cice_1.txt' using ($1/3600):13 with lines, 'Tin_cice_1.txt' using ($1/3600):14 with lines, 'Tin_cice_1.txt' using ($1/3600):15 with lines

! plot 'Sin_cice_1.txt' using ($1/3600):2 with lines, 'Sin_cice_1.txt' using ($1/3600):3 with lines, 'Sin_cice_1.txt' using ($1/3600):4 with lines, 'Sin_cice_1.txt' using ($1/3600):5 with lines, 'Sin_cice_1.txt' using ($1/3600):6 with lines, 'Sin_cice_1.txt' using ($1/3600):7 with lines, 'Sin_cice_1.txt' using ($1/3600):8 with lines, 'Sin_cice_1.txt' using ($1/3600):9 with lines, 'Sin_cice_1.txt' using ($1/3600):10 with lines, 'Sin_cice_1.txt' using ($1/3600):11 with lines, 'Sin_cice_1.txt' using ($1/3600):12 with lines, 'Sin_cice_1.txt' using ($1/3600):13 with lines, 'Sin_cice_1.txt' using ($1/3600):14 with lines, 'Sin_cice_1.txt' using ($1/3600):15 with lines

! plot 'phi_cice_1.txt' using ($1/3600):2 with lines, 'phi_cice_1.txt' using ($1/3600):3 with lines, 'phi_cice_1.txt' using ($1/3600):4 with lines, 'phi_cice_1.txt' using ($1/3600):5 with lines, 'phi_cice_1.txt' using ($1/3600):6 with lines, 'phi_cice_1.txt' using ($1/3600):7 with lines, 'phi_cice_1.txt' using ($1/3600):8 with lines, 'phi_cice_1.txt' using ($1/3600):9 with lines, 'phi_cice_1.txt' using ($1/3600):10 with lines, 'phi_cice_1.txt' using ($1/3600):11 with lines, 'phi_cice_1.txt' using ($1/3600):12 with lines, 'phi_cice_1.txt' using ($1/3600):13 with lines, 'phi_cice_1.txt' using ($1/3600):14 with lines, 'phi_cice_1.txt' using ($1/3600):15 with lines

! plot 'Sbr_cice_1.txt' using ($1/3600):2 with lines, 'Sbr_cice_1.txt' using ($1/3600):3 with lines, 'Sbr_cice_1.txt' using ($1/3600):4 with lines, 'Sbr_cice_1.txt' using ($1/3600):5 with lines, 'Sbr_cice_1.txt' using ($1/3600):6 with lines, 'Sbr_cice_1.txt' using ($1/3600):7 with lines, 'Sbr_cice_1.txt' using ($1/3600):8 with lines, 'Sbr_cice_1.txt' using ($1/3600):9 with lines, 'Sbr_cice_1.txt' using ($1/3600):10 with lines, 'Sbr_cice_1.txt' using ($1/3600):11 with lines, 'Sbr_cice_1.txt' using ($1/3600):12 with lines, 'Sbr_cice_1.txt' using ($1/3600):13 with lines, 'Sbr_cice_1.txt' using ($1/3600):14 with lines, 'Sbr_cice_1.txt' using ($1/3600):15 with lines

  end subroutine restart_simulation_notz_fieldwork_times

!=======================================================================

  subroutine notz_diagnostics()
    
    call delete_diagnostic_file("history/Temp_gnuplot.txt")
    call delete_diagnostic_file("history/Sin_gnuplot.txt")
    call delete_diagnostic_file("history/Sbr_gnuplot.txt")
    call delete_diagnostic_file("history/phi_gnuplot.txt")
    
  end subroutine notz_diagnostics

!=======================================================================
  
  subroutine delete_diagnostic_file(filename)
    
    character(len=*), intent(in) :: filename
    
    open(11,file=filename)
    close(11,status='delete')
    
  end subroutine delete_diagnostic_file

!=======================================================================

  subroutine restart_simulation_notz_fieldwork(time, nx_block, ny_block, vicen, vsnon, trcrn, phi_init)

    use ice_state, only: nt_Tsfc, nt_sice, nt_qice, nt_qsno
    
    real(kind=dbl_kind), intent(in) :: time
    integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block  ! block dimensions
    real(kind=dbl_kind), dimension(nx_block,ny_block,ncat), intent(inout) :: &
         vicen, vsnon
    real(kind=dbl_kind), dimension(nx_block,ny_block,max_ntrcr,ncat), intent(inout) :: &
         trcrn
    real(kind=dbl_kind), intent(in) :: phi_init

    real(kind=dbl_kind), parameter :: saltmax = 35.0_dbl_kind
    real(kind=dbl_kind), parameter :: dSin0_frazil = 3.0
    
    integer, parameter :: i = 4
    integer, parameter :: j = 4
    integer, parameter :: n = 1
    
    real(kind=dbl_kind) :: Ti, salin
    integer :: k
    
    write(*,*) "Notz restart: ", istep
    
    salin = saltmax - dSin0_frazil

    ! ice volume
    vicen(i,j,n) = 0.01_dbl_kind
    
    ! snow volume
    vsnon(i,j,n) = c0
    
    ! surface temperature
    trcrn(i,j,nt_Tsfc,n) = liquidus_temperature_mush(salin / phi_init)
    
    ! snow temperature
    do k = 1, nslyr

       Ti = min(c0, trcrn(i,j,nt_Tsfc,n))
       trcrn(i,j,nt_qsno+k-1,n) = -rhos*(Lfresh - cp_ice*Ti)
       
    enddo ! k
    
    ! ice temperature
    do k = 1, nilyr
       
       Ti = liquidus_temperature_mush(salin / phi_init)
       trcrn(i,j,nt_qice+k-1,n) = &
            enthalpy_mush(Ti, salin)
       
    enddo ! k
    
    ! salinity
    do k = 1, nilyr
       
       trcrn(i,j,nt_sice+k-1,n) = salin
       
    enddo ! k
    
  end subroutine restart_simulation_notz_fieldwork
  
!=======================================================================
  
  subroutine load_notz_fieldwork_data()
    
    character(len=200), parameter :: filenotz = &
         '/Users/akt/Work/Development/NotzFieldwork/weather_notz_4m.txt'
    
    integer :: nt
    
    open(11,file=filenotz)
    
    do nt = 1, ndata_notz
       
       read(11,*) time_notz(nt), tair_notz(nt), relh_notz(nt), wind_notz(nt)
       
    enddo ! nt
    
    close(11)
    
    ! convert values
    time_notz = time_notz
    tair_notz = tair_notz + Tffresh        ! convert to Kelvin from deg C
    relh_notz = relh_notz / 100.0_dbl_kind ! convert from percent to fraction
    
    ! load era data
    call load_era_data()
    
  end subroutine load_notz_fieldwork_data

!=======================================================================

  function calc_notz_tair(time, tair_alternate) result(tair)
    
    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: tair_alternate
    real(kind=dbl_kind) :: tair
    
    tair = calc_notz_data(time, tair_notz, tair_alternate)
    
    !open(13,file='history/Tair.txt',position='append')
    !write(13,*) istep, istep1, time, tair
    !close(13)
    
  end function calc_notz_tair
  
!=======================================================================
  
  function calc_notz_relh(time, relh_alternate) result(relh)
    
    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: relh_alternate
    real(kind=dbl_kind) :: relh
    
    relh = calc_notz_data(time, relh_notz, relh_alternate)
    
  end function calc_notz_relh
  
!=======================================================================
  
  function calc_notz_wind(time, wind_alternate) result(wind)
    
    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: wind_alternate
    real(kind=dbl_kind) :: wind

    wind = calc_notz_data(time, wind_notz, wind_alternate)
    
  end function calc_notz_wind

!=======================================================================

  function calc_notz_Qa(time, Tair, Qa_alternate) result(Qa)
    
    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: Tair
    real(kind=dbl_kind), intent(in) :: Qa_alternate
    real(kind=dbl_kind) :: Qa
    
    real(kind=dbl_kind) :: relh
    
    relh = calc_notz_relh(time, -999.0_dbl_kind)
    
    if (relh == -999.0_dbl_kind) then
       
       Qa = Qa_alternate
       
    else
       
       Qa = specific_from_relative_humidity(relh, Tair)
       
    endif
    
  end function calc_notz_Qa

!=======================================================================

  function calc_notz_data(time, data_notz, data_alternate) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), dimension(1:ndata_notz), intent(in) :: data_notz
    real(kind=dbl_kind), intent(in) :: data_alternate
    real(kind=dbl_kind) :: data
    
    integer, parameter :: nt_start = 1
    
    real(kind=dbl_kind) :: x0, x1, x2, y1, y2
    
    integer :: nt
    
    do nt = nt_start, ndata_notz-1
       
       x0 = time
       
       if (x0 > time_notz(nt) .and. x0 <= time_notz(nt+1)) then
          
          x1 = time_notz(nt)
          x2 = time_notz(nt+1)
          
          y1 = data_notz(nt)
          y2 = data_notz(nt+1)
          
          data = ((x2 - x0) * y1 + (x0 - x1) * y2) / (x2 - x1)
          
          return
          
       endif
       
    enddo ! nt
    
    data = data_alternate
    
  end function calc_notz_data

!=======================================================================

  function sst_notz_fieldwork() result(sst)
    
    real(kind=dbl_kind) :: sst
    
    sst = -1.4_dbl_kind
    
    if (lrestart_notz1 .and. .not. lrestart_notz2) then
       
       ! first period
       sst = -1.4570137414141411_dbl_kind
       
    else if (lrestart_notz1 .and. lrestart_notz2) then

       ! second period
       sst = -1.3922273647727270_dbl_kind
       
    endif
    
  end function sst_notz_fieldwork

!=======================================================================

  function specific_from_relative_humidity(RH, Tair) result(Qa)

    real(kind=dbl_kind), intent(in) :: RH
    real(kind=dbl_kind), intent(in) :: Tair
    real(kind=dbl_kind) :: Qa
    
    real(kind=dbl_kind) :: Qa_sat
    
    Qa_sat = saturated_specific_humidity(Tair)
    
    Qa = (RH * Qa_sat) / (c1 - Qa_sat + RH * Qa_sat)

  end function specific_from_relative_humidity

!=======================================================================

  function saturated_specific_humidity(Tair) result(Qa_sat)
    
    real(kind=dbl_kind), intent(in) :: Tair
    real(kind=dbl_kind) :: Qa_sat
    
    real(kind=dbl_kind) :: T
    real(kind=dbl_kind) :: vp
    
    real(kind=dbl_kind), parameter :: a = 0.7859_dbl_kind
    real(kind=dbl_kind), parameter :: b = 0.03477_dbl_kind
    real(kind=dbl_kind), parameter :: c = 0.00412_dbl_kind
    real(kind=dbl_kind), parameter :: d = 0.00422_dbl_kind
    
    real(kind=dbl_kind), parameter :: Pa = 1.0e5_dbl_kind ! atmospheric pressure
    
    ! temp
    T = Tair - Tffresh
    
    ! vapour pressure
    vp = (a + b * T) / (c1 + c * T) + d * T ! log10(vp)
    
    vp = 10.0_dbl_kind**vp ! vp in mb
    
    vp = vp * 100.0_dbl_kind ! vp in Pa
    
    ! saturated specific humidity
    Qa_sat = (0.622_dbl_kind * vp) / (Pa - 0.378_dbl_kind * vp)
    
  end function saturated_specific_humidity
  
!=======================================================================

  subroutine load_era_data()
    
    character(len=200), parameter :: file_era_tcc = &
         '/Users/akt/Work/Development/AOMIP/ERA-Interim/notz_field_ERA_tcc.txt'
    
    character(len=200), parameter :: file_era_ssrd = &
         '/Users/akt/Work/Development/AOMIP/ERA-Interim/notz_field_ERA_ssrd.txt'    
    
    character(len=200), parameter :: file_era_strd = &
         '/Users/akt/Work/Development/AOMIP/ERA-Interim/notz_field_ERA_strd.txt'
    
    integer :: n
    
    real(kind=dbl_kind), parameter :: dt = 3.0_dbl_kind * 3600.0_dbl_kind
    
    ! tcc
    open(11,file=file_era_tcc)
    do n = 1, ndata_era
       read(11,*) time_era(n), tcc_era(n)
    enddo ! n
    close(11)
    
    ! ssrd
    open(11,file=file_era_ssrd)
    do n = 1, ndata_era
       read(11,*) time_era(n), ssrd_era(n)
    enddo ! n
    close(11)
    
    ! strd
    open(11,file=file_era_strd)
    do n = 1, ndata_era
       read(11,*) time_era(n), strd_era(n)
    enddo ! n
    close(11)
    
    ! scale ERA data
    ssrd_era = ssrd_era / dt
    strd_era = strd_era / dt
    
  end subroutine load_era_data

!=======================================================================

  function calc_era_tcc(time, tcc_alternate) result(tcc)
    
    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: tcc_alternate
    real(kind=dbl_kind) :: tcc
    
    tcc = calc_era_data(time, tcc_era, tcc_alternate)
    
  end function calc_era_tcc
  
!=======================================================================
  
  function calc_era_ssrd(time, ssrd_alternate) result(ssrd)
    
    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: ssrd_alternate
    real(kind=dbl_kind) :: ssrd
    
    ssrd = calc_era_data(time, ssrd_era, ssrd_alternate)
    
  end function calc_era_ssrd
  
!=======================================================================
  
  function calc_era_strd(time, strd_alternate) result(strd)
    
    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: strd_alternate
    real(kind=dbl_kind) :: strd
    
    strd = calc_era_data(time, strd_era, strd_alternate)
    
  end function calc_era_strd

!=======================================================================

  function calc_era_data(time, data_era, data_alternate) result(data)
    
    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), dimension(1:ndata_era), intent(in) :: data_era
    real(kind=dbl_kind), intent(in) :: data_alternate
    real(kind=dbl_kind) :: data
    
    integer, parameter :: nt_start = 1
    
    real(kind=dbl_kind) :: x0, x1, x2, y1, y2
    
    integer :: nt
    
    do nt = nt_start, ndata_era-1
       
       x0 = time
       
       if (x0 > time_era(nt) .and. x0 <= time_era(nt+1)) then
          
          x1 = time_era(nt)
          x2 = time_era(nt+1)
          
          y1 = data_era(nt)
          y2 = data_era(nt+1)
          
          data = ((x2 - x0) * y1 + (x0 - x1) * y2) / (x2 - x1)
          
          return
          
       endif
       
    enddo ! nt
    
    data = data_alternate
    
  end function calc_era_data
      
!=======================================================================

#endif
#if defined notz_experiment

  subroutine init_notz_experiment(time, nx_block, ny_block, vicen, vsnon, trcrn, phi_init, Tf)
    
    use ice_state, only: nt_Tsfc, nt_sice, nt_qice, nt_qsno
    
    real(kind=dbl_kind), intent(in) :: time
    integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block  ! block dimensions
    real(kind=dbl_kind), dimension(nx_block,ny_block,ncat), intent(inout) :: &
         vicen, vsnon
    real(kind=dbl_kind), dimension(nx_block,ny_block,max_ntrcr,ncat), intent(inout) :: &
         trcrn
    real(kind=dbl_kind), intent(in) :: phi_init
    real(kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: Tf

    integer :: k
    
    real(kind=dbl_kind) :: slope, Ti, salin
    
    integer, parameter :: i = 4
    integer, parameter :: j = 4
    integer, parameter :: n = 1

    salin = 34.0_dbl_kind
    
    trcrn(i,j,nt_Tsfc,n) = -10.0_dbl_kind
    
    do k = 1, nilyr
       
       Tf(i,j) = liquidus_temperature_mush(salin)
       
       slope = Tf(i,j) - trcrn(i,j,nt_Tsfc,n)
       Ti = trcrn(i,j,nt_Tsfc,n) &
            + slope*(real(k,kind=dbl_kind)-p5) &
            /real(nilyr,kind=dbl_kind)
       
       trcrn(i,j,nt_qice+k-1,n) = &
            enthalpy_mush(Ti, salin)
       
    enddo ! k

  end subroutine init_notz_experiment
#endif
!=======================================================================
#if defined flushing_notz
  subroutine init_flushing_notz(nx_block, ny_block, hpond, apond)

    integer, intent(in) :: nx_block, ny_block
    real(kind=dbl_kind), dimension(:,:), intent(inout) :: hpond
    real(kind=dbl_kind), dimension(:,:), intent(inout) :: apond

    integer :: i, j

    do j = 1, ny_block
       do i = 1, nx_block

          hpond(i,j) = 0.1_dbl_kind
          apond(i,j) = c1

       enddo ! i
    enddo ! j

  end subroutine init_flushing_notz
#endif
!=======================================================================
#if defined era-interim
#if defined snowice_maksym

  subroutine restart_simulation_snowice_maksym(time, nx_block, ny_block, vicen, vsnon, trcrn, phi_init, sss)

    use ice_state, only: nt_Tsfc, nt_sice, nt_qice, nt_qsno

    real(kind=dbl_kind), intent(in) :: time
    integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block  ! block dimensions
    real(kind=dbl_kind), dimension(nx_block,ny_block,ncat), intent(inout) :: &
         vicen, vsnon
    real(kind=dbl_kind), dimension(nx_block,ny_block,max_ntrcr,ncat), intent(inout) :: &
         trcrn
    real(kind=dbl_kind), intent(in) :: phi_init
    real(kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: sss

    integer, parameter :: i = 4
    integer, parameter :: j = 4
    integer, parameter :: n = 1

    real(kind=dbl_kind) :: Ti
    integer :: k

    if (time >= time_snowice_maksym .and. .not. lrestart_snowice_maksym) then

       write(*,*) "Maksym restart", istep, time, time_snowice_maksym, lrestart_snowice_maksym

       ! ice volume
       vicen(i,j,n) = 0.01_dbl_kind
       
       ! snow volume
       vsnon(i,j,n) = c0
       
       ! surface temperature
       trcrn(i,j,nt_Tsfc,n) = liquidus_temperature_mush(sss(i,j) / phi_init)
       write(*,*) liquidus_temperature_mush(sss(i,j) / phi_init)

       ! snow temperature
       do k = 1, nslyr
          
          Ti = min(c0, trcrn(i,j,nt_Tsfc,n))
          trcrn(i,j,nt_qsno+k-1,n) = -rhos*(Lfresh - cp_ice*Ti)
          
       enddo ! k
       
       ! ice temperature
       do k = 1, nilyr

          Ti = liquidus_temperature_mush(sss(i,j) / phi_init)
          trcrn(i,j,nt_qice+k-1,n) = &
               enthalpy_mush(Ti, sss(i,j))

          write(*,*) k, Ti, trcrn(i,j,nt_qice+k-1,n), sss(i,j)
          
       enddo ! k
       
       ! salinity
       do k = 1, nilyr
         
          trcrn(i,j,nt_sice+k-1,n) = sss(i,j)
          
       enddo ! k

       lrestart_snowice_maksym = .true.

       call restart_diagnostics(c0, 0.01_dbl_kind)

       !stop

    endif

  end subroutine restart_simulation_snowice_maksym
#endif
!=======================================================================

  subroutine restart_diagnostics(hsn, hin)

    real(kind=dbl_kind), intent(in) :: hsn
    real(kind=dbl_kind), intent(in) :: hin

    !call delete_diagnostic_file("./history/hin_gnuplot.txt")

    z_interface_snow = hsn
    z_interface_snic = c0
    z_interface_ices = c0
    z_interface_iceb = -hin

    !call delete_diagnostic_file("./history/Temp_gnuplot.txt")
    !call delete_diagnostic_file("./history/Sin_gnuplot.txt")

  end subroutine restart_diagnostics

!=======================================================================
  
  subroutine delete_diagnostic_file(filename)
    
    character(len=*), intent(in) :: filename
    
    open(11,file=filename)
    close(11,status='delete')
    
  end subroutine delete_diagnostic_file

!=======================================================================

  subroutine load_era_data()

#if defined snowice_maksym
    character(len=200), parameter :: t2m_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Possession/possession_t2m.txt"
    character(len=200), parameter :: d2m_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Possession/possession_d2m.txt"
    character(len=200), parameter :: u10_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Possession/possession_u10.txt"
    character(len=200), parameter :: v10_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Possession/possession_v10.txt"
    character(len=200), parameter :: tcc_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Possession/possession_tcc.txt"
    character(len=200), parameter :: tp_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Possession/possession_tp.txt"
    character(len=200), parameter :: sf_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Possession/possession_sf.txt"
#elif defined pond_barrow
    character(len=200), parameter :: t2m_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Barrow/barrow_t2m.txt"
    character(len=200), parameter :: d2m_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Barrow/barrow_d2m.txt"
    character(len=200), parameter :: u10_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Barrow/barrow_u10.txt"
    character(len=200), parameter :: v10_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Barrow/barrow_v10.txt"
    character(len=200), parameter :: tcc_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Barrow/barrow_tcc.txt"
    character(len=200), parameter :: tp_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Barrow/barrow_tp.txt"
    character(len=200), parameter :: sf_filename = &
         "/Users/akt/Work/Development/ERA-forcing/ERA-Interim/Barrow/barrow_sf.txt"
#endif

    call load_single_era_data(t2m_filename, ntimes1_era, times1_era, t2m_era)
    call load_single_era_data(d2m_filename, ntimes1_era, times1_era, d2m_era)
    call load_single_era_data(u10_filename, ntimes1_era, times1_era, u10_era)
    call load_single_era_data(v10_filename, ntimes1_era, times1_era, v10_era)
    call load_single_era_data(tcc_filename, ntimes1_era, times1_era, tcc_era)
    call load_single_era_data(tp_filename,  ntimes2_era, times2_era, tp_era)
    call load_single_era_data(sf_filename,  ntimes2_era, times2_era, sf_era)

    ! convert data
    call dew_point_to_specific_humidity(d2m_era)
    call convert_precip_era(tp_era)
    call convert_precip_era(sf_era)

  end subroutine load_era_data

!=======================================================================

  subroutine convert_precip_era(precip)

    real(kind=dbl_kind), dimension(:), intent(inout) :: precip

    real(kind=dbl_kind), parameter :: secs_in_period = 10800.0_dbl_kind
    real(kind=dbl_kind), parameter :: m_to_mm = 1000.0_dbl_kind

    precip = m_to_mm * (precip / secs_in_period)

  end subroutine convert_precip_era

!=======================================================================

  subroutine dew_point_to_specific_humidity(dew)

    real(kind=dbl_kind), dimension(:), intent(inout) :: dew

    integer :: i

    do i = 1, size(dew)

       dew(i) = (qqqocn/1.3_dbl_kind) * exp(-TTTocn/dew(i))
       dew(i) = max(dew(i),c0)

    enddo ! i

  end subroutine dew_point_to_specific_humidity

!=======================================================================

  subroutine load_single_era_data(filename, ntimes, times, data)
    
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ntimes
    real(kind=dbl_kind), dimension(ntimes), intent(out) :: times
    real(kind=dbl_kind), dimension(ntimes), intent(out) :: data

    integer :: ntime

    open(11,file=filename)

    do ntime = 1, ntimes

       read(11, *) times(ntime), data(ntime)

    enddo ! ntime

    close(11)

  end subroutine load_single_era_data

!=======================================================================

  function calc_era_t2m(time, data_alternate) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: data_alternate
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes1_era, times1_era, t2m_era, data_alternate)

  end function calc_era_t2m

!=======================================================================

  function calc_era_d2m(time, data_alternate) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: data_alternate
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes1_era, times1_era, d2m_era, data_alternate)

  end function calc_era_d2m

!=======================================================================

  function calc_era_u10(time, data_alternate) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: data_alternate
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes1_era, times1_era, u10_era, data_alternate)

  end function calc_era_u10

!=======================================================================

  function calc_era_v10(time, data_alternate) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: data_alternate
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes1_era, times1_era, v10_era, data_alternate)

  end function calc_era_v10

  !=======================================================================

  function calc_era_tcc(time, data_alternate) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: data_alternate
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes1_era, times1_era, tcc_era, data_alternate)

  end function calc_era_tcc

  !=======================================================================

  function calc_era_tp(time, data_alternate) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: data_alternate
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes2_era, times2_era, tp_era, data_alternate)

  end function calc_era_tp

  !=======================================================================

  function calc_era_sf(time, data_alternate) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind), intent(in) :: data_alternate
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes2_era, times2_era, sf_era, data_alternate)

  end function calc_era_sf

!=======================================================================

  function calc_era_data(time, ntimes, times, data_era, data_alternate) result(data)
    
    real(kind=dbl_kind), intent(in) :: time
    integer, intent(in) :: ntimes
    real(kind=dbl_kind), dimension(1:ntimes), intent(in) :: times
    real(kind=dbl_kind), dimension(1:ntimes), intent(in) :: data_era
    real(kind=dbl_kind), intent(in) :: data_alternate
    real(kind=dbl_kind) :: data
    
    integer, parameter :: nt_start = 1
    
    real(kind=dbl_kind) :: x0, x1, x2, y1, y2
    
    integer :: nt
    
    do nt = nt_start, ntimes-1
       
       x0 = time
       
       if (x0 > times(nt) .and. x0 <= times(nt+1)) then
          
          x1 = times(nt)
          x2 = times(nt+1)
          
          y1 = data_era(nt)
          y2 = data_era(nt+1)
          
          data = ((x2 - x0) * y1 + (x0 - x1) * y2) / (x2 - x1)
          
          return
          
       endif
       
    enddo ! nt
    
    data = data_alternate
    
  end function calc_era_data

!=======================================================================

  subroutine load_ocean_data()

#if defined snowice_maksym
    character(len=200), parameter :: filein = &
         "/Users/akt/Work/Development/gx1v3/gx1v3/forcing/ocean_possession.txt"
#elif defined pond_barrow
    character(len=200), parameter :: filein = &
         "/Users/akt/Work/Development/gx1v3/gx1v3/forcing/ocean_barrow.txt"
#endif

    integer :: ntime

    open(11,file=filein)

    do ntime = 1, ntimes_ocean
       
       read(11,*) times_ocean(ntime), &
                  S_ocean(ntime), &
                  T_ocean(ntime), &
                  U_ocean(ntime), &
                  V_ocean(ntime), &
                  dhdx_ocean(ntime), &
                  dhdy_ocean(ntime), &
                  hblt_ocean(ntime), &
                  qdp_ocean(ntime)

       !write(*,*) "oce: ", times_ocean(ntime), &
       !     S_ocean(ntime), &
       !     T_ocean(ntime), &
       !     U_ocean(ntime), &
       !     V_ocean(ntime), &
       !     dhdx_ocean(ntime), &
       !     dhdy_ocean(ntime), &
       !     hblt_ocean(ntime), &
       !     qdp_ocean(ntime)

    enddo ! ntime

    close(11)

  end subroutine load_ocean_data

!=======================================================================

  function calc_ocean_S(time) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes_ocean, times_ocean, S_ocean, c0)

  end function calc_ocean_S

!=======================================================================

  function calc_ocean_T(time) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes_ocean, times_ocean, T_ocean, c0)

  end function calc_ocean_T

!=======================================================================

  function calc_ocean_U(time) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes_ocean, times_ocean, U_ocean, c0)

  end function calc_ocean_U

!=======================================================================

  function calc_ocean_V(time) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes_ocean, times_ocean, V_ocean, c0)

  end function calc_ocean_V

!=======================================================================

  function calc_ocean_dhdx(time) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes_ocean, times_ocean, dhdx_ocean, c0)

  end function calc_ocean_dhdx

!=======================================================================

  function calc_ocean_dhdy(time) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes_ocean, times_ocean, dhdy_ocean, c0)

  end function calc_ocean_dhdy

!=======================================================================

  function calc_ocean_hblt(time) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes_ocean, times_ocean, hblt_ocean, c0)

  end function calc_ocean_hblt

!=======================================================================

  function calc_ocean_qdp(time) result(data)

    real(kind=dbl_kind), intent(in) :: time
    real(kind=dbl_kind) :: data

    data = calc_era_data(time, ntimes_ocean, times_ocean, qdp_ocean, c0)

  end function calc_ocean_qdp

!=======================================================================

#endif
!=======================================================================
! Oned diagnostic
!=======================================================================

    subroutine thermo_vertical_diag(nx_block, ny_block, indxi, indxj, icells, n, istep, istep1, dt, &
         aicen, vicen, vsnon, hin, hsn, Tsf, qin, qsn, Sin, Tin, Tsn, Tbot, fbot, &
         meltt, melts, meltb, congel, evapn, snoice, worki, works, fsnow, &
         fsurfn, fsensn, flatn, flw, flwoutn, fswsfc, fswint, shcoef, lhcoef, hpond, apond, potT, sss)

      integer (kind=int_kind), intent(in) :: &
           nx_block, ny_block, & ! block dimensions
           icells, n                ! number of cells with ice present
      
      integer (kind=int_kind), dimension (nx_block*ny_block), &
           intent(in) :: &
           indxi, indxj     ! compressed indices for cells with ice
      
      integer, intent(in) :: istep, istep1

      real(kind=dbl_kind), intent(in) :: dt

      ! variables
      real(kind=dbl_kind), dimension(icells), intent(in) :: hin, hsn, Tsf, worki, works

      real(kind=dbl_kind), dimension(icells,nilyr), intent(in) :: qin, Sin, Tin

      real(kind=dbl_kind), dimension(icells,nslyr), intent(in) :: qsn, Tsn

      real(kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: fbot, Tbot, meltt, melts, meltb, congel, evapn, snoice, &
           aicen, vicen, vsnon

      real(kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: fsnow, fsurfn, fsensn, flatn, flw, flwoutn, &
           fswsfc, fswint, shcoef, lhcoef, hpond, apond, potT, sss

#if diag == 1

#ifdef oned
      integer, parameter :: my_taskex = 0
      integer, parameter :: nex = 1
      integer, parameter :: iex = 4
      integer, parameter :: jex = 4
#else
      integer, parameter :: my_taskex = 23
      integer, parameter :: nex = 1
      integer, parameter :: iex = 5
      integer, parameter :: jex = 155
#endif

#ifdef notz_fieldwork
      !integer, parameter :: istep1_every = 10 ! notz fieldwork
      integer, parameter :: istep1_every = 1 ! notz fieldwork
#else
      integer, parameter :: istep1_every = 1
#endif

      integer :: i, j, k, ij

      real(kind=dbl_kind) :: dhi, dhs
      real(kind=dbl_kind), dimension(nilyr) :: data_nilyr
      real(kind=dbl_kind) :: z_ice_draft, phi
      real(kind=dbl_kind) :: perm_harm

      real(kind=dbl_kind), dimension(nilyr) :: Tin_local
      real(kind=dbl_kind), dimension(nslyr) :: Tsn_local

#ifdef notz_fieldwork
      real(kind=dbl_kind) :: time_field
      integer :: field_period
      character(len=200) :: filename
#endif

      if (my_task == my_taskex) then

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         
         if (n == nex .and. i == iex .and. j == jex) then

            if (ktherm == 2) then

               do k = 1, nilyr
                  Tin_local(k) = temperature_mush(qin(ij,k),Sin(ij,k))
               enddo ! k

               if (hsn(ij) > c0) then
                  do k = 1, nslyr
                     Tsn_local(k) = temperature_snow(qsn(ij,k))
                  enddo ! k
               else
                  Tsn_local(k) = c0
               endif

            else

               Tin_local = Tin(ij,:)
               Tsn_local = Tsn(ij,:)

            endif


            ! precalc - happens every timestep

            if (istep == 1) then

               z_interface_snow = hsn(ij)
               z_interface_snic = c0
               z_interface_ices = c0
               z_interface_iceb = -hin(ij)

            endif

            dhi = hin(ij) - worki(ij)
            dhs = hsn(ij) - works(ij)

            z_interface_snow = z_interface_snow + dhs + dhi + meltb(i,j) - congel(i,j)
            z_interface_snic = z_interface_snic + dhi + meltb(i,j) - congel(i,j)
            z_interface_ices = z_interface_ices + dhi + meltb(i,j) - congel(i,j) - snoice(i,j)
            z_interface_iceb = z_interface_iceb + meltb(i,j) - congel(i,j)

            z_ice_draft = c0
            do k = 1, nilyr
               
               phi = liquid_fraction(Tin_local(k), Sin(ij,k))
               
               z_ice_draft = z_ice_draft + (hin(ij) / real(nilyr, dbl_kind)) * (phi * rhow + (c1 - phi) * rhoi)
               
            enddo ! k
            
            !z_ice_draft = (z_ice_draft + hpond * apond * rhow + hsn * rhos) / rhow
            z_ice_draft = (z_ice_draft + hsn(ij) * rhos) / rhow
            !z_ice_draft = (hin(ij) * rhoi + hsn(ij) * rhos) / rhow

            z_interface_snow_float = - z_ice_draft + z_interface_snow - z_interface_iceb
            z_interface_snic_float = - z_ice_draft + z_interface_snic - z_interface_iceb
            z_interface_ices_float = - z_ice_draft + z_interface_ices - z_interface_iceb
            z_interface_iceb_float = - z_ice_draft

            ! output that must be reduced in output number
            if (mod(istep1,istep1_every) == 0 .or. istep1 == 0) then ! experiment

               ! hin
               open(44,file='./history/hin_gnuplot.txt',position='append')
               write(44,*) real(istep1,dbl_kind), hin(ij), hsn(ij), aicen(i,j), vicen(i,j), vsnon(i,j), &
                    z_interface_snow, z_interface_snic, z_interface_ices, z_interface_iceb, &
                    z_interface_snow_float, z_interface_snic_float, z_interface_ices_float, z_interface_iceb_float
               close(44)
               
#if defined notz_fieldwork
               if (lrestart_notz1 .and. .not. lrestart_notz2) then
                  open(44,file='history/hin_gnuplot_f1.txt',position='append')
                  write(44,*) istep1*dt - 5741712.0_dbl_kind, hin(ij), hsn(ij)
                  close(44)
               else if (lrestart_notz1 .and. lrestart_notz2) then
                  open(44,file='history/hin_gnuplot_f2.txt',position='append')
                  write(44,*) istep1*dt - 6270048.0_dbl_kind, hin(ij), hsn(ij)
                  close(44)
               endif
#endif

               ! pm3d gnuplots
               !open(50,file='./history/Temp_gnuplot.txt',position='append')
               !do k = 1, nilyr
               !   write(50,*) real(istep1_prev,dbl_kind), real(istep1,dbl_kind), &
               !        -hin(ij) * ((real(k,dbl_kind) - c1) / real(nilyr,dbl_kind)), &
               !        -hin(ij) * ((real(k,dbl_kind)) / real(nilyr,dbl_kind)), &
               !        Tin_local(k)
               !enddo ! k
               !do k = 1, nslyr
               !   write(50,*) real(istep1_prev,dbl_kind), real(istep1,dbl_kind), &
               !        hsn(ij) * ((real(nslyr-k+1,dbl_kind)) / real(nslyr,dbl_kind)), &
               !        hsn(ij) * ((real(nslyr-k,dbl_kind)) / real(nslyr,dbl_kind)), &
               !        Tsn_local(k)
               !enddo ! k
               !close(50)
               
               
               !open(50,file='./history/Sin_gnuplot.txt',position='append')
               !do k = 1, nilyr
               !   write(50,*) real(istep1_prev,dbl_kind), real(istep1,dbl_kind), &
               !        -hin(ij) * ((real(k,dbl_kind) - c1) / real(nilyr,dbl_kind)), &
               !        -hin(ij) * ((real(k,dbl_kind)) / real(nilyr,dbl_kind)), &
               !        Sin(ij,k)
               !enddo ! k
               !close(50)
               
               
               !open(50,file='./history/Sbr_gnuplot.txt',position='append')
               !do k = 1, nilyr
               !   write(50,*) real(istep1_prev,dbl_kind), real(istep1,dbl_kind), &
               !        -hin(ij) * ((real(k,dbl_kind) - c1) / real(nilyr,dbl_kind)), &
               !        -hin(ij) * ((real(k,dbl_kind)) / real(nilyr,dbl_kind)), &
               !        max(liquidus_brine_salinity_mush(Tin_local(k)), Sin(ij,k))
               !enddo ! k
               !close(50)
               
               !open(50,file='./history/phi_gnuplot.txt',position='append')
               !do k = 1, nilyr
               !   write(50,*) real(istep1_prev,dbl_kind), real(istep1,dbl_kind), &
               !        -hin(ij) * ((real(k,dbl_kind) - c1) / real(nilyr,dbl_kind)), &
               !        -hin(ij) * ((real(k,dbl_kind)) / real(nilyr,dbl_kind)), &
               !        min(Sin(ij,k) / max(liquidus_brine_salinity_mush(Tin_local(k)), Sin(ij,k)), c1)
               !enddo ! k
               !close(50)

               ! global temperature gradient
               !open(44,file='./history/global_grads.txt',position='append')
               !write(44,*) istep1, Tsf(ij) - Tbot(i,j), (Tsf(ij) - Tbot(i,j)) / hin(ij), &
               !     liquidus_brine_salinity_mush(Tsf(ij)) - liquidus_brine_salinity_mush(Tbot(i,j)), &
               !     (liquidus_brine_salinity_mush(Tsf(ij)) - liquidus_brine_salinity_mush(Tbot(i,j))) / hin(ij)
               !close(44)   

               ! details output
               open(44,file='./history/qin_det.txt',position='append')
               write(44,*) istep1, qin(ij,:)
               close(44)
               
               open(44,file='./history/Tin_det.txt',position='append')
               write(44,*) istep1, Tin_local(:), Tbot(i,j)
               close(44)

               open(44,file='./history/qsn_det.txt',position='append')
               write(44,*) istep1, qsn(ij,:)
               close(44)
               
               open(44,file='./history/Tsn_det.txt',position='append')
               write(44,*) istep1, Tsf(ij), Tsn_local(:)
               close(44)
               
               do k = 1, nilyr
                  data_nilyr(k) = liquidus_brine_salinity_mush(Tin_local(k))
               enddo ! k
               open(44,file='./history/Sbr_det.txt',position='append')
               write(44,*) istep1, data_nilyr(:)
               close(44)

               do k = 1, nilyr
                  data_nilyr(k) = liquid_fraction(Tin_local(k), Sin(ij,k))
               enddo ! k
               open(44,file='./history/phi_det.txt',position='append')
               write(44,*) istep1, data_nilyr(:)
               close(44)

               do k = 1, nilyr
                  data_nilyr(k) = liquidus_temperature_mush(Sin(ij,k))
               enddo ! k
               open(44,file='./history/Tmlt_det.txt',position='append')
               write(44,*) istep1, data_nilyr(:), Tin_local(:) - data_nilyr(:)
               close(44)

               open(44,file='./history/Sin_det.txt',position='append')
               write(44,*) istep1, Sin(ij,:)
               close(44)

               open(44,file='./history/meltgrow.txt',position='append')
               write(44,*) istep1, meltt(i,j), melts(i,j), meltb(i,j), congel(i,j), evapn(i,j), snoice(i,j), fsnow(i,j), fbot(i,j)
               close(44)

               !perm_harm = c0
               !do k = 1, nilyr
               !   phi = liquid_fraction(Tin_local(k), Sin(ij,k))
               !   data_nilyr(k) = 3.0e-8 * phi**3
               !   perm_harm = perm_harm + c1 / max(data_nilyr(k),1.0e-30_dbl_kind)
               !enddo ! k
               !perm_harm = real(nilyr,dbl_kind) / perm_harm
               !open(44,file='./history/perm_flood_det.txt',position='append')
               !write(44,*) istep1, data_nilyr(:), perm_harm
               !close(44)
  
               open(44,file='./history/surfflux.txt',position='append')
               write(44,*) istep1, fsurfn(i,j), fsensn(i,j), flatn(i,j), flw(i,j)*emissivity, &
                    flwoutn(i,j), fswsfc(i,j), fswint(i,j), shcoef(i,j), lhcoef(i,j)
               close(44)

               ! size of pipe K
              ! do k = 1, nilyr
              !    data_nilyr(k) = heat_conductivity(Tin_local(k), Sin(ij,k))
              ! enddo ! k
               !open(44,file='./history/cond_compare.txt',position='append')
               !write(44,*) istep1, maxval(g_k_pipe(1:nilyr) / data_nilyr(:)), &
               !     maxval(g_k_pipe) / maxval(data_nilyr), maxval(g_k_pipe), maxval(data_nilyr)
               !close(44)

               ! flushing
               open(44,file='./history/flushing.txt',position='append')
               write(44,*) istep1, hpond(i,j), apond(i,j)
               close(44)

               open(44,file='./history/potT.txt',position='append')
               write(44,*) istep1, potT(i,j)
               close(44)

               open(44,file='./history/sss.txt',position='append')
               write(44,*) istep1, sss(i,j)
               close(44)

!plot 'results/cloud_orig/surfflux.txt' using 1:2 with lines title "fsurfn", 'results/cloud_orig/surfflux.txt' using 1:3 with lines title "fsensn", 'results/cloud_orig/surfflux.txt' using 1:4 with lines title "flatn", 'results/cloud_orig/surfflux.txt' using 1:5 with lines title "flw", 'results/cloud_orig/surfflux.txt' using 1:6 with lines title "flwoutn", 'results/cloud_orig/surfflux.txt' using 1:7 with lines title "fswsfc"

!plot 'results/cloud_new/surfflux.txt' using 1:2 with lines title "fsurfn", 'results/cloud_new/surfflux.txt' using 1:3 with lines title "fsensn", 'results/cloud_new/surfflux.txt' using 1:4 with lines title "flatn", 'results/cloud_new/surfflux.txt' using 1:5 with lines title "flw", 'results/cloud_new/surfflux.txt' using 1:6 with lines title "flwoutn", 'results/cloud_new/surfflux.txt' using 1:7 with lines title "fswsfc", 'results/cloud_new/surfflux.txt' using 1:($3+$4+$5+$6+$7) with lines title "sum"

               ! update the prev istep1
               istep1_prev = istep1

            end if

#ifdef notz_experiment
            if (istep1 * dt == 4.0_dbl_kind * 3600.0_dbl_kind) then
               open(11,file='history/sin_profile_04.txt')
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 34.0_dbl_kind
               write(11,*) -c1, 34.0_dbl_kind
               close(11)
            endif
            if (istep1 * dt == 12.0_dbl_kind * 3600.0_dbl_kind) then
               open(11,file='history/sin_profile_12.txt')
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 34.0_dbl_kind
               write(11,*) -c1, 34.0_dbl_kind
               close(11)
            endif
            if (istep1 * dt == 24.0_dbl_kind * 3600.0_dbl_kind) then
               open(11,file='history/sin_profile_24.txt')
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 34.0_dbl_kind
               write(11,*) -c1, 34.0_dbl_kind
               close(11)
            endif
            if (istep1 * dt == 36.0_dbl_kind * 3600.0_dbl_kind) then
               open(11,file='history/sin_profile_36.txt')
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 34.0_dbl_kind
               write(11,*) -c1, 34.0_dbl_kind
               close(11)
            endif
            if (istep1 * dt == 48.0_dbl_kind * 3600.0_dbl_kind) then
               open(11,file='history/sin_profile_48.txt')
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 34.0_dbl_kind
               write(11,*) -c1, 34.0_dbl_kind
               close(11)
            endif
#endif

#ifdef notz_fieldwork
            if (lrestart_notz1 .and. .not. lrestart_notz2) then
               !time_field = istep1 * dt - time1_notz
               time_field = (istep1 - istep_notz1) * dt
               field_period = 1
            else if (lrestart_notz1 .and. lrestart_notz2) then
               !time_field = istep1 * dt - time2_notz
               time_field = (istep1 - istep_notz2) * dt
               field_period = 2
            else 
               time_field = -999.0_dbl_kind
               field_period = 0
            endif

            !write(*,*) istep1, istep1 * dt, time_field, 24.0_dbl_kind * 3600.0_dbl_kind, field_period

            if (time_field == 6.0_dbl_kind * 3600.0_dbl_kind) then
               write(filename,fmt='(a,i1,a,i3.3,a)') 'history/sin_profile_f',field_period,'_',6,'.txt'
               open(11,file=filename)
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 35.0_dbl_kind
               write(11,*) -c1, 35.0_dbl_kind
               close(11)
            endif
            if (time_field == 12.0_dbl_kind * 3600.0_dbl_kind) then
               write(filename,fmt='(a,i1,a,i3.3,a)') 'history/sin_profile_f',field_period,'_',12,'.txt'
               open(11,file=filename)
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 35.0_dbl_kind
               write(11,*) -c1, 35.0_dbl_kind
               close(11)
            endif
            if (time_field == 24.0_dbl_kind * 3600.0_dbl_kind) then
               write(filename,fmt='(a,i1,a,i3.3,a)') 'history/sin_profile_f',field_period,'_',24,'.txt'
               open(11,file=filename)
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 35.0_dbl_kind
               write(11,*) -c1, 35.0_dbl_kind
               close(11)
            endif
            if (time_field == 48.0_dbl_kind * 3600.0_dbl_kind) then
               write(filename,fmt='(a,i1,a,i3.3,a)') 'history/sin_profile_f',field_period,'_',48,'.txt'
               open(11,file=filename)
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 35.0_dbl_kind
               write(11,*) -c1, 35.0_dbl_kind
               close(11)
            endif
            if (time_field == 72.0_dbl_kind * 3600.0_dbl_kind) then
               write(filename,fmt='(a,i1,a,i3.3,a)') 'history/sin_profile_f',field_period,'_',72,'.txt'
               open(11,file=filename)
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 35.0_dbl_kind
               write(11,*) -c1, 35.0_dbl_kind
               close(11)
            endif
            if (time_field == 96.0_dbl_kind * 3600.0_dbl_kind) then
               write(filename,fmt='(a,i1,a,i3.3,a)') 'history/sin_profile_f',field_period,'_',96,'.txt'
               open(11,file=filename)
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 35.0_dbl_kind
               write(11,*) -c1, 35.0_dbl_kind
               close(11)
            endif
            if (time_field == 120.0_dbl_kind * 3600.0_dbl_kind) then
               write(filename,fmt='(a,i1,a,i3.3,a)') 'history/sin_profile_f',field_period,'_',120,'.txt'
               open(11,file=filename)
               do k = 1, nilyr
                  write(11,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               write(11,*) -hin(ij), 35.0_dbl_kind
               write(11,*) -c1, 35.0_dbl_kind
               close(11)
            endif
#endif 

#if defined flushing_notz

            if (istep1 == npt) then

               open(44,file='history/Sin_flush.txt',position='append')
               do k = 1, nilyr
                  write(44,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Sin(ij,k)
               enddo ! k
               close(44)

               open(44,file='history/phi_flush.txt',position='append')
               do k = 1, nilyr
                  write(44,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), &
                       liquid_fraction(Tin_local(k), Sin(ij,k))
               enddo ! k
               close(44)

               open(44,file='history/Tin_flush.txt',position='append')
               write(44,*) c0, Tsf(ij)
               do k = 1, nilyr
                  write(44,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), Tin_local(k)
               enddo ! k
               write(44,*) -hin(ij), Tbot(i,j)
               close(44)

               open(44,file='history/qin_flush.txt',position='append')
               do k = 1, nilyr
                  write(44,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), qin(ij,k)
               enddo ! k
               close(44)

               open(44,file='history/Sbr_flush.txt',position='append')
               do k = 1, nilyr
                  write(44,*) -hin(ij) * ((real(k,dbl_kind) - p5) / (real(nilyr,dbl_kind))), &
                       liquidus_brine_salinity_mush(Tin_local(k))
               enddo ! k
               close(44)

            endif

#endif

         endif

      enddo ! ij

#endif

      endif
            
    end subroutine thermo_vertical_diag

!=======================================================================

    subroutine diagnose_itd(nx_block, ny_block, aicen, vicen, vsnon, trcrn)

      use ice_state, only: nt_sice, nt_qice

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block

      real(kind=dbl_kind), dimension(nx_block,ny_block,ncat), intent(in) :: &
           aicen, vicen, vsnon

      real(kind=dbl_kind), dimension(nx_block,ny_block,max_ntrcr,ncat), intent(in) :: &
           trcrn

      integer, parameter :: my_taskex = 23
      integer, parameter :: nex = 1
      integer, parameter :: iex = 5
      integer, parameter :: jex = 155

      if (my_task == my_taskex) then
         
         open(44,file='history/itd_aicen.txt',position="append")
         write(44,*) istep1, aicen(iex,jex,:)
         close(44)

         open(44,file='history/itd_vicen.txt',position="append")
         write(44,*) istep1, vicen(iex,jex,:)
         close(44)

         open(44,file='history/itd_vsnon.txt',position="append")
         write(44,*) istep1, vsnon(iex,jex,:)
         close(44)

         open(44,file='history/itd_qin.txt',position="append")
         write(44,*) istep1, trcrn(iex,jex,nt_qice,:)
         close(44)

         open(44,file='history/itd_sin.txt',position="append")
         write(44,*) istep1, trcrn(iex,jex,nt_sice,:)
         close(44)

      endif

    end subroutine diagnose_itd

!=======================================================================

end module ice_therm_oned



