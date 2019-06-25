c $Id: ice_ocean.F,v 1.2 2001/02/23 21:06:55 schramm Exp $
c=======================================================================
!---! ocean boundary interface (stability based flux calculations)
!---! for uncoupled runs and for coupled runs
!---!
!---! author John Weatherly, C.M. Bitz,
!---!        Elizabeth C. Hunke, Bruce P. Briegleb
c=======================================================================

      module ice_ocean

      use ice_flux
      use ice_calendar
      use ice_constants

      implicit none

      real (kind=dbl_kind), parameter ::
     &   hmix = 20._dbl_kind     ! ocean mixed layer depth (m)
     &,  cphm = cp_ocn*rhow*hmix

      real (kind=dbl_kind) ::
     &   sst_mixedlayer(ilo:ihi,jlo:jhi)     ! save value of sst

c=======================================================================

      contains

c=======================================================================

      subroutine mixed_layer(i,j)

      use ice_grid
      use ice_atmo

      integer (kind=int_kind), intent(in) :: i,j
      real (kind=dbl_kind) ::
     &   delt   ! potential T difference   (K)
     &,  delq   ! humidity difference      (kg/kg)
     &,  rdn_in  ! initial value for rdn in stability routine
     &,  vmag   ! surface wind magnitude   (m/s)
     &,  strx   ! dummy argument
     &,  stry   ! dummy argument
     &,  Trf   ! dummy argument
     &,  dssqdt   ! dummy argument

      real (kind=dbl_kind), parameter ::
     &   umin =  c1          ! minimum wind speed (m/s)

        if (istep.gt.1) sst(i,j) = sst_mixedlayer(i,j)

        if (tmask(i,j)) then

         vmag   = max(umin, wind(i,j))      ! CSM Ocean wind-speed formula 
         rdn_in = sqrt(0.0027_dbl_kind/vmag 
     &     + .000142_dbl_kind + .0000764_dbl_kind*vmag) ! for rdn_in

         call stability( i,j,1, sst(i,j), rdn_in,
     $         strx   ,stry   ,Trf, dssqdt,  delt, delq)

         call ocnheat(i,j, delt, delq)  ! update sst, determine frzmlt

        endif

      end subroutine mixed_layer

c=======================================================================

      subroutine ocnheat(i,j, delt, delq)

!---!-------------------------------------------------------------------
!---! update sst and freezing/melting potential (atm/ocn forcing)
!---! sst changes due to ice fluxes are updated in ice_therm_driver.F
!---!-------------------------------------------------------------------

      use ice_state
      use ice_albedo

      integer (kind=int_kind) :: i,j

      real (kind=dbl_kind), intent(in) ::
     &   delt   ! potential T difference   (K)
     &,  delq   ! humidity difference      (kg/kg)

      real (kind=dbl_kind) ::
     &   TsfK   ! surface temperature in Kelvin (K)
     &,  fsh    ! sensible         heat flux  (W/m**2)
     &,  flh    ! latent           heat flux  (W/m**2)
     &,  fswabs ! srfc absrbd short-wave heat flux (W/m**2)
     &,  flwup  ! long-wave upward heat flux  (W/m**2)

      real (kind=dbl_kind) ::
     &   ft     ! fraction reduction of positive qdp
     &,  qdp    ! deep ocean heat flux

      ! specify as constant for now
      qdp = -c2 ! negative upward

      ! ocean surface temperature in Kelvin
      TsfK = sst(i,j)+Tffresh

      ! shortwave radiative flux
      fswabs = (c1 - albocn) * (swvdr(i,j)
     & + swvdf(i,j) + swidr(i,j) + swidf(i,j))

      ! longwave radiative flux
      flwup  = -emissivity*stefan_boltzmann * TsfK**4

      ! downward latent and sensible heat fluxes
      flh = lhcoef(i,j,1) * delq  ! using lhcoef(i,j,1), shcoef(i,j,1)
      fsh = shcoef(i,j,1) * delt  ! as temporary arrays for ocean

      ! first, compute sst change due to exchange with atm/ice above
      sst(i,j) = sst(i,j) +((fsh + flh + flwup + Flw(i,j) + fswabs)
     &            *(c1-aice(i,j)))*dt/cphm

      ! computed T change due to exchange with deep layers:
      sst(i,j) = sst(i,j) - qdp*dt/cphm

      ! compute potential to freeze or melt ice
      frzmlt(i,j) = (Tf(i,j)-sst(i,j))*cphm/dt          
      frzmlt(i,j) = min(max(frzmlt(i,j),-c1000),c1000)

      ! if sst is below freezing, reset sst to Tf
      if(sst(i,j) .le. Tf(i,j)) sst(i,j) = Tf(i,j)

      ! save mixed layer sst
      sst_mixedlayer(i,j) = sst(i,j)

      end subroutine ocnheat

c=======================================================================

      end module ice_ocean

c=======================================================================
