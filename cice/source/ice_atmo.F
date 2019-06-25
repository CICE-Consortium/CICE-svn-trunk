c $Id: $
c=======================================================================
!---! atmospheric boundary interface (stability based flux calculations)
!---!
!---! author Elizabeth C. Hunke, LANL
c=======================================================================

      module ice_atmo

      use ice_domain
      use ice_constants
      use ice_flux
      use ice_state

      implicit none

c=======================================================================

      contains

c=======================================================================

      subroutine stability( i,j,nc, Tsf, rdn_in,
     $         strx   ,stry   ,Trf, dssqdt,  delt, delq)

!---!-------------------------------------------------------------------
!---! Compute coefficients for atm/ice fluxes, stress, and reference 
!---! temperature
!---! NOTE: 
!---!  o all fluxes are positive downward
!---!  o net heat flux = fswabs + flwup + (1-emissivity)flwdn + fsh + flh
!---!  o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!---!  o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
!---!
!---! ASSUME:
!---!  o The saturation humidity of air at T(K): qsat(T)  (kg/m**3)
!---!
!---! code originally based on CSM1
!---! see Large and Pond (1981, 1982) ???
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: i,j,nc

      real (kind=dbl_kind), intent(in) ::
     &   Tsf    ! surface temperature of ice or ocean
     &,  rdn_in ! initial value for rdn, dependent on ice/ocean surface

      real (kind=dbl_kind), intent(out) ::
     &   strx   ! x surface stress (N)
     &,  stry   ! y surface stress (N)
     &,  Trf    ! reference height temperature  (K)
     &,  dssqdt ! derivative of ssq wrt Ti (kg/kg/K)
     &,  delt   ! potential T difference   (K)
     &,  delq   ! humidity difference      (kg/kg)
 
      integer (kind=int_kind) :: k      ! iteration index
      real (kind=dbl_kind) :: 
     &   TsfK   ! surface temperature in Kelvin (K)
     &,  thva   ! virtual temperature      (K)
     &,  stable ! stability factor
     &,  rdn    ! sqrt of neutral exchange coefficient (momentum)
     &,  rhn    ! sqrt of neutral exchange coefficient (heat)
     &,  ren    ! sqrt of neutral exchange coefficient (water)
     &,  hol    ! H (at zlvl  ) over L
     &,  xqq    ! temporary variable
     &,  psimh  ! stability function at zlvl   (momentum)
     &,  psixh  ! stability function at zlvl   (heat and water)
     &,  alz    ! ln(zlvl  /z10)
     &,  tau    ! stress at zlvl
     &,  bn     ! exchange coef funct for interpolation
     &,  bh     ! exchange coef funct for interpolation
     &,  fac    ! interpolation factor
     &,  ln0    ! log factor for interpolation
     &,  ln3    ! log factor for interpolation
     &,  ustar  ! ustar (m/s)
     &,  tstar  ! tstar
     &,  qstar  ! qstar
     &,  rd     ! sqrt of exchange coefficient (momentum)
     &,  re     ! sqrt of exchange coefficient (water)            
     &,  rh     ! sqrt of exchange coefficient (heat)
     &,  vmag   ! surface wind magnitude   (m/s)
     &,  ssq    ! sat surface humidity     (kg/kg)
     &,  cp     ! specific heat of moist air
     &,  psimhs ! stable profile 

      real (kind=dbl_kind), parameter ::
     &   cpvir = cpwv/cp_air - c1  ! Defined as cpwv/cp_air - 1.
     &,  zTrf  = c2          ! reference height for air temperature (m)
     &,  umin  = c1          ! minimum wind speed (m/s)
     &,  zvir  = 0.606_dbl_kind    ! rho_vapor/rho_air - 1.0
     &,  qqq   = 640380._dbl_kind  ! for qsat, dqsatdt
     &,  TTT   = 5107.4_dbl_kind   ! for qsat, dqsatdt

      ! local functions
      real (kind=dbl_kind) :: 
     &   Tk      ! temperature (K)
     &,  qsat    ! the saturation humididty of air (kg/m**3)
     &,  dqsatdt ! derivative of qsat wrt surface temperature
     &,  xd      ! dummy argument  
     &,  psimhu  ! unstable part of psimh
     &,  psixhu  ! unstable part of psimx

      qsat(Tk)    = qqq / exp(TTT/Tk)

      dqsatdt(Tk) = (TTT / Tk**2) * qqq / exp(TTT/Tk)

      psimhu(xd)  = log((c1+xd*(c2+xd))*(c1+xd*xd)/c8)
     $              - c2*atan(xd) + 1.571_dbl_kind

      psixhu(xd)  =  c2 * log((c1 + xd*xd)/c2)

      ! define some needed variables
      TsfK = Tsf+Tffresh                         !  surface temp (K)
      vmag   = max(umin, wind(i,j))
      thva   = potT(i,j) * (c1 + zvir * Qa(i,j)) ! virtual pot temp (K)
      ssq    = qsat   (TsfK) / rhoa(i,j)         ! sat surf hum (kg/kg)
      dssqdt = dqsatdt(TsfK) / rhoa(i,j)         ! deriv of ssq wrt Ti 
      delt   = potT(i,j) - TsfK                  ! pot temp diff (K)
      delq   = Qa(i,j) - ssq                     ! spec hum dif (kg/kg)
      alz    = log(zlvl(i,j)/zref) 
      cp     = cp_air*(c1 + cpvir*ssq)

      !------------------------------------------------------------
      ! first estimate of Z/L and ustar, tstar and qstar
      !------------------------------------------------------------

      ! neutral coefficients, z/L = 0.0 
      rdn = rdn_in
      rhn = rdn
      ren = rdn

      ! ustar,tstar,qstar
      ustar = rdn * vmag
      tstar = rhn * delt  
      qstar = ren * delq  

      !------------------------------------------------------------
      ! iterate to converge on Z/L, ustar, tstar and qstar
      !------------------------------------------------------------

c      do k=1,2             !!! ech thinks this should be at least 5
      do k=1,5

        ! compute stability & evaluate all stability functions 
        hol    = vonkar * gravit * zlvl(i,j)
     $           * (tstar/thva+qstar/(c1/zvir+Qa(i,j))) / ustar**2
        hol    = sign( min(abs(hol),c10), hol )
        stable = p5 + sign(p5 , hol)
        xqq    = max(sqrt(abs(c1 - c16*hol)) , c1)
        xqq    = sqrt(xqq)
c        psimh  = -c5*hol*stable + (c1-stable)*psimhu(xqq)
c        psixh  = -c5*hol*stable + (c1-stable)*psixhu(xqq)

        psimhs = -(0.7_dbl_kind*hol + 0.75_dbl_kind*(hol-14.3_dbl_kind)
     &    * exp(-0.35_dbl_kind*hol) + 10.7_dbl_kind) ! Jordan et al 1999

        psimh  = psimhs*stable + (c1-stable)*psimhu(xqq)
        psixh  = psimhs*stable + (c1-stable)*psixhu(xqq)

        ! shift all coeffs to measurement height and stability
        rd = rdn / (c1+rdn/vonkar*(alz-psimh))
        rh = rhn / (c1+rhn/vonkar*(alz-psixh))
        re = ren / (c1+ren/vonkar*(alz-psixh))

        ! update ustar, tstar, qstar using updated, shifted coeffs 
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 

      enddo    ! end iteration

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------
c      shcoef(i,j,nc) = rhoa(i,j)*ustar*cp  *rh
      ! add windless coefficient for sensible heat flux
      ! as in Jordan et al (JGR, 1999)
      shcoef(i,j,nc) = rhoa(i,j)*ustar*cp  *rh + c1
      lhcoef(i,j,nc) = rhoa(i,j)*ustar*Lsub*re

      !------------------------------------------------------------
      ! momentum flux
      !------------------------------------------------------------
      ! tau = rhoa(i,j) * ustar * ustar 
      ! strx = tau * uatm(i,j) / vmag 
      ! stry = tau * vatm(i,j) / vmag 
      !------------------------------------------------------------

      tau = rhoa(i,j) * ustar * rd    ! not the stress at zlvl(i,j)
      strx = tau * uatm(i,j) 
      stry = tau * vatm(i,j) 

      !------------------------------------------------------------
      ! reference temperature interpolation to height zTrf
      !------------------------------------------------------------
      ! Assume that 
      ! cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
      ! 1/sqrt(cn(i,j))=1/rdn and sqrt(cm(i,j))/ch(i,j)=1/rh 
      !------------------------------------------------------------
      bn = vonkar/rdn
      bh = vonkar/rh

      ! Interpolation factor for stable and unstable cases
      ln0 = log(c1 + (zTrf/zlvl(i,j))*(exp(bn) - c1))
      ln3 = log(c1 + (zTrf/zlvl(i,j))*(exp(bn - bh) - c1))
      fac = (ln0 - zTrf/zlvl(i,j)*(bn - bh))/bh * stable
     $    + (ln0 - ln3)/bh * (c1-stable)
      fac = min(max(fac,c0),c1)

      Trf = TsfK + (Tair(i,j) - TsfK)*fac

      end subroutine stability

c=======================================================================

      subroutine ice_sfc_flux(i,j,nc,Tsf,flwdabs,flwup,
     &    fsh,flh,dflhdT,dfshdT,dflwdT)

!---!-------------------------------------------------------------------
!---! compute ice-atm surface fluxes
!---!-------------------------------------------------------------------

      use ice_kinds_mod
      use ice_albedo

      integer (kind=int_kind), intent(in) :: i,j,nc

      real (kind=dbl_kind), intent(in) ::
     &     Tsf     ! ice/snow surface temperature   (C)

      real (kind=dbl_kind), intent(out) ::
     &     flwdabs   ! down long-wave  absorbed heat flx   (W/m**2)
     &,    flwup     ! emitted long-wave upward heat flux  (W/m**2)
     &,    fsh      ! sensible         heat flux  (W/m**2)
     &,    flh      ! latent           heat flux  (W/m**2)
     &,    dflhdT     ! d(flh)/d(T)      (W/m**2/K)
     &,    dfshdT     ! d(fsh)/d(T)      (W/m**2/K)
     &,    dflwdT    ! d(flwup)/d(T)     (W/m**2/K)

      real (kind=dbl_kind) ::
     &     TsfK     ! ice/snow surface temperature   (K)
     &,    Qsfc        ! saturated surface specific humidity

      real (kind=dbl_kind), parameter ::
     &   qqq   = 640380._dbl_kind  ! for qsat, dqsatdt
     &,  TTT   = 5107.4_dbl_kind   ! for qsat, dqsatdt

      ! local function
      real (kind=dbl_kind) :: 
     &   Tk      ! temperature (K)
     &,  dqsatdt ! derivative of qsat wrt surface temperature
     &,  qsat    ! the saturation humididty of air (kg/m**3)

      qsat(Tk)    = qqq / exp(TTT/Tk)
      dqsatdt(Tk) = (TTT / Tk**2) * qqq / exp(TTT/Tk)

      ! ice surface temperature in Kelvin
      TsfK = Tsf+Tffresh
      ! sat surf hum (kg/kg)
      Qsfc    = qsat(TsfK) / rhoa(i,j)

      ! longwave radiative flux
      flwdabs = emissivity*Flw(i,j)
      flwup   = -emissivity*stefan_boltzmann * TsfK**4

      ! downward latent and sensible heat fluxes
      fsh  = shcoef(i,j,nc)*(potT(i,j) - TsfK)
      flh  = lhcoef(i,j,nc)*(Qa(i,j) - Qsfc)

      ! derivatives wrt surface temp
      dflwdT = - emissivity*stefan_boltzmann * c4*TsfK**3 
      dfshdT = - shcoef(i,j,nc)
      dflhdT = - lhcoef(i,j,nc)*dqsatdt(TsfK)

      end subroutine ice_sfc_flux

c=======================================================================

      end module ice_atmo

c=======================================================================
