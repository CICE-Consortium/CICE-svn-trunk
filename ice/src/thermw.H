c $Id: thermw.H,v 2.0 1999/12/09 22:25:36 eclare Exp $
c-----------------------------------------------------------------------
c.. Parameters and common blocks for the thermodynamics model
c-----------------------------------------------------------------------

c.. General
      real sigma          ! Stefan-Boltzmann constant (W/m^2 K^4)
      real himin          ! minimum thickness for ice (m)
      real Terrmax        ! temperature convergence precision (C) 
      parameter (sigma = 5.67e-8)
      parameter (himin = 0.01)
      parameter (Terrmax = 2.0e-3)

c.. Air
      real cpair          ! specific heat of air (J/kg C)
      real p0             ! surface pressure (mb)
      parameter(cpair = 1004.)
      parameter(p0 = 1013.)  

c.. Ice and Snow
      real kfresh         ! thermal cond. of pure ice (W/m C)
      real ksnow          ! thermal conductivity of snow (W/m C)
      real cfresh         ! specific heat of pure ice (J/kg C)
      real Lfresh         ! latent heat of melting of pure ice (J/kg)
      real i0             ! fraction of penetrating solar rad
      real ksolar         ! bulk solar abs coeff of sea ice (m-1)
      real Cheat          ! bulk heat transfer coeff for thick ice 
      real Tf0kel         ! Freezing temp of fresh ice in Kelvin = 273.15
      parameter (kfresh = 2.03)
      parameter (ksnow = 0.30)
      parameter (cfresh = 2106.)
      parameter (Lfresh = 3.34e5)
      parameter (i0 = 0.3)
      parameter (ksolar = 1.5)
      parameter (Cheat = 1.75e-3)
      parameter (Tf0kel = 273.15)

c.. Water
      real cwater         ! specific heat of water (J/kg C)
      real Lvap           ! latent heat of vaporizn at 0 C (J/kg)
      real Rvap           ! gas constant for H2O vapor (J/kg K)

      parameter (cwater = 4218.)
      parameter(Lvap = 2.5e+6)
      parameter(Rvap = 461.)  

c.. Combinations used for efficiency
      real kfresh2        ! 2*kfresh
      real kfresh4        ! 4*kfresh
      real rhoi_cfresh    ! rhoi*cfresh
      real rhos_cfresh    ! rhos*cfresh
      real rhoi_Lfresh    ! rhoi*Lfresh
      real Qcoef          ! another constant for latent heat flux
      real lvrrv          ! (Lvap + Lfresh) / Rvap

      parameter (kfresh2 = 2.*kfresh)
      parameter (kfresh4 = 4.*kfresh)
      parameter (rhoi_cfresh = rhoi*cfresh)
      parameter (rhos_cfresh = rhos*cfresh)
      parameter (rhoi_Lfresh = rhoi*Lfresh)
      parameter (Qcoef = 6.11*0.622/p0)
      parameter (lvrrv = (Lvap + Lfresh) / Rvap)

c.. Enthalpy
      real qsnow            ! snow enthalpy (J/kg)
      parameter (qsnow = Lfresh)

