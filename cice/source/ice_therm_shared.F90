!=========================================================================
!BOP
!
! !MODULE: ice_therm_shared
!
! !DESCRIPTION:
!
! Shared thermo variables, subroutines
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! authors: Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      module ice_therm_shared
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size, only: ncat, nilyr, nslyr, max_ntrcr
!
!EOP
!
      implicit none
      save

      private
      public :: calculate_Tin_from_qin

      integer (kind=int_kind), public :: &
         ktherm          ! type of thermodynamics
                         ! 0 = 0-layer approximation
                         ! 1 = Bitz and Lipscomb 1999
                         ! 2 = mushy layer theory

      !real (kind=dbl_kind), dimension(nilyr+1) :: &
      !   salin        ! salinity (ppt)   

      real (kind=dbl_kind), dimension(nilyr+1), public :: &
         Tmlt            ! melting temp, -depressT * salinity
                         ! nilyr + 1 index is for bottom surface

      real (kind=dbl_kind), parameter, public :: &
         ferrmax = 1.0e-3_dbl_kind    ! max allowed energy flux error (W m-2)
                                      ! recommend ferrmax < 0.01 W m-2

      character (char_len), public :: &
         conduct         ! 'MU71' or 'bubbly'

      logical (kind=log_kind), public :: &
         l_brine         ! if true, treat brine pocket effects

      logical (kind=log_kind), public :: &
         heat_capacity, &! if true, ice has nonzero heat capacity
                         ! if false, use zero-layer thermodynamics
         calc_Tsfc,     &! if true, calculate surface temperature
                         ! if false, Tsfc is computed elsewhere and
                         ! atmos-ice fluxes are provided to CICE
         read_Sin,      &! if true, update salinity profile from file
         solve_Sin       ! if true, update salinity profile from solve_S_dt

      real (kind=dbl_kind), parameter, public :: &
         hfrazilmin = 0.05_dbl_kind ! min thickness of new frazil ice (m)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: calculate_Tin_from_qin  - calculate internal ice temperatures
!
! !DESCRIPTION:
!
!  Compute the internal ice temperatures from enthalpy using
!  quadratic formula
!
! !REVISION HISTORY:
!
! !INTERFACE:
!
      function calculate_Tin_from_qin (qin, Tmltk) &
               result(Tin)
!
! !USES:
!
      use ice_constants
!
! !INPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         qin   , &              ! enthalpy
         Tmltk                  ! melting temperature at one level
!
! !OUTPUT PARAMETERS
!
     real (kind=dbl_kind) :: &
         Tin                 ! internal temperature
!
!EOP
!
      real (kind=dbl_kind) :: &
         aa1,bb1,cc1         ! quadratic solvers

      if (l_brine) then
         aa1 = cp_ice
         bb1 = (cp_ocn-cp_ice)*Tmltk - qin/rhoi - Lfresh 
         cc1 = Lfresh * Tmltk
         Tin =  min((-bb1 - sqrt(bb1*bb1 - c4*aa1*cc1)) /  &
                         (c2*aa1),Tmltk)

      else                ! fresh ice
         Tin = (Lfresh + qin/rhoi) / cp_ice
      endif
 
      end function calculate_Tin_from_qin

!=======================================================================

      end module ice_therm_shared

!=======================================================================
