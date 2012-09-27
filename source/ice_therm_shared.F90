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
      use ice_domain_size, only: ncat, nilyr, nslyr, ntilyr, ntslyr, max_ntrcr
!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) :: &
         ktherm          ! type of thermodynamics
                         ! 0 = 0-layer approximation
                         ! 1 = Bitz and Lipscomb 1999
                         ! 2 = mushy layer theory

      real (kind=dbl_kind), dimension(nilyr+1) :: &
         salin       , & ! salinity (ppt)   
         Tmlt            ! melting temp, -depressT * salinity
                         ! nilyr + 1 index is for bottom surface

      real (kind=dbl_kind), parameter :: &
         ferrmax = 1.0e-3_dbl_kind    ! max allowed energy flux error (W m-2)
                                      ! recommend ferrmax < 0.01 W m-2

      character (char_len) :: &
         conduct         ! 'MU71' or 'bubbly'

      logical (kind=log_kind) :: &
         l_brine         ! if true, treat brine pocket effects

      logical (kind=log_kind) :: &
         heat_capacity, &! if true, ice has nonzero heat capacity
                         ! if false, use zero-layer thermodynamics
         calc_Tsfc       ! if true, calculate surface temperature
                         ! if false, Tsfc is computed elsewhere and
                         ! atmos-ice fluxes are provided to CICE

!=======================================================================

      end module ice_therm_shared

!=======================================================================
